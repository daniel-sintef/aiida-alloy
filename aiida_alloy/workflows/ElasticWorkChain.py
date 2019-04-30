import copy
import pymatgen as mg
from pymatgen.analysis.elasticity import DeformedStructureSet
from pymatgen.analysis.elasticity.tensors import symmetry_reduce
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity.elastic import ElasticTensor
import numpy as np
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.data.base import Str, Bool, List
from aiida.orm.data.array import ArrayData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.parameter import ParameterData
from aiida.work.workchain import WorkChain, ToContext, if_, append_
from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain
from aiida.work.workfunctions import workfunction

@workfunction
def get_conventional_structure(structuredata):
    '''
    Standardize an AiiDA StructureData object via pymatgen Structure
        using spglib
    :param structuredata: original StructureData
    '''
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    mg_structure = structuredata.get_pymatgen()
    sga = SpacegroupAnalyzer(mg_structure)
    standard_structure = sga.get_conventional_standard_structure()
    standard_structuredata = StructureData(pymatgen_structure=standard_structure)

    return standard_structuredata

def get_qerelax_stress(workchain):
    '''
    Get the stress output of a PwRelaxWorkchain
    '''
    outputs_dict = workchain.get_outputs_dict()
    output_parameters = outputs_dict[u'output_parameters']

    output_parameters_dict = output_parameters.get_dict()
    stress = np.array(output_parameters_dict[u'stress'])
    return stress

class ElasticWorkChain(WorkChain):

    @classmethod
    def define(cls, spec):
        super(ElasticWorkChain, cls).define(spec)

        spec.input('structure', valid_type=StructureData)
        spec.input('symmetric_strains_only', valid_type=Bool, default=Bool(True))
        spec.input('strain_magnitudes', valid_type=List,
                                        default=List(list=[-0.01,-0.005,0.005,0.01]))
        spec.expose_inputs(PwRelaxWorkChain, namespace='initial_relax',
                           exclude=('structure', 'clean_workdir')) 
        spec.expose_inputs(PwRelaxWorkChain, namespace='elastic_relax',
                           exclude=('structure', 'clean_workdir'))

        spec.outline(
           cls.get_conventional_structure,
           cls.relax_conventional_structure,
           cls.get_relaxed_structure_stress,
           cls.get_deformed_structures,
           cls.compute_deformed_structures,
           cls.gather_computed_stresses,
           cls.fit_elastic_tensor, #NOTE: may wish to add a check of elastic constant quality
           cls.set_outputs
        )

        spec.output('relaxed_conventional_structure', valid_type=StructureData)
        spec.output('elastic_outputs', valid_type=ArrayData)
        spec.output('symmetry_mapping', valid_type=ParameterData)

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX',
                       message='one of the PwRelaxWorkChain subprocesses failed')
    
    def get_conventional_structure(self):
        '''
        Standardize the input structure and set the initial structure
        '''
        
        structure = self.inputs.structure
        conventional_structure = get_conventional_structure(structure)
        self.ctx.conventional_structure = conventional_structure

    def relax_conventional_structure(self):
        '''
        Run a relax/vc-relax calculation to find the ground state structure

        NOTE: I would really like to modify this to be flexible (LAMMPS or QE)
        '''
        inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain,
                                                   namespace='initial_relax'))
        inputs.structure = self.ctx.conventional_structure

        future = self.submit(PwRelaxWorkChain, **inputs)
        self.report('Launching PwRelaxWorkChain<{}>'.format(future.pk))
        self.to_context(initial_relax_workchain=future)

    def get_relaxed_structure_stress(self):
        '''
        Verify that the PwRelaxWorkChain finished successfully

        NOTE: I would really like to modify this to be flexible (LAMMPS or QE)
        '''
        workchain = self.ctx.initial_relax_workchain
        self.report('Getting relaxed structure stresses')

        if not workchain.is_finished_ok:
            self.report('PwRelaxWorkChain<{}> failed with exit status {}'
                        .format(workchain.pk, workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX
        else:
            outputs_dict = workchain.get_outputs_dict()
            output_parameters_dict = outputs_dict[u'output_parameters'].get_dict()
            stress = np.array(output_parameters_dict['stress'])

            self.ctx.ground_state_structure = workchain.out.output_structure
            self.ctx.ground_state_stress = Stress(stress)
        self.report('Finished getting relaxed structure stresses')

    def get_deformed_structures(self):
        """
        Determines the set of strains to be applied and generates the deformed structures
        """
        self.report('Getting deformed structures')
        strain_magnitudes = self.inputs.strain_magnitudes
        structure_mat = self.ctx.ground_state_structure.get_pymatgen_structure()
        deformed_mat_set = DeformedStructureSet(structure_mat, 
                                                norm_strains=strain_magnitudes, 
                                                shear_strains=strain_magnitudes)

        self.ctx.symmetry_operations_dict = {}
        self.ctx.deformations = deformed_mat_set.deformations
        if self.inputs.symmetric_strains_only:
            symmetry_operations_dict = symmetry_reduce(self.ctx.deformations,
                                                                structure_mat)
            self.ctx.deformations = [x for x in symmetry_operations_dict]
        
        self.ctx.deformed_structures = []
        self.ctx.strains = []
        for i in range(len(self.ctx.deformations)):
            mg_deformed_structure = self.ctx.deformations[i].apply_to_structure(structure_mat)
            deformed_structure = StructureData(pymatgen_structure=mg_deformed_structure)
            self.ctx.deformed_structures.append(deformed_structure)
            self.ctx.strains.append(self.ctx.deformations[i].green_lagrange_strain)
        self.report('Finished getting deformed structures')



    def compute_deformed_structures(self):
        '''
        Run relax workflows for each deformed structure

        NOTE: I would really like to modify this to be flexible (LAMMPS or QE)
        '''
        self.report('Computing deformed structures')
        deformed_structures = self.ctx.deformed_structures

        for deformed_structure in deformed_structures:
            inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain,
                                                       namespace='elastic_relax'))
            inputs.structure = deformed_structure

            future = self.submit(PwRelaxWorkChain, **inputs)
            self.report('launching PwRelaxWorkChain<{}>'.format(future.pk))
            self.to_context(deformed_workchains=append_(future))

    def gather_computed_stresses(self):
        '''
        Retrieve final stress from the defomed structure relax workflows
        '''
        self.report('Gathering computed stresses')
        deformed_workchains = self.ctx.deformed_workchains

        computed_stresses = []
        for workchain in deformed_workchains:
            if not workchain.is_finished_ok:
                self.report('PwRelaxWorkChain failed with exit status {}'
                            .format(workchain.exit_status))
                return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX
            else:
                computed_stress = get_qerelax_stress(workchain)
                computed_stresses.append(computed_stress)
        
        pymatgen_stresses = []
        for computed_stress in computed_stresses:
            pymatgen_stresses.append(Stress(computed_stress))
        
        self.ctx.stresses = pymatgen_stresses 
    
    def fit_elastic_tensor(self):
        """
        Fit the elastic tensor to the computed stresses & strains
        """
        self.report('Fitting elastic tensor')
        symm_equivalent_strains = copy.deepcopy(self.ctx.strains)
        symm_equivalent_stresses = copy.deepcopy(self.ctx.stresses)

        # Add in the symmetrically-equivalent stresses & strains for the fitting
        for i in range(len(self.ctx.deformations)):
            deformation = self.ctx.deformations[i]
            symmetry_operations = []
            if self.inputs.symmetric_strains_only:
                structure_mat = self.ctx.ground_state_structure.get_pymatgen_structure()
                symmetry_operations_dict = symmetry_reduce(self.ctx.deformations,
                                                                structure_mat)
                symmetry_operations = [x for x in symmetry_operations_dict[deformation]]
            for symm_op in symmetry_operations:
                symm_equivalent_strains.append(self.ctx.strains[i].transform(symm_op))
                symm_equivalent_stresses.append(self.ctx.stresses[i].transform(symm_op))
        
        # Now fit the elastic constants
        compliance_tensor = ElasticTensor.from_independent_strains(
                            stresses=symm_equivalent_stresses,
                            strains=symm_equivalent_strains,
                            eq_stress=self.ctx.ground_state_stress
        )
        compliance_tensor = -1.0 * compliance_tensor # pymatgen has opposite sign convention

        self.ctx.symm_equivalent_strains = symm_equivalent_strains
        self.ctx.symm_equivalent_stresses = symm_equivalent_stresses
        self.ctx.elastic_tensor = compliance_tensor.voigt


    def set_outputs(self):
        self.report('Setting Outputs')
        elastic_outputs = ArrayData()

        #An ugly ugly function to make the symmetry_operations_dict storable
        def make_symmopdict_aiidafriendly(symmopdict):
            aiida_symmopdict = dict((str(k).replace('.',','), [x.as_dict() for x in v])
                                    for k, v in symmopdict.iteritems()) 
            return aiida_symmopdict

        structure_mat = self.ctx.ground_state_structure.get_pymatgen_structure()
        symmetry_operations_dict = symmetry_reduce(self.ctx.deformations,
                                                        structure_mat)
        symmetry_mapping = make_symmopdict_aiidafriendly(symmetry_operations_dict)
        symmetry_mapping = ParameterData(dict=symmetry_mapping)

        elastic_outputs.set_array('strains', np.array(self.ctx.strains))
        elastic_outputs.set_array('stresses', np.array(self.ctx.stresses))
        elastic_outputs.set_array('elastic_tensor',
                                 np.array(self.ctx.elastic_tensor))
        elastic_outputs.set_array("symm_equivalent_strains",
                                           np.array(self.ctx.symm_equivalent_strains))
        elastic_outputs.set_array("symm_equivalent_stresses",
                                           np.array(self.ctx.symm_equivalent_stresses))

        self.out('ground_state_structure', self.ctx.ground_state_structure)
        self.out('elastic_outputs', elastic_outputs)
        self.out('symmetry_mapping', symmetry_mapping)

