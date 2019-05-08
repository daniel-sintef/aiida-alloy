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

def get_qerelax_stress(workchain):
    '''
    Get the stress output of a PwRelaxWorkchain
    '''
    outputs_dict = workchain.get_outputs_dict()
    output_parameters = outputs_dict[u'output_parameters']

    output_parameters_dict = output_parameters.get_dict()
    stress = np.array(output_parameters_dict[u'stress'])
    return stress

def _get_deformed_structures(equilibrium_structure, strain_magnitudes, 
                             symmetric_strains_only=True):
    deformed_mat_set = DeformedStructureSet(equilibrium_structure, 
                                            norm_strains=strain_magnitudes, 
                                            shear_strains=strain_magnitudes)

    symmetry_operations_dict = {}
    deformations = deformed_mat_set.deformations
    if symmetric_strains_only:
        symmetry_operations_dict = symmetry_reduce(deformations,equilibrium_structure)
        deformations = [x for x in symmetry_operations_dict]

    deformed_structures = []
    for i in range(len(deformations)):
        mg_deformed_structure = deformations[i].apply_to_structure(equilibrium_structure)
        deformed_structure = StructureData(pymatgen_structure=mg_deformed_structure)
        deformed_structures.append(deformed_structure)

    return deformations, deformed_structures

def _fit_elastic_tensor(stresses, strains, strain_magnitudes, 
                        equilibrium_structure, equilibrium_stress,  
                        symmetric_strains_only=True):

    symm_eql_stresses = copy.deepcopy(stresses)
    symm_eql_strains = copy.deepcopy(strains)

    if symmetric_strains_only:
        all_deformations = _get_deformed_structures(equilibrium_structure, 
                                                    strain_magnitudes,
                                                    symmetric_strains_only=False)[0]
        symmetry_operations_dict = symmetry_reduce(all_deformations, equilibrium_structure)
        deformations = [deformation for deformation in symmetry_operations_dict]
        for i in range(len(deformations)):
            deformation = deformations[i]
            symmetry_operations = [x for x in symmetry_operations_dict[deformation]]
            for symm_op in symmetry_operations:
                symm_eql_strains.append(strains[i].transform(symm_op))
                symm_eql_stresses.append(stresses[i].transform(symm_op))

    # Fit the elastic constants
    compliance_tensor = ElasticTensor.from_independent_strains(
                        stresses=symm_eql_stresses,
                        strains=symm_eql_strains,
                        eq_stress=equilibrium_stress
    )
    compliance_tensor = -1.0 * compliance_tensor # pymatgen has opposite sign convention
    return symm_eql_stresses, symm_eql_strains, compliance_tensor

class ElasticWorkChain(WorkChain):

    @classmethod
    def define(cls, spec):
        super(ElasticWorkChain, cls).define(spec)

        spec.input('structure', valid_type=StructureData)
        spec.input('symmetric_strains_only', valid_type=Bool, default=Bool(True))
        spec.input('skip_input_relax', valid_type=Bool, default=Bool(False))
        spec.input('strain_magnitudes', valid_type=List,
                                        default=List(list=[-0.01,-0.005,0.005,0.01]))
        spec.expose_inputs(PwRelaxWorkChain, namespace='initial_relax',
                           exclude=('structure', 'clean_workdir')) 
        spec.expose_inputs(PwRelaxWorkChain, namespace='elastic_relax',
                           exclude=('structure', 'clean_workdir'))

        spec.outline(
           cls.relax_input_structure,
           cls.get_relaxed_structure_stress,
           cls.get_deformed_structures,
           cls.compute_deformed_structures,
           cls.gather_computed_stresses,
           cls.fit_elastic_tensor, #NOTE: may wish to add a check of elastic constant quality
           cls.set_outputs
        )

        spec.output('equilibrium_structure', valid_type=StructureData)
        spec.output('elastic_outputs', valid_type=ArrayData)
        spec.output('symmetry_mapping', valid_type=ParameterData)

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX',
                       message='one of the PwRelaxWorkChain subprocesses failed')
    
    def relax_input_structure(self):
        '''
        Run a relax/vc-relax calculation to find the ground state structure

        NOTE: I would really like to modify this to be flexible (LAMMPS or QE)
        '''
        inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain,
                                                   namespace='initial_relax'))
        inputs.structure = self.inputs.structure

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

            self.ctx.equilibrium_structure = workchain.out.output_structure
            self.ctx.equilibrium_stress = Stress(stress)

            ## Uncomment to control the input structure while debugging 
            #self.report('WARNING: wrong eql structure!')
            #self.ctx.equilibrium_structure = self.inputs.structure 

    def get_deformed_structures(self):
        """
        Determines the set of strains to be applied and generates the deformed structures
        """
        self.report('Getting deformed structures')
        strain_magnitudes = self.inputs.strain_magnitudes
        equilibrium_structure = self.ctx.equilibrium_structure.get_pymatgen_structure()
        symmetric_strains_only = self.inputs.symmetric_strains_only 

        deformations, deformed_structures = _get_deformed_structures(equilibrium_structure,
                                            strain_magnitudes, 
                                            symmetric_strains_only=symmetric_strains_only)

        strains = []
        for i in range(len(deformations)):
            strains.append(deformations[i].green_lagrange_strain)
        
        self.ctx.deformations = deformations
        self.ctx.deformed_structures = deformed_structures
        self.ctx.strains = strains

    def compute_deformed_structures(self):
        '''
        Run relax workflows for each deformed structure

        NOTE: I would really like to modify this to be flexible (LAMMPS or QE)
        '''
        self.report('Computing deformed structures')
        deformed_structures = self.ctx.deformed_structures

        for key_index in range(len(deformed_structures)):
            inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain,
                                                       namespace='elastic_relax'))
            inputs.structure = deformed_structures[key_index]


            future = self.submit(PwRelaxWorkChain, **inputs)

            deformation_key = 'deformation_{}'.format(key_index)
            self.report('launching PwRelaxWorkChain<{}>'.format(future.pk))
            self.to_context(**{deformation_key: future})

    def gather_computed_stresses(self):
        '''
        Retrieve final stress from the defomed structure relax workflows
        '''
        self.report('Gathering computed stresses')
        deformed_structures = self.ctx.deformed_structures

        computed_stresses = []
        for key_index in range(len(deformed_structures)):
            deformation_key = 'deformation_{}'.format(key_index)
            deformation_key = 'deformation_{}'.format(key_index)
            workchain = self.ctx[deformation_key]
            if not workchain.is_finished_ok:
                self.report('PwRelaxWorkChain failed with exit status {}'
                            .format(workchain.exit_status))
                return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX
            else:
                computed_stress = get_qerelax_stress(workchain)
                computed_stresses.append(Stress(computed_stress))
        
        self.ctx.stresses = computed_stresses 
    
    def fit_elastic_tensor(self):
        """
        Fit the elastic tensor to the computed stresses & strains
        """
        self.report('Fitting elastic tensor')
        stresses = self.ctx.stresses
        strains = self.ctx.strains
        strain_magnitudes = self.inputs.strain_magnitudes
        equilibrium_structure = self.ctx.equilibrium_structure.get_pymatgen_structure()
        equilibrium_stress = self.ctx.equilibrium_stress
        symmetric_strains_only = self.inputs.symmetric_strains_only        

        symm_eql_stresses, symm_eql_strains, compliance_tensor = _fit_elastic_tensor(
            stresses, strains, strain_magnitudes,  equilibrium_structure, equilibrium_stress, 
            symmetric_strains_only=symmetric_strains_only) 


        self.ctx.symm_eql_strains = symm_eql_strains
        self.ctx.symm_eql_stresses = symm_eql_stresses
        self.ctx.elastic_tensor = compliance_tensor.voigt

    def set_outputs(self):
        self.report('Setting Outputs')
        elastic_outputs = ArrayData()

        #An ugly ugly function to make the symmetry_operations_dict storable
        def make_symmopdict_aiidafriendly(symmopdict):
            aiida_symmopdict = dict((str(k).replace('.',','), [x.as_dict() for x in v])
                                    for k, v in symmopdict.iteritems()) 
            return aiida_symmopdict

        equilibrium_structure = self.ctx.equilibrium_structure.get_pymatgen_structure()
        symmetry_operations_dict = symmetry_reduce(self.ctx.deformations,
                                                        equilibrium_structure)
        symmetry_mapping = make_symmopdict_aiidafriendly(symmetry_operations_dict)
        symmetry_mapping = ParameterData(dict=symmetry_mapping)

        elastic_outputs.set_array('strains', np.array(self.ctx.strains))
        elastic_outputs.set_array('stresses', np.array(self.ctx.stresses))
        elastic_outputs.set_array('elastic_tensor',
                                 np.array(self.ctx.elastic_tensor))
        elastic_outputs.set_array("symm_eql_strains",
                                           np.array(self.ctx.symm_eql_strains))
        elastic_outputs.set_array("symm_eql_stresses",
                                           np.array(self.ctx.symm_eql_stresses))

        self.out('equilibrium_structure', self.ctx.equilibrium_structure)
        self.out('elastic_outputs', elastic_outputs)
        self.out('symmetry_mapping', symmetry_mapping)

