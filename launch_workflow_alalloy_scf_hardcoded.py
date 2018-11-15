#!/usr/bin/env python
import click
from aiida.cmdline.params import options, arguments
from aiida.cmdline.params.types import DataParamType
from aiida.cmdline.utils import decorators
import sys

import aiida
aiida.load_dbenv()
from aiida.work.workfunctions import workfunction


def retrieve_alluncalculated_structures(structure_group_name,
                                        workchain_group_name=None):
    from aiida.orm.group import Group
    from aiida.orm.data.structure import StructureData
    from aiida.orm.calculation import WorkCalculation
    from aiida.orm.querybuilder import QueryBuilder

    sqb = QueryBuilder()
    sqb.append(Group, filters={'name': structure_group_name}, tag='g')
    sqb.append(StructureData, project='id', tag='s', member_of='g')
    sqb.append(WorkCalculation, tag='job', output_of='s')

    filters = {}
    if workchain_group_name:
        filters = {'name': workchain_group_name}
    sqb.append(Group, group_of='job', filters=filters)

    ids_dealt_with = [_ for _, in sqb.distinct().all()] or [-1]  # prevent empty list

    # # Now the main query:
    qb = QueryBuilder()
    qb.append(Group, filters={'name': structure_group_name}, tag='g')
    qb.append(StructureData, project='*', tag='s', member_of='g',
              filters={'id': {'!in': ids_dealt_with}})  # filter out calculated '!in' for not in

    res = [x[0] for x in qb.all()]
    return res


def get_numelectrons_structure_upffamily(structure, pseudos):

    def parse_numelectrons_upfpath(upfpath):
        for line in open(upfpath):
            if "valence" in line.lower() and "z" in line.lower():
                if len(line.split("=")) == 2:
                    num_e = int(float((line.split("=")[-1].strip().strip('"'))))
                elif len(line.split()) == 3:
                    num_e = int(float(line.split()[0]))
                else:
                    raise Exception("Could not parse {}".format(upfpath))
        return num_e

    def build_upf_numelectrons_dict(structure_ase, pseudos):
        element_nume_dict = {}
        for element in set(structure_ase.get_chemical_symbols()):
            upfpath = pseudos[element].get_file_abs_path()
            element_nume_dict[element] = parse_numelectrons_upfpath(upfpath)
        return element_nume_dict

    structure_ase = structure.get_ase()

    element_nume_dict = build_upf_numelectrons_dict(structure_ase, pseudos)

    num_e = 0
    for element in structure_ase.get_chemical_symbols():
        num_e += element_nume_dict[element]

    return num_e


def get_kmeshfrom_kptper_recipang(aiida_structure, kptper_recipang):
    import numpy as np

    ase_structure = aiida_structure.get_ase()
    reci_cell = ase_structure.get_reciprocal_cell()
    kmesh = [np.ceil(kptper_recipang * np.linalg.norm(reci_cell[i]))
             for i in range(len(reci_cell))]
    return kmesh


def get_nummachines(structure, a=2*(10**-4), b=2):
    # NOTE: used very adhoc guess for nodes, assuming quadratic scaling
    # NOTE: perhaps num_electrons may give a better estimate
    try:
        num_atoms = structure.get_extras()['num_atoms']
    except KeyError:
        num_atoms = len(structure.get_ase())
    numnodes = a*num_atoms**2. + b
    numnodes = max(round(numnodes/2)*2, b)  # force even # of nodes
    return numnodes


def get_nk(num_machines, code):
    def nk_nump_evenlydivisible(nk, nump):
        nk = float(nk)
        nump = float(nump)

        if round(nump/nk) == nump/nk:
            return True
        else:
            return False

    nk = str(max(4, int(num_machines/2)))  # adhoc guess

    # check if our choice is valid
    computer = code.get_computer()
    ppm = computer.get_default_mpiprocs_per_machine()
    nump = num_machines * ppm
    if not nk_nump_evenlydivisible(nk, nump):
        raise Exception("Error number processors: {} "
                        "is not divisible by nk: {}".format(nump, nk))

    return nk


@workfunction
def wf_getkpoints(aiida_structure, kptper_recipang):
    from aiida.orm.data.array.kpoints import KpointsData

    def get_kmeshfrom_kptper_recipang(aiida_structure, kptper_recipang):
        import numpy as np
        kptper_recipang = kptper_recipang.value

        ase_structure = aiida_structure.get_ase()
        reci_cell = ase_structure.get_reciprocal_cell()
        kmesh = [np.ceil(kptper_recipang * np.linalg.norm(reci_cell[i]))
                 for i in range(len(reci_cell))]
        return kmesh

    kpoints_mesh = get_kmeshfrom_kptper_recipang(aiida_structure, kptper_recipang)
    kpoints = KpointsData(kpoints_mesh=kpoints_mesh)

    return kpoints


@workfunction
def wf_setupparams(base_parameter, structure,
                   pseudo_familyname, nume2bnd_ratio):
        from aiida.orm.data.upf import get_pseudos_from_structure
        from aiida.orm.data.parameter import ParameterData

        pseudos = get_pseudos_from_structure(structure, pseudo_familyname.value)
        nelec = get_numelectrons_structure_upffamily(structure, pseudos)
        nbnd = nelec * nume2bnd_ratio.value
        nbnd = max(nbnd, 20) # minimum of 20 bands to avoid certain crashes
        parameter_dict = base_parameter.get_dict()
        parameter_dict['SYSTEM']['nbnd'] = nbnd
        parameters = ParameterData(dict=parameter_dict)

        return parameters


#@decorators.with_dbenv()
def launch():
    from aiida.orm.group import Group
    from aiida.orm.utils import load_node, WorkflowFactory
    from aiida.orm.data.base import Float, Int, Str
    from aiida.orm.data.parameter import ParameterData
    from aiida.work.launch import submit


    ################ HARDCODED PARAMETERS REPLACE WITH CLICK! ############
    code_node = 10402 #daint
    #code_node = 15716 #fidis
    code = load_node(code_node)
    #code = "pw-v6.3@daint"

    structure_group_name = "allatoms_fiverandom1"
    structure_group = Group.get_from_string(structure_group_name)
    # NOTE: in the real version we are likely to be passed groups, not names

    workchain_group_name = "allatoms_fiverandom1_calc1"
    workchain_group = Group.get_or_create(name=workchain_group_name)[0]
    # NOTE: in the real version we are likely to be passed groups, not names

    base_parameter_node = 15705
    base_parameter = load_node(base_parameter_node)

    pseudo_familyname = "SSSPefv1.1_alalloy"

    kptper_recipang = 80

    nume2bnd_ratio = 0.75

    run_debug = False

    max_wallclock_seconds=6*60*60 # Try to scale nodes s.t. we definitely finish in time
    ######################################################################

    # Load all the structures in the structure group, not-yet run in workchain_group_name
    uncalculated_structures = retrieve_alluncalculated_structures(
                                structure_group_name,
                                workchain_group_name=workchain_group_name
    )

    if len(uncalculated_structures) == 0:
        print("All structures in {} already have associated workchains in "
              "the group {}".format(structure_group_name, workchain_group_name))
        sys.exit()

    for structure in uncalculated_structures:

        # determine number of bands & setup the parameters
        parameters = wf_setupparams(base_parameter,
                                    structure,
                                    Str(pseudo_familyname),
                                    Float(nume2bnd_ratio))

        # determine kpoint mesh & setup kpoints
        kpoints = wf_getkpoints(structure, Int(kptper_recipang))

        # determine parallelization & resources (setup the settings & options)
        num_machines = get_nummachines(structure)
        options_dict = {
            'max_wallclock_seconds': max_wallclock_seconds,
            'resources': {'num_machines': num_machines}
        }
        if run_debug:
            num_machines = 2
            options_dict['resources']['num_machines'] = num_machines
            options_dict['max_wallclock_seconds'] = int(0.5*60*60)
            options_dict['queue_name'] = 'debug'
        workchain_options = ParameterData(dict=options_dict)

        nk = get_nk(num_machines, code)
        settings_dict = {'cmdline': ['-nk', nk]}
        settings = ParameterData(dict=settings_dict)

        # setup inputs & submit workchain
        inputs = {
            'code': code,
            'structure': structure,
            'pseudo_family': Str(pseudo_familyname),
            'kpoints': kpoints,
            'parameters': parameters,
            'options': workchain_options,
            'settings': settings
        }

        PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
        node = submit(PwBaseWorkChain, **inputs)
        if run_debug:
            sys.exit()
        workchain_group.add_nodes([node])
        print "WorkChain: {} submitted".format(node)


if __name__ == "__main__":
    launch()

