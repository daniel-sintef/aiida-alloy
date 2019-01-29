#!/usr/bin/env python
import click
import time
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

def retrieve_numactive_calculations():
    from aiida.orm.calculation import JobCalculation
    from aiida.orm.querybuilder import QueryBuilder
    qb = QueryBuilder()
    qb.append(JobCalculation,
              filters={'attributes.process_state':
                       {'!in': ['finished', 'excepted', 'killed']}}
    )
    return len(qb.all())


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


@click.command()
@click.option('-c', '--code_node', required=True,
              help="node of code to use")
@click.option('-sg', '--structure_group_name', required=True,
              help='input group of structures to submit workchains on')
@click.option('-wg', '--workchain_group_name', required=True,
              help='output group of workchains')
@click.option('-bp', '--base_parameter_node', required=True,
              help='node of base ParameterData to setup calculations')
@click.option('-pfn', '--pseudo_familyname', required=True,
              help='name of pseudopotential family to use')
@click.option('-kra', '--kptper_recipang', required=True,
              help='number of kpoints to use per reciprocal angstrom')
@click.option('-ber', '--nume2bnd_ratio', required=True,
              help='band to electron ratio')
@click.option('-cm', '--calc_method', default='scf',
              help='The calculation to perform, supported types are: scf, relax, vc-relax')
@click.option('-mws', '--max_wallclock_seconds', default=6*60*60,
              help='maximum wallclock time per job in seconds')
@click.option('-mac', '--max_active_calculations', default=300,
              help='maximum number of active calculations')
@click.option('-sli', '--sleep_interval', default=10*60,
              help='time to wait (sleep) between calculation submissions')
@click.option('-rdb', '--run_debug', is_flag=True,
              help='run the script in debug mode. Submits one structure only'
                   ' and does not attach the output to the workchain_group')
@click.option('-cwd', '--keep_workdir', is_flag=True,
              help='Keep the workdir files after running')
def launch(code_node, structure_group_name, workchain_group_name,
           base_parameter_node, pseudo_familyname, kptper_recipang,
           nume2bnd_ratio, calc_method, max_wallclock_seconds, max_active_calculations,
           sleep_interval, run_debug, keep_workdir):
    from aiida.orm.group import Group
    from aiida.orm.utils import load_node, WorkflowFactory
    from aiida.orm.data.base import Bool, Float, Int, Str
    from aiida.orm.data.parameter import ParameterData
    from aiida.work.launch import submit

    valid_calc_methods = ['scf', 'relax', 'vc-relax']
    if calc_method not in valid_calc_methods:
        raise Exception("Invalid calc_method: {}".format(calc_method))

    # setup parameters
    code = load_node(code_node)
    structure_group = Group.get_from_string(structure_group_name)
    workchain_group = Group.get_or_create(name=workchain_group_name)[0]
    base_parameter = load_node(base_parameter_node)
    # announce if running in debug mode
    if run_debug:
        print "Running in debug mode!"

    # Load all the structures in the structure group, not-yet run in workchain_group_name
    uncalculated_structures = retrieve_alluncalculated_structures(
                                structure_group_name,
                                workchain_group_name=workchain_group_name
    )
    if len(uncalculated_structures) == 0:
        print("All structures in {} already have associated workchains in "
              "the group {}".format(structure_group_name, workchain_group_name))
        sys.exit()

    # determine number of calculations to submit
    running_calculations = retrieve_numactive_calculations()
    calcs_to_submit = max_active_calculations - running_calculations


    # submit calculations
    for structure in uncalculated_structures:

        # ensure no more than the max number of calcs are submitted
        while (calcs_to_submit <= 0):
            running_calculations = retrieve_numactive_calculations()
            calcs_to_submit = max_active_calculations - running_calculations
            if calcs_to_submit <= 0:  # in case jobs finished during submission
                print("{} calcs running,"
                      "max num calcs {} waiting....".format(
                          running_calculations, max_active_calculations))
                time.sleep(sleep_interval)

        # start timer to inspect job submission times
        from timeit import default_timer as timer
        start = timer()

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
            'resources': {'num_machines': num_machines},
        }
        if run_debug:
            num_machines = 2
            options_dict['resources']['num_machines'] = num_machines
            options_dict['max_wallclock_seconds'] = int(30*60)
            options_dict['queue_name'] = 'debug'
        workchain_options = ParameterData(dict=options_dict)

        nk = get_nk(num_machines, code)
        settings_dict = {
            'cmdline': ['-nk', nk],
            'no_bands': True
            }
        settings = ParameterData(dict=settings_dict)

        # setup inputs & submit workchain
        clean_workdir = not keep_workdir
        inputs = {
                  'structure': structure,
                  'settings': settings,
                  'clean_workdir': Bool(clean_workdir)
                  }
        base_inputs = {
                'code': code,
                'pseudo_family': Str(pseudo_familyname),
                'kpoints': kpoints,
                'parameters': parameters,
                'options': workchain_options,
                'settings': settings,
                }
        if calc_method == 'scf':
            PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
            inputs.update(base_inputs)
        elif calc_method == 'relax':
            PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.relax')
            inputs['base'] = base_inputs
            inputs['relaxation_scheme'] = Str('relax')
            inputs['final_scf'] = Bool(False)
            inputs['meta_convergence'] = Bool(False)
        elif calc_method == 'vc-relax':
            PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.relax')
            inputs['base'] = base_inputs
            inputs['relaxation_scheme'] = Str('vc-relax')
            inputs['final_scf'] = Bool(True)
            inputs['meta_convergence'] = Bool(False)
        else:
            raise Exception("Invalid calc_method: {}".format(calc_method))

        node = submit(PwBaseWorkChain, **inputs)
        end = timer()
        time_elapsed = end - start
        print "WorkChain: {} submitted, took {}s".format(node, time_elapsed)

        if run_debug:
            sys.exit()

        workchain_group.add_nodes([node])

        calcs_to_submit -= 1



if __name__ == "__main__":
    launch()

