#!/usr/bin/env python

import aiida
aiida.try_load_dbenv()
from aiida.common import constants
from aiida_create_solutesupercell_structures import *
from aiida.orm import Node
from aiida.orm import QueryBuilder
from aiida.orm import Calculation, Group, WorkCalculation
from aiida.orm import StructureData
from aiida.orm.nodes.data.array.trajectory import TrajectoryData
from aiida.orm.utils import load_node, WorkflowFactory
import click
from ase import Atoms
from ase.io import write as ase_write
import aiida_utils
import sys
import os
import numpy as np

#Define unit conversions
ANGSTROM_TO_BOHRRADIUS = 1./constants.bohr_to_ang
EV_PER_ANGSTROM_TO_HARTREE_PER_BOHRRADIUS = constants.bohr_to_ang/constants.hartree_to_ev
EV_TO_HARTREE = 1/constants.hartree_to_ev

def get_allnodes_fromgroup(group_name):
    qb = QueryBuilder()
    qb.append(Group, filters={'name': group_name}, tag='g')
    qb.append(Node, tag='job', with_group='g')
    all_nodes = [x[0] for x in qb.all()]
    return all_nodes

def get_outputcalcs(node):
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": node.uuid}, tag="worknode")
    q.append(WorkCalculation, tag="worknode2", output_of="worknode", project=["id", "ctime",  "*"])
    q.order_by({"worknode2": "ctime"})
    child_nodes = [x[2] for x in q.all()]
    return child_nodes


def write_runner_commentline(fileout, uuid, extra_comments={}):
    fileout.write("begin\ncomment ")
    fileout.write("uuid: {} ".format(uuid))
    for label in extra_comments:
        fileout.write("{}: {} ".format(label, extra_comments[label]))
    fileout.write("\n")
    return

def write_runner_cell(fileout, cell):
    cell = cell * ANGSTROM_TO_BOHRRADIUS
    for idx_cell, i in enumerate(cell):
        fileout.write("lattice %.10f %.10f %.10f\n" %
                      (cell[idx_cell][0], cell[idx_cell][1], cell[idx_cell][2]))
    return

def write_runner_atomlines(fileout, atomiccoord_array, elements, atomicforce_array=None):
    """
    Assumes input units of eV and angstrom
    """
    if atomicforce_array is None:
        atomicforce_array = np.zeros(atomiccoord_array.shape)

    atomiccoord_array = atomiccoord_array * ANGSTROM_TO_BOHRRADIUS
    atomicforce_array = atomicforce_array * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHRRADIUS

    for i in range(len(atomiccoord_array)):
        atomiccoord = atomiccoord_array[i]
        element = elements[i]
        atomicforce = atomicforce_array[i]

        fileout.write("atom   %.6f    %.6f   %.6f "
                      "%s  0.0   0.0  "
                      "%.10f  %.10f  %.10f\n" %
                      (atomiccoord[0], atomiccoord[1], atomiccoord[2],
                       element,
                       atomicforce[0], atomicforce[1], atomicforce[2]))
    return

def write_runner_finalline(fileout, energy=0, charge=0):
    fileout.write("energy %.15f\n" % (energy*EV_TO_HARTREE))
    fileout.write("charge %.15f\nend\n" % charge)
    return


def write_structure_torunner(fileout, structure_node, extra_comments={}):
    # get structure path, if applicable
    ase_structure = structure_node.get_ase()

    cell = ase_structure.get_cell()
    positions = ase_structure.get_positions()
    elements = ase_structure.get_chemical_symbols()

    write_runner_commentline(fileout, structure_node.uuid, extra_comments=extra_comments)
    write_runner_cell(fileout, cell)
    write_runner_atomlines(fileout, positions, elements)
    write_runner_finalline(fileout)
    return

def get_timesorted_basenodes(relaxworknode):
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": relaxworknode.uuid}, tag="relaxworknode")
    q.append(WorkCalculation, output_of="relaxworknode",
             project=["id", "ctime",  "*"],  tag="calc")
    q.order_by({"calc": "ctime"})
    timesorted_scf = [x[2] for x in q.all()]
    return timesorted_scf

def get_timesorted_scfs(worknode, relax_worknode=False):
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": worknode.uuid}, tag="worknode")
    output_tag = "worknode"
    if relax_worknode:
        output_tag = "worknode2"
        q.append(WorkCalculation, tag=output_tag, output_of="worknode")
    q.append(Calculation, output_of=output_tag, project=["id", "ctime",  "*"],  tag="calc")
    q.order_by({"calc": "ctime"})
    timesorted_scf = [x[2] for x in q.all()]
    return timesorted_scf

def write_pwbase_torunner(fileout, pwbasenode, extra_comments={}):
    scf_node = get_timesorted_scfs(pwbasenode)[-1]

    ase_structure = scf_node.inp.structure.get_ase()
    cell = ase_structure.get_cell()
    positions = ase_structure.get_positions()
    elements = ase_structure.get_chemical_symbols()

    atomicforce_array = scf_node.out.output_array.get_array('forces')[-1]
    try:
        atomicforce_array = scf_node.out.output_array.get_array('forces')[-1]
    except KeyError:
        print('Error obtaining forces for: {} skipping..'.format(pwbasenode))
        return
    energy = scf_node.out.output_parameters.get_attr('energy')

    write_runner_commentline(fileout, pwbasenode.uuid, extra_comments=extra_comments)
    write_runner_cell(fileout, cell)
    write_runner_atomlines(fileout, positions, elements, atomicforce_array=atomicforce_array)
    write_runner_finalline(fileout, energy=energy)
    return

def get_relaxnode_calcoutputs(node):
    from aiida_quantumespresso.calculations.pw import PwCalculation
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": node.uuid}, tag="worknode")
    q.append(WorkCalculation, tag="worknode2", output_of="worknode")
    q.append(PwCalculation, tag="pwcalc", output_of="worknode2", project=["id", "ctime",  "*"])
    q.order_by({"pwcalc": "ctime"})
    child_nodes = [x[2] for x in q.all()]
    return child_nodes

def get_timesorted_values(relax_node, arrayname, np_concatenate=True,
                          check_outputparams=False):
    # check_outputparams because some values are not added to trajectory if numsteps =1 (e.g energy)
    child_nodes = get_relaxnode_calcoutputs(relax_node)
    output_array = []
    num_steps = 0
    addextra_vcinfo = False # some vc-relax calcs omit final data in the trajectory
    for node in child_nodes:
        # vc-relax calcuations require some double-counting
        if node.inp.parameters.get_dict()['CONTROL']['calculation'] == 'vc-relax':
            addextra_vcinfo = True
        try:
            parser_warning = bool(node.out.output_parameters.get_dict()['parser_warnings'])
        except AttributeError:
            parser_warning = True
        if parser_warning != 0:
            print("Skipping failed child {} of {}".format(node, relax_node))
            continue
        try:
            num_steps = len(node.out.output_trajectory.get_array('forces'))
        except AttributeError:
            print("No trjactories in child {} of {}".format(node, relax_node))
            continue
        if num_steps == 1 and check_outputparams:
            output_array.append([node.out.output_parameters.get_dict()[arrayname]])
        else:
            output_array.append(node.out.output_trajectory.get_array(arrayname))

    # vc-relax nodes are really fussy when it comes to appending to trajectory data
    if relax_node.exit_status != 0:
       addextra_vcinfo = False
    if num_steps == 1:
       addextra_vcinfo = False
    if addextra_vcinfo and arrayname in ['steps', 'cells', 'positions']:
        output_array.append([output_array[-1][-1]])

    if np_concatenate:
        try:
            output_array = np.concatenate(output_array)
        except ValueError:
            return []
    return output_array

def write_pwrelax_torunner(fileout, relax_node, write_only_relaxed,
                           energy_tol, verbose,  extra_comments={}):
    timesorted_steps = get_timesorted_values(relax_node, 'steps')
    num_steps = len(timesorted_steps)
    if num_steps == 0:
        print("WARNING: {} is empty, skipping!".format(relax_node))
        return
    timesorted_cells = get_timesorted_values(relax_node, 'cells')
    timesorted_positions = get_timesorted_values(relax_node, 'positions')
    timesorted_energy = get_timesorted_values(relax_node, 'energy', check_outputparams=True)
    timesorted_forces = get_timesorted_values(relax_node, 'forces')
    timesorted_elements = get_timesorted_values(relax_node, 'symbols', np_concatenate=False)
    elements = timesorted_elements[0] #assume unchanging element positions

    #Sometimes the final forces are not parsed. Trim out the last energy in that case
    if relax_node.exit_status == 401:
        timesorted_energy = timesorted_energy[0:len(timesorted_forces)]

    if len(timesorted_cells) != num_steps:
       raise Exception("Cells does not have the proper number of steps!")
    if len(timesorted_positions) != num_steps:
       raise Exception("Positions does not have the proper number of steps!")
    if len(timesorted_energy) != num_steps:
       print(timesorted_forces, timesorted_energy, timesorted_steps)
       print(len(timesorted_energy), len(timesorted_steps), len(timesorted_forces))
       raise Exception("Energies does not have the proper number of steps!")
    if len(timesorted_forces) != num_steps:
       raise Exception("Forces does not have the proper number of steps!")

    extra_comments={"trajectory_step":None}

    if verbose:
       print('relaxation steps:',len(timesorted_steps))

    final_scf = bool(relax_node.inp.final_scf)
    if not write_only_relaxed:
        trajectory_looprange = list(range(len(timesorted_cells)))
    elif not final_scf:
        trajectory_looprange = [-1]
    else:
        trajectory_looprange = []
    final_loop = trajectory_looprange[-1]
    old_energy = timesorted_energy[0]
    for i in trajectory_looprange:
        del_e = np.abs(timesorted_energy[i] - old_energy)
        if i != 0 and i != final_loop and del_e  < energy_tol:
            print("Skipping step {} E_new: {} E_old:{} delE: {}".format(
                 i, timesorted_energy[i], old_energy, del_e))
            continue
        else:
            old_energy = timesorted_energy[i]
        extra_comments["trajectory_step"] = i
        write_runner_commentline(fileout, relax_node.uuid, extra_comments=extra_comments)
        write_runner_cell(fileout, timesorted_cells[i])
        write_runner_atomlines(fileout,
           timesorted_positions[i],
           elements,
           atomicforce_array=timesorted_forces[i])
        write_runner_finalline(fileout, energy=timesorted_energy[i])

    if final_scf:
        print('final')
        final_basenode = get_timesorted_basenodes(relax_node)[-1]
        extra_comments["trajectory_step"] = "final_scf"
        extra_comments["parent_uuid"] = relax_node.uuid
        write_pwbase_torunner(fileout, final_basenode, extra_comments=extra_comments)

    return


# show default values in click
orig_init = click.core.Option.__init__
def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True
click.core.Option.__init__ = new_init
# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-gn', '--group_name', default=None,
              type=str, help="Group to export identified by name")
@click.option('-f', '--filename', required=False, default=False,
         type=str, help="filename for outputfile, default = work_group+\".input.data\"")
@click.option('-wor', '--write_only_relaxed', required=False, default=False,
         is_flag=True, help="only write the final relaxed structure")
@click.option('-et', '--energy_tol', required=False, default=0.5,
         help="Only dumps relaxation steps of minimum energy_tol(eV) apart")
@click.option('-sreadme', '--supress_readme', is_flag=True,
         help="supresses the generation of a readme file")
@click.option('-v', '--verbose', is_flag=True,
         type=str, help="Enables verbosity")
@click.option('-oe', '--output_elements', required=False,
              help="only output structures containing these elements")


def createjob(group_name, filename, write_only_relaxed, energy_tol,
              supress_readme, verbose, output_elements):
    ''' e.g.
    ./aiida_export_group_to_runner.py -gn Al6xxxDB_structuregroup
    '''
    energy_tol = energy_tol*EV_TO_HARTREE
    all_nodes = get_allnodes_fromgroup(group_name)
    if not supress_readme:
        aiida_utils.create_READMEtxt()

    if output_elements:
        output_elements = prep_elementlist(output_elements)

    add_to_filename = "__all_steps"
    if write_only_relaxed == True:
        add_to_filename = "__only_relaxed"

    if filename == False:
        file = "aiida_exported_group_"+group_name+add_to_filename+".input.data"
        fileout = open(file, "w")
    else:
        file = filename
        fileout = open(filename, "w")

    if verbose:
        print('file           :', file)
        print('structure_group:', group_name)

    for node in all_nodes:
        exit_status = node.exit_status
        if int(exit_status) in [401]:
            print("WARNING {} had a non-critical non-zero exit!".format(node, exit_status))
        if int(exit_status) in [104]:
            print("WARNING {} had a critical non-zero exit, skipping!".format(node, exit_status))
            continue

        if verbose:
            print("Writing node: {}".format(node.uuid))
        if output_elements:
            input_ase = node.inp.structure.get_ase()
            only_output_elements = all([x in output_elements
                                        for x in input_ase.get_chemical_symbols()])
            if not only_output_elements:
                print("Skipping: {}/{}".format(input_ase, node))
                continue

        if isinstance(node, StructureData):
            print('using write_structure_torunner')
            write_structure_torunner(fileout, node)
        elif isinstance(node, WorkCalculation):
            process_label = node.get_attrs()['_process_label']
            if process_label == "PwBaseWorkChain":
                print('using write_pwbase_torunner')
                write_pwbase_torunner(fileout, node)
            elif process_label == "PwRelaxWorkChain":
                print('using write_pwrelax_torunner')
                write_pwrelax_torunner(fileout, node,write_only_relaxed,
                                       energy_tol,verbose)
            elif process_label == "ElasticWorkChain":
                print("recursively adding Elastic nodes")
                elastic_children = get_outputcalcs(node)
                all_nodes += elastic_children
            else:
                print("Could not identify node, skipping")
        else:
            print("Could not identify node, skipping")

    fileout.close()
    return


if __name__ == "__main__":
    createjob()
