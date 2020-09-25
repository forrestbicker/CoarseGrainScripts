#! /Library/Frameworks/Python.framework/Versions/3.7/bin/python3
# ============================================== #
# ====== B: Coarse Grain Parameterization ====== #
# ============================================== #
# Written by Forrest Bicker
# August 2019
#

# O()

# ================ Dependencies ================ #
import numpy as np
import multiprocessing
import time

import MDAnalysis as mda
import BINAnalysis as boandi

from util import colorify
from json import load
from util import MesType
import math


# ================ Input Files ================  #
topology = 'outputs/CoarseGrain/alanin_CoarseGrain.pdb'
trajectory = 'outputs/CoarseGrain/alanin_CoarseGrain.dcd'
# '/Users/forrestbicker/Documents/Code/Python/csi_projects/protines/twoBlobs-CG.pdb'

max_frame = 1000
stride = 1
block_count = 3


# ============= Pattern Generation ============= #
residue_list = ['ALA']  # list of ammino acids to be CoarseGrained

# This is a dictionary of nested dictionaries containing specification of the structure of all desired measurements. A new mold will likely have to be added for each coarse grained simulation you wish to analyze, however, since the format acts as a mold rather an explicit mapping for each measurement, molds are highly reusable with little to no modification. Each new mold entry should be in the following format:
amino_acid_blueprint = {  # TODO: autodetect number and pattern of segments given a name using config
    4: {  # 4 segments
        MesType.BOND: [['30', '20'], ['20', '10'], ['10', 'B0'], ['B0', 'B1']],
        MesType.ANGLE: [['30', '20', '10'], ['20', '10', 'B0'], ['10', 'B0', 'B1'], ['B0', 'B1', '11']],
        MesType.DIHEDRAL: [['30', '20', '10', 'B0'], ['20', '10', 'B0', 'B1'], ['10', 'B0', 'B1', '11'], ['B0', 'B1', '11', '21']]
    },
    3: {  # 3 segments
        MesType.BOND: [['20', '10'], ['10', 'B0'], ['B0', 'B1']],
        MesType.ANGLE: [['20', '10', 'B0'], ['10', 'B0', 'B1'], ['B0', 'B1', '11']],
        MesType.DIHEDRAL: [['20', '10', 'B0', 'B1'], ['10', 'B0', 'B1', '11'], ['B0', 'B1', '11', '21']]
    },
    2: {  # 2 segments
        MesType.BOND: [['10', 'B0'], ['B0', 'B1']],
        MesType.ANGLE: [['10', 'B0', 'B1'], ['B0', 'B1', '11']],
        MesType.DIHEDRAL: [['10', 'B0', 'B1', '11']]
    },
    1: {  # 1 segment
        MesType.BOND: [['B0', 'B1']],
        MesType.ANGLE: [['B0', 'B1', 'B2']],
        MesType.DIHEDRAL: [['B0', 'B1', 'B2', 'B3']]
    }
}

with open('mapping_dict.json', "r") as f:
    mapping_dict = load(f)

with open('abrev_dict.json', "r") as f:
    abrev_dict = load(f)

# ========= Multiprocessing Functions =========  #
# master multiprocessing function
def measure_all_connections(u, block_count, max_frame, stride):
    atms_dict = {}
    for resname_key in residue_list:  # three nested loops to acsess the beads
        component_count = len(mapping_dict[resname_key].keys())
        resname_dict = amino_acid_blueprint[component_count]
        for mes_type, component_list in resname_dict.items():
            sel = u.atoms.select_atoms(f'resname {resname_key}')
            sel_resids = sel.residues.resids
            max_resid = max(sel_resids)
            for resid in sel_resids:  # loops thru each resid
                for name_list in component_list:  # generating info for specific measurement
                    # generates selection paramets
                    last_name = name_list[-1]
                    last_resid = resid + int(last_name[1:])
                    if last_resid <= max_resid:
                        params = gen_params(abrev_dict[resname_key], name_list, mes_type, max_resid, resid)
                        mes_name = '_'.join(params[5:].split())
                        atms = u.atoms.select_atoms(params)
                        if mes_type.value == len(atms):  # ensures measurement is valid
                            atms_dict[mes_name] = atms


    if block_count > 1:
        num_frames = len(u.trajectory[:max_frame])
        # determines size of each block
        block_size = int(np.ceil(num_frames / float(block_count)))

        # starting frame for each block
        starts = [n * block_size for n in range(block_count)]
        # ending frame for each block
        stops = [(n + 1) * block_size for n in range(block_count)]
        # stride length for each block (constant)
        strides = [stride for n in range(block_count)]
        block_id = [n for n in range(block_count)]  # numeric ID for each block
        atms_dicts = [atms_dict for n in range(block_count)]

        # zips block information into one iterable
        arg_iter = zip(starts, stops, strides, block_id, atms_dicts)

        with multiprocessing.Pool() as pool:  # pools blocks and runs in parallel
            # calls get_containers(params) for each block in arg_iter
            output_dict_list = pool.map(get_containers, arg_iter)

    elif block_count == 1:
        output_dict_list = get_containers([0, max_frame, stride, 0, atms_dict])

    else:
        raise ValueError(
            f'block_count must be greater than or equal to 1, but is {block_count}')

    return(output_dict_list)  # returns output from every block


# ============ Measurment Functions ============ #
def get_containers(arglist):
    # retrives block information from argument list
    start, stop, step, block_id, atms_dict = arglist

    print(f'Initiating Block {block_id} for frames {start} through {stop} with stride {step}')

    value_dict = {}
    for frame in u.trajectory:  # TODO: get traj splicing to work
        f = frame.frame
        if start <= f < stop:  # alterantive to slicing trajectory, because slicing breaks MDAnalysis in strange ways
            if f % step == 0:
                for mes_name, atms in atms_dict.items():
                    mes_type = MesType(mes_name.count('_') + 1)
                    if mes_type != MesType.NULL:  # ensures that measurment exists
                        value = measure(mes_type, atms)
                        if value != None:
                            if mes_name not in value_dict.keys():
                                value_dict[mes_name] = boandi.Container(mes_name, mes_type)
                            value_dict[mes_name].add_values([value])
        elif f >= stop:
            break
    print(colorify('32', f'Block {block_id} completed!'))
    return(value_dict)

def gen_params(resname_key, name_list, mes_type, max_resid, resid):
    params = 'name'
    for name in name_list:
        name = resname_key + name
        start_resid = int(name[2:])
        # ensures the function only works on mes_types that actually exist
        name = name[:2] + str(resid + start_resid)
        params += (f' {name}')
    return(params)


def measure(mes_type, atms):
    if mes_type.value == len(atms):
        if mes_type == MesType.BOND:  # bond
            return(atms.bond.length())
        elif mes_type == MesType.ANGLE:  # angle
            return(math.radians(atms.angle.angle()))
        elif mes_type == MesType.DIHEDRAL:  # dihedral
            return(math.radians(atms.dihedral.value()))
    else:
        raise Exception(f'{len(atms)} is an invalid number of atoms for measurement type {mes_type.name}')


# ============= Object Measurment =============  #
print('Generating Universe...')
u = mda.Universe(topology, trajectory)
print('Universe Generated!\nBegining Measurments:')


# logging information to screen
for resname_key in residue_list:  # three nested loops to acsess the beads
    component_count = len(mapping_dict[resname_key].keys())
    resname_dict = amino_acid_blueprint[component_count]
    try:
        for mes_type, connection_map in resname_dict.items():
            mes_count = len(connection_map)
            res_count = len(u.atoms.select_atoms(f'resname {resname_key}').residues)
            frame_count = len(u.trajectory[:max_frame]) + 1
            print(f'- Measuring {mes_count} {resname_key} {mes_type.name}s in {res_count} residues over {frame_count} frames')
    except IndexError:
        pass


s_time = time.time()
output_dict_list = measure_all_connections(u, block_count, max_frame, stride)
exec_time = time.time() - s_time
master_container_dict = output_dict_list.pop()

# NOTE: it is more efficent to split MEASURMEENTS across processors, not FRAMES
for output_dict in output_dict_list:
    for container in output_dict.values():
        master_container_dict[container.mes_name].add_values(container.values)  # combining outputs


# ========= Measruement Data Generation ========== #
print('\nExporting {} measurement datasets to file...'.format(
    len(master_container_dict)))
for container in master_container_dict.values():  # loops thru each measurement
    filename = f'outputs/measurement_data/{container.mes_name}.dat'
    with open(filename, 'w+') as instance_output:  # writes measurment list data to file
        # writes integer denoting mes_type to file
        unit = "angstrom" if container.mes_type.name == MesType.BOND else "radians"
        instance_output.write(f'mes_type: {container.mes_type.name} unit: {unit}\n')

        value_str = ""
        for value in container.values:
            value_str += value + '\n'
        instance_output.write(value_str)


# ============ File Object Cleanup ============  #
print('Outputs Written!\nTask Complete!')
print(f'Finished execution in {exec_time}s')