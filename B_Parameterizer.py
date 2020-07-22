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


# ================ Input Files ================  #
topology = 'outputs/CoarseGrain/alanin_CoarseGrain.pdb'
trajectory = 'outputs/CoarseGrain/alanin_CoarseGrain.dcd.dcd'
# '/Users/forrestbicker/Documents/Code/Python/csi_projects/protines/twoBlobs-CG.pdb'

max_frame = 1000
stride = 1
block_count = 4


# ============= Pattern Generation ============= #
residue_list = ['ALA']  # list of ammino acids to be CoarseGrained

amino_acid_blueprint = { # TODO: autodetect number and pattern of segments given a name using config
    3: {  # 3 segments
        'Bond': [['21', '11'], ['11', 'B1'], ['B1', 'B2']],
        'Angle': [['21', '11', 'B1'], ['11', 'B1', 'B2'], ['B1', 'B2', '12']],
        'Dihedral': [['21', '11', 'B1', 'B2'], ['K11', 'B1', 'B2', '12'], ['B1', 'B2', '12', '22']]
    },
    2: {  # 2 segments
        'Bond': [['11', 'B1'], ['B1', 'B2']],
        'Angle': [['11', 'B1', 'B2'], ['B1', 'B2', '12']],
        'Dihedral': [['11', 'B1', 'B2', '12']]
    },
    1: {  # 1 segment
        'Bond': [['B1', 'B2']],
        'Angle': [['B1', 'B2', 'B3']],
        'Dihedral': [['B1', 'B2', 'B3', 'B4']]
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
            for resid in sel_resids:  # loops thru each resid
                for name_list in component_list:
                    # generates selection paramets
                    params = gen_params(abrev_dict[resname_key], name_list, mes_type, sel_resids, resid)
                    mes_name = '_'.join(params[5:].split())
                    atms = u.atoms.select_atoms(params)
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
    for frame in u.trajectory: # TODO: get traj splicing to work 
        f = frame.frame
        if start <= f < stop:  # alterantive to slicing trajectory, because slicing breaks MDAnalysis in strange ways
            if f % step == 0:
                for mes_name, atms in atms_dict.items():
                    mes_type = MesType(mes_name.count('_') + 1)
                    if mes_type != MesType.NULL:  # ensures that measurment exists
                        value = measure(mes_type, atms)
                        if value != None:
                            value_dict.setdefault(mes_name, []).append(value)  # safe append to dict
        elif f >= stop:
            print(colorify('32', f'Block {block_id} completed!'))
            return(value_dict)


def gen_params(resname_key, name_list, mes_type, sel_resids, resid):
    mes_type_list = ['Bond', 'Angle', 'Dihedral'] # TODO: create enum for mestypes
    params = 'name'
    for name in name_list:
        name = resname_key + name
        i = int(name[2:])
        # ensures the function only works on mes_types that actually exist
        bool_list = [mes_type == item and resid + i < max(sel_resids) for item in enumerate(mes_type_list)]
        if bool_list:
            name = name[:2] + str(resid + i)
        else:
            pass  # ignores non-existent measurements
        params += (f' {name}')
    return(params)


def measure(mes_type, atms):
    if mes_type == len(atms):
        if mes_type == MesType.BOND:  # bond
            return(atms.bond.length())
        elif mes_type == MesType.ANGLE:  # angle
            return(atms.angle.angle())
        elif mes_type == MesType.DIHEDRAL:  # dihedral
            return(atms.dihedral.value())


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
            print(f'- Measuring {mes_count} {resname_key} {mes_type}s in {res_count} residues over {frame_count} frames')
    except IndexError:
        pass


master_container_dict = {}
s_time = time.time()
output_dict_list = measure_all_connections(u, block_count, max_frame, stride)
output_dict_list = filter(None, output_dict_list)
exec_time = time.time() - s_time


for output_dict in output_dict_list:  # combines each block's output
    for mes_name, values in output_dict.items():

        if mes_name not in master_container_dict.keys():  # ensures container exists in dict
            master_container_dict[mes_name] = boandi.Container(mes_name)

        container = master_container_dict[mes_name]
        container.add_values(values)  # combining outputs


# ========= Measruement Data Generation ========== #
print('\nExporting {} measurement datasets to file...'.format(
    len(master_container_dict)))
for container in master_container_dict.values():  # loops thru each measurement
    filename = f'outputs/measurement_data/{container.name}.dat'
    with open(filename, 'w+') as instance_output:  # writes measurment list data to file
        # writes integer denoting mes_type to file
        instance_output.write(f'mes_type: {container.mes_type}\n')
        if container.mes_type == 0:
            str_values = [str(value) for value in container.values]
        else:
            str_values = [str(value * 3.14159 / 180) for value in container.values]
        instance_output.write('\n'.join(str_values))


# ============ File Object Cleanup ============  #
print('Outputs Written!\nTask Complete!')
print(f'Finished execution in {exec_time}s')
