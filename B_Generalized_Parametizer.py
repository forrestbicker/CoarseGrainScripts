#! /Library/Frameworks/Python.framework/Versions/3.7/bin/python3
# ============================================== #
# ====== B: Generalized Coarse Grain Parameterization ====== #
# ============================================== #
# Written by Forrest Bicker
# September 2020
#

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

max_frame = 1000
stride = 1
block_count = 3
# ================ Input Files ================  #
topology = '/Users/forrestbicker/Documents/GitHub/CoarseGrainScrips/outputs/CoarseGrain/1kx5_CG.pdb'
trajectory = ''
simulation_name = os.path.basename(topology).split(".")[0]

print('Generating Universe...')
if trajectory != "":
    u = mda.Universe(topology, trajectory)
else:
    u = mda.Universe(topology)
print('Universe Generated!')

print('Begining Angle and Dihedral Calculation...')
u.atoms.guess_bonds({'D': 1, 'RB': 1})

print('Begining Measurements!')

def gen_name(measurement):
    mes_name = ''
    for atom in list(measurement.atoms):
        mes_name += atom.name + '-'
    return(mes_name[:-1])

def measure(measurement):
    if measurement.btype == 'bond':  # bond
        return(measurement.length())
    elif measurement.btype == 'angle':  # angle
        return(math.radians(measurement.angle()))
    elif measurement.btype == 'dihedral':  # dihedral
        return(math.radians(measurement.value()))

measurements = list(u.bonds) + list(u.angles) + list(u.dihedrals)
measurement_dict = {}
measurement_names = []

for measurement in measurements:
    name = gen_name(measurement)
    measurement_dict[name] = []
    measurement_names.append(name)

if trajectory != "":
    for t in u.trajectory[0:max_frame:stride]:
        for i, measurement in enumerate(measurements):
            measurement_dict[measurement_names[i]].append(measure(measurement))
else:
    for i, measurement in enumerate(measurements):
        measurement_dict[measurement_names[i]].append(measure(measurement))


print('\nExporting {} measurement datasets to file...'.format(
    len(measurement_names)))
for i, name in enumerate(measurement_names):  # loops thru each measurement
    filename = f'outputs/measurement_data/{simulation_name}/{name}.dat'
    with open(filename, 'w+') as instance_output:  # writes measurment list data to file
        # writes integer denoting mes_type to file
        unit = "angstrom" if i < len(u.bonds) else "radians"
        instance_output.write(f'mes_type: {"a"} unit: {unit}\n')

        value_str = ""
        for value in measurement_dict[name]:
            value_str += str(value) + '\n'
        instance_output.write(value_str)
# print('Done!')
