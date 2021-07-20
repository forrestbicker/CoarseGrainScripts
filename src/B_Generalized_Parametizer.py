#! /Library/Frameworks/Python.framework/Versions/3.7/bin/python3
# ============================================== #
# ====== B: Generalized Coarse Grain Parameterization ====== #
# ============================================== #
# Written by Forrest Bicker
# September 2020
#

# ================ Dependencies ================ #
# import numpy as np
# import multiprocessing
# import time
import os

import MDAnalysis as mda
# import BINAnalysis as boandi

# from util import colorify
# from json import load
# from util import MesType
import config
import math

from util import progress

max_frame = 1000
stride = 1
block_count = 3


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

def parametize(u, simulation_name='simulation_name'):


    print('Begining Angle and Dihedral Calculation...')
    # u.atoms.guess_bonds(vdwradii=config.vdw_radi)

    print('Begining Measurements!')

    measurements = list(u.bonds) + list(u.angles) + list(u.dihedrals)
    measurement_dict = {}
    measurement_names = []

    for measurement in measurements:
        name = gen_name(measurement)
        measurement_dict[name] = []
        measurement_names.append(name)

    number_of_frames = len(u.trajectory)
    if number_of_frames > 1:
        for t in u.trajectory[0:max_frame:stride]:
            for i, measurement in enumerate(measurements):
                measurement_dict[measurement_names[i]].append(measure(measurement))
    else:
        for i, measurement in enumerate(measurements):
            measurement_dict[measurement_names[i]].append(measure(measurement))


    print('\nExporting {} measurement datasets to file...'.format(
        len(measurement_names)))


    if not os.path.isdir(f'outputs/measurement_data/{simulation_name}'):
        os.mkdir(f'outputs/measurement_data/{simulation_name}')

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
    
