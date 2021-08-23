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

max_frame = 10000
stride = 1
block_count = 3


def gen_name(measurement):
    names = [atom.name for atom in measurement.atoms]
    return '-'.join(names)

def measure(measurement):
    if measurement.btype == 'bond':  # bond
        return(measurement.length())
    elif measurement.btype == 'angle':  # angle
        return(math.radians(measurement.angle()))
    elif measurement.btype == 'dihedral':  # dihedral
        return(math.radians(measurement.value()))

def parametize(u, simulation_name='simulation_name', export=False):

    print(f'Found {len(u.bonds)} bonds, {len(u.angles)} angles, and {len(u.dihedrals)} dihedrals!')
    print(f'Begining Bond, Angle, and Dihedral measurements over {len(u.trajectory)} frames...')

    # u.atoms.guess_bonds(vdwradii=config.vdw_radi)

    # print('Begining Measurements!')

    measurements = list(u.bonds) + list(u.angles) + list(u.dihedrals)
    measurement_blueprint_dict = {}
    for measurement in measurements:
        atoms = measurement.atoms
        atom_names = tuple([atom.name for atom in atoms])
        name = gen_name(measurement)

        # get just the first 2 letters of the atom names, which encode the type of atom and the segment of the residue
        # the ordered combination of all the atom types define the type of measurement.
        # a measurement is of a certian type if the fist 2 letters of its constituent atoms in order match the type. it is okay if the order is reversed, but nor permuted.
        measurement_blueprint = tuple([name[:2] for name in atom_names])

        if tuple(reversed(measurement_blueprint)) in measurement_blueprint_dict:
            measurement_blueprint = tuple(reversed(measurement_blueprint))

        if measurement_blueprint not in measurement_blueprint_dict:
            measurement_blueprint_dict[measurement_blueprint] = {}

        # create a dictioary that contains dictonaries of all the measurements of that type, e.g.
        # {
        #     (AB, A1): { 
        #         'AB1-A12': AtomGroup
        #     }
        # }
        measurement_blueprint_dict[measurement_blueprint][name] = {
            'measurement': measurement,
            'values': [],
        }

    measurement_dict = {}
    measurement_names = []

    progress(0)
    for t in u.trajectory[0:max_frame:stride]:
        for measurement_blueprint in measurement_blueprint_dict:
            for measurement_name in measurement_blueprint_dict[measurement_blueprint]:
                measurement_data = measurement_blueprint_dict[measurement_blueprint][measurement_name]
                measurement = measurement_data['measurement']
                measurement_data['values'].append(measure(measurement))
        progress(t.frame / len(u.trajectory))
    progress(1)
    print('\n')

    if export:
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
    print('Finished measuring bonds angles and dihedreals!')

    return measurement_blueprint_dict
