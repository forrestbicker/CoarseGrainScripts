#! /Library/Frameworks/Python.framework/Versions/3.7/bin/python3
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

