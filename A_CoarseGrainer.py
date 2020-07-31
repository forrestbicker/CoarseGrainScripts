#! /Library/Frameworks/Python.framework/Versions/3.7/bin/python3
# ============================================== #
# ============= A: Coarse Grainer =============  #
# ============================================== #
# Written by Forrest Bicker
# August 2019
#

# O(fr)
# f = number of frames
# r = number of residue segments

# ================= Requiremet ================= #
from collections import defaultdict
import os

import MDAnalysis as mda
from util import progress

from json import load


# ================ Input Files ================  #
topology = 'inputs/alanin.pdb'
trajectory = 'inputs/alanin.dcd'
simulation_name = 'alanin'

# ================= User Input ================= #
residue_list = ['ALA']  # list of ammino acids to be CoarseGrained

# ============== Misc Initiation ==============  #
with open('mapping_dict.json', "r") as f:
    mapping_dict = load(f)

abrev_dict = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M',
    'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}


# ================= Execution =================  #
print('Generating Universe...')
u = mda.Universe(topology, trajectory)
print('Universe Generated!')

print('Genarating Coarse Gained Molecules...')

number_of_frames = len(u.trajectory)


bead_data = []
for resname in residue_list:  # loops tru each residue to be coarse grained
    # extracts the residue name in amino_acid_dict format
    resname_root = resname[-3:]
    resname_sel = f'resname {resname}'
    # selects all resname-specific atoms
    resname_atoms = u.atoms.select_atoms(resname_sel)
    # identifys all resname-specific residues
    residues = resname_atoms.residues
    for residue in residues:  # loops thu each matching residue id
        resid = residue.resid
        try:
            segments = mapping_dict[resname_root].keys()
        except KeyError:
            print('{resname_root} was not found in amino_acid_dict. Please add its parameters to the dictionary. (See README section A3. for help)')
            raise
        for segment in segments:  # loops thru each segment of each residue
            name_params = ' '.join(mapping_dict[resname_root][segment])
            params = f'name {name_params}'
            # selects all atoms in a given residue segment
            atms = residue.atoms.select_atoms(params)
            dummy = atms[0]
            # positions a dummy atom at the center of mass
            dummy.position = atms.center_of_mass()
            # names dummy atom in propper format
            dummy.name = '{}{}{}'.format(
                abrev_dict[resname[-3:]], segment[0], resid)

            bead_data.append((dummy, atms))


fools = mda.AtomGroup([pair[0] for pair in bead_data])

print('Writing Output Files...')

progress(0)
with mda.Writer(f'outputs/CoarseGrain/{simulation_name}_CG.dcd', fools.n_atoms) as w:
    for frame in u.trajectory:  # loops tru each frame
        f = frame.frame

        # positions a dummy atoms at cluster center of mass
        for dummy, atms in bead_data:
            dummy.position = atms.center_of_mass()

        w.write(fools)
        progress(f / number_of_frames)
progress(1)
print('\nGenerated All Coarse Grained Molecules!')
print('Trajectory written!')

fools.write(f'outputs/CoarseGrain/{simulation_name}_CG.pdb')
print('Topology written!')

print('Task complete!')
