#! /Library/Frameworks/Python.framework/Versions/3.7/bin/python3
# ============================================== #
# ============= A: Coarse Grainer =============  #
# ============================================== #
# Written by Forrest Bicker
# August 2019
#


# ================= Requiremet ================= #
from collections import defaultdict
import os

import MDAnalysis as mda
from commands import cd
from commands import progress


# ================ Input Files ================  #
topology = '/Users/forrestbicker/Desktop/alanin.pdb'
trajectory = '/Users/forrestbicker/Desktop/alanin.dcd'
simulation_name = 'alaninSIM'


# ================= User Input ================= #
residue_list = ['ALA']  # list of ammino acids to be CoarseGrained


# ============== Misc Initiation ==============  #
amino_acid_dict = {  # dictionary dictating how to split ammino acid segments
    'ALA': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2', 'CB', 'HB1', 'HB2', 'HB3'],
    },
    'ARG': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2'],
        '2': ['CD', 'HD1', 'HD2', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH22', 'HH21'],
    },
    'ASN': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'OD1', 'ND2', 'HD21', 'HD22'],
    },
    'ASP': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'OD1', 'OD2', 'HD22'],
    },
    'ASX': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'OD1', 'OD2', 'ND2', 'HD21', 'HD22'],
    },
    'CYS': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'SG', 'HG1'],
    },
    'GLU': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'OE2'],
    },
    'GLN': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'NE2', 'HE21', 'HE22'],
    },
    'GLX': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2', 'HA1', 'HA2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'OE2'],
    },
    'GLY': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2', 'H', 'HA1', 'HA2'],
    },
    'HIS': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'CD2', 'HD2', 'NE2', 'HE2']
    },
    'ILE': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB', 'CG1', 'HG11', 'HG12', 'CG2', 'HG21', 'HG22', 'HG23', 'CD', 'HD1', 'HD2', 'HD3'],
    },
    'LEU': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23'],
    },
    'LYS': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2'],
        '2': ['CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3'],
    },
    'MET': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'SD', 'CE', 'HE1', 'HE2', 'HE3'],
    },
    'PHE': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'CD1', 'CD2', 'HD1', 'HD2'],
        '2': ['CE1', 'HE1', 'CE2', 'HE2', 'CZ', 'HZ'],
    },
    'PRO': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2', 'CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2'],
    },
    'SER': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'OG', 'HG1'],
    },
    'THR': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB', 'OG1', 'HG1', 'CG2', 'HG21', 'HG22', 'HG23'],
    },
    'TRP': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'CD2'],
        '2': ['CE3', 'HE3', 'CZ3', 'HZ3', 'CH2', 'HH2'],
        '3': ['CZ2', 'HZ2', 'CE2', 'NE1', 'HE1', 'CD1', 'HD1'],
    },
    'TYR': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CD2', 'HD2'],
        '2': ['CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2', 'HE2'],
    },
    'VAL': {
        'BB': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
        '1': ['CB', 'HB', 'CG1', 'HG11', 'HG12', 'CG2', 'HG21', 'HG22', 'HG23'],
    },
}

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

dummies = []
residue_atoms = []
for resname in residue_list:  # loops tru each residue to be coarse grained
    # extracts the residue name in amino_acid_dict format
    resname_root = resname[-3:]
    resname_sel = f'resname {resname}'
    # selects all resname-specific atoms
    resname_atoms = u.atoms.select_atoms(resname_sel)
    # identifys all resname-specific residues
    resids = resname_atoms.residues.resids
    for resid in resids:  # loops thu each matching residue id
        try:
            segments = amino_acid_dict[resname_root].keys()
        except KeyError:
            print('{} was not found in amino_acid_dict. Please add its parameters to the dictionary. (See README section A3. for help)'.format(
                resname_root))
            raise
        for segment in segments:  # loops thru each segment of each residue
            name_params = ' '.join(amino_acid_dict[resname_root][segment])
            params = 'resname {} and resid {} and (name {})'.format(
                resname, resid, name_params)
            # selects all atoms in a given residue segment
            atms = u.atoms.select_atoms(params)
            dummy = atms[0]
            # positions a dummy atom at the center of mass
            dummy.position = atms.center_of_mass()
            # names dummy atom in propper format
            dummy.name = '{}{}{}'.format(
                abrev_dict[resname[-3:]], segment[0], resid)

            dummies.append(dummy)
            residue_atoms.append(atms)

progress(0)
for frame in u.trajectory:  # loops tru each frame
    f = frame.frame
    for i, dummy in enumerate(dummies):
        dummy.position = residue_atoms[i].center_of_mass()
    progress(f / number_of_frames)
progress(1)
print('\nGenerated All Coarse Grained Molecules!')


# =================== Output =================== #
print('Writing Output Files...')

cd('CoarseGrain')
fools = mda.AtomGroup(dummies)
fools.write('{}_CoarseGrain.pdb'.format(simulation_name))
print('Topology written!')
with mda.Writer('{}_CoarseGrain.dcd'.format(simulation_name), fools.n_atoms) as w:
    for frame in u.trajectory:
        w.write(fools)
print('Trajectory written!\nTask Complete')
