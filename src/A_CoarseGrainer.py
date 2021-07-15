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
# from collections import defaultdict
import os
# import math

import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import Atomindices
from util import generate_universe, get_file_name, progress

from json import load

def coarse_grain(universe, residue_list, simulation_name='simulation_name'):
    # ============== Misc Initiation ==============  #

    with open('src/mapping_dict.json', "r") as f:
        mapping_dict = load(f)

    with open('src/abrev_dict.json', "r") as f:
        abrev_dict = load(f)

    u = universe


    # ================= Execution =================  #
    

    print('Genarating Coarse Gained Molecules...')

    print(u.bonds)
    print('Calculating Bond connections...')
    resnames = ' '.join(residue_list)
    u.select_atoms(f'resname {resnames}').guess_bonds(vdwradii={
        'MN': 1,
        'CT3': 1.7,
        'OC': 1.5,
        'CT1': 1.7,
        'HB1': 1,
        'HA3': 1,
        'CC': 1.7,
        'NH1': 1.625,
        })
    print(u.bonds)

    bead_data = []
    cg_beads = []
    dummy_parents = {}
    for resname in residue_list:  # loops tru each residue to be coarse grained
        if resname == "PHOSPHATE" or resname == "RIBOSE":
            resname_atoms = u.atoms.select_atoms('resname DA DT DG DC DU')
        else:
            resname_atoms = u.atoms.select_atoms(f'resname {resname}')  # selects all resname-specific atoms
        residues = resname_atoms.residues  # identifys all resname-specific residues
        for residue in residues:  # loops thu each matching residue id
            resid = residue.resid  # store int id
            try:
                segments = mapping_dict[resname].keys()
                for segment in segments:  # loops thru each segment of each residue
                    params = 'name ' + ' '.join(mapping_dict[resname][segment])  # generates param
                    # selects all atoms in a given residue segment
                    atms = residue.atoms.select_atoms(params)
                    dummy = atms[0]
                    # names dummy atom in propper format
                    dummy.name = str(abrev_dict[resname]) + str(segment[0]) + str(resid)
                    dummy.type = str(segment[0])

                    bead_data.append((dummy, atms))
                    cg_beads.append(dummy)

                    for atm in atms:
                        dummy_parents[atm.ix] = dummy
                    
            except KeyError:
                print(f'{resname} was not found in abrev_dict, skipping coarse grain. Please add its parameters to the dictionary. (See README section A3. for help)')

    cg_beads = mda.AtomGroup(cg_beads)

    new_bonds = []
    # for residue in residue_list:
    #     for mapping in mapping_dict[residue]["Bonds"]:
    #         first_code = mapping[0]  # segment, resid offset, resname
    #         seccond_code = mapping[1] if isinstance(mapping[1], list) else [mapping[1], 0, residue]

    #         type_params = list(mapping_dict[residue]["Mapping"].keys())[first_code]

    #         first_atoms = cg_beads.select_atoms(f'resname {residue} and type {type_params}')
    #         for first_atom in first_atoms:
    #             type_params = list(mapping_dict[residue]["Mapping"].keys())[seccond_code[0]]  # segment
    #             seccond_atom_resid = int(first_atom.resid) + int(seccond_code[1])
    #             try:
    #                 seccond_atom = cg_beads.atoms.select_atoms(f'resname {seccond_code[2]} and type {type_params} and resid {seccond_atom_resid}')
    #             except IndexError:
    #                 pass

    #             if isinstance(seccond_atom, mda.core.groups.AtomGroup):
    #                 closest = seccond_atom[0]
    #                 closest_dist = mda.AtomGroup([first_atom, seccond_atom[0]]).bond.length()
    #                 for atom in seccond_atom:
    #                     dist = mda.AtomGroup([first_atom, atom]).bond.length()
    #                     if dist < closest_dist:
    #                         closest = atom
    #                         closest_dist = dist
    #                 seccond_atom = closest

    #             new_bonds.append([first_atom.index, seccond_atom.index])

    # new_bonds = []
    for dummy, atms in bead_data:
        for bond in dummy.bonds:
            for atom in bond.atoms:
                if atom != dummy:
                    if atom not in atms:
                        try:
                            new_bonds.append([dummy.ix, dummy_parents[atom.ix].ix]) # type is used to store the cluster dummy
                        except KeyError: # raises if connected atom is annother dummy
                            new_bonds.append([dummy.ix, atom.ix])

    # #     # purge existing reminant bonds

    for bond in u.bonds:
        u.delete_bonds([bond])

    u.add_TopologyAttr('bonds', new_bonds)

    print('Writing Output Files...')

    number_of_frames = len(u.trajectory)
    if number_of_frames > 1:
        progress(0)
        with mda.Writer(f'outputs/CoarseGrain/{simulation_name}_CG.dcd', cg_beads.n_atoms, multiframe=True, bonds='all') as w:
            for frame in u.trajectory:  # loops tru each frame
                f = frame.frame

                # positions a dummy atoms at cluster center of mass
                for dummy, atms in bead_data:
                    dummy.position = atms.center_of_mass()

                w.write(cg_beads)
                progress(f / number_of_frames)
        progress(1)

        print('\nGenerated All Coarse Grained Molecules!')
        print(f'Trajectory written to {simulation_name}_CG.dcd!')
    else:
        for dummy, atms in bead_data:
            dummy.position = atms.center_of_mass()

    # for dummy, atms in bead_data:
    #         dummy.type = ''

    print(u.bonds)

    out_file = f'outputs/CoarseGrain/{simulation_name}_CG.pdb'
    with open(out_file, 'w+') as _:
        cg_beads.write(out_file, bonds='all')
    print(f'Topology written to {simulation_name}_CG.pdb!')
    print(f'Reduced {len(u.atoms)} atoms to {len(cg_beads)} beads!')

    print('Task complete!')

    return u
