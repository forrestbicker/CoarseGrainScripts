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
from MDAnalysis.core.groups import AtomGroup

from MDAnalysis.core.universe import Universe
from MDAnalysis.topology.guessers import guess_angles, guess_dihedrals
import config
# import math

import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import Atomindices
from util import generate_universe, get_file_name, progress

from json import load
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction

def coarse_grain(universe, residue_list, simulation_name='simulation_name', export=False):
    # ============== Misc Initiation ==============  #

    with open('src/mapping_dict.json', "r") as f:
        mapping_dict = load(f)

    with open('src/abrev_dict.json', "r") as f:
        abrev_dict = load(f)

    u = universe


    # ================= Execution =================  #
    
    print('Calculating Bond connections...')
    resnames = ' '.join(residue_list)
    original_bond_count = len(u.bonds)
    u.select_atoms(f'resname {resnames}').guess_bonds(vdwradii=config.vdw_radi)

    print(f'Original file contained {original_bond_count} bonds. {len(u.bonds) - original_bond_count} additional bonds infered.')

    print(f'Begining Coarse-Graining process...')

    bead_data = []
    cg_beads = []
    dummy_parents = {}
    for resname in residue_list:  # loops tru each residue to be coarse grained
        if resname == "PHOSPHATE" or resname == "RIBOSE":
            resname_atoms = u.atoms.select_atoms('resname DA DT DG DC DU')
        else:
            resname_atoms = u.atoms.select_atoms(f'resname {resname}')  # selects all resname-specific atoms

        if len(resname) == 4 and resname[0] == 'D': # for D-varants
            resname_key = resname[1:]
        else:
            resname_key = resname

        residues = resname_atoms.residues  # identifys all resname-specific residues
        for residue in residues:  # loops thu each matching residue id
            resid = residue.resid  # store int id
            try:
                segments = mapping_dict[resname_key].keys()
                for segment in segments:  # loops thru each segment of each residue
                    params = 'name ' + ' '.join(mapping_dict[resname_key][segment])  # generates param
                    # selects all atoms in a given residue segment
                    atms = residue.atoms.select_atoms(params)
                    dummy = atms[0]
                    # names dummy atom in propper format
                    dummy.name = str(abrev_dict[resname_key]) + str(segment[0]) + str(resid)
                    dummy.type = str(segment[0])

                    bead_data.append((dummy, atms)) 
                    cg_beads.append(dummy)

                    for atm in atms:
                        dummy_parents[atm.ix] = dummy
                    
            except KeyError:
                print(f'{resname_key} was not found in mapping/abrev_dict, skipping coarse grain. Please add its parameters to the dictionary. (See README section A3. for help)')

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
        # connect all parents with connected children
        for atom in atms:
            for bond in atom.bonds:
                for bonded_atom in bond.atoms:
                    if bonded_atom not in atms: # make more efficent if atms were a set
                        # by the end of all these loops and ifs, every bonded_atom that gets to this point is an atom connected to the edge of the cluster of atoms assigned to the coarse grain dummy bead in question
                        try:
                            new_bonds.append([cg_beads.index(dummy), cg_beads.index(dummy_parents[bonded_atom.ix])]) # type is used to store the cluster dummy
                        except KeyError: # raises if atom does not belong to a coarse grain bead
                            pass
                        #     try:
                        #         new_bonds.append([cg_beads.index(dummy), cg_beads.index(bonded_atom)]) # adds the bond between the dummies
                        #     except ValueError:  # if the other atom is just an atom withouot a coarse grain bead parent, ignore it
                        #        pass


    cg_beads = mda.AtomGroup(cg_beads)

    # TODO: EXPORT NEW_U INSTEAD OF OLD U
    # TODO: EXPORT NEW_U TO HAVE APPROPRIATE FRAMES
    # TODO: SHIFT THE DEFINITION OF CENTERS IN THE UNIVERSE EVEN IF NOT EXPORTING
    # TODO: AUTOTUNE THE CURVE TO FIND THE RIGHT STEP

    # #     # purge existing reminant bonds
    # u.delete_bonds(u.bonds)
    # u.delete_angles(u.angles)
    # u.delete_dihedrals(u.dihedrals)

    progress(0)
    number_of_frames = len(u.trajectory)
    for frame in u.trajectory:  # loops tru each frame
        f = frame.frame

        # positions a dummy atoms at cluster center of mass
        for dummy, atms in bead_data:
            dummy.position = AtomGroup(atms).center_of_mass()
        progress(f / number_of_frames)
    progress(1)
    print()

    print(f'Building new coarse-grained universe...')
    coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(), cg_beads).run().results
    new_u = mda.Merge(cg_beads)
    new_u.load_new(coordinates, format=MemoryReader)
    new_u.add_TopologyAttr('bonds', new_bonds)
    new_u.add_TopologyAttr('angles', guess_angles(new_u.bonds))
    new_u.add_TopologyAttr('dihedrals', guess_dihedrals(new_u.angles))
    print(f'Built universe with {len(new_u.atoms)} coarse-grained beads, {len(new_u.bonds)} bonds, {len(new_u.angles)} angles, and {len(new_u.dihedrals)} dihedrals')


    if export:
        print('Writing Output Files...')

        is_multiframe = number_of_frames > 1
        with mda.Writer(f'outputs/CoarseGrain/{simulation_name}_CG.dcd', new_u.atoms.n_atoms, multiframe=is_multiframe, bonds='all') as w:
            for frame in new_u.trajectory[1:]:  # loops tru each frame
                w.write(new_u.atoms)


        print('Generated All Coarse Grained Molecules!')
        print(f'Trajectory written to {simulation_name}_CG.dcd!')


        # for dummy, atms in bead_data:
        #         dummy.type = ''

        out_file = f'outputs/CoarseGrain/{simulation_name}_CG.pdb'
        with open(out_file, 'w+') as _:
            new_u.atoms.write(out_file, bonds='all')
        print(f'Topology written to {simulation_name}_CG.pdb!')

    print(f'Reduced {len(u.atoms)} atoms to {len(new_u.atoms)} beads!')

    print('Coarse Graining Task complete!')

    return new_u
