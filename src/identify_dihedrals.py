import sys
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_bonds, guess_angles, guess_dihedrals

if len(sys.argv) != 2:
    print("Usage: python identify_dihedrals.py <topology>")
    sys.exit(1)
else:
    u = mda.Universe(sys.argv[1])

    n_bonds = len(u.bonds)
    print(f'{len(u.bonds)} bonds specified in original file')

    u.add_TopologyAttr('bonds', guess_bonds(u.atoms, u.atoms.positions))  # add additional bonds using same bond guessing algoritm as VMD
    print(f'{len(u.bonds) - n_bonds} additional bonds guessed')
    u.add_TopologyAttr('angles', guess_angles(u.bonds))  # add angles in order to be able to calculate dihedrals
    u.add_TopologyAttr('dihedrals', guess_dihedrals(u.angles))  # add dihedrals

    print(f'{len(u.dihedrals)} dihedrals found. Outputting atom IDs:')

    # print id of all atoms in every dihedral
    for i, dihedral in enumerate(u.dihedrals):
        ids = [str(a.id) for a in dihedral]
        print('  ' + ' '.join(ids))
