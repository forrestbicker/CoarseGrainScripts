# Coarse Grain Scripts v1.0.0

## Script A - Coarse Grainer

###### **A1. Description**
Converts an atomistic simulation to a coarse grained one using SDK coarse graining mapping

###### **A2. Input Files**
+ `topology`: Atomistic topology file
+ `trajectory`: Atomistic trajectory file

###### **A3. Input Parameters**
+ `residue_list`: List of three-letter amino acid abbreviations. The coarse grained file will only contain beads from amino acids included in this list.
+ `mapping_dict`: Dictionary containing Coarse Grain mappings which dictate the atoms that each coarse grained bead encapsulates. This version includes mappings for all 20 natural amino acids, however if you wish to coarse grain a molecule whose mapping is not contained in `mapping_dict`, you will have to add its mapping to the dictionary. If you so desire to use a mapping not included in the `mapping_dict`, you will need to add it yourself. This can be done by appending a new dictionary to `mapping_dict` in the following format:

```
'NAME': {
     'segment ID': [component_atoms],
     'segment ID': [component_atoms],
}
```

    + _NAME_: The three-letter IUPAC abbreviation for the amino acid
        + **WARNING**: The input topology and trajectory _must_ specify the NAME of each atom using the three-letter IUPAC NAME of its containing amino acid which matches how amino acids are named in `mapping_dict`. This is especially important to take note of when adding a mapping to the `mapping_dict`, as the three-letter IUPAC NAME you input must correspond to atom names in the topology and trajectory. Additionally, all NAMEs must be exactly three letters long or errors may occur.
    + _component_atoms_: A list of all atoms which compose a given bead.
        + **WARNING**: List contents will be directly passed to MDAnalysis to select the individual atoms, and as such must be in a format MDAnalysis can understand which matches the topology and trajectory contents. To ensure proper functionality, check how atoms are named in the source topology and trajectory if unsure.
        + **WARNING**: Be sure to pay attention of changes in atomic structure at the amino acid termini, as neglecting to do so will lead to incorrect bead placement.
        + **NOTE**: All entries in the  _component_atoms_ are selected with OR logic, meaning every atom list in _component_atoms_ does not necessarily have to appear in every residue. This means it is theoretically possible to accommodate for ambiguous amino acids such as GLX under this framework.
        + **NOTE**: The _NAME_ in `mapping_dict` represents both L-Chiral and D-Chiral amino acids.
    + _segment ID_: A one-character identification for each segment
        + **NOTE**: Will be used to name atoms in the outgoing topology and trajectory.

e.g.

```
'LYS': {
     'B': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2'],
     '1': ['CB', 'HB1', 'HB2', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2'],
     '2': ['CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3'],
}
```

###### **A4. Output**
+ **Topology**: Coarse Grained topology file
+ **Trajectory**: Coarse Grained trajectory file
+ **NOTE**: Coarse grained beads output in the coarse grained topology and trajectory are named in the format **amino acid** + **segment ID** + **residue ID**

     - **amino acid**: One-letter IUPAC code corresponding with the identity of the amino acid the bead is a part of. If the amino acid is D-chiral, the one-character IUPAC code is prefixed with the letter `D`, e.g. DG represents a D-chiral Glutamic Acid bead. (L-chirality is implicit in a pure one-letter code.)

     - **segment ID**: One-character code indicating which part of the amino acid mapping the bead corresponds to. The code is derived directly from the `mapping_dict`. Possible segment IDs are `B`, `1`, and `2`—although modification of the `mapping_dict` can yield new codes.

     - **residue ID**: Integer corresponding with the residue ID of the amino acid.

     - e.g. KB4 denotes the backbone of the fourth L-chiral Lysine residue


## Script B - Parameterizer

###### **B1. Description**
Measures all bond lengths, angle measures, and dihedral angles between coarse grain beads. Can be configured to run computations in parallel on multiple CPUs.

###### **B2. Input Files**
+ `topology`: Coarse grained topology file (generated from script A)
+ `trajectory`: Coarse grained trajectory file (generated from script A)

###### **B3. Input Parameters**
+ `block_count`: Determines the number of blocks which the computation will be broken into. A value of `1` will make the computation run on a single CPU. Any value greater than 1 will cause the computation to be spread across all available CPUs on the device, with each CPU computing one block at a time, and moving onto a new block if once it has completed the prior computation. Setting this value to the exact number of available CPUs will yield the best performance.

+ `max_frame`: Determines the final frame of `trajectory` to be analyzed. A value of `-1` will analyze all frames.
+ `stride`: Determines the stride to be used to analyze `trajectory`. A value of `1` will analyze all frames.

+ `amino_acid_molds`: This is a dictionary of nested dictionaries containing specification of the structure of all desired measurements. A new mold will likely have to be added for each coarse grained simulation you wish to analyze, however, since the format acts as a mold rather an explicit mapping for each measurement, molds are highly reusable with little to no modification. Each new mold entry should be in the following format:

```
'NAME': {
    'Bond': [bond_list],
    'Angle': [angle_list],
    'Dihedral': [dihedral_list]
},
```

    + _NAME_: Three-letter IUPAC abbreviation identifying the amino acid.
        + **NOTE**: Unlike in script A's `mapping_dict`, the _NAME_ in `amino_acid_molds` requires chirality specification: If the amino acid is D-chiral, the three-character IUPAC code must be prefixed with the letter `D`. (L-chirality is implicit in a pure three-character IUPAC code)
        + **NOTE**: The _NAME_ in `amino_acid_molds` actually supports multiple, space-delimited amino acid NAMEs, allowing you to notate and measure connections between different kinds of amino acids.

    + _bond_list_/_angle_list_/_dihedral_list_: A list of lists (either pairs, triplets, or quadruplets) of beads which define a bond/angle/dihedral. Beads are referenced in accordance to how their NAME was defined in script A: following the format **amino acid** + **segment ID** + **residue ID** (See Section A5.1.).
        + As an example, to measure the bond between the backbone bead, backbone bead and 1st segment bead of Lysine residue number 1, the bond list would be `['KB1', 'K11']`, where each item in the list denotes a simulated bead.
        + For convenience, rather than _only_ measuring the Lysine backbone-segment1 bond in residue 1, `['KB1', 'K11']` will measure the backbone-segment1 bond across _all_ Lysine residues. It does this using the roots of each bead NAME pattern (in this case, `'KB', 'K1'`) as a mold, to which incrementing resid values are successively appended. This means you will only have to specify a handful of _bond_list_/_angle_list_/_dihedral_list_ indicating general **patterns** which will then be applied across the entire simulation, rather than explicitly specifying every bond, angle, and dihedral connection.
        + **NOTE**: You are safely able to provide patterns that span across multiple residues, for example, there is full support for measuring the `['EB1', 'EB2', 'E12']` angle in a poly-E chain.
            + **WARNING**: If you wish to measure inter-residue bonds/angles/dihedrals, the RESIDs of the topology and trajectory must be ordered in a consistent fashion. If not, the script may incorrectly map the
        + **WARNING**: Cyclical chains of amino acids support has failed. Support will be re-implemented sometime in the future.

e.g.

```
'GLU DGLU': {
    'Bond': [['E11', 'EB1'], ['EB1', 'EB2']],
    'Angle': [['E11', 'EB1', 'EB2'], ['EB1', 'EB2', 'E12']],
    'Dihedral': [['E11', 'EB1', 'EB2', 'E12']]
}
```

###### **B4. Output**
+ **measurement data** files: Outputs a dat file containing the measured length/angle for all bonds/angles/dihedrals in `amino_acid_molds` across every observed frame. Each file is named by joining the names of its component beads.
    + **NOTE**: Units for length are Armstrongs; units for angles are degrees

## Script C - Curve Fitter
###### **C1. Description**
Plots a series of measurement values in relation to their Boltzmann inversion on an xy-plane, and fits a curve to said points.

###### **C2. Input Files**
+ `value_file`: A dat file containing any number of newline-delimited values corresponding to the measurements of  (generated from script B)

###### **C3. Input Parameters**
+ `view_range`: Defines the range of measurement values surrounding the global minima to display and fit a curve to. The value of this number in the units of the measurement in question (Armstrongs/degrees) will determine range of the displayed data  Additionally, there are two special values which `view_range` can be set to for useful effects: -1 to view _all_ data points, and 0 to display one x-axis standard deviation of data points.

###### **C4. Output**

## Dependencies
+ `MatPlotLib`
+ `SciPy`
+ `MDAnalysis`
+ `Boandi`*
    + Histogram, bin, and measurement assistance package
+ `commands`*
    + mini-library of three general-use commands implemented in these scripts

\* Indicates a custom-made pseudo-package which is

###### **C4. Output**
