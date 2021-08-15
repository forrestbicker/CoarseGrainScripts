# Coarse Grain Scripts
A Generalized Algorithm for Coarse-Graining Molecular Dynamics Simulations with Boltzmann Parametrization

[![Latest Version](https://img.shields.io/github/v/release/forrestbicker/CoarseGrainScripts.svg)](https://github.com/forrestbicker/CoarseGrainScripts/releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)


###### Forrest Bicker
###### College Of Staten Island - Loverde Labratory
<p align="center">
  <img src="https://github.com/forrestbicker/CoarseGrainScripts/blob/d30803f818a1d56d78ea24317423301c154cce7d/assets/beforeAndAfter.png" alt="CoarseGrainedComparison">
</p>
<h6 align="center">Coiled DNA before and after Coarse-Graining</h6>

## Table of Contents  
* [Installation](#installation)
  * [Dependencies](#dependencies)  
  * [Setup](#setup)  
* [Usage](#usage)
  * [Coarse-Grainer](#coarse-grainer)
  * [Measurer](#measurer)
  * [Curve-Fitter](#curve-fitter)
* [Acknowledgements](#acknowledgements)  

## Installation
#### Dependencies
+ [`Python 3.8.10`](https://www.python.org/downloads/)
+ [`MDAnalysis 1.1.1`](https://www.mdanalysis.org/pages/installation_quick_start/)
+ [`SciPy 1.6.2`](https://scipy.org/install.html)
+ [`Dash 1.19.0`](https://dash.plotly.com/installation)
+ [`Plotly 4.14.3`](https://plotly.com/python/getting-started/)

Although other versions may work, the code was designed for and tested with these ones so it's reccomended you use them.

#### Setup
To use this program, you must first install the dependencies listed above to your local environment. It's reccomended to use a virtual environment such as conda to install and manage the versioning of these packages. Installation guides for each package can be found in the above section by clicking the hyperlink on the package name. Then, once you clone this repository to your local machine you should be able to execute the code found in `src/main.py`

## Usage
This program preforms three main tasks: converting an atomistic universe to a coarse-grained universe, measuring the bonds angles and dihedrals in a coarse-grained universe, and fitting curves to said measurements to help parameterize the coarse-grained simulation. Each of these tasks is abstracted into a single, modular function call that can be called from the default `main.py` file or imported into other projects if you want to apply pre- or post-processing to the inputs and outputs.

### Coarse-Grainer

###### **A1. Description**
Converts an atomistic molecular dynamics simulation consisting of a pre-calculated topology and a trajectory to a coarse grained universe. Pre-existing SDK coarse-grain mappings are provided for the fundamental amino acids and DNA nucleic acids + backbone. However, the mapping blueprint is designed to be highly customizable for indivudal expansion meaning any simulation can be corase grained in accordance to user-specified input.

###### **A2. Input Files**
+ `topology`: Atomistic topology file
+ `trajectory`: Atomistic trajectory file (if you wish to coarse-grain a static frame, use an empty string for trajectory)
+ `mapping_dict`: Dictionary containing Coarse Grain mappings which dictate the atoms that each coarse grained bead contains. This version includes mappings for all 20 natural amino acids, DNA nucleic acids, and DNA backbone components, however if you wish to coarse grain a molecule whose mapping is not contained in `mapping_dict`, you will have to add its mapping to the dictionary. If you so desire to use a mapping not included in the `mapping_dict`, you will need to add it yourself by following the format bellow. This can be done by appending a new dictionary to `mapping_dict.json` in the following format:

###### **A3. Input Parameters**
+ `residue_list`: List of three-letter amino acid abbreviations. The coarse grained file will only contain beads from amino acids included in this list.

```
'NAME': {
     'segment ID': [
         'type': type,
         'charge': charge
         'atoms': component_atoms,
     ],
}
```

   + _NAME_: The three-letter IUPAC abbreviation for the amino acid; DA, DT, DG, DC, PHOSPHATE, RIBOSE for DNA components
       + **WARNING**: The input topology and trajectory _must_ specify the NAME of each atom using the three-letter IUPAC NAME of its containing amino acid which matches how amino acids are named in `mapping_dict`. This is especially important to take note of when adding a mapping to the `mapping_dict`, as the three-letter IUPAC NAME you input must correspond to atom names in the topology and trajectory.
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
     'B': [
         'type': 'GBB',
         'charge': 0,
         'atoms': ['C', 'CA', 'O', 'N', 'HN', 'HA', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'OT1', 'OT2']
     ],
     '1': [
         'type': 'LY1',
         'charge': 0,
         'atoms': ['CB', 'HB1', 'HB2', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2']
     ],
     '2': [
         'type': 'LY2',
         'charge': 0.118,
         'atoms': ['CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3']
     ],
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


### Measurer

###### **B1. Description**
Measures all bonds, angles, and dihedrals between coarse grain beads. Can be configured to run computations in parallel on multiple CPUs.

###### **B2. Input Files**
+ `topology`: Coarse grained topology file (MUST specify bond connections between CG beads, this is automatically garunteed when directly using output from script A)
+ `trajectory`: Coarse grained trajectory file

###### **B3. Input Parameters**

###### **B4. Output**
+ **measurement data** files: Outputs a dat file containing the measured length/angle for all bonds/angles/dihedrals in `amino_acid_molds` across every observed frame. Each file is named by joining the names of its component beads.
    + **NOTE**: Units for length are Armstrongs; units for angles are degrees


### Curve-Fitter

###### **C1. Description**
Plots a series of measurement values in relation to their Boltzmann inversion on an xy-plane, and fits a curve to said points.

###### **C2. Input Files**
+ `value_file`: A dat file containing any number of newline-delimited values corresponding to the measurements of  (generated from script B)

###### **C3. Input Parameters**
+ `view_range`: Defines the range of measurement values surrounding the global minima to display and fit a curve to. The value of this number in the units of the measurement in question (Armstrongs/degrees) will determine range of the displayed data  Additionally, a negative number will evaluate to the corresponding number standard deviations of the dataset (e.g. -2 equals two standard deviations). To view _all_ data points set `view_range` to 0.
+ `step`: Defines the size of each bin of measurement values to be plotted. Should generally be somewhere between 0.001 and 0.1, although I have unfortunately been unable to find an effective way of programming this, so the value will need to be manually optimized depending on the `value_file`. To do this, run the program with a random value and repeatedly adjust it until a distinct curve starts to appear. Be careful not to make the value too small, or the graph will take a long time to render and look very strange when it does.

###### **C4. Output**
+ **Scatterplot**: A scatterplot containing all the bins within the `view_range` of the global minima, where the x-axis is the average measurement value of the contents of a bin, and the y-axis is its Boltzmann inversion.
+ `view_range`: Defines the range of measurement values surrounding the global minima to display and fit a curve to. The value of this number in the units of the measurement in question (Armstrongs/degrees) will determine range of the displayed data  Additionally, there are two special values which `view_range` can be set to for useful effects: -1 to view _all_ data points, and 0 to display one x-axis standard deviation of data points.
+ `step`: Defines the size of each bin of measurement values to be plotted. Should generally be somewhere between 0.001 and 0.1, although I have unfortunatley been unable to find an efficent way of prograing this, so the value will need to be manually optimized depending on the `value_file`. To do this, run the program with a random value and repeatedly adjust it until a distinct curve starts to appear. Be careful not to make the value too small, or the graph will take a long time to render and look very strange when it does.

###### **C4. Output**
+ **Scatterplot**: A scatterplot containg all the bins within the `view_range` of the global minima, where the x-axis is the average measurement value of the contents of a bin, and the y-axis is its Boltzmann inversion.

## Acknowledgements
+ Dr. Sharon Loverde, Principal Investigator
+ Phu Tang, project mentor
+ Tania Rajpersaud, provided debugging help and testing data
+ William Hu, project tester and contributor
