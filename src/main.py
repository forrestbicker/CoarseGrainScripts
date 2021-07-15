from A_CoarseGrainer import coarse_grain
from B_Generalized_Parametizer import parametize
from util import generate_universe, get_file_name


# ================= User Input ================= #
topology = 'inputs/trimmed_traj/pkA/pkA_run1.psf'
trajectory = 'inputs/trimmed_traj/pkA/pkA_run1.dcd'

# residue_list = ['GLY', 'DGLU', 'DLYS']  # list of ammino acids to be CoarseGrained
residue_list = ['ALA']
# residue_list = ['DA', 'DT', 'DG', 'DC', 'PHOSPHATE', 'RIBOSE']


# ================= Execution ================= #
u = generate_universe('inputs\\alanin.pdb')
sim_name = get_file_name(topology)
u_cg = coarse_grain(u, residue_list, simulation_name=sim_name)
parametize(u_cg)
