from A_CoarseGrainer import coarse_grain
from B_Generalized_Parametizer import parametize
from C_CurveFitter import spline_fit
from app import startup_manual_refining
from util import MesType, generate_universe, get_file_name


# ================= User Input ================= #
topology = 'inputs/pkA_run3.psf'
trajectory = 'inputs/pkA_run3.dcd'

residue_list = ['ALA', 'DGLU', 'DLYS']  # list of ammino acids to be CoarseGrained
# residue_list = ['ALA']
# residue_list = ['DA', 'DT', 'DG', 'DC', 'PHOSPHATE', 'RIBOSE']


# ================= Execution ================= #
u = generate_universe(topology, trajectory)
sim_name = get_file_name(topology)
u_cg = coarse_grain(u, residue_list, simulation_name=sim_name, export=False)
measurement_dict = parametize(u_cg, export=True)
startup_manual_refining(measurement_dict, u_cg)

# for measurement_blueprint_name in measurement_dict:
#     measurement_blueprint = measurement_dict[measurement_blueprint_name]
#     aggregate_values = []
#     for measurement_name, measurement_data in measurement_blueprint.items():
#         aggregate_values += measurement_data['values']
#     measurement_type = MesType(len(measurement_blueprint_name))
#     spline_fit(aggregate_values, measurement_type, 0.01)

# with open('outputs/{sim_name}_parametized.txt', 'w+') as outfile:
#     for key, value in measurement_dict.items():
#         atoms = key.count('-') + 1
#         measurement_type = MesType(atoms)

#         k, x0 = fit_curve(value, measurement_type, 0.01)
        
#         outfile.write(f"{measurement_type.name} {key.replace('-', ' ')} {k} {x0}")
    
