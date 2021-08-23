# ================ Dependencies ================ #
from enum import Enum
from MDAnalysis import Universe
from os import path
import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go

# ===================== cd ===================== #
# def cd(dir):  # cleanly changed cwd
#     native = os.path.dirname(os.path.abspath(__file__))
#     path = os.path.join(native,'outputs',dir)
#     try:  # makes output folder
#         os.mkdir(path)
#     except FileExistsError:
#         pass
#     finally:  # sets cwd to output folder
#         os.chdir(path)


# ================== colorify ================== #
def colorify(code, str):
    return(f'\033[{code}m{str}\033[0m')


# ================== progress ================== #
def progress(float, width=25):  # renders ACSII progress bar from precent completeion
    percent = int(float * 100)
    print('\r[{:{}}] {}%'.format('#' * (width * percent // 100), width, percent), end='', flush=True)

class MesType(Enum):
    NULL = 0
    BOND = 2
    ANGLE = 3
    DIHEDRAL = 4


def generate_universe(topology, trajectory=None):
    print('Generating Universe...')
    if trajectory is None or trajectory == '':
        u = Universe(topology)
    else:
        u = Universe(topology, trajectory)

    x, y, z = u.dimensions[:3]
    print(f'Universe with dimensions x: {x}, y: {y}, z: {z} loaded!')
    n_waters = u.select_atoms('resname WAT').n_residues
    print(f'{n_waters} water molecules detected!')

    return u

def get_file_name(file):
    return path.basename(file).split(".")[0]

def dfs(atom, visited):
    if not visited[atom.ix]:
        dfs(atom, visited)
        visited[atom.ix] = True


def func_to_xy(x, y, func, *argv):
    x_points = np.linspace(min(x), max(x), 2048)
    y_points = [None] * 2048
    for i in range(len(x_points)):
        y_points[i] = func(x_points[i], *argv)

    return([x_points, y_points])


# construct plotly scatter plot from histogram
def generate_figure(x_data, y_data, name, vertex=(1, 1)):
    def f(x, k, x0, c):
        return k * (x - x0)**2 + c

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x_data,
            y=y_data,
            name=name,
            mode='markers'
        )
    )

    k, x0, c = curve_fit(f, x_data, y_data, maxfev=100000, p0=[1, vertex[0], vertex[1]])[0]
    x_curve, y_curve = func_to_xy(x_data, y_data, f, k, x0, c)
    fig.add_trace(go.Scatter(x=x_curve, y=y_curve, mode='lines'))

    return fig, k, x0, c

def generate_figure_dihedral(x_data, y_data, name, vertex=(1, 1)):
    def f(x, k, n, d):
        x = np.rad2deg(x)
        return k * (1 + np.cos(n * x - d))

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x_data,
            y=y_data,
            name=name,
            mode='markers'
        )
    )

    k, n, d = curve_fit(f, x_data, y_data, maxfev=100000)[0]
    x_curve, y_curve = func_to_xy(x_data, y_data, f, k, n, d)
    fig.add_trace(go.Scatter(x=x_curve, y=y_curve, mode='lines', line_shape='spline'))

    return fig, k, n, d, 0
