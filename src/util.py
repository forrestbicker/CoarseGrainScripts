# ================ Dependencies ================ #
from enum import Enum
from MDAnalysis import Universe
from os import path
    
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

    print('Universe Generated!')

    return u

def get_file_name(file):
    return path.basename(file).split(".")[0]