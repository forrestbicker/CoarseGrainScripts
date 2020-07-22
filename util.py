# ================ Dependencies ================ #
from enum import Enum


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

