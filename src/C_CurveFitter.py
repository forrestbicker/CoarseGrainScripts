#! /Library/Frameworks/Python.framework/Versions/3.7/bin/python3
# ============================================== #
# =========== C: Coarse Grain Fitter =========== #
# ============================================== #
# Written by Forrest Bicker
# August 2019
#



# ================ Requiremets ================  #
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import os
# import math

import BINAnalysis as boandi
from util import func_to_xy


# ================= Input File ================= #
value_file = 'outputs/measurement_data/EB2_EB3.dat'


# ============= Specifications ============== #
view_range = 0  # negative integer for a specific number of standard deviations; 0 for full view; positve integer for a specific value
# step = 0.01  # step needs manual adjustment (See README section C2. for help)


# ============= Fitting Functions ============== #
def f(x, k, x0, c):
    return(k * (x - x0)**2 + c)




# ================= Execution ================== #
def spline_fit(values, mes_type, step, export=False, plot=True, view_range=0):
    # file_name = os.path.splitext(os.path.basename(value_file))[0]
    # with open(value_file, 'r') as file:
    #     dataset = file.read()
    #     dataset_list = dataset.split('\n')
    #     mes_type = dataset_list.pop(0)[-1]
    #     values = [float(value) for value in dataset_list]

    # === Calculating bin information === #
    if view_range <= -1:
        view_range = abs(view_range) * np.std(values)
    elif view_range == 0:
        view_range = len(values)


    # initiates histogram object
    histogram = boandi.Histogram(mes_type, values, step)
    # for value in values:
    #     histogram.add_instance(value)  # places each value into its appropriate bin
    histogram.clear_empty_bins()


    x = histogram.get_floors()
    y = histogram.get_boltzes()


    # === Fitting curve === #
    spl = UnivariateSpline(x, y)
    xs = np.linspace(min(x), max(x), 1000)
    # k, x0, c = curve_fit(f, x, y, maxfev=1000000, p0=[10, 5, 5])[0]


    if plot:
        lo_cut = histogram.get_biggest(1)[0].floor - (view_range / 2)
        up_cut = histogram.get_biggest(1)[0].floor + (view_range / 2)
        x = [bin.floor for bin in histogram if lo_cut <= bin.floor <= up_cut]
        y = [bin.boltz() for bin in histogram if lo_cut <= bin.floor <= up_cut]
        # === Plotting points === #
        x2 = [bin.floor for bin in histogram.get_biggest(1)]
        y2 = [bin.boltz() for bin in histogram.get_biggest(1)]
        plt.scatter(x2, y2, s=0.5, c='#eb5600')  # plots biggest bin in red
        plt.scatter(x, y, s=0.5)  # plots points within view_range of biggest bin
        plt.plot(xs, spl(xs), 'g', lw=3)


        # === Adusting display settings === #
        fmt = '{} {}\nBin: {:.2f} - {:.2f}\n{:.3f}(x-{:.3f})+{:.3f}'
        # plt.suptitle(fmt.format(histogram.mes_type, histogram.name,
                                # min_bin, max_bin, k, x0, c), fontname='Courier New')
        plt.ylabel('Boltzmann Inversion')
        plt.xlabel(histogram.mes_type.name.title() + ' Measurment')
        plt.savefig(f'outputs/fit/hist_{histogram.name}.png')
        plt.show()


    # === Histogram Data Writing === #
    # edge_flor = edges[biggest_bin_ix]
    # edge_ceil = edges[biggest_bin_ix+1]
    # percent_within = 100*(max(hist)/sum(hist))
    # hist_output.write('{:9} {:23} {:22.16f} {:22.16f}  {:.1f}%\n'.format(mes_type,histogram.name,edge_flor,edge_ceil,percent_within))
    #
    # # === Angle Data Writing === #
    # cd('angles')
    # angle_output = open('{}.dat'.format(histogram.name),'w')
    #
    # for bin in bins.values():
    #     for angle in sorted(bin.contents,key=float):
    #         output = '{} {}\n'.format(angle.val,angle.prob)
    #         angle_output.write(output)



    print('Outputs Written!\nTask Complete!')
    return spl
