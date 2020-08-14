#! /Library/Frameworks/Python.framework/Versions/3.7/bin/python3
# ============================================== #
# ============== Boandi Analysis =============== #
# ============================================== #
# Written by Forrest Bicker
# August 2019
#



# ================ Requiremets ================  #
import numpy as np
from numpy import log
from numpy import sin



# ================= Histogram =================  #
class Histogram:  # a Histogram containing many bin objects
    def __init__(self, name, mes_type, values, bin_width):
        self.name = name
        self.mes_type = mes_type
        self.bins = []
        self.count = 0
        self.counter = -1

        self.bin_width = bin_width
        self.bin_min = min(values) - self.bin_width
        self.bin_max = max(values) + self.bin_width

        # generates bin edges and creates respective empty bin objects
        self.edges = [edge for edge in np.arange(
            self.bin_min, self.bin_max, self.bin_width)]
        for i, floor in enumerate(self.edges[:-1]):
            ceil = self.edges[i + 1]
            self.bins.append(Bin(self, floor, ceil))
            self.count += 1

        for value in values:
            self.add_instance(value)


    # creates a BoAnDi instance and places it in its appropriate bin
    def add_instance(self, value):
        displacement = value - self.bin_min
        ix = int(displacement / self.bin_width)
        bin = self.bins[ix]
        bin.add_instance(value)


    # returns list of non-empty bin objects
    def clear_empty_bins(self):
        self.bins = [bin for bin in self.bins if bin.count != 0]
        self.count = len(self.bins)


    def get_biggest(self, num):
        return([bin for bin in sorted(self, key=lambda bin: bin.boltz(), reverse=False)][:num])


    def __iter__(self):
        self.counter = 0
        return(self)


    # iterates over all contained bin objects
    def __next__(self):
        self.counter += 1
        if self.counter >= self.count:
            raise StopIteration
        else:
            return(self.bins[self.counter])



# ==================== Bin ====================  #
class Bin:  # a bin containing many boandi objects
    def __init__(self, p_histogram, floor, ceil):
        self.floor = floor
        self.ceil = ceil
        self.contents = []
        self.count = 0
        self.p_histogram = p_histogram
        self.mes_type = p_histogram.mes_type


    def __int__(self):
        return(self.count)


    def add_instance(self, value):
        self.count += 1
        self.contents.append(BoAnDi(self, value))


    def boltz(self):
        kB = 1.38e-23 * .1 * 6.022e23 / 418.4
        T = 298
        x = abs(self.floor)
        Px = self.count / self.p_histogram.count
        if self.mes_type == '0':
            return(-kB * T * log(Px / (4 * math.pi * x ** 2)))
        else:
            return(-kB * T * log(Px / sin(x)))


# ================= BoAnDi ================= #
class BoAnDi:  # a single instance of an bond/angle/dihedral
    def __init__(self, p_bin, value):
        self.value = value
        self.p_bin = p_bin
        self.mes_type = p_bin.mes_type


    def __float__(self):
        return(self.value)


# ================= Container =================  #
class Container:  # A glorified list to be converted to a histogram
    def __init__(self, mes_name, mes_type):
        self.mes_name = mes_name    
        self.mes_type = mes_type
        self.values = []

    def add_values(self, value_list):
        for value in value_list:
            self.values.append(value)
