# Coarse Grain Scripts v0.9.0
Programmed by Forrest Bicker

## v0.9.0

Initial Stable Public Beta release

## v0.9.1

### A

+ Added dynamic simulation detection to remove an input parameter

### C

+ Prevented Armstrongs from mistakenly being converted to Radians
+ Added support for `view_range` to be set to any amount of standard deviations

### Boandi

+ Corrected Container().get_biggest() to accurately return the global minima, rather than the biggest bin

## v0.9.2

### A

+ Fixed crash that caused the script to crash if the input simulation `type` data did not match `name` data
+ Also significantly decreases runtime for script A

## v0.9.3

### B

+ Reduced script computations to run in ~1/3rd the time

### Boandi

+ Fixed crash that prevented the boltzmann inversion of bond measurements from being calculated correctly
