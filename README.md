# ScarpOffsetAnalysis
Matlab code to calculate the vertical separation, dip slip, heave, and throw from a topographic profile across a fault scarp

Written by Christine Regalla 
last modified  3/23/16 and 12/30/24

This script will calculate vertical separation, heave, throw, and fault slip from a topographic profile across a fault scarp using a Monte Carlo Simulation. This code reads a tab delimited text file of topographic profile data, which must be formatted in two columns of x (distance) and z (elevation) data points. See example text file included (test_profile.txt). Once text profile is entered, the profile can be saved as a .mat file and reloaded in subsequent runs. Midpoint and regressions through lower and upper surfaces can be chosen graphically or entered numerically. Code accounts for uncertainty in the geometry of linear regressions through the upper and lower offset surfaces, the dip of the fault, and the intersection point of the fault with the ground surface. Requires Matlab Statistics Toolbox. 
