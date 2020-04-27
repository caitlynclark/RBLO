This directory contains the code to process FAST.Farm data generated over 36 inflow wind conditions to complete RBLO based on planet bearing L10 life. This directory contains three files:

1] make_surrogates.py contains the module that reads in the FAST.Farm data, converts the binary data files and concatenates the seed data into readable dataframes, calculates L10 based on that FAST.Farm data, and then builds and pickles an interpolation surface.

2] example_make_surrogates.py contains code to make surrogates for turbine 1 (upwind) and 2 (downwind) over 36 inflow conditions

3] make_surrogates.sh is the shell script to execute the surrogate-making operation on NREL's high-performance computer. 

4] /Surrogate_Models/ contains the resulting 2D interpolation functions for 36 inflow conditions using FAST.Farm data and their visualized surfaces
