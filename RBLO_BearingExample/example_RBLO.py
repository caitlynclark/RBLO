#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:09:54 2019

@author: cclark2
"""

#Internal Modules
import random

#External Modules
import floris.tools.rblo as opt


seeds = 1                           # for multi-starts

# Define wind rose data
wd = [270, 315, 0]                  # wind direction vector
ws = [8]                            # wind speed (make sure it aligns with interpolant surfaces)
freq = [(1/3), (1/3), (1/3)]        # frequency of wind direction

# Set turbine and array parameters 
D = 126                             # diameter of rotor (m)
N = 10                              # number of turbines

min_dist = 5*D                      # minimum separation distance

x_min = 0
x_max = 1 #design variables are set between 0 and 2*D*N, farm boundaries
y_min = 0
y_max = 1

# Set array boundaries
bndx_min = 0
bndx_max = 1890
bndy_min = 0
bndy_max = 1890

for jj in range(seeds):
        
    # Random layouts    
    layout_x = [random.uniform(0, x_max) for ii in range(N)] 
    layout_y = [random.uniform(0, y_max) for ii in range(N)] 
    
    # Set optimization options
    opt_options = {'maxiter': 400, 'disp': True, 'iprint': 2, 'ftol': 1e-8}
    
    # Instantiate the layout optimization object
    layout_opt = opt.RBLO(
                            N = N, 
                            D = D, 
                            layout_x = layout_x, 
                            layout_y = layout_y, 
                            x_min = x_min, 
                            x_max = x_max, 
                            y_min = y_min, 
                            y_max = y_max,
                            bndx_min = bndx_min,
                            bndy_min = bndy_min,
                            bndx_max = bndx_max,                        
                            bndy_max = bndy_max,
                            wd=wd,
                            ws=ws,
                            freq=freq,
                            min_dist = min_dist,
                            cons = None,
                            opt_options=opt_options,
    )
        
    # Run optimization
    layout_results = layout_opt.optimize()
    print('Layout Results: ', layout_results)

    # Plot the new layout vs the old layout
    layout_opt.plot_layout_opt_results(jj)
