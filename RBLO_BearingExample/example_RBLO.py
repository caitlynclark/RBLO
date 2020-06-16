#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:09:54 2019

@author: cclark2
"""

#Internal Modules
import random
from itertools import combinations

#External Modules
import floris.tools.rblo as opt


'''Parameters to test
#wd = [0, 90, 135, 180]
#ws = [8, 12, 18]
#shear = [0.1, 0.2, 0.3]
#TI = [5, 10, 15, 20]
'''

seeds = 2                           # for multi-starts

# Define wind rose data
wd = [0]                            # wind direction 
wd_freq = [1]                       # wind direction frequency

ws = [8]                            # wind speed
ws_freq = [1]                       # wind speed frequency

shear = 0.1 
TI = 10 

# Set turbine and array parameters 
D = 126                             # diameter of rotor (m)
N = 5                               # number of turbines

min_dist = 5*D                      # minimum separation distance

x_min = 0
x_max = 1                           
y_min = 0
y_max = 1

# Set array boundaries
bndx_min = 0
bndx_max = 1890
bndy_min = 0
bndy_max = 1890

average = []

#for ii in range(2, 16):

#    N = ii
objective = []
cost = []
power = []
replacements = []

count = 0

for jj in range(seeds):

    while count <= seeds:        
    #    # Random layouts    
        layout_x = [random.uniform(0, x_max) for ii in range(N)] 
        layout_y = [random.uniform(0, y_max) for ii in range(N)] 
            
        # Set optimization options
        opt_options = {'maxiter': 600, 'disp': True, 'iprint': 2, 'ftol': 1e-8}
        
        # Instantiate the layout optimization object
        layout_opt = opt.RBLO(  N = N, 
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
                                wd = wd,
                                wd_freq = wd_freq,
                                ws = ws,
                                ws_freq = ws_freq,
                                TI = TI,
                                shear = shear,
                                min_dist = min_dist,
                                cons = None,
                                opt_options=opt_options,
        )
                
        
        # Run optimization
        opt_locs = layout_opt.optimize()
        locsx_D = [x/D for x in opt_locs[0]]
        locsy_D = [x/D for x in opt_locs[1]]
    
        
        loc_combos = list(combinations(zip(locsx_D, locsy_D), 2))
        
        distance = [layout_opt._calc_distance(ii[0][0], ii[0][1], ii[1][0], ii[1][1]) for ii in loc_combos]
        if all(ii >= 5 for ii in distance): 
                    
            # Plot the new layout vs the old layout
            objective_temp, cost_temp, power_temp, replacements_temp = layout_opt.plot_layout_opt_results(count)
            objective.append(objective_temp)
            cost.append(cost_temp)
            power.append(power_temp)
            replacements.append(replacements_temp)
            count += 1
    
        else:
            print('Constraint violated')
            
    
print('Objective: ', objective)
print('Cost: ', cost)
print('Power: ', power)
print('Replacements: ', replacements)
print('Wind Direction: ', wd, 'No. Turbines: ', N, 'TI: ', TI, 'Shear: ', shear, 'Velocity: ', ws)
