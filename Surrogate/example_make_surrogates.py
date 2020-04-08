#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 09:55:09 2020

@author: cclark2
"""

#Internal Modules
import numpy as np
import math
from datetime import datetime
import multiprocessing 
from itertools import repeat
import dill
import mpi4py.MPI as MPI
from mpi4py.futures import MPIPoolExecutor 

MPI.pickle.__init__(dill.dumps, dill.loads)
num_cores = multiprocessing.cpu_count()


#External Modules
import make_surrogates

if __name__ == '__main__':
    inflows = ['Inflow' + '%02d' % (i) for i in range(36)]
    cases = ['Case' + '%03d' % (i) for i in range(625)]
    turbines = ['T1', 'T2']

    seeds = ['Seed_0', 'Seed_1', 'Seed_2', 'Seed_3', 'Seed_4', 'Seed_5']
    path_to_surrogate = '/scratch/cclark2/Surrogate_Models/'    
    
    D_rotor = 126.0 #rotor diameter of 5MW Ref turbine
    x_cross = np.arange(-3, 3.01, 0.25) * D_rotor #offset distances from FastFarm to train interpolant
    y_down  = np.arange(3, 15.1, 0.5) * D_rotor #distances downstream from FastFarm to train interpolant
    
    # Define turbine and drivetrain characteristics
    FF_timestep = 0.025                         # FAST.Farm timestep
    m_c = 5800                                  # mass of carrier, kg cc
    d_c = 863E-3                                # center distance, m
    #pa = 20/180 * math.pi                      # pressure angle, radians
    m_s = 18E3                                  # mass of shaft, kg
    m_p = 1500                                  # mass of planet, kg
    N = 3                                       # number of planets
    g = 9.81 * math.cos(5/180*math.pi)          # 5° is the bedplate tilting angle rho cc
    beta = [0.0,0.0,0.0]                        # mounting angle variation from standard division
    L_c = 4.25                                  # distance from hub to carrier center, m
    L_s = 3.22                                  # length of shaft, m
    L_p = 5                                     # distance between main shaft and pl bearing system center, m
    N_r = 56                                    # number of teeth in ring gear
    N_p = 17                                    # number of teeth in planet gear
    rho = 5*(math.pi/180)                       # bedplate tilting angle, radians
    C = 4730000                                 # capacity from bearing, N
    e = 10/3                                    # constant for roller bearings
    N_r = 56                                    # number of teeth in ring gear
    N_p = 17                                    # number of teeth in planet gear
    
    
    # Instantiate the drivetrain object
    surrogate = make_surrogates.Surrogate(
                            FF_timestep = FF_timestep,
                            m_c = m_c,
                            d_c = d_c,
                            m_s = m_s,
                            m_p = m_p,
                            N = N,
                            g = g,
                            beta = beta,
                            L_c = L_c,
                            L_s = L_s,
                            L_p = L_p,
                            rho = rho,
                            C = C,
                            e = e,
                            N_r = N_r,
                            N_p = N_p
    )
    
    startTime = datetime.now()
    
    for turbine in turbines:
        print("Starting " + str(turbine) + ": ", datetime.now() - startTime)
    
        if turbine == 'T1':
                case = cases[0]
                
                L10_list = []
                power_list = []
                
                outfile = "BINARY"
                with MPIPoolExecutor(max_workers=num_cores) as executor: 
                    for L10 in executor.map(surrogate.calc_L10, repeat(case), inflows, repeat(seeds), repeat(turbine), repeat(outfile)):
                        L10_list.append(L10)
                        
                    for power in executor.map(surrogate.calc_power, repeat(case), inflows, repeat(seeds), repeat(turbine), repeat(outfile)):
                        power_list.append(power)
                
                with open(str(path_to_surrogate) + str(turbine) + '_L10_List', "wb") as dill_file:
                    dill.dump(L10_list, dill_file)    
                print('Made ' + str(turbine) + '_L10 List (for all inflows) pickle', datetime.now() - startTime)   
                
                with open(str(path_to_surrogate) + str(turbine) + '_Power_List', "wb") as dill_file:
                    dill.dump(power_list, dill_file)
                print('Made ' + str(turbine) + '_Power List (for all inflows) pickle', datetime.now() - startTime)    

        else:            
            for inflow in inflows:
                print("Starting ", str(inflow), datetime.now() - startTime)

                turbine == 'T2'
                f_interp = surrogate.make_surrogate_L10(inflow, cases, seeds, turbine, x_cross, y_down, outfile = "BINARY")
                with open(str(path_to_surrogate) + str(inflow) + '_' + str(turbine) + '_L10', "wb") as dill_file:
                    dill.dump(f_interp, dill_file)
                print('Made ' + str(inflow) + '_' + str(turbine) + '_L10 pickle', datetime.now() - startTime)   
        
            
                f_interp = surrogate.make_surrogate_power(inflow, cases, seeds, turbine, x_cross, y_down, outfile = "BINARY")
                with open(str(path_to_surrogate) + str(inflow) + '_' + str(turbine) + '_Power', "wb") as dill_file:
                    dill.dump(f_interp, dill_file)
                print('Made ' + str(inflow) + '_' + str(turbine) + '_Power pickle', datetime.now() - startTime)   

                print("Time to Finish: ", inflow, datetime.now() - startTime)

