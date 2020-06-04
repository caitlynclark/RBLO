#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 08:48:39 2020

@author: cclark2
"""

import numpy as np
import math
from scipy.interpolate import interp2d
from scipy.optimize import fsolve
import pandas as pd
import os
import struct
import multiprocessing #import Pool
from itertools import repeat
import dill
import mpi4py.MPI as MPI
from mpi4py.futures import MPIPoolExecutor #executor map

MPI.pickle.__init__(dill.dumps, dill.loads)
num_cores = multiprocessing.cpu_count()

#%% READ BINARY FILES AND CONCATENATE SEEDS

class Surrogate():
    """
    Base drivetrain class that calculates forces and L10 lifetime for planet bearings. 
    """

    def __init__(self, FF_timestep, m_c, d_c, m_s, m_p, N, g, beta, L_c, L_s, L_p, rho, C, e, N_r, N_p, omega):
        
        '''Instantiate LayoutOptimization object and parameter values.'''

        self.FF_timestep = FF_timestep      # FAST.Farm timestep for outputs
        self.m_c = m_c                      # carrier mass
        self.d_c = d_c                      # center distance 
        self.m_s = m_s                      # shaft mass
        self.m_p = m_p                      # planet bearing mass
        self.N = N                          # number of planet bearings
        self.g = g                          # gravitational force
        self.beta = beta                    # mounting angle
        self.L_c = L_c                      # distance from main bearing to the carrier's center of gravity
        self.L_s = L_s                      # distance from main bearing to the main shaft's center of gravity
        self.L_p = L_p                      # distance from main bearing to the planet bearing's center of gravity
        self.rho = rho                      # bedplate tilting angle (if don't want to include, set to 0 degrees)
        self.C = C                          # bearing basic dynamic load rating or capacity, N (the load that a bearing can carry for 1 million inner-race revolutions with a 90% probability of survival)
        self.e = e                          # constant for roller bearings
        self.N_r = N_r                      # ring gear teeth (#)
        self.N_p = N_p                      # planet gear teeth (#) 
        self.omega = omega                  # bearing mount


    def fread(self, fid, n, type):
        fmt, nbytes = {'uint8': ('B', 1), 'int16':('h', 2), 'int32':('i', 4), 'float32':('f', 4), 'float64':('d', 8)}[type]
    
        return struct.unpack(fmt * n, fid.read(nbytes * n))


    def load_binary_output(self, filename):
        
        '''Ported from ReadFASTbinary.m by Mads M Pedersen, DTU Wind
        Info about ReadFASTbinary.m:
        % Author: Bonnie Jonkman, National Renewable Energy Laboratory
        % (c) 2012, National Renewable Energy Laboratory
        %
        %  Edited for FAST v7.02.00b-bjj  22-Oct-2012
        '''
        
    
        FileFmtID_WithTime = 1                                          # File identifiers used in FAST
        LenName = 10                                                    # number of characters per channel name
        LenUnit = 10                                                    # number of characters per unit name
    
        with open(filename, 'rb') as fid:
            FileID = self.fread(fid, 1, 'int16')                        # FAST output file format, INT(2)
    
            NumOutChans = self.fread(fid, 1, 'int32')[0]                # The number of output channels, INT(4)
            NT = self.fread(fid, 1, 'int32')[0]                         # The number of time steps, INT(4)
    
            if FileID == FileFmtID_WithTime:
                TimeScl = self.fread(fid, 1, 'float64')                 # The time slopes for scaling, REAL(8)
                TimeOff = self.fread(fid, 1, 'float64')                 # The time offsets for scaling, REAL(8)
            else:
                TimeOut1 = self.fread(fid, 1, 'float64')                # The first time in the time series, REAL(8)
                TimeIncr = self.fread(fid, 1, 'float64')                # The time increment, REAL(8)
    
            ColScl = self.fread(fid, NumOutChans, 'float32')            # The channel slopes for scaling, REAL(4)
            ColOff = self.fread(fid, NumOutChans, 'float32')            # The channel offsets for scaling, REAL(4)
    
            LenDesc = self.fread(fid, 1, 'int32')[0]                    # The number of characters in the description string, INT(4)
            DescStrASCII = self.fread(fid, LenDesc, 'uint8')            # DescStr converted to ASCII
            DescStr = "".join(map(chr, DescStrASCII)).strip()
    
            ChanName = []  # initialize the ChanName cell array
            for iChan in range(NumOutChans + 1):
                ChanNameASCII = self.fread(fid, LenName, 'uint8')       # ChanName converted to numeric ASCII
                ChanName.append("".join(map(chr, ChanNameASCII)).strip())
    
            ChanUnit = []  # initialize the ChanUnit cell array
            for iChan in range(NumOutChans + 1):
                ChanUnitASCII = self.fread(fid, LenUnit, 'uint8')       # ChanUnit converted to numeric ASCII
                ChanUnit.append("".join(map(chr, ChanUnitASCII)).strip()[1:-1])
    
    
            # Get the channel time series    
            nPts = NT * NumOutChans                                     # number of data points in the file
    
            if FileID == FileFmtID_WithTime:
                PackedTime = self.fread(fid, NT, 'int32')                    # read the time data
                cnt = len(PackedTime)
                if cnt < NT:
                    raise Exception('Could not read entire %s file: read %d of %d time values' % (filename, cnt, NT))
            PackedData = self.fread(fid, nPts, 'int16')                      # read the channel data
            cnt = len(PackedData)
            if cnt < nPts:
                raise Exception('Could not read entire %s file: read %d of %d values' % (filename, cnt, nPts))
    
        # Scale the packed binary to real data    
        data = np.array(PackedData).reshape(NT, NumOutChans)
        data = (data - ColOff) / ColScl
    
        if FileID == FileFmtID_WithTime:
            time = (np.array(PackedTime) - TimeOff) / TimeScl;
        else:
            time = TimeOut1 + TimeIncr * np.arange(NT)
    
        data = np.concatenate([time.reshape(NT, 1), data], 1)
    
        info = {'name': os.path.splitext(os.path.basename(filename))[0],
                'description': DescStr,
                'attribute_names': ChanName,
                'attribute_units': ChanUnit}
    
        return data, ChanName #data, info
    
    
    def concatenate_seeds(self, inflow, case, seeds, turbine, outfile):
        
        '''Concatenate seeds data from FAST.Farm into a single dataframe.'''
                
        if outfile == "BINARY":
            data = []
            for seed in seeds:
                file = '/projects/windse/kshaler/SystemsEngineering/GriddedDatabase/FFarm/NewCases/{0}/{1}/{2}/FFarm_mod.{3}.outb'.format(inflow, case, seed, turbine) 
                temp_data, channel = self.load_binary_output(file)
                print(str(seed) + 'temp_data size:' + print(temp_data.size))
                data.append(temp_data)
            concatenated_data = np.concatenate(data)
            frame = pd.DataFrame(concatenated_data, columns = channel)
    
    
        elif outfile == "ASCII":
            result_files = ['/projects/windse/kshaler/SystemsEngineering/GriddedDatabase/FFarm/NewCases/{0}/{1}/{2}/FFarm_mod.{3}.out'.format(inflow, case, seed, turbine) for seed in seeds]
            df_list = [pd.read_csv(file, delim_whitespace=True, header = [0,1], skiprows=6, error_bad_lines=False) for file in result_files]
            frame = pd.concat(df_list, axis=0, ignore_index = True)
    
        return channel, frame
            

    def pl(self, R, alpha, torque, m_y, m_z):
    
        '''Calculation of forces via Eq 9 in Guo et al. 2015, Wind Energy. Inconsistencies between this code and the equations in the paper have been confirmed with Guo in January 2020 and again June 2020'''

        temp = np.zeros((2 * self.N,))      # sun planet mesh loads

        for i in range(self.N): # parts of Eq. 9 that require summation over the bearings
            temp[0] = temp[0] - (R[i] * np.cos(self.omega * i + alpha + self.beta[i]) - (self.m_p * self.g * np.sin(self.omega * i + alpha + self.beta[i]) ** 2) * np.cos(self.rho))
            temp[1] = temp[1] + R[i]
            temp[2] = temp[2] + R[i] * np.sin(self.omega * i + alpha + self.beta[i]) + self.m_p * self.g * np.sin(self.omega * i + alpha + self.beta[i]) * np.cos(self.omega * i + alpha + self.beta[i]) * np.cos(self.rho)

        z = np.zeros((self.N,))
        z[0] = temp[0] + (self.m_c * self.g * np.cos(self.rho))*(self.L_c/self.L_p) + ((0.5 * self.m_s * self.g * self.L_s * np.cos(self.rho)) / self.L_p) - (m_y / self.L_p)        
        z[1] = temp[1] * self.d_c - torque 
        z[2] = -temp[2] * self.L_p - m_z  

        return z #this value makes up the Jacobian       
  

    def calc_L10(self, case, inflow, seeds, turbine, outfile): 
        
        '''Calculates L10 for the planet bearings given the planet bearing force and speed time history. L10 life is the time that 90% of a group of bearings will exceed without failing by rolling-element fatigue.'''

        print(str(case))
        channel, frame = self.concatenate_seeds(inflow, case, seeds, turbine, outfile)
        print(str(case) + ' frame.size: ' +  str(frame.size))

        try:
            alpha = frame['LSSTipPxa'].apply(lambda x: x * (math.pi/180)).values
        except:	
       	    print(str(case) + ' FAILED for LSSTipPxa')  
            alpha = self.concatenate_seeds(inflow, 'Case000', seeds, turbine, outfile)
 
        planet_speed = frame['RotSpeed'].apply(lambda x: x * abs(1-(self.N_r/self.N_p))).values #translate rotor speed to planet speed (rpm)
        torque = frame['RotTorq'].apply(lambda x: x * 1E3).values # in N-m
        
        try:
            m_y = frame['LSSGagMys'].apply(lambda x: x * 1E3).values # in N-m
        except:
            print(str(case) + ' FAILED for LSSGagMys')
            m_y = self.concatenate_seeds(inflow, 'Case000', seeds, turbine, outfile)    
       
        try:
            m_z = frame['LSSGagMzs'].apply(lambda x: x * 1E3).values # in N-m
        except:
            print(str(case) + ' FAILED for LSSGagMzs')
            m_z = self.concatenate_seeds(inflow, 'Case000', seeds, turbine, outfile)    
   
        planet_forces = np.zeros((len(torque), self.N)) 
        R0 = np.asarray([1E3,1E3,1E3])    #initial guess for ring planet mesh forces

        for j in range(len(torque)):
            planet_forces[j,:] = fsolve(self.pl, R0, args = (alpha[j], torque[j], m_y[j], m_z[j]), xtol = 0.01, maxfev = 200) # planet forces (tangential) (requires inputs/outputs in N-m)
            R0 = planet_forces[j,:] #updates initial guess for faster computation
        
        #Define tangential forces for a single planet bearing (pb = 0)
        f_t = planet_forces[:, 0]
        
        #Calculate radial forces for a single planet bearing (come back and make this for all bearings so that planet forces can be calc'd for all three in future)
        f_r = -self.m_p * self.g * np.sin(self.omega + alpha + self.beta[0]) * np.cos(self.rho)
        
        #Combine tangential and radial forces via sum of squares
        F = [((ii**2 + jj**2)**(0.5)) for ii, jj in zip(f_t, f_r)]

        #Calculate L10 using planet forces (tangential and radial) and planet speed
        T = self.FF_timestep / (len(planet_speed) * self.FF_timestep - self.FF_timestep) # fraction of total running time at a given load and speed
        L10 = [T/((10**6/(60*i))*(self.C/abs(j))**self.e) for i,j in zip(planet_speed, F)] # life at each load and speed combination, in hours
        L10_total = 1/sum(L10) # total life in hours over varying loads and speeds, Eq. 14 from "Rolling Bearing Life Prediction, Theory, and Application," by Zaretsky (2016)
        print(str(inflow), (case), 'L10_total: ', L10_total)         
            
        return L10_total #returns L10, or the vector of life calculations at each point in the time series


    def calc_power(self, case, inflow, seeds, turbine, outfile):
        
        channel, frame = self.concatenate_seeds(inflow, case, seeds, turbine, outfile)
        power = frame['GenPwr'].sum(axis = 0, skipna = True)/(len(frame['GenPwr'])/40)*3600
        
        return power


    def make_surrogate_L10(self, inflow, cases, seeds, turbine, x_cross, y_down, outfile):
        grid = []         
        with MPIPoolExecutor(max_workers=num_cores) as executor: 
            for L10 in executor.map(self.calc_L10, cases, repeat(inflow), repeat(seeds), repeat(turbine), repeat(outfile)):
                grid.append(L10)
        grid = np.reshape(grid, (len(x_cross), len(y_down)))
        f_interp = interp2d(x_cross, y_down, grid[:,:])

        return f_interp


    def make_surrogate_power(self, inflow, cases, seeds, turbine, x_cross, y_down, outfile):
        
        grid = []
        with MPIPoolExecutor(max_workers=num_cores) as executor: 
            for power in executor.map(self.calc_power, cases, repeat(inflow), repeat(seeds), repeat(turbine), repeat(outfile)):
                grid.append(power)
        grid = np.reshape(grid, (len(x_cross), len(y_down)))
        f_interp = interp2d(x_cross, y_down, grid[:,:])
        
        return f_interp 
    
    
