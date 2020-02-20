#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:49:46 2019

@author: cclark2
"""

#%%LOAD MODULES

#Internal Modules
import pandas as pd
import math
from datetime import datetime
from scipy import stats

#External Modules
import drivetrain

# Define turbine and drivetrain characteristics

FF_timestep = 0.025                         #FAST.Farm timestep

m_c = 5800                                    # mass of carrier, kg cc
d_c = 863E-3                                 # center distance, m
#pa = 20/180 * math.pi                       # pressure angle, radians
m_s = 18E3                                   # mass of shaft, kg
m_p = 1500                                   # mass of planet, kg
N = 3                                       # number of planets
g = 9.81*math.cos(5/180*math.pi)            # 5° is the bedplate tilting angle rho cc
beta = [0.0,0.0,0.0]                        # mounting angle variation from standard division
L_c = 4.25                                   # distance from hub to carrier center, m
L_s = 3.22                                   # length of shaft, m
L_p = 5                                      # distance between main shaft and pl bearing system center, m
N_r = 56                                    # number of teeth in ring gear
N_p = 17                                    # number of teeth in planet gear
rho = 5                                     # bedplate tilting angle

C = 4730000                                 #capacity from bearing, N
e = 10/3                                    #constant for roller bearings

#Define load channel inputs
frame = pd.read_csv('/Users/cclark2/Desktop/Case_A_8/Seed_1/Case_A_8_dist630.0_off0.0/Case_A_8_dist630.0_off0.0.T2.out', delim_whitespace=True, header = [0,1], skiprows=6, error_bad_lines=False)
planet_speed = frame['RotSpeed'].apply(lambda x: x * (FF_timestep/60) * abs(1-(N_r/N_p))).values #translate rotor > planet speed in rotations per minute (rpm)
torque = frame['RotTorq'].apply(lambda x: x * 1E3).values # in N-m
m_y = frame['LSSGagMys'].apply(lambda x: x * 1E3).values # in N-m
m_z = frame['LSSGagMzs'].apply(lambda x: x * 1E3).values # in N-m

#t = np.arange(0, 21, 0.01) #s
#omega = 8.52385848

# -------AMP--------
# planet_speed = 5.515058825 * np.sin(omega*t) + 19.55473415 #mean
# planet_speed = 11.03011765 * np.sin(omega*t) + 19.55473415 #max
# planet_speed = 1 * np.sin(omega*t) + 19.55473415 #min

# torque = 695923.3922075948 * np.sin(omega*t) + 1391846.78441519 #mean
# torque = 3849000 * np.sin(omega*t) + 1391846.78441519 #max
# torque = 250000 * np.sin(omega*t) + 1391846.78441519 #min
# torque = 1 * np.sin(omega*t) + 1391846.78441519 #min

# m_y = 3999500 * np.sin(omega*t) - 2301860.71486219 #mean
# m_y = 7999000 * np.sin(omega*t) - 2301860.71486219 #max
# m_y = 1 * np.sin(omega*t) - 2301860.71486219 #min

# m_z = 1924500 * np.sin(omega*t) + 77878.60127798 #mean
# m_z = 8659000 * np.sin(omega*t) + 77878.60127798 #max
# m_z = 1 * np.sin(omega*t) + 77878.60127798 #min

# -------MEAN-------- 
# planet_speed = 5.515058825 * np.sin(omega*t) + 19.55473415 #mean
# planet_speed = 5.515058825 * np.sin(omega*t) + 27.23117647 #max
# planet_speed = 5.515058825 * np.sin(omega*t) + 16.20105882 #min

# torque = 695923.3922075948 * np.sin(omega*t) + 1391846.78441519 #mean
# torque = 695923.3922075948 * np.sin(omega*t) + 3849000 #max
# torque = 695923.3922075948 * np.sin(omega*t) + 4.423e-10 #min
# torque = 695923.3922075948 * np.sin(omega*t) + 500000 #min

# m_y = 3999500 * np.sin(omega*t) - 2301860.71486219 #mean
# m_y = 3999500 * np.sin(omega*t) + 1789000 #max
# m_y = 3999500 * np.sin(omega*t) - 6210000 #min

# m_z = 1924500 * np.sin(omega*t) + 77878.60127798 #mean
# m_z = 1924500 * np.sin(omega*t) + 4149000 #max
# m_z = 1924500 * np.sin(omega*t) - 4510000 #min

# -------FREQUENCY--------
# omega = 7.062, 8.52385848, 11.87

#planet_speed = 5.515058825 * np.sin(8.52385848*t) + 19.55473415 #mean
# planet_speed = 5.515058825 * np.sin(11.87*t) + 19.55473415 #max
# planet_speed = 5.515058825 * np.sin(7.062*t) + 19.55473415 #min

#torque = 695923.3922075948 * np.sin(8.52385848*t) + 1391846.78441519 #mean
# torque = 695923.3922075948 * np.sin(11.87*t) + 1391846.78441519 #max
# torque = 695923.3922075948 * np.sin(7.062*t) + 1391846.78441519 #min

# m_y = 3999500 * np.sin(8.52385848*t) - 2301860.71486219 #mean
# m_y = 3999500 * np.sin(11.87*t) - 2301860.71486219 #max
# m_y = 3999500 * np.sin(7.062*t) - 2301860.71486219 #min

# m_z = 1924500 * np.sin(8.52385848*t) + 77878.60127798 #mean
# m_z = 1924500 * np.sin(11.87*t) + 77878.60127798 #max
# m_z = 1924500 * np.sin(7.062*t) + 77878.60127798 #min

#------------CONSTANT------------------
# planet_speed = np.full((m_z.shape), 19.55633078)
# torque = np.full((planet_speed.shape), 1392201.74747816) 
#m_y = np.full((planet_speed.shape), -2302301.76495669)
#m_z = np.full((planet_speed.shape), 77878.60127798)

#-----------------------------------------------------------------------------------------
# planet_speed = np.full((m_z.shape), 19.55633078)
# planet_speed = np.arange(16.0, 28.0, 0.1)
# planet_speed = frame['RotSpeed'].apply(lambda x: x * abs(1-(N_r/N_p))).values #translate rotor > planet speed in rotations per minute (rpm) (also carrier speed (assumption) (w_c)) w_p = (1-(N_r/N_p))*w_c
# torque = np.full((planet_speed.shape), 1392201.74747816) 
# torque = np.arange(0, 3849000, 32075)
# torque = np.zeros(planet_speed.shape) 
# torque = frame['RotTorq'].apply(lambda x: x * 1E3).values # in N-m
# m_y = np.full((80001, 1), -2302301.76495669)
# m_y = np.full((planet_speed.shape), -2302301.76495669)
# m_y = np.arange(-6210000, 6210000, 103500)
# m_y = np.arange(-6210000, 1789001, 66658.3416667)
# m_y = np.zeros(planet_speed.shape) 
# m_y = frame['LSSGagMys'].apply(lambda x: x * 1E3).values # in N-m
# m_z = np.full((80001, 1), 77900.17874652)
# m_z = np.arange(-4510000, 4149001, 72158.3416667)
# m_z = np.full((planet_speed.shape), 77900.17874652)
# m_z = np.zeros(planet_speed.shape)
# m_z = frame['LSSGagMzs'].apply(lambda x: x * 1E3).values # in N-m
# wind_vel = frame['Wind1VelX']

#%%

# Instantiate the drivetrain object
drivetrain_model = drivetrain.Drivetrain(
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
)

# Calculate planet forces
startTime = datetime.now()
planet_forces, planet_speed = drivetrain_model.calc_planet_forces(planet_speed, torque, m_y, m_z)
print('Planet Forces Calculated: ', datetime.now() - startTime)

# Calculate L10
startTime = datetime.now()
#L10, eq_load, eq_speed, eq_load_case, eq_speed_case = drivetrain_model.calc_L10(planet_forces, planet_speed)
L10 = drivetrain_model.calc_L10(planet_forces, planet_speed)

print('L10 Calculated: ', datetime.now() - startTime)

#%%
# Get statistics

print('Torque Statistics:', stats.describe(torque))
print('M_y Statistics:', stats.describe(m_y))
print('M_z Statistics:', stats.describe(m_z))
# print('Wind Velocity Statistics:', stats.describe(wind_vel))
print('Planet Force Statistics:', stats.describe(planet_forces))
print('Planet Speed Statistics:', stats.describe(planet_speed))

# Plot

drivetrain_model.plot_loads(torque, m_y, m_z, planet_forces, 'Torque', 'M_y', 'M_z', 'Planet Forces', 'Time (s)', 'Loads (N-m)', 'SineAverageTorqueMean.png')

#drivetrain_model.plot_timeseries(L10, 0, 1e30, 'Timestep (0.025s)', 'L10 (hrs)', 'RelationshipMzPlanetSpeed.png')    

# drivetrain_model.plot_timeseries(planet_speed, 'Timestep (0.025s)', 'Planet Speed (rpm)', 'RelationshipMzPlanetSpeed.png')    
# drivetrain_model.plot_timeseries(planet_forces, 'Timestep (0.025s)', 'Planet Force (N-m)', 'RelationshipMzPlanetForce.png') 
# drivetrain_model.plot_loads(torque, m_y, m_z, 'Torque', 'M_y', 'M_z', 'Timestep (0.025s)', 'Main Shaft Loads (N-m)', 'RelationshipMzLoads.png')

# drivetrain_model.plot_scatter(planet_speed, planet_forces, min(planet_speed), 28, 0, 1000000, 'Planet Speed (rpm)', 'Planet Forces (N-m)', 'RelationshipPlanetSpeedPlanetForce.png')




#%%
# Correlation matrix (generate and visualize)

# frame = pd.read_csv('/Users/cclark2/Desktop/Case_A_8/Seed_1/Case_A_8_dist630.0_off0.0/Case_A_8_dist630.0_off0.0.T2.V4.txt', delim_whitespace=True, error_bad_lines=False) #delim_whitespace=True, header = [0,1], skiprows=4, error_bad_lines=False)

# corr = frame.corr()
# ax = sns.heatmap(
#     corr, 
#     vmin=-1, vmax=1, center=0,
#     cmap=sns.diverging_palette(20, 220, n=200),
#     square=True
# )

# ax.grid(False, 'major')
# ax.grid(True, 'minor')
# ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
# ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)
# ax.set_xticklabels(
#     ax.get_xticklabels(),
#     rotation=45,
#     horizontalalignment='right'   
# );