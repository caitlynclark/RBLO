#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:49:46 2019
@author: cclark2
"""


#Internal Modules
import pandas as pd
import math
from datetime import datetime
from scipy import stats

#External Modules
import drivetrain

# Define turbine and drivetrain characteristics
FF_timestep = 0.025                         # FAST.Farm timestep
m_c = 5800                                  # mass of carrier, kg cc
d_c = 863E-3                                # center distance, m
#pa = 20/180 * math.pi                      # pressure angle, radians
m_s = 18E3                                  # mass of shaft, kg
m_p = 1500                                  # mass of planet, kg
N = 3                                       # number of planets
g = 9.81 * math.cos(5/180*math.pi)          # 5° is the bedplate tilting angle rho cc
beta = [0.0, 120.0, 240.0]                        # mounting angle variation from standard division
L_c = 4.25                                  # distance from hub to carrier center, m
L_s = 3.22                                  # length of shaft, m
L_p = 5                                     # distance between main shaft and pl bearing system center, m
N_r = 56                                    # number of teeth in ring gear
N_p = 17                                    # number of teeth in planet gear
rho = 5*(math.pi/180)                       # bedplate tilting angle, radians
C = 4730000                                 # capacity from bearing, N
e = 10/3                                    # constant for roller bearings
omega = 2 * math.pi/N                       # angle from planetary gear center to bearing center


#Define load channel inputs
frame = pd.read_csv('/Users/cclark2/Desktop/Case_A_8/Case_A_8_dist630.0_off0.0_T2.csv', delim_whitespace=False, header = [0], skiprows=[1], error_bad_lines=False)
alpha = frame['Azimuth'].apply(lambda x: x * (math.pi/180)).values
planet_speed = frame['RotSpeed'].apply(lambda x: x * abs(1-(N_r/N_p))).values #translate rotor speed to planet speed (rpm)
torque = frame['RotTorq'].apply(lambda x: x * 1E3).values # in N-m
m_y = frame['LSSGagMys'].apply(lambda x: x * 1E3).values # in N-m
m_z = frame['LSSGagMzs'].apply(lambda x: x * 1E3).values # in N-m

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
                        omega = omega
)

# Calculate planet forces
startTime = datetime.now()
F, f_t, f_r, planet_speed = drivetrain_model.calc_planet_forces(planet_speed, alpha, torque, m_y, m_z)
#planet_forces, planet_speed = drivetrain_model.calc_planet_forces(planet_speed, alpha, torque, m_y, m_z)
print('Planet Forces Calculated: ', datetime.now() - startTime)

# Calculate L10
startTime = datetime.now()
L10, L10_total = drivetrain_model.calc_L10(F, planet_speed) #L10 is a vector representing the life at each timestep, while L10 is the total life of the bearing given varying loads and speeds

#L10, L10_total = drivetrain_model.calc_L10(planet_forces, planet_speed) #L10 is a vector representing the life at each timestep, while L10 is the total life of the bearing given varying loads and speeds
print('L10 Calculated: ', L10_total, datetime.now() - startTime)

# Plot
drivetrain_model.plot_loads(torque, m_y, m_z, F, 'Torque', 'M_y', 'M_z', 'Planet Forces', 'Time (s)', 'Loads (N-m)')
