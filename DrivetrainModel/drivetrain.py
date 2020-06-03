#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 10:10:11 2019
@author: cclark2
"""

import numpy as np
import matplotlib.pyplot as plt


class Drivetrain():
    """
    Base drivetrain class that calculates forces and L10 lifetime for planet bearings. 
    """

    def __init__(self, FF_timestep, m_c, d_c, m_s, m_p, N, g, beta, L_c, L_s, L_p, rho, C, e, omega):
        
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
        self.omega = omega                  # angle from planetary gear center to bearing center


    def pl(self, R, alpha, torque, m_y, m_z):
    
        '''Calculation of forces via Eq 9 in Guo et al. 2015, Wind Energy. Inconsistencies between this code and the equations in the paper have been confirmed with Guo in January 2020 and again June 2020'''

        temp = np.zeros((2 * self.N,))      # sun planet mesh loads

        for i in range(self.N): # parts of Eq. 9 that require summation over the bearings
            temp[0] = temp[0] - (R[i] * np.cos(self.omega * i + alpha + self.beta[i]) - (self.m_p * self.g * np.sin(self.omega * i + alpha + self.beta[i]) ** 2) * np.cos(self.rho))
            temp[1] = temp[1] + R[i]
            temp[2] = temp[2] + R[i] * np.sin(self.omega * i + alpha + self.beta[i]) + self.m_c * self.g * np.sin(self.omega * i + alpha + self.beta[i]) * np.cos(self.omega * i + alpha + self.beta[i]) * np.cos(self.rho)

        z = np.zeros((self.N,))
        z[0] = temp[0] + (self.m_c * self.g * np.cos(self.rho))*(self.L_c/self.L_p) + ((0.5 * self.m_s * self.g * self.L_s * np.cos(self.rho)) / self.L_p) - (m_y / self.L_p)        
        z[1] = temp[1] * self.d_c - torque 
        z[2] = -temp[2] * self.L_p - m_z  

        return z #this value makes up the Jacobian       
    

    def fsolve_TE3(self, fun, x0, tol, max_fsolve_iter, alpha, torque, m_y, m_z):
        
        '''Function that solves for the Jacobian matrix. Can be replaced by other built-in solvers.'''
        
        fsolve_iter = 0 
        error = 10 
        n = len(x0) 
        x0 = np.asarray(x0).ravel().conj().T 

        E = []
        H = [] 
        Jacobian = np.zeros((n,n)) 
        while fsolve_iter < max_fsolve_iter and error > tol: 
            fsolve_iter += 1 
            
            for i in range(n):
                E.append([x0[i] * 0.95, x0[i] * 1.05]) 
                H.append(E[i][1] - E[i][0])

                xplus = x0.copy() 
                xminus = x0.copy()
                
                xplus[i]  = E[i][1] 
                xminus[i] = E[i][0] 

                J = (self.pl(xplus, alpha, torque, m_y, m_z) - self.pl(xminus, alpha, torque, m_y, m_z)) / H[i] 
                Jacobian[:,i] = J 
                
            f = self.pl(x0, alpha, torque, m_y, m_z).conj().T 
            error = np.sqrt(np.sum(f ** 2)) 
            inv_result = np.linalg.solve(Jacobian, f)
            x1 = x0 - inv_result
            x0 = x1

        #max iteration error
        if fsolve_iter == max_fsolve_iter:
            print('Max fsolve_iterations reached');
            x1 = np.nan + np.zeros_like(x1)
        
        #imaginary number error
        if np.abs(np.sum(np.imag(x0))) > 0:
            print('Imaginary');
            x1 = np.nan + np.zeros_like(x1)
            
        return x1 #planet forces (f_t, i)


    def calc_planet_forces(self, planet_speed, alpha, torque, m_y, m_z):
        
        '''Calculate planet bearing forces: calls the fsolve_TE3 solver which builds and solves the Jacobian matrix made of planetary forces calculated in the pl function'''

        planet_forces = np.zeros((len(torque), self.N)) 
        R0 = np.asarray([1E3,1E3,1E3])    #initial guess for ring planet mesh forces

        for j in range(len(torque)):
            planet_forces[j,:] = self.fsolve_TE3(self.pl, R0, 0.01, 200, alpha[j], torque[j], m_y[j], m_z[j]) # planet forces (tangential) (requires inputs/outputs in N-m)
            R0 = planet_forces[j,:] #updates initial guess for faster computation
        
        #Define tangential forces for a single planet bearing (pb = 0)
        f_t = planet_forces[:, 0]
        
        #Calculate radial forces for a single planet bearing (come back and make this for all bearings so that planet forces can be calc'd for all three in future)
        f_r = -self.m_p * self.g * np.sin(self.omega + alpha + self.beta[0]) * np.cos(self.rho)
        
        #Combine tangential and radial forces via sum of squares
        F = [((ii**2 + jj**2)**(0.5)) for ii, jj in zip(f_t, f_r)]
                
        return F, f_t, f_r, planet_speed #only return one bearing's forces, if you are assuming even load distribution 


    def calc_L10(self, planet_forces, planet_speed): 
        
        '''Calculates L10 for the planet bearings given the planet bearing force and speed time history. L10 life is the time that 90% of a group of bearings will exceed without failing by rolling-element fatigue.'''

        T = self.FF_timestep / (len(planet_speed) * self.FF_timestep - self.FF_timestep) # fraction of total running time at a given load and speed
        L10 = [T/((10**6/(60*i))*(self.C/abs(j))**self.e) for i,j in zip(planet_speed, planet_forces)] # life at each load and speed combination, in hours
        # 10**6 million race revolutions conversion factor
        # 60 min/hr conversion factor
        # i: planet speed 
        # C: bearing basic dynamic load rating or capacity, N (the load that a bearing can carry for 1 million inner-race revolutions with a 90% probability of survival)
        # e: load-life exponent (determined by Lundberg and Palmgren to be 3 for ball bearings and 10/3 for cylindrical roller bearings)
        L10_total = 1/sum(L10) # total life in hours over varying loads and speeds, Eq. 14 from "Rolling Bearing Life Prediction, Theory, and Application," by Zaretsky (2016)
                     
        return L10, L10_total #returns L10, or the vector of life calculations at each point in the time series
  

    def plot_loads(self, x1, x2, x3, x4, x1_label, x2_label, x3_label, x4_label, xlabel, ylabel):
        
        '''Plot torque and non-torque loads'''

        plt.plot(range(len(x1)), x1, alpha=0.5, label = str(x1_label)) 
        plt.plot(range(len(x2)), x2, alpha=0.5, label = str(x2_label))  
        plt.plot(range(len(x3)), x3, alpha=0.5, label = str(x3_label)) 
        plt.plot(range(len(x4)), x4, alpha=0.5, label = str(x4_label))
        
        plt.tight_layout()
        plt.xlabel(str(xlabel))
        plt.ylabel(str(ylabel))
        plt.legend(loc='lower right')
#        plt.savefig(str(filename))
        plt.show()
        
