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

    def __init__(self, FF_timestep, m_c, d_c, m_s, m_p, N, g, beta, L_c, L_s, L_p, rho, C, e):
        
        '''Instantiate LayoutOptimization object and parameter values.'''

        self.FF_timestep = FF_timestep
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
        self.C = C                          # capacity from bearing, N
        self.e = e                          # constant for roller bearings
        

    def pl(self, R, torque, m_y, m_z):
    
        '''Calculation of forces via Eq 9 in Guo et al. 2015, Wind Energy'''

        omega = 2 * np.pi/self.N            # angle from planetary gear center to bearing center
        temp = np.zeros((2 * self.N,))      # sun planet mesh loads

        for i in range(self.N): 
            temp[0] = temp[0] - R[i] * np.cos(omega * i + self.beta[i]) + (self.m_p * self.g * np.sin(omega * i + self.beta[i]) ** 2) * np.cos(self.rho)
            temp[1] = temp[1] + R[i]
            temp[2] = temp[2] + R[i] * np.sin(omega * i + self.beta[i]) + self.m_c * self.g * np.sin(omega * i + self.beta[i]) * np.cos(omega * i + self.beta[i]) * np.cos(self.rho)

        z = np.zeros((self.N,))
        z[0] = temp[0] + (self.m_c * self.g * np.cos(self.rho)) + ((0.5 * self.m_s * self.g * self.L_s * np.cos(self.rho)) / self.L_c) + (m_y / self.L_c)        
        z[1] = temp[1] * self.d_c - torque 
        z[2] = -temp[2] * self.L_p - m_z #z not included, but could, here (currently yaw moment is equal to 0) 

        return z #this value makes up the Jacobian


    def fsolve_TE3(self, fun, x0, tol, max_fsolve_iter, torque, m_y, m_z):
        
        '''Solver function. Can be replaced by other solvers.'''
        
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

                J = (self.pl(xplus, torque, m_y, m_z) - self.pl(xminus, torque, m_y, m_z)) / H[i] 
                Jacobian[:,i] = J 
                
            f = self.pl(x0, torque, m_y, m_z).conj().T 
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
            
        return x1 #planet forces


    def calc_planet_forces(self, planet_speed, torque, m_y, m_z):
        
        '''Calculate planet bearing forces'''

        planet_forces = np.zeros((len(torque), self.N)) 
        R0 = np.asarray([1E3,1E3,1E3])    #initial guess for ring planet mesh forces

        for j in range(len(torque)):
            planet_forces[j,:] = self.fsolve_TE3(self.pl, R0, 0.01, 200, torque[j], m_y[j], m_z[j]) # planet forces (tangential) (requires inputs/outputs in N-m)
            R0 = planet_forces[j,:] #updates initial guess for faster computation
        return planet_forces, planet_speed


    def calc_L10(self, planet_forces, planet_speed): 
        
        '''Calculate L10 for the planet bearings given the planet bearing force and speed time history'''
        L10 = []
        T = self.FF_timestep / (len(planet_speed) * self.FF_timestep - self.FF_timestep) # contribution
        L10 = [float(1/(T/((10**6)/(60*j)*(self.C/abs(i))**self.e))) for i,j in zip(planet_forces, planet_speed)] #in hours

#        L10 = sum(L10)
                     
        return L10
    

    def plot_timeseries(self, x, ylim_min, ylim_max, xlabel, ylabel, filename):    
        '''Plot timeseries of torque, planet forces, or planet speed'''

        plt.scatter(range(len(x)), x, s = 1)
        # plt.title(str(title))
        plt.ylim(ylim_min, ylim_max)
        plt.xlabel(str(xlabel))
        plt.ylabel(str(ylabel))
        plt.savefig(str(filename))
        plt.show()


    def plot_scatter(self, x, y, xlim_min, xlim_max, ylim_min, ylim_max, xlabel, ylabel, filename):    
        '''Plot timeseries of torque, planet forces, or planet speed'''

        plt.scatter(x, y, s = 1)
        # plt.title(str(title))
        plt.xlabel(str(xlabel))
        plt.xlim(xlim_min, xlim_max)
        plt.ylim(ylim_min, ylim_max)
        plt.xticks(rotation=90)
        plt.ylabel(str(ylabel))
        plt.savefig(str(filename))
        plt.show()


    def plot_loads(self, x1, x2, x3, x4, x1_label, x2_label, x3_label, x4_label, xlabel, ylabel, filename):
        '''Plot torque and non-torque loads'''

        plt.plot(range(len(x1)), x1, alpha=0.5, label = str(x1_label)) 
        plt.plot(range(len(x2)), x2, alpha=0.5, label = str(x2_label))  
        plt.plot(range(len(x3)), x3, alpha=0.5, label = str(x3_label)) 
        plt.plot(range(len(x4)), x4, alpha=0.5, label = str(x4_label))
        
        plt.tight_layout()
        plt.xlabel(str(xlabel))
        plt.ylabel(str(ylabel))
        plt.legend(loc='lower right')
        plt.savefig(str(filename))
        plt.show()
