#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 13:39:26 2019

@author: cclark2
"""

'''This code optimizes the layout of a wind farm based on calculations of the L10 and power.'''

#%% 1) Import modules and define functions and path

import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import math
import dill as pickle
from scipy.optimize import minimize
 

#%% 2) Load L10 and power interpolant surfaces

path = '/Users/cclark2/floris/examples/'

with open(str(path) + 'Case_A_8_L10_T1', "rb") as dill_file:
     L10_T1_surface = pickle.load(dill_file)
     
with open(str(path) + 'Case_A_8_L10_T2', "rb") as dill_file:
     L10_T2_surface = pickle.load(dill_file)

with open(str(path) + 'Case_A_8_T1_GenPwr', "rb") as dill_file:
     Power_T1_surface = pickle.load(dill_file)

with open(str(path) + 'Case_A_8_T2_GenPwr', "rb") as dill_file:
     Power_T2_surface = pickle.load(dill_file)

#%%

class Optimization():
    """
    Base optimization class.

    Args:
        properties

    Returns:
        Optimization: An instantiated Optimization object.
    """

    def __init__(self, wd, freq, ws):
        """
        Instantiate Optimization object and its parameters.
        """

        self.wd = wd
        self.freq = freq
        self.ws = ws

    def _norm(self, val, x1, x2):
        return (val - x1)/(x2 - x1)
    
    def _unnorm(self, val, x1, x2):
        return np.array(val)*(x2 - x1) + x1

       
class RBLO(Optimization):
    """
    Sub class of the :py:class`floris.tools.Optimization`
    object class that performs reliability-based layout optimization.
    """

    def __init__(self, N, D, layout_x, layout_y, x_min, x_max, 
                           y_min, 
                           y_max, 
                           bndx_min, 
                           bndy_min, 
                           bndx_max, 
                           bndy_max, 
                           wd,
                           ws,
                           freq,
                           x0=None,
                           bnds=None,
                           min_dist=None,
                           opt_method='SLSQP',
                           cons = None,
                           opt_options=None
):


        super().__init__(wd, freq, ws)

        self.N = N
        self.D = D
        self.layout_x = layout_x
        self.layout_y = layout_y
        self.x0 = layout_x + layout_y
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.bndx_min = bndx_min
        self.bndy_min = bndy_min
        self.bndx_max = bndx_max
        self.bndy_max = bndy_max
        self.min_dist = min_dist 
        self.opt_method = 'SLSQP'
        self.constraints = cons

        if opt_options is None:
            self.opt_options = {'maxiter': 100, 'disp': True, \
                         'iprint': 2, 'ftol': 1e-8}
        else:
            self.opt_options = opt_options
            
        self._set_opt_bounds()
        self._generate_constraints()
        self.cost_initial, self.power_initial, self.replacements = self.calc_cost_power(self.layout_x + self.layout_y)


    def rotate_array(self, xy, radians):
        """Use numpy to build a rotation matrix and take the dot product about the origin."""
        x, y = xy
        c, s = np.cos(radians), np.sin(radians)
        j = np.matrix([[c, s], [-s, c]])
        m = np.dot(j, [x, y])
    
        return float(m.T[0]), float(m.T[1])
    
    
    def calc_L10_Power(self, layout_x, layout_y):    
        '''Assign L10 values to all turbines'''
    
        Turbine_L10s = np.zeros((len(self.wd), self.N))
        Turbine_PowerVals = np.zeros((len(self.wd), self.N))
    
        for aa, wind_direction in enumerate(self.wd):    
            
            # Define the wind array and rotate towards the wind direction
            turbs = []
            theta = math.radians(wind_direction) #rotate grid 90 deg
            for ii in list(zip(layout_x, layout_y)):
                turbs.append((self.rotate_array(ii, theta)))
            
            layout_x = [ii[0] for ii in turbs]
            layout_y = [ii[1] for ii in turbs]
            
            turb_label = range(0, len(layout_x))    
    
            # Assign L10 values to all turbines
            # GB: DEFAULT VALUE FOR L10 SHOULD BE T1 VALUE LIKE IT IS FOR POWER
            L10_matrix = np.full( (len(turb_label), len(turb_label)), float(L10_T1_surface[0])) #np.full( (len(turb_label), len(turb_label)), float(L10_T1_surface[0]) )
            Power_matrix = np.full( (len(turb_label), len(turb_label)), float(Power_T1_surface[0]) )
    
            turbs_to_analyse = []
            
            # Delete pairs that fall outside space constraints
            # GB: THERE IS NO LOGIC HERE TO DECIDE WHICH IS T1 AND WHICH IS T2 BASED ON THE INFLOW DIRECTION
            # I DON'T THINK THIS COMBINATIONS IS GIVING YOU WHAT YOU WANT.  WITH THIS APPROACH EVERYTHING IS A T2 AND NOTHING IS A T1
            for ii in list(combinations(turb_label, 2)):
                turb1 = ii[0] 
                turb2 = ii[1]

                # GB: ARE THESE DISTANCES CONSISTENT WITH THE DOMAIN LIMITS OF THE SURROGATE MODEL?
                if abs(turbs[turb2][1] - turbs[turb1][1]) <= 1260 and abs(turbs[turb2][0] - turbs[turb1][0]) <= 252:
                    turbs_to_analyse.append(ii)
                    
            # Assign L10 & power vals; turbs are ordered in proximity to the incoming wind direction
            # GB: AGAIN, YOU ARE ASSIGNING T1 AND T2 BUT THERE IS NO LOGIC TO DETERMINE THEM, JUST INDEX ORDERING IN YOUR ARRAY
            for ii in turbs_to_analyse:
                turb1 = ii[0]
                turb2 = ii[1]
                                
                L10_matrix[turb1, turb2] = L10_T2_surface(layout_y[turb2]-layout_y[turb1],layout_x[turb2]-layout_x[turb1])
                Power_matrix[turb1, turb2] = Power_T2_surface(layout_y[turb2]-layout_y[turb1],layout_x[turb2]-layout_x[turb1]) #kW/hr
                        
            # Find the minimum of each row (this is the minimum )
            # GB: SHOULD BE SIMPLY SOMETHING LIKE THIS:
            #------
            Turbine_L10s[aa,:] = L10_matrix.min(axis=0)
            Turbine_PowerVals[aa,:] = Power_matrix.min(axis=0)
            #------
            for column in L10_matrix.T:
                Turbine_L10s[aa] = [min(column) for column in L10_matrix.T]
                # ARGMIN OF COLUMN MEANS ARGMIN OF THE L10_MATRIX COLUMN, WHICH IS NOT WHAT YOU WANT FOR POWER
                Turbine_PowerVals[aa] = [Power_matrix.T[ii][np.argmin(column)] for ii, column in enumerate(L10_matrix.T)]

        # Determine L10 for each turbine based on weighted averaging the wind directions (in hours)
        # GB: THIS NEEDS TO USE THE LINEAR DAMAGE SUMMATION AS WE DISCUSSED.
        L10_AllWindConds = []
        for column in Turbine_L10s.T:
            L10_AllWindConds.append(float(sum([a*b for a,b in zip(column, self.freq)])) ) #1E5
        
        # Determine power for each turbine based on weighted averaging the wind directions (in hours)
        Power_AllWindConds = []
        for column in Turbine_PowerVals.T:
            Power_AllWindConds.append(float(sum([a*b for a,b in zip(column, self.freq)])) )   #kW/hr

        return L10_AllWindConds, Power_AllWindConds                     # in hours
    
    
    def calc_cost_power(self, layout):
        '''Failure costs over AEP'''
        
        layout_x = layout[0:int(len(layout)/2)]
        layout_y = layout[int(len(layout)/2):] 
        layout_x = [self._unnorm(valx, self.bndx_min, self.bndx_max) \
            for valx in layout_x]
        layout_y = [self._unnorm(valy, self.bndy_min, self.bndy_max) \
            for valy in layout_y]
        
        hours_lifetime = 20*365*24*(len(layout_x))                      # hours in 20 year lifespan
        downtime = 231                                                  # hours per replacement
        L10, Power = self.calc_L10_Power(layout_y, layout_x)            # L10 in hours, Power in kW/hr

        # GB: I DON'T SEE THE FLOOR FUNCTION WE DISCUSSED
        # GB: THIS MEANS EVERY TURBINE GEARBOX IS REPLACED INDEPENDENTLY.
        #     THAT IS A FINE ASSUMPTION, BUT THE PAPER WAS UNCLEAR ABOUT THIS DETAIL
        replacements = (sum([hours_lifetime/ii for ii in L10]))         # for all turbines in farm took away the "int" in front of the parenthesis bc potlly causing Exit mode 5
        replacement_cost = 715920                                       # cost of GB replacement (ref: Yeng, Sheng, and Court: $628,000 USD2011, $715920 USD2019)
        power_hours = hours_lifetime-(downtime*replacements)    
        cost = replacement_cost*replacements                            # for all turbines, failure cost
        # GB: THIS IS AVERAGE POWER PER TURBINE OVER ALL INFLOWS.
        power = (sum(Power)/len(Power)) * power_hours                   # kW/hr*hr for all turbines
        return cost, power, replacements


    def objective(self, layout):
        '''Failure costs over AEP'''
        
        cost, power, replacements = self.calc_cost_power(layout)        # L10 in hours, Power in kW/hr
#        return (1E12/power)                                            # power-only objective
        return ((0.9*cost)/(0.1*power))*1E03                            # multi-objective for Pareto
#        return (cost/1E10)                                             # cost-only objective


    def _space_constraint(self, x_in, min_dist):
        x = x_in[0:self.N] 
        y = x_in[self.N:]
        

        dist = [np.sqrt((abs(x[i])-abs(x[j]))**2 + (abs(y[i])-abs(y[j]))**2) \
                for i in range(self.N) \
                for j in range(self.N) if i != j]
        return np.min(dist) - self._norm(min_dist, self.bndx_min, self.bndx_max)


    def _generate_constraints(self):

        tmp1 = {'type': 'ineq','fun' : lambda x,*args: \
                self._space_constraint(x, self.min_dist), \
                'args':(self.min_dist,)}

        self.cons = [tmp1] 


    def _optimize(self):
        self.residual_plant = minimize(self.objective,
                                self.x0,
                                method=self.opt_method,
                                bounds=self.bnds,
                                constraints=self.cons,
                                options=self.opt_options) #

        opt_results = self.residual_plant.x #need the .x for later parsing
#        opt_results_all = self.residual_plant #need the .x for later parsing
#        print('opt_results_all: ', opt_results_all)
        return opt_results


    def _set_opt_bounds(self):
        bnds_x = [(self.x_min, self.x_max) for i in range(self.N)]
        bnds_y = [(self.y_min, self.y_max) for i in range(self.N)]
        self.bnds = bnds_x + bnds_y


    def optimize(self):
        """
        Find optimized layout of wind turbines for power production given
        fixed atmospheric conditions (wind speed, direction, etc.).
        
        Returns:
            opt_locs (iterable): optimized locations of each turbine.
        """
        print('=====================================================')
        print('Optimizing turbine layout...')
        print('Number of parameters to optimize = ', len(self.x0))
        print('=====================================================')

        opt_locs_norm = self._optimize()
        
        print('Optimization complete!')

        opt_locs = [[self._unnorm(valx, self.bndx_min, self.bndx_max) \
            for valx in opt_locs_norm[0:self.N]], \
            [self._unnorm(valy, self.bndy_min, self.bndy_max) \
            for valy in opt_locs_norm[self.N:2*self.N]]]
            
        return opt_locs #_norm



#%%
    def plot_layout_opt_results(self, jj):
        """
        Method to plot the old and new locations of the layout optimization.
        """

        locsx_old = [self._unnorm(valx, self.bndx_min, self.bndx_max) \
                     for valx in self.x0[0:self.N]]
        locsy_old = [self._unnorm(valy, self.bndy_min, self.bndy_max) \
                     for valy in self.x0[self.N:2*self.N]]
        locsx = [self._unnorm(valx, self.bndx_min, self.bndx_max) \
                 for valx in self.residual_plant.x[0:self.N]]
        locsy = [self._unnorm(valy, self.bndy_min, self.bndy_max) \
                 for valy in self.residual_plant.x[self.N:2*self.N]]
                
        cost, power, replacements = self.calc_cost_power(self.residual_plant.x)   #L10 in hours, Power in kW/hr
        print('Replacements: ', replacements, 'Cost: ', cost, 'Power: ', power)

        plt.figure(figsize=(8,7))
        fontsize= 16
        plt.plot(locsx_old, locsy_old, '3b', markersize=16)
        plt.plot(locsx, locsy, '2r', markersize=16)
        plt.xlabel('x (m)', fontsize=fontsize)
        plt.ylabel('y (m)', fontsize=fontsize)
        plt.grid()
        plt.axis([-100, self.bndx_max + 100, -100, self.bndy_max + 100])
        plt.tick_params(which='both', labelsize=fontsize)
        plt.legend(['Old locations', 'New locations'], loc='lower center', \
            bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=fontsize)
        plt.show()

