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
from matplotlib.patches import Rectangle
from itertools import combinations
import math
import dill as pickle
from scipy.optimize import minimize
 

#%% 2) Load L10 and power interpolant surfaces

#path = '/Users/cclark2/Desktop/RBLO/Surrogate_Models/surrogate_models/'
path = '/Users/cclark2/floris/examples/surrogate_models/'

#%%

class Optimization():
    """
    Base optimization class.

    Args:
        properties

    Returns:
        Optimization: An instantiated Optimization object.
    """

    def __init__(self, wd, wd_freq, ws, ws_freq):
        """
        Instantiate Optimization object and its parameters.
        """

        self.wd = wd
        self.wd_freq = wd_freq
        self.ws = ws
        self.ws_freq = ws_freq

    def _norm(self, val, x1, x2):
        return (val - x1)/(x2 - x1)
    
    def _unnorm(self, val, x1, x2):
        return np.array(val)*(x2 - x1) + x1

    def _calc_distance(self, x1, y1, x2, y2):
        return math.sqrt(math.pow(x2 - x1, 2) +
                         math.pow(y2 - y1, 2) * 1.0)
       
class RBLO(Optimization):
    """
    Sub class of the :py:class`floris.tools.Optimization`
    object class that performs reliability-based layout optimization.
    """

    def __init__(self, N, D, layout_x, 
                           layout_y, 
                           x_min, 
                           x_max, 
                           y_min, 
                           y_max, 
                           bndx_min, 
                           bndy_min, 
                           bndx_max, 
                           bndy_max, 
                           wd,
                           wd_freq,
                           ws,
                           ws_freq,
                           TI,
                           shear,
                           x0=None,
                           bnds=None,
                           min_dist=None,
                           opt_method='SLSQP',
                           cons = None,
                           opt_options=None
):


        super().__init__(wd, wd_freq, ws, ws_freq)

        self.N = N
        self.D = D
        self.layout_x = layout_x
        self.layout_y = layout_y
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.bndx_min = bndx_min
        self.bndy_min = bndy_min
        self.bndx_max = bndx_max
        self.bndy_max = bndy_max
        self.wd = wd
        self.wd_freq = wd_freq
        self.ws = ws
        self.ws_freq = ws_freq
        self.TI = TI
        self.shear = shear
        self.x0 = layout_x + layout_y
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
        self.cost_initial, self.power_initial, self.replacements, self.L10 = self.calc_cost_power(self.layout_x + self.layout_y)


    def rotate_array(self, xy, radians):
        """Use numpy to build a rotation matrix and take the dot product about the origin."""
        x, y = xy
        c, s = np.cos(radians), np.sin(radians)
        j = np.matrix([[c, s], [-s, c]])
        m = np.dot(j, [x, y])
    
        return float(m.T[0]), float(m.T[1])
    
        
    def calc_L10_Power(self, layout_x, layout_y, Vhub):    
        '''Assign L10 values to all turbines'''
        
#        with open(str(path) + 'T1_L10_List', "rb") as dill_file:
        with open(str(path) + 'Vhub_' + str(Vhub) + '_TI_' + str(self.TI) + '_Shear_' + str(self.shear) + '_L10_T1', "rb") as dill_file:
             L10_T1 = pickle.load(dill_file)
             
#        with open(str(path) + 'T1_Power_List', "rb") as dill_file:
        with open(str(path) + 'Vhub_' + str(Vhub) + '_TI_' + str(self.TI) + '_Shear_' + str(self.shear) + '_Power_T1', "rb") as dill_file:
             Power_T1 = pickle.load(dill_file)
    
#        with open(str(path) + 'Inflow00_T2_L10_2D_cubic', "rb") as dill_file: 
        with open(str(path) + 'Vhub_' + str(Vhub) + '_TI_' + str(self.TI) + '_Shear_' + str(self.shear) + '_L10_T2', "rb") as dill_file: 
             L10_T2_surface = pickle.load(dill_file)

#        with open(str(path) + 'Inflow00_T2_Power_2D_cubic', "rb") as dill_file:
        with open(str(path) + 'Vhub_' + str(Vhub) + '_TI_' + str(self.TI) + '_Shear_' + str(self.shear) + '_Power_T2', "rb") as dill_file:
             Power_T2_surface = pickle.load(dill_file)

        
        Turbine_L10s = [] #np.zeros((len(wd), N))
        Turbine_PowerVals = [] #np.zeros((len(wd), N))
    
        for aa, wind_direction in enumerate(self.wd):    

            # Define the wind array and rotate towards the wind direction, this rotates the array cw, which is equivalent to rotating the array ccw
            theta = -(math.radians(wind_direction)) #convert degrees to radians 
            layout = [self.rotate_array(ii, theta) for ii in list(zip(layout_x, layout_y))]
            
            # Assign L10 values to all turbines
            L10_matrix = np.full( (len(layout), len(layout)), 1E15) # float(L10_T1) 
            Power_matrix = np.full( (len(layout), len(layout)), float(Power_T1) )

            # Sort the turbines from upwind to downwind (with y = 0 in front)
            layout.sort(key = lambda x: x[1])
            
            # If two turbines are close enough together, assign the downwind turbine a reliability and power value
            for turb in list(combinations(enumerate(layout), 2)):
#                print('turb: ', turb)
                if abs(turb[1][1][0] - turb[0][1][0]) <= 378 and abs(turb[1][1][1] - turb[0][1][1]) <= 1890:
                    print('x: ', abs(turb[1][1][0] - turb[0][1][0]))
                    print('y: ', abs(turb[1][1][1] - turb[0][1][1]))
                    L10_matrix[turb[0][0], turb[1][0]] = L10_T2_surface(turb[1][1][0] - turb[0][1][0], turb[1][1][1] - turb[0][1][1])
#                    print(L10_matrix[turb[0][0], turb[1][0]])
#                    print(L10_matrix)
                    Power_matrix[turb[0][0], turb[1][0]] = Power_T2_surface(turb[1][1][0] - turb[0][1][0], turb[1][1][1] - turb[0][1][1]) #kW/hr
#                    print(Power_matrix[turb[0][0], turb[1][0]])
#                    print(Power_matrix)

            # Find the minimum L10 of each row (this identifies which wake limits the life of the bearing)
            print('L10_matrix.T: ', L10_matrix.T)
#            print('Power_matrix.T: ', Power_matrix.T)

#            Turbine_L10s = [min(row) for row in L10_matrix.T] #this is the minimum L10 for each turbine
#            print('Turbine L10s: ', Turbine_L10s)
#            Turbine_L10s.append( [max(row) for row in L10_matrix.T] ) #this is the minimum L10 for each turbine
            for i, j in zip(L10_matrix.T, Power_matrix.T):
                if all(x ==1E15 for x in i):
                    Turbine_L10s.append(float(L10_T1))
                    Turbine_PowerVals.append(float(Power_T1))
                else:
                    Turbine_L10s.append(min(i))
                    minpos = np.argmin(i)
                    Turbine_PowerVals.append(j[minpos])
#            print(Turbine_L10s)
            print(Turbine_PowerVals)
                    
#            Turbine_L10s = [max(row) for row in L10_matrix.T]  #this is the minimum L10 for each turbine

#            Turbine_L10s = [x if x < 1E15 else float(L10_T1) for x in Turbine_L10s]
#            print('Turbine L10s: ', Turbine_L10s)


#            print('Turbine L10s: ', Turbine_L10s)

            
#            result = [np.where(row == np.amax(row)) for row in L10_matrix.T] # a list of each min L10 value for each column

#            print('L10_matrix.T: ', L10_matrix.T)
#            print('result1: ', result)
#            result = [ii[0][0] for ii in result]    
#            print('result1: ', result)
#            Turbine_PowerVals = [column[ii] for column, ii in zip(Power_matrix.T, result)]   #this is the corresponding power value for each turbine                    

#            Turbine_PowerVals.append( [column[ii] for column, ii in zip(Power_matrix.T, result)] )   #this is the corresponding power value for each turbine                    
#            print('Turbine Power: ', Turbine_PowerVals)
#        # Use linear damage accumulation to determine L10 for each turbine for the wind directions (in hours)
#
#        L10_OverAllWindDirs = [1/sum([b/a for a,b in zip(column, self.wd_freq)]) for column in Turbine_L10s] #1E5 #self.freq        
#        print('L10_OverAllWindDirs: ', L10_OverAllWindDirs)
#        # Use weighted sums to determine power for each turbine for the wind directions (in hours)
#        Power_OverAllWindDirs = [sum([a*b for a,b in zip(column, self.wd_freq)]) for column in Turbine_PowerVals]   #kW/hr #self.freq
#        print('Power_OverAllWindDirs: ', Power_OverAllWindDirs)
#
#        return L10_OverAllWindDirs, Power_OverAllWindDirs   # in hours, will be the same length of N
        return Turbine_L10s, Turbine_PowerVals
    
    
    def calc_cost_power(self, x0):
        '''Failure costs over AEP'''
                
        layout_x = [self._unnorm(valx, self.bndx_min, self.bndx_max) \
            for valx in x0[0:self.N]]
        
        layout_y = [self._unnorm(valy, self.bndy_min, self.bndy_max) \
            for valy in x0[self.N:2*self.N]]
        
        if len(self.ws) > 1:
            L10_ws = []      
            Power_ws = []
                                    
            for Vhub in self.ws:
                L10_temp, Power_temp = self.calc_L10_Power(layout_y, layout_x, Vhub)  # L10 in hours, Power in kW/hr, len(N)
                L10_ws.append(L10_temp)
                Power_ws.append(Power_temp)    
    
            L10 = []
            Power = []
    
            for aa in range(self.N):
    
                Turbine_L10 = [ii[aa] for ii in L10_ws] #get list of L10s per turbine over ws
                L10.append(1/sum([b/a for a,b in zip(Turbine_L10, self.ws_freq)]))
                
                Turbine_Power = [ii[aa] for ii in Power_ws]
                Power.append(sum([a*b for a,b in zip(Turbine_Power, self.ws_freq)]))

        else:
#            print('YAS')
            L10, Power = self.calc_L10_Power(layout_x, layout_y, self.ws[0])    # L10 in hours, Power in kW/hr #outputs L10 values for each turbine (len(N))

        hours_lifetime = 20*365*24                                              # hours in 20 year lifespan
        downtime = 231                                                          # hours per replacement
        replacement_cost = 721000                                               # cost of GB replacement (ref: Yeng, Sheng, and Court: $628,000 USD2011, $715920 USD2019, 720623.01 USD2020)
                                                      
#        replacements = [math.floor(hours_lifetime/ii) for ii in L10]            # for all turbines in farm took away the "int" in front of the parenthesis bc potlly causing Exit mode 5
        replacements = [hours_lifetime/ii for ii in L10]            # for all turbines in farm took away the "int" in front of the parenthesis bc potlly causing Exit mode 5
#        print([hours_lifetime/ii for ii in L10])
        print('L10: ', L10)
#        print('Power: ', Power)
#        print('replacements: ', replacements)
        power_hours = [hours_lifetime-ii*downtime for ii in replacements]       # hours turbine is producing power over its life
#        print('power_hours: ', power_hours)
        
        cost = replacement_cost*(sum(replacements))                             # for all turbines, failure cost   
#        print('Cost: ', cost)
        power = sum([a*b for a, b in zip(power_hours, Power)])                  # kW/hr*hr for all turbines
#        print('PPPPPower: ', power)

#        L10_sum = sum(L10) #get rid: experiment

        return cost, power, replacements, L10                                        # of array


    def objective(self, layout):
        '''Failure costs over AEP'''
        
        cost, power, replacements, L10 = self.calc_cost_power(layout)                # L10 in hours, Power in kW/hr
#        print('Obj Cost: ', cost)
#        print('Obj Power: ', power)
#        print('Obj Rep: ', replacements)
#        print('Obj L10: ', L10)
#        print('Obj sum(L10): ', )
        print('obj: ', (1/(sum(L10)))*1E4)
        return (1/(sum(L10)))*1E4
#        return (sum(replacements))* 0.01
#        return (cost/power)*1E7 
#        return cost/1E7 
#        return (1/power)*1E15 


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
                                options=self.opt_options)

        opt_results = self.residual_plant.x #need the .x for later parsing
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
        locsx_old_D = [valx/self.D for valx in locsx_old]
        locsy_old = [self._unnorm(valy, self.bndy_min, self.bndy_max) \
                     for valy in self.x0[self.N:2*self.N]]
        locsy_old_D = [valy/self.D for valy in locsy_old]
        locsx = [self._unnorm(valx, self.bndx_min, self.bndx_max) \
                 for valx in self.residual_plant.x[0:self.N]]
        locsx_D = [valx/self.D for valx in locsx]

        locsy = [self._unnorm(valy, self.bndy_min, self.bndy_max) \
                 for valy in self.residual_plant.x[self.N:2*self.N]]
        locsy_D = [valx/self.D for valx in locsy]

        cost, power, replacements, L10 = self.calc_cost_power(self.residual_plant.x)   #L10 in hours, Power in kW/hr
        objective = self.objective(self.residual_plant.x) #((0.3*cost)/(0.9*power))*1E08


##        locsx_old = [self._unnorm(valx, self.bndx_min, self.bndx_max) \
##                     for valx in self.x0[0:self.N]]
##        locsx_old_D = [valx/self.D for valx in locsx_old]
##        locsy_old = [self._unnorm(valy, self.bndy_min, self.bndy_max) \
##                     for valy in self.x0[self.N:2*self.N]]
##        locsy_old_D = [valy/self.D for valy in locsy_old]
##        locsx = [self._unnorm(valx, self.bndx_min, self.bndx_max) \
##                 for valx in layout_x]
#        locsx_D = [valx*15 for valx in layout_x]
#
##        locsy = [self._unnorm(valy, self.bndy_min, self.bndy_max) \
##                 for valy in layout_y]
#        locsy_D = [valy*15 for valy in layout_y]
#
#        layout = layout_x + layout_y

#        cost, power, replacements = self.calc_cost_power(layout)   #L10 in hours, Power in kW/hr
#        objective = self.objective(layout) #((0.3*cost)/(0.9*power))*1E08


#        print('Replacements: ', replacements, 'Cost: ', cost, 'Power: ', power, 'Objective: ', objective)

#        print(locsx_D, locsy_D)
        plt.figure(figsize=(8,7))
#        plt.figure(figsize=(3,6))
        fontsize= 16
        plt.plot(locsx_old_D, locsy_old_D, '3', color = 'darkcyan', markersize=16)
        plt.plot(locsx_D, locsy_D, '2', color = 'red', markersize=16)
        plt.xlabel('x (rotor diameters)', fontsize=fontsize)
        plt.ylabel('y (rotor diameters)', fontsize=fontsize)
        plt.grid()
        plt.axis([-0.25, 1.25, -2.0, 16])

        # Add the rectangle
#        rect = Rectangle( (0,0), 1, 15, linestyle = 'dashed', facecolor = 'None', clip_on=False)
#        plt.axes.Axes.add_patch(rect)
#        plt.axis([-1.5, 16, -1.5, 16])
        plt.tick_params(which='both', labelsize=fontsize)
        plt.legend(['Original locations', 'Optimized locations'], loc='lower center', \
#        plt.legend(['Optimized locations'], loc='lower center', \
            bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=fontsize)
#        plt.text(3.5, -0.75, 'Objective = ' + str(objective))
#        plt.text(2.5, -1.25, 'Cost = ' + str(cost) + '   Power = ' + str(power))
        plt.text(0.3, -1.00, 'Objective = Cost/Power') #+ str(cost/power))
        plt.text(0.1, -1.50, 'Cost = $' + str(round(cost)) + '   Power = ' + '{:.5e}'.format(power) + ' kW')

        plt.show()

        plt.close()
        return objective, cost, power, replacements 
