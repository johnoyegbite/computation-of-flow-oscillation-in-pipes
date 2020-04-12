# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 01:00:48 2019

@author: USER
"""
"""
LAMINAR RESISTANCE:
    Displacement z of one meniscus of the column as a function of time,
    starting with the meniscus at z = 0 when t = 0 and rising with velocity V0.
    
    Solution to d^2(z)/dt^2 + (32v/D^2)dz/dt + (2g/L)z = 0 is given by:
        
    z = (V0/n) * e^(-mt) * sinh(nt)
    
    where m = 16v/D^2
          b = 2g/L
          n = sqrt(m^2 - b)
          where v is the kinematic viscosity
                D is the diameter of pipe
                L is the lenght of the pipe
                g is the acceleration due to gravity, 9.81m/s^2
                V0 is the velocity at t=0, z=0 
    Case 1: when m^2 < b
        z = (V0/n') * e^(-mt) * sin(n't)
            where n' = sqrt(b - m^2)
    Case 2: when m^2 = b
        z = V0*t*e^(-mt)
    Case 3: when m^2 > b
        z = (V0/sqrt(m^2 - n^2)) * ((m-n)/(m +n))^(m/2n)
    
""" 

from math import e as E, sin, sqrt

from Designs.general import numerical_method

# This analytical method is for laminar flow only
class Z_Analytical(numerical_method.NumericalMethod):
    def __init__(self, D, L, v, e, g=9.81):
        numerical_method.NumericalMethod.__init__(self, D, L, v, e, g)
        self.order = 1.0
        
    def get_m(self):
        # See definition above
        if self.A1 and self.A2:
            m = ((16 * self.v) / (self.D ** 2)) * (self.L/self.Le)
            return m
        return (16 * self.v) / (self.D ** 2)
    
    def get_z_analytical(self, t):
        # Take care of all the Cases that can be seen in the  definition above
        m = self.get_m()
        a = pow(m, 2)
        if self.A1 and self.A2:
            b = self.b
        else:
            b = (2 * self.g)/self.L
        discriminant = a - b
        z = None
        if discriminant < 0: # discriminant is negative
            n = sqrt(-discriminant) # make it positive by negating it
            z = (self.V0/n) * pow(E, -m * t) * sin(n*t)
        elif discriminant == 0:
            z = self.V0 * t * pow(E, -m * t)
        else:
            n = sqrt(discriminant)
            first_diff = m - n
            second_diff = m + n
            z = (self.V0 / (sqrt(first_diff*second_diff))) \
                * (first_diff/second_diff) ** (m/(2*n))
        return z

    def exact_solution(self):
        # Count number of iterations using step size or 
        # step height h.
        n = int((self.t - self.t0) / self.h)
        
        # Inital values
        t0 = self.t0
        z0 = self.z0
        h = self.h
        
        # Here n = 0
        tn = t0
        zn = z0

        # Store the initial in terms of n: (tn, zn).
        z_s = {}
        z_s_1 = {}
        z_s_2 = {}

        z_s = {0: (t0, z0)} # here n = 0
        if self.A1 and self.A2: # if we have two reservoirs
            z1 = (z0 * self.A) / self.A1
            z_s_1 = {0: (t0, z1)}
            z2 = (z0 * self.A) / self.A2
            z_s_2 = {0: (t0, z2)}
        
        i = 1
        # Store all the value of z 
        # from t (corresponding to i = 1), 
        # to t (corresponding to i = n).
        # Also store the values of the reservoirs if defined
        while i <= n:
            # Update the next value of t
            tn += h
            zn = self.get_z_analytical(tn)
            
            z_s[i] = (tn, zn)
            if self.A1 and self.A2: # if we have two reservoirs
                z1 = (zn * self.A) / self.A1
                z_s_1[i] = (tn, z1)
                
                z2 = (zn * self.A) / self.A2
                z_s_2[i] = (tn, z2)
            i += 1
            
        if self.A1 and self.A2: # if we have two reservoirs
            return {self.pipe_label: z_s, 
                    self.reservoir_1_label: z_s_1, 
                    self.reservoir_2_label: z_s_2
                   }
        return {self.pipe_label: z_s}

    def plot_single_graph(self, fig="", _title="", t_label="", z_label="", 
                          t_vals=[], z_vals=[], _label="", 
                          color_plus_lineType=['r--', 'b--'], line_thickness=1.0):
        '''t_vals: a tuple containinig two list object
            z_vals: a tuple containinig two list object
        '''
        import matplotlib.pyplot as plt

        plt.figure(fig)
        plt.clf()
        plt.xlabel(t_label)
        plt.ylabel(z_label)
        plt.plot(t_vals, z_vals, color_plus_lineType[0], 
                 label = _label, linewidth = line_thickness)
        plt.legend(loc = 'upper right')
        plt.title(_title)
        plt.show()

    def plot_single(self, method_vals, _type):
        tvals = []
        surge_vals = []
        time, surge = 0, 1
        surge_dict = method_vals[_type]
        for key in ((surge_dict)):
            tvals.append(surge_dict[key][time])
            surge_vals.append(surge_dict[key][surge])
            
        self.plot_single_graph(fig=f"Exact method ({_type})", 
                               _title=f"Graph of {_type} Surge Vs. Time",
                               t_label="time, t",
                               z_label="surge, z", 
                               t_vals=tvals, 
                               z_vals=surge_vals,
                              _label=f"{_type} (viscosity: {self.v})", 
                              color_plus_lineType=['b--', 'r--'])

    def plot_multiple_graphs(self, fig="", _title="", t_label="", 
                             z_label="", t_vals=[], z_vals=[], labels=(), 
                          color_plus_lineType=['b--', 'g^', 'ro'], 
                          line_thickness=.5):
        '''t_vals: a tuple containinig two list object
            z_vals: a tuple containinig two list object
        '''
        import matplotlib.pyplot as plt
        plt.figure(fig)
        plt.clf()
        plt.xlabel(t_label)
        plt.ylabel(z_label)
        y0, y1, y2 = 0, 1, 2
        plt.plot(t_vals, z_vals[y0], color_plus_lineType[y0], 
                 label = labels[y0], linewidth = line_thickness)
        plt.plot(t_vals, z_vals[y1], color_plus_lineType[y1], 
                 label = labels[y1], linewidth = line_thickness)
        plt.plot(t_vals, z_vals[y2], color_plus_lineType[y2], 
                 label = labels[y2], linewidth = line_thickness)
        plt.legend(loc = 'upper right')
        plt.title(_title)
        plt.show()
        
    def plot_multiple(self, method_vals, _type):
        tvals = []
        pipe_vals, reservoir_1_vals, reservoir_2_vals = [], [], []
        time, surge = 0, 1
        pipe_dict = method_vals[self.pipe_label]
        reservoir_1_dict = method_vals[self.reservoir_1_label]
        reservoir_2_dict = method_vals[self.reservoir_2_label]
        for key in ((pipe_dict)):
            tvals.append(pipe_dict[key][time])
            pipe_vals.append(pipe_dict[key][surge])
            reservoir_1_vals.append(reservoir_1_dict[key][surge])
            reservoir_2_vals.append(reservoir_2_dict[key][surge])
            
        self.plot_multiple_graphs(fig=f"Exact method", 
                               _title=f"Graph of Pipe-Reservoirs Surge Vs. Time",
                               t_label="time, t", 
                               z_label="surge, z", 
                               t_vals=tvals, 
                               z_vals=(pipe_vals, 
                                       reservoir_1_vals, 
                                       reservoir_2_vals),
                              labels=(self.pipe_label, 
                                      self.reservoir_1_label,
                                      self.reservoir_2_label), 
                              color_plus_lineType=['b--', 'g^', 'ro'])
                             
    def display(self, _type):
        if self.A1 and self.A2: # If we have two reservoirs
            
            method_vals = self.exact_solution()
            # Draw table for all the values
            self.drawTable(method_vals, self.pipe_label, 
                           self.reservoir_1_label, self.reservoir_2_label)
            # Plot all values
            self.plot_single({self.pipe_label: method_vals[self.pipe_label]}, 
                              self.pipe_label)
#            self.plot_single({self.reservoir_1_label: method_vals[self.reservoir_1_label]}, 
#                              self.reservoir_1_label)
#            self.plot_single({self.reservoir_2_label: method_vals[self.reservoir_2_label]},
#                              self.reservoir_2_label)
#            
#            self.plot_multiple(method_vals, _type + " (all-values)")
            
            # Turning point. Since the shape of the graph would be sinusoidal
#            turning_points = self.get_all_turning_points_for_pipe_reservoir1_reservoir2(method_vals)
            # Draw table for all the turning points
#            self.drawTable(turning_points, self.pipe_label, 
#                           self.reservoir_1_label, self.reservoir_2_label)
            # Plot values for all the turning points
#            self.plot_multiple(turning_points, _type+" (subsequent: minimum-maximum)")

        else:
            method_vals = self.exact_solution()
            # Draw table for all the values
            self.drawTable(method_vals, self.pipe_label)
            # Plot all values
#            self.plot_single(method_vals, self.pipe_label)
            
#            turning_points = self.get_all_turning_points_for_pipe(method_vals)
             # Draw table for all the turning points.
#            self.drawTable(turning_points, self.pipe_label)
             # Plot values for all the turning points.
#            self.plot_single(turning_points, _type+" (subsequent: minimum-maximum)")
    
    
    # To get just one value of z, given a t value.
    def get_single(self, t):
        z = self.get_z_analytical(t)
        return (t, z)


if __name__ == "__main__":
    import time
    t1 = time.time()

    # pipe properties
#    D = (1/4)/12.0 # diameter of pipe; m
#    L = (30/12.0) # length of pipe; m
#    v = 1e-4 # kinematic viscosity; m2/s
#    e = 0.007093 # Roughness height; m
#
#    # initial conditions
#    t0 = 0 # time
#    z0 = 0 # head 
#    dzdt = 6 # velocity m/s
##    t = 1.25 # final time
#    t = 4 # final time
#    h = (1/16.0) # time step
#
#    g = 32.2
    
###### Streeter reservoir example (Ex. 10.15)
        # pipe properties
#    D = 3 # diameter of pipe; ft
#    L = 2000 # length of pipe; ft
#    v = 1e-4 # kinematic viscosity; f^2/s
#    e = 0.007093 # Roughness height; ft
#        # initial conditions
#    t0 = 0 # time; secs
#    z0 = 1132 # head; ft
#    dzdt = 0 # velocity; ft/s
#    t = 3600 # final time; secs
#    h = 0.5 # time step; secs
#    g = 32.2 # acceleration due to gravity; ft/secs^2
#     
#    ex = Z_Analytical(D, L, v, e, g)
#    ex.set_initial_conditions(t0, z0, dzdt, h, t)
    
#    D1 = 15.958 # diameter of first reservoir; ft
#    D2 = 19.544 # diameter of second reservoir; ft
#    f = 0.024 # constant friction factor
#    length_minor_losses  = 438 # minor losses; ft
#    ex.set_reservoir_param(D1, D2, length_minor_losses, f)
    
##### CEG 848, UNILAG, M.Sc(Civil Engr) Exam, Second Semester 2017/2018 Session, 
    # Reservoir Question 1b.
        # pipe properties
#    D = 0.8 # diameter of pipe; m
#    L = 2000 # length of pipe; m
#    v = 1e-4 # kinematic viscosity; m^2/s
#    e = 0.007093 # Roughness height; m
#        # initial conditions
#    t0 = 0 # time; secs
#    z0 = 7.8125 # head; m
#    dzdt = 0 # velocity; m/s
#    t = 500 # final time; secs
#    h = 0.5 # time step; secs
#
#    g = 9.81 # acceleration due to gravity; m/secs^2
##    
##    
#    ex = Z_Analytical(D, L, v, e, g)
#    ex.set_initial_conditions(t0, z0, dzdt, h, t)
##
#    D1 = 10 # diameter of first reservoir; m
#    D2 = 8 # diameter of second reservoir; m
#    f = 0.012 # constant friction factor
#    length_minor_losses  = 160 # minor losses == 8% of pipe length; m
#    ex.set_reservoir_param(D1, D2, length_minor_losses, f)
        
##### Formulated Example
        #<b> pipe properties </b> <br>
    D = 3 # diameter of pipe; ft <br>
    L = 2000 # length of pipe; ft <br>
    v = 0.0001 # kinematic viscosity; f^2/s <br>
    e = 0.006 # Roughness height; ft <br>
        #<b> initial conditions </b>  <br>
    t0 = 0 # time; secs <br>
    z0 = 0 # head; ft <br>
    dzdt = 3 # velocity; ft/s <br>
    t = 10 # final time; secs <br>
    h = 0.5 # time step; secs <br>
    g = 32.2 # acceleration due to gravity; ft/secs^2

    ex = Z_Analytical(D, L, v, e, g)
    ex.set_initial_conditions(t0, z0, dzdt, h, t)
      
    print("Exact...")
    print()
    ex.display("Exact")

    
    t2 = time.time()
    print()
    print((t2 - t1))

    
      
  