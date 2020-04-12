# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 04:47:38 2019

@author: USER
"""

from Designs.general import numerical_method

class TrapezoidalRule(numerical_method.NumericalMethod):
    # Assign given values of "Trapezoidal Rule"
    def __init__(self, D, L, v, e, g=9.81):
        numerical_method.NumericalMethod.__init__(self, D, L, v, e, g)
    
    def trapezoidal_2nd_order(self, tn, zn, Vn):
        # Trapezoidal / Modified Euler numerical method.
        zn_bar = zn + self.h * self.fn(tn, zn, Vn)
        Vn_bar = Vn + self.h * self.gn1(tn, zn, Vn)
        tn_bar = tn + self.h
        zn_new = zn + (self.h/2.0) * (self.fn(tn_bar, zn_bar, Vn_bar) + self.fn(tn, zn, Vn))
        Vn_new = Vn + (self.h/2.0) * (self.gn1(tn_bar, zn_bar, Vn_bar) + self.gn1(tn, zn, Vn))

        return (zn_new, Vn_new)
      
    # Finds value of z for a given t and dzdt using step size h 
    # and initial value dzdt, z0 at t0. 
    def trapezoidal_solution(self):
        # Count number of iterations using step size or 
        # step height h.
        n = int((self.t - self.t0) / self.h)

        # Inital values
        t0 = self.t0
        z0 = self.z0
        V0 = self.V0
        h = self.h

        # here n = 0
        tn = t0
        zn = z0
        Vn = V0
        
        z_s = {}
        z_s_1 = {}
        z_s_2 = {}

        z_s = {0: (t0, z0)} # here n = 0
        if self.A1 and self.A2: # if we have two reservoirs
            # z*A = z1*A1 = z2*A2 (from continuity equation)
            z1 = (z0 * self.A) / self.A1
            z_s_1 = {0: (t0, z1)}
            z2 = (z0 * self.A) / self.A2
            z_s_2 = {0: (t0, z2)}
        
        # Iterate for number of iterations
        i = 1
        while i <= n:
            zn, Vn = self.trapezoidal_2nd_order(tn, zn, Vn)
            
            # Update the next value of t
            tn += h
            
            z_s[i] = (tn, zn)
            if self.A1 and self.A2: # if we have two reservoirs
                # z*A = z1*A1 = z2*A2 (from continuity equation)
                z1 = (zn * self.A) / self.A1
                z_s_1[i] = (tn, z1)
                
                z2 = (zn * self.A) / self.A2
                z_s_2[i] = (tn, z2)
#               
            i += 1
        if self.A1 and self.A2: # if we have two reservoirs
            return {self.pipe_label: z_s, 
                    self.reservoir_1_label: z_s_1, 
                    self.reservoir_2_label: z_s_2
                   }

        return {self.pipe_label: z_s}
                              
    def display(self, _type):
        if self.A1 and self.A2: # if we have two reservoirs
            method_vals = self.trapezoidal_solution()
            # Draw table for all the values
            self.drawTable(method_vals, self.pipe_label, 
                           self.reservoir_1_label, self.reservoir_2_label)
            # Plot all values
            self.plot_single({self.pipe_label: method_vals[self.pipe_label]}, 
                              self.pipe_label, method_label=self.trapezoidal_label)
#            self.plot_single({self.reservoir_1_label: method_vals[self.reservoir_1_label]}, 
#                              self.reservoir_1_label, method_label=self.trapezoidal_label)
#            self.plot_single({self.reservoir_2_label: method_vals[self.reservoir_2_label]},
#                              self.reservoir_2_label, method_label=self.trapezoidal_label)
#            self.plot_multiple(method_vals, _type + " (all-values)", 
#                               method_label=self.trapezoidal_label)
#            
#            turning_points = self.get_all_turning_points_for_pipe_reservoir1_reservoir2(method_vals)
            # Draw table for all the turning points
#            self.drawTable(turning_points, self.pipe_label, 
#                           self.reservoir_1_label, self.reservoir_2_label)
             # Plot values for all the turning points
#            self.plot_multiple(turning_points, _type+" (subsequent: minimum-maximum)",
#                               method_label=self.trapezoidal_label)
        else:
            method_vals = self.trapezoidal_solution()
            # Draw table for all the values
            self.drawTable(method_vals, self.pipe_label)

            # Plot all values
#            self.plot_single(method_vals, self.pipe_label, 
#                             method_label=self.trapezoidal_label)
            
#            turning_points = self.get_all_turning_points_for_pipe(method_vals)
#            # Draw table for all the turning points
#            self.drawTable(turning_points, self.pipe_label)
#            # Plot values for all the turning points
#            self.plot_single(turning_points, _type+" (subsequent: minimum-maximum)",
#                             method_label=self.trapezoidal_label)

if __name__ == "__main__":
    print("Trapezoidal Rule...")
    print()
    from time import time
    t1 = time()
######
     # pipe properties
#    D = (1/4)/12.0 # diameter of pipe; m
#    L = (30/12) # length of pipe; m
#    v = 1e-4 # kinematic viscosity; m2/s
#    e = 0.007093 # Roughness height; m
#    
#    # initial conditions
#    t0 = 0 # time
#    z0 = 0 # head 
#    dzdt = 6 # velocity m/s
#    t = 4 # final time
#    h = (1/16.0) # time step
#    g = 32.2
#
#    tr = TrapezoidalRule(D, L, v, e, g)
#    tr.set_initial_conditions(t0, z0, dzdt, h, t)

#### Streeter reservoir example (Ex. 10.15)
#        # pipe properties
#    D = 3 # diameter of pipe; ft
#    L = 2000 # length of pipe; ft
#    v = 1e-4 # kinematic viscosity; f^2/s
#    e = 0.007093 # Roughness height; ft
#        # initial conditions
#    t0 = 0 # time; secs
#    z0 = 1132 # head; ft
#    dzdt = 0 # velocity; ft/s
#    t = 200 # final time; secs
#    h = (0.5) # time step; secs
#    g = 32.2 # acceleration due to gravity; ft/secs^2
#     
#    tr = TrapezoidalRule(D, L, v, e, g)
#    tr.set_initial_conditions(t0, z0, dzdt, h, t)
#
#    D1 = 15.958 # diameter of first reservoir; ft
#    D2 = 19.544 # diameter of second reservoir; ft
#    f = 0.024 # constant friction factor
#    length_minor_losses  = 438 # minor losses; ft
#    tr.set_reservoir_param(D1, D2, length_minor_losses, f)
    
##### CEG 848, UNILAG, M.Sc(Civil Engr) Exam, Second Semester 2017/2018 Session, 
    # Reservoir Question 1b.
        # pipe properties
    D = 0.8 # diameter of pipe; m
    L = 2000 # length of pipe; m
    v = 1e-4 # kinematic viscosity; m^2/s
    e = 0.007093 # Roughness height; m
        # initial conditions
    t0 = 0 # time; secs
    z0 = 7.8125 # head; m
    dzdt = 0 # velocity; m/s
    t = 250 # final time; secs
    h = 5 # time step; secs
    g = 9.81 # acceleration due to gravity; m/secs^2 

    tr = TrapezoidalRule(D, L, v, e, g)
    tr.set_initial_conditions(t0, z0, dzdt, h, t)
#
    D1 = 10 # diameter of first reservoir; m
    D2 = 8 # diameter of second reservoir; m
#    f = 0.012 # constant friction factor
    length_minor_losses  = 160 # minor losses == 8% of pipe length; m
#    tr.set_reservoir_param(D1, D2, length_minor_losses, f)
        
##### Formulated Example
        #<b> pipe properties </b> <br>
#    D = 3 # diameter of pipe; ft <br>
#    L = 2000 # length of pipe; ft <br>
#    v = 0.0001 # kinematic viscosity; f^2/s <br>
#    e = 0.006 # Roughness height; ft <br>
#        #<b> initial conditions </b>  <br>
#    t0 = 0 # time; secs <br>
#    z0 = 0 # head; ft <br>
#    dzdt = 3 # velocity; ft/s <br>
#    t = 50 # final time; secs <br>
#    h = 1/16.0 # time step; secs <br>
#    g = 32.2 # acceleration due to gravity; ft/secs^2
#
#    tr = TrapezoidalRule(D, L, v, e, g)
#    tr.set_initial_conditions(t0, z0, dzdt, h, t)
######
    
    tr.display("Trapezoidal")
  
    t2 = time()
    print()
    print((t2 - t1))
#

