# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 04:40:26 2019

@author: USER
"""

from Designs.general import numerical_method


class AdamsMoulton(numerical_method.NumericalMethod):
  
    def __init__(self, D, L, v, e, g=9.81):
        numerical_method.NumericalMethod.__init__(self, D, L, v, e, g)
        self.precomputed_fns_bashforth = {}
                    
    '''# Adams Bashforth 4th order as the Predictor to Adams Moulton
    # tn = tn;
    # tn_1 = t_(n-1) (i.e, t subscript (n - 1));
    # tn_2 = t_(n-2) (i.e, t subscript (n - 2));
    # tn_3 = t_(n-3) (i.e, t subscript (n - 3));
    # say n = 3; we have t3, t2, t1, t0 in the order shown above.
    # same for z and V.'''
    def bash_forth_4th(self, tn, zn, Vn, tn_1, zn_1, Vn_1, 
                       tn_2, zn_2, Vn_2, tn_3, zn_3, Vn_3):
        # zn+1 = zn + (h/24)*(55*f(tn, zn, Vn) - 59*f(tn-1, zn-1, Vn-1)
        #        + 37*f(tn-2, zn-2, Vn-2) - 9*f(tn-3, zn-3, Vn-3))
        if (not self.A1) or (not self.A2):
            ff = self.getF(Vn)
            f = ff["f"] # get friction factor
            R = ff["R"] # get Reynold's number\
            self.results[tn] = (tn, zn, Vn, f, R) # track all the required values
        
        fn, fn_1, fn_2, fn_3 = 0, 0, 0, 0
        gn, gn_1, gn_2, gn_3 = 0, 0, 0, 0
        
        # Dynamic Programming
        # Checks if a function with the same t, z, V is already evaluated
        # if it is look it up, else compute the function
        if (tn, zn, Vn) in self.precomputed_fns_bashforth:
            fn = self.precomputed_fns_bashforth[(tn, zn, Vn)][0]
            gn = self.precomputed_fns_bashforth[(tn, zn, Vn)][1]
        else:
            fn = self.fn(tn, zn, Vn)
            gn = self.gn1(tn, zn, Vn)
            self.precomputed_fns_bashforth[(tn, zn, Vn)] = fn, gn
            
        if (tn_1, zn_1, Vn_1) in self.precomputed_fns_bashforth:
            fn_1 = self.precomputed_fns_bashforth[(tn_1, zn_1, Vn_1)][0]
            gn_1 = self.precomputed_fns_bashforth[(tn_1, zn_1, Vn_1)][1]
        else:
            fn_1 = self.fn(tn_1, zn_1, Vn_1)
            gn_1 = self.gn1(tn_1, zn_1, Vn_1)
            self.precomputed_fns_bashforth[(tn_1, zn_1, Vn_1)] = fn_1, gn_1
            
        if (tn_2, zn_2, Vn_2) in self.precomputed_fns_bashforth:
            fn_2 = self.precomputed_fns_bashforth[(tn_2, zn_2, Vn_2)][0]
            gn_2 = self.precomputed_fns_bashforth[(tn_2, zn_2, Vn_2)][1]
        else:
            fn_2 = self.fn(tn_2, zn_2, Vn_2)
            gn_2 = self.gn1(tn_2, zn_2, Vn_2)
            self.precomputed_fns_bashforth[(tn_2, zn_2, Vn_2)] = fn_2, gn_2
            
        if (tn_3, zn_3, Vn_3) in self.precomputed_fns_bashforth:
            fn_3 = self.precomputed_fns_bashforth[(tn_3, zn_3, Vn_3)][0]
            gn_3 = self.precomputed_fns_bashforth[(tn_3, zn_3, Vn_3)][1]
        else:
            fn_3 = self.fn(tn_3, zn_3, Vn_3)
            gn_3 = self.gn1(tn_3, zn_3, Vn_3)
            self.precomputed_fns_bashforth[(tn_3, zn_3, Vn_3)] = fn_3, gn_3
        
        zn_new =  zn + ((1/24.0) * self.h) * \
                    (
                      55*fn
                      - 59*fn_1
                      + 37*fn_2
                      - 9*fn_3
                    )
        Vn_new =  Vn + ((1/24.0) * self.h) * \
                    (
                      55*gn
                      - 59*gn_1
                      + 37*gn_2
                      - 9*gn_3
                    )
                    
        return (zn_new, Vn_new)
    
    def moulton_4th(self, tn, zn, Vn, tn_1, zn_1, Vn_1, 
                       tn_2, zn_2, Vn_2, tn_3, zn_3, Vn_3):
        if (not self.A1) or (not self.A2):
            ff = self.getF(Vn)
            f = ff["f"] # get friction factor
            R = ff["R"] # get Reynold's number\
            self.results[tn] = (tn, zn, Vn, f, R) # track all the required values
        # Use Adams Bashforth as the Predictor
        zn_bar, Vn_bar = self.bash_forth_4th(tn, zn, Vn, tn_1, zn_1, Vn_1, 
                       tn_2, zn_2, Vn_2, tn_3, zn_3, Vn_3)
        
        tn_bar = tn + self.h
        
        # Use Adams Moulton as the Corrector
        fn, fn_1, fn_2, fn_bar = 0, 0, 0, 0
        gn, gn_1, gn_2, gn_bar = 0, 0, 0, 0
        
        # Dynamic Programming
        # Checks if a function with the same t, z, V is already evaluated
        # if it is look it up, else compute the function
        
        # 6 values (fn, fn_1, fn_2, gn, gn_1, gn_2) would have been computed 
        # from the call of Bashforth as the predictor.
        fn = self.precomputed_fns_bashforth[(tn, zn, Vn)][0]
        gn = self.precomputed_fns_bashforth[(tn, zn, Vn)][1]
        
        fn_1 = self.precomputed_fns_bashforth[(tn_1, zn_1, Vn_1)][0]
        gn_1 = self.precomputed_fns_bashforth[(tn_1, zn_1, Vn_1)][1]
        
        fn_2 = self.precomputed_fns_bashforth[(tn_2, zn_2, Vn_2)][0]
        gn_2 = self.precomputed_fns_bashforth[(tn_2, zn_2, Vn_2)][1]
            
        if (tn_bar, zn_bar, Vn_bar) in self.precomputed_fns_bashforth:
            fn_bar = self.precomputed_fns_bashforth[(tn_bar, zn_bar, Vn_bar)][0]
            gn_bar = self.precomputed_fns_bashforth[(tn_bar, zn_bar, Vn_bar)][1]
        else:
            fn_bar = self.fn(tn_bar, zn_bar, Vn_bar)
            gn_bar = self.gn1(tn_bar, zn_bar, Vn_bar)
            self.precomputed_fns_bashforth[(tn_bar, zn_bar, Vn_bar)] = fn_bar, gn_bar
     
        zn_new =  zn + ((1/24.0) * self.h) * \
                    (
                      9*fn_bar
                      + 19*fn 
                      - 5*fn_1
                      + fn_2
                    )
        Vn_new =  Vn + ((1/24.0) * self.h) * \
                    (
                      9*gn_bar
                      + 19*gn
                      - 5*gn_1 
                      + gn_2
                    )
        return (zn_new, Vn_new)
    
    def moulton_solution(self):
        # Inital values
        t0 = self.t0
        z0 = self.z0
        V0 = self.V0
        h = self.h

        # here n = 0
        tn = t0
        zn = z0
        Vn = V0
          
        # Get all starting value from runge-kutta:
        # 1. Get the first 3 unknown values from runge-kutta:
        # n-3, n-2, n-1
        rn = 3
        
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
        
        rk_starting_values = []
        
        # Iterate for number of iterations
        i = 1
        while i <= rn:
            # Update the next value of z and V
            zn, Vn = self.runge_kutta_4th(tn, zn, Vn)

            # Update the next value of t
            tn += h
            rk_starting_values.append((tn, zn, Vn))

            z_s[i] = (tn, zn)
            if self.A1 and self.A2: # if we have two reservoirs
                z1 = (zn * self.A) / self.A1
                z_s_1[i] = (tn, z1)
                
                z2 = (zn * self.A) / self.A2
                z_s_2[i] = (tn, z2)
           
            i += 1

        # note that "xn_a" implies "x subscript n-a"
        # => that "tn_2" implies "t subscript n-2"
        tn_3, zn_3, Vn_3 = t0, z0, V0
        tn_2, zn_2, Vn_2 = rk_starting_values[0] # t1, z1, V1
        tn_1, zn_1, Vn_1 = rk_starting_values[1] # t2, z2, V2
        tn, zn, Vn = rk_starting_values[2] # t3, z3, V3
        
        # note that the moulton continuing values does not contain t0, z0 and V0
        # hence it has just the '3' values for the Runge-Kutta starting values.
        am_continuing_values = []
        am_continuing_values.extend(rk_starting_values)
        
        # Count number of iterations using step size or 
        # step height h.
        n = int((self.t - self.t0) / self.h)
        # note that 'i' ends at '4' from the frist 'while loop iteration'
        # and we want to begin 'i' also at '4'.
        while i <= n:
            # Update the next value of z and V
            zn, Vn = self.moulton_4th(tn, zn, Vn, tn_1, zn_1, Vn_1,
                                         tn_2, zn_2, Vn_2, tn_3, zn_3, Vn_3)
            
            # Update the next value of t
            tn += h

            am_continuing_values.append((tn, zn, Vn))

            # the first continuing values only has the Runge-kutta '3' starting
            # values.
            tn_3, zn_3, Vn_3 = am_continuing_values[0]
            tn_2, zn_2, Vn_2 = am_continuing_values[1]
            tn_1, zn_1, Vn_1 = am_continuing_values[2]
            tn, zn, Vn = am_continuing_values[3]
            
            # now, am_continuing_values has four (4) values but
            # the last three (3) values would be picked for the next iteration.
            am_continuing_values = am_continuing_values[1:]
            
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
    
    def display(self, _type):
        if self.A1 and self.A2: # if we have two reservoirs
            method_vals = self.moulton_solution()
            # Draw table for all the values
            self.drawTable(method_vals, self.pipe_label, 
                           self.reservoir_1_label, self.reservoir_2_label)
            # Plot all values
            self.plot_single({self.pipe_label: method_vals[self.pipe_label]}, 
                              self.pipe_label, method_label=self.moulton_label)
#            self.plot_single({self.reservoir_1_label: method_vals[self.reservoir_1_label]}, 
#                              self.reservoir_1_label, method_label=self.moulton_label)
#            self.plot_single({self.reservoir_2_label: method_vals[self.reservoir_2_label]},
#                              self.reservoir_2_label, method_label=self.moulton_label)
            self.plot_multiple(method_vals, _type + " (all-values)", 
                               method_label=self.moulton_label)
#            
#            turning_points = self.get_all_turning_points_for_pipe_reservoir1_reservoir2(method_vals)
            # Draw table for all the turning points
#            self.drawTable(turning_points, self.pipe_label, 
#                           self.reservoir_1_label, self.reservoir_2_label)
             # Plot values for all the turning points
#            self.plot_multiple(turning_points, _type+" (subsequent: minimum-maximum)",
#                               method_label=self.moulton_label)
        else:
            method_vals = self.moulton_solution()
            # Draw table for all the values
            self.drawTable(method_vals, self.pipe_label)

            # Plot all values
#            self.plot_single(method_vals, self.pipe_label, 
#                             method_label=self.moulton_label)
            
#            turning_points = self.get_all_turning_points_for_pipe(method_vals)
#            # Draw table for all the turning points
#            self.drawTable(turning_points, self.pipe_label)
#            # Plot values for all the turning points
#            self.plot_single(turning_points, _type+" (subsequent: minimum-maximum)",
#                             method_label=self.moulton_label)
        

if __name__ == "__main__":
    print("Adams Moulton...")
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
#    t = 2 # final time
#    h = (1/16.0) # time step
#    g = 32.2

#    am = AdamsMoulton(D, L, v, e, g)
#    am.set_initial_conditions(t0, z0, dzdt, h, t)
    
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
#    am = AdamsMoulton(D, L, v, e, g)
#    am.set_initial_conditions(t0, z0, dzdt, h, t)
#
#    D1 = 15.958 # diameter of first reservoir; ft
#    D2 = 19.544 # diameter of second reservoir; ft
#    f = 0.024 # constant friction factor
#    length_minor_losses  = 438 # minor losses; ft
#    am.set_reservoir_param(D1, D2, length_minor_losses, f)
        
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
 
    am = AdamsMoulton(D, L, v, e, g)
    am.set_initial_conditions(t0, z0, dzdt, h, t)

    D1 = 10 # diameter of first reservoir; m
    D2 = 8 # diameter of second reservoir; m
#    f = 0.012 # constant friction factor
    length_minor_losses  = 160 # minor losses == 8% of pipe length; m
#    am.set_reservoir_param(D1, D2, length_minor_losses, f)

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
#    t = 100 # final time; secs <br>
#    h = 1/16.0 # time step; secs <br>
#    g = 32.2 # acceleration due to gravity; ft/secs^2
#
#    am = AdamsMoulton(D, L, v, e, g)
#    am.set_initial_conditions(t0, z0, dzdt, h, t)
######
    
    am.display("Moulton")
    t2 = time()
    print()
    print((t2 - t1))