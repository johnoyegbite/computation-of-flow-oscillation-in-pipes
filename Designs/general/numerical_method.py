# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

from Designs.general import friction_factor
from prettytable import PrettyTable
"""This module is built to solve a differential equation of the form:
        d^2(z)/dt^2 + a * |dz/dt| * (dz/dt) + b*z = 0.
        where a = f/2D
              b = 2g/L
              
        if there are two reservoirs
        a = fk/2D where k = Le/L
        b = (gA/L)*((1/A1)+(1/A2)) 
        where A = area of pipe
              A1 = area of first reservoir
              A2 = area of second reservoir
    
    reducing this 2nd order ODE to a system of 1st order ODE:
    we have: dz/dt = f(t, z, V) = V
             d^2(z)/dt^2 = dV/dt = g(t, z, V) = - a*|V|*V - b*z    
"""   


class NumericalMethod(object):
    # Assign given values of all numerical methods
    def __init__(self, D, L, v, e, g=9.81):
        """where v is the kinematic viscosity; unit in m^2/s.
                 D is the diameter of pipe; unit in meters.
                 L is the length of the pipe; unit in meters.
                 e is the roughness height; unit in m.
                 g is the acceleration due to gravity, 9.81m/s^2.
                 
                 A1 is the area of first reservoir
                 A2 is the area of second reservoir
                 D1 is the diameter of first reservoir
                 D2 is the diameter of second reservoir
                 f is the friction factor
        """
        self.D = D
        self.L = L
        self.v = v
        self.e = e
        self.g = g

        self.b = (2*self.g)/self.L
        self.A1 = None # no reservoir 1
        self.A2 = None # no reservoir 2
        self.D1 = None
        self.D2 = None
        self.length_minor_losses = None
        self.Le = L
        self.f = None # 
        
        self.pipe_label = "Pipe"
        self.reservoir_1_label = "Reservoir 1"
        self.reservoir_2_label = "Reservoir 2"
        
        self.exact_label = "Exact"
        self.trapezoidal_label = "Trapezoidal"
        self.runge_kutta_label = "Runge-Kutta"
        self.bashforth_label = "Adams Bashforth"
        self.moulton_label = "Adams Moulton"
        
        self.results = {}
        
    def set_initial_conditions(self, t0, z0, dzdt, h, t):
        """where t0 is the initial curr_time; unit in secs.
                 z0 is the initial head of fluid at t0; unit in m.
                 dzdt is the initial velocity at t0; unit in m/s.
                 h is the curr_time step; unit in secs.
                 t is the final curr_time; unit in secs.
        """
        self.t0 = t0
        self.z0 = z0
        self.t = t
        self.h = h if h!=0 else 1
        self.V0 = dzdt
        
    def set_reservoir_param(self, D1, D2, length_minor_losses, f=None):
        from math import pi
        A1 = (pi * (D1**2)) / 4.0 # Cross sectional area of a circular shape
        self.A1 = A1
        A2 = (pi * (D2**2)) / 4.0 # Cross sectional area of a circular shape
        self.A2 = A2
        A = (pi * (self.D**2)) / 4.0 # Cross sectional area of a circular shape
        self.A = A
        self.length_minor_losses = length_minor_losses
        Le = self.L + self.length_minor_losses # To get equivalent loss in pipe lenght
        self.Le = Le
        self.f = f
        # Formula for coefficient b when we have two reservoirs.
        self.b = ( (self.g * self.A) / self.L ) * ( (1/self.A1) + (1/self.A2) )
        
    def interpolate(self, boundary_1, boundary_2,  param_2):
        """Assumes boundary_1 is a list of two Real numbers
                   boundary_2 is a list of two Real numbers
                   param_2_middle is a Real numbers.
        Returns: a Real number that correspond to param_2_middle, 
                which is the interpolation of boundary_1 and boundary_2
        """
        lower, upper = 0, 1
        lower_1, upper_1 = boundary_1[lower], boundary_1[upper]
        lower_2, upper_2 = boundary_2[lower], boundary_2[upper]
        param_1 = (((param_2 - lower_2)*(upper_1 - lower_1))/(upper_2 - lower_2)) + lower_1
        return param_1

    def getF(self, V):
        # Returns a dictionary that contains friction factor and Reynolds number.
        ff = friction_factor.FrictionFactor(V, self.v, self.e, self.D)
        return ff.getF_and_R()

    # This represent dz/dt = V; where V is the velocity.
    def fn(self, t, z, V): 
        return V
    
    # This control the value of 'a' in the 2nd order ODE 
    def get_k(self):
        return self.Le / self.L

    # Please don't use 'g' variable for the function name because 'g' 
    # has been used for the acceleration due to gravity; g = 9.81.
    # This represent d^2(z)/dt^2 = dv/dt = - a|V|V - bz; where V is the velocity
    def gn1(self, t, z, V): 
        # function to get friction factor and Reynold's number given a velocity.
        # note that new friction factor is computed for each call of this
        # function
        ff = self.getF(V) if self.f == None else {"f": self.f, "R": None}
        f = ff["f"] # get friction factor
        R = ff["R"] # get Reynold's number

        k = self.get_k()
        a = (f*k)/(2*self.D) # k would be 1 if there is no reservoir since Le == L.
                            # which would reduce a to a = f/2D  
        return -(a * V*abs(V)) - (self.b * z)
    
    def k(self, t, z, V):
        return self.fn(t, z, V)
  
    def m(self, t, z, V):
        return self.gn1(t, z, V)

    # Runge-Kutta constant computation for 2nd order
    def runge_kutta_4th(self, tn, zn, Vn):
        if (not self.A1) or (not self.A2):
            ff = self.getF(Vn)
            f = ff["f"] # get friction factor
            R = ff["R"] # get Reynold's number
            self.results[tn] = (tn, zn, Vn, f, R) # track all the required values
        h = self.h
        
        k1 = h * self.k(tn, zn, Vn)
        m1 = h * self.m(tn, zn, Vn)
        
        k2 = h * self.k(tn + h/2.0, zn + k1/2.0, Vn + m1/2.0)
        m2 = h * self.m(tn + h/2.0, zn + k1/2.0, Vn + m1/2.0)
        
        k3 = h * self.k(tn + h/2.0, zn + k2/2.0, Vn + m2/2.0)
        m3 = h * self.m(tn + h/2.0, zn + k2/2.0, Vn + m2/2.0)
        
        k4 = h * self.k(tn + h, zn + k3, Vn + m3)
        m4 = h * self.m(tn + h, zn + k3, Vn + m3)
        
        zn = zn + (1/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        Vn = Vn + (1/6.0) * (m1 + 2*m2 + 2*m3 + m4)

        return (zn, Vn)
    
    def get_turning_point(self, index, elem_dict):
        """Assumes index: an int
                   elem_dict: a dictionary, 
                       where the key is an integer and 
                                 value is a tuple of (time, head).
                           
        Turning point definition:
            A number is a turning point if 
            1. its absolute value is greater than or less than both absolute 
                value of previous number and the next number
            2. and all three element are of the same sign.
        
        """
        is_turning_point = False
        head_index = 1
        time_index = 0
        curr_head = elem_dict[index][head_index]
        curr_time = elem_dict[index][time_index]

        prev_index = index - 1
        next_index = index + 1

        # If element is at the start of the dictionary,
        # just check if the next element is a turning point neighbour. 
        # See doc string for turning point definition.
        if index == 0:
            next_head = elem_dict[next_index][head_index]
            if abs(curr_head) >= abs(next_head):
                if curr_head < 0 and next_head < 0:
                    is_turning_point = True
                if curr_head > 0 and next_head > 0:
                    is_turning_point = True
        # If element is at the end of the dictionary,
        # just check if the previous element is a turning point neighbour.
        # See doc string for turning point definition.
        elif index == (len(elem_dict) - 1):
            prev_head = elem_dict[prev_index][head_index]
            if abs(curr_head) >= abs(prev_head):
                if curr_head < 0 and prev_head < 0:
                    is_turning_point = True
                if curr_head > 0 and prev_head > 0:
                    is_turning_point = True
        # Check if the both the previous and next element are turning point
        # neighbour.
        # See doc string for turning point definition.
        else:
            prev_head = elem_dict[prev_index][head_index]
            next_head = elem_dict[next_index][head_index]
            
            if abs(curr_head) > abs(prev_head) and abs(curr_head) > abs(next_head):
                # all element have to have the same sign to be a turning point
                if curr_head < 0 and prev_head < 0 and next_head < 0:
                    is_turning_point = True
                if curr_head > 0 and prev_head > 0 and next_head > 0:
                    is_turning_point = True
            if abs(curr_head) < abs(prev_head) and abs(curr_head) < abs(next_head):
                if curr_head < 0 and prev_head < 0 and next_head < 0:
                    is_turning_point = True
                if curr_head > 0 and prev_head > 0 and next_head > 0:
                    is_turning_point = True
        if is_turning_point:
            return (curr_time, curr_head)
                    
    def get_all_turning_points_for_pipe_reservoir1_reservoir2(self, elems_dict):
        '''Assumes elems_dict is a dictionary that contains:
            key: a string (Pipe, Reservoir 1 or Reservoir 2)
            values: dictionary:
                        key: an integer
                        value: tuple of (time, head).
        '''
        pipe_turning_points = {}
        reservoir_1_turning_points = {}
        reservoir_2_turning_points = {}
        # we must have at least 3 element to get a turning point
        first_dict = elems_dict[self.pipe_label]
        if len(first_dict) > 2:
            pipe_dict = elems_dict[self.pipe_label] 
            reservoir_1_dict = elems_dict[self.reservoir_1_label]
            reservoir_2_dict  = elems_dict[self.reservoir_2_label]
            for index in range(len(first_dict)):
                pipe_turning_point = self.get_turning_point(index, pipe_dict)
                if pipe_turning_point != None:
                    pipe_turning_points[index] = pipe_turning_point
                    
                reservoir_1_turning_point = self.get_turning_point(index, reservoir_1_dict)
                if reservoir_1_turning_point != None:
                    reservoir_1_turning_points[index] = reservoir_1_turning_point

                reservoir_2_turning_point = self.get_turning_point(index, reservoir_2_dict)
                if reservoir_2_turning_point != None:
                    reservoir_2_turning_points[index] = reservoir_2_turning_point
        # return turnining points for the pipe and the two reservoirs
        return {self.pipe_label: pipe_turning_points, 
                self.reservoir_1_label: reservoir_1_turning_points, 
                self.reservoir_2_label: reservoir_2_turning_points
                }
    
    def get_all_turning_points_for_pipe(self, elem_dict):
        '''Assumes elems_dict is a dictionary that contains:
            key: a string (Pipe, Reservoir 1 or Reservoir 2)
            value: dictionary:
                        key: an integer
                        value: tuple of (time, head).
        '''
        pipe_dict = elem_dict[self.pipe_label]
        pipe_turning_points = {}
        # we must have at least 3 element to get a turning point
        if len(pipe_dict) > 2:
            for index in range(len(pipe_dict)):
                turning_point = self.get_turning_point(index, pipe_dict)
                if turning_point != None:
                    pipe_turning_points[index] = turning_point
        return {self.pipe_label: pipe_turning_points}
    
    def plot_single_graph(self, fig="", _title="", t_label="", z_label="", 
                          t_vals=[], z_vals=[], _label="", 
                          color_plus_lineType=['r--', 'b--'], line_thickness=1.0):
        plt.figure(fig)
        plt.clf()
        plt.xlabel(t_label)
        plt.ylabel(z_label)
        plt.plot(t_vals, z_vals, color_plus_lineType[0], 
                 label = _label, linewidth = line_thickness)
        plt.legend(loc = 'upper right')
        plt.title(_title)
        plt.show()

    def plot_single(self, method_vals, _type, method_label="Numerical"):
        tvals = []
        surge_vals = []
        time, head = 0, 1
        surge_dict = method_vals[_type]

        for key in ((surge_dict)):
            tvals.append(surge_dict[key][time])
            surge_vals.append(surge_dict[key][head])
            
        self.plot_single_graph(fig=f"{method_label} method ({_type})", 
                               _title=f"Graph of {_type} head Vs. Time",
                               t_label="time, t",
                               z_label="head, z", 
                               t_vals=tvals, 
                               z_vals=surge_vals,
                              _label=f"{_type} (viscosity: {self.v})", 
                              color_plus_lineType=['b--', 'r--'])

    def plot_multiple_graphs(self, fig="", _title="", t_label="", 
                             z_label="", t_vals=[], z_vals=[], labels=(), 
                          color_plus_lineType=['b--', 'g^', 'ro'], 
                          line_thickness=.5):
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
    
    def plot_multiple(self, method_vals, _type, method_label="Numerical"):
        '''Assumes _type: string
                  method_vals is a dictionary that contains:
                    key: a string (Pipe, Reservoir 1 or Reservoir 2)
                    value: dictionary;
                                key: an integer
                                value: tuple of (time, head).
        '''
        tvals = []
        pipe_vals, reservoir_1_vals, reservoir_2_vals = [], [], []
        time, head = 0, 1
        pipe_dict = method_vals[self.pipe_label]
        reservoir_1_dict = method_vals[self.reservoir_1_label]
        reservoir_2_dict = method_vals[self.reservoir_2_label]

        for key in ((pipe_dict)):
            tvals.append(pipe_dict[key][time])
            pipe_vals.append(pipe_dict[key][head])
            reservoir_1_vals.append(reservoir_1_dict[key][head])
            reservoir_2_vals.append(reservoir_2_dict[key][head])
            
        self.plot_multiple_graphs(fig=f"{method_label} method", 
                               _title=f"Graph of Pipe-Reservoirs head Vs. Time",
                               t_label="time, t", 
                               z_label="head, z", 
                               t_vals=tvals, 
                               z_vals=(pipe_vals, 
                                       reservoir_1_vals, 
                                       reservoir_2_vals),
                              labels=(self.pipe_label, 
                                      self.reservoir_1_label,
                                      self.reservoir_2_label), 
                              color_plus_lineType=['b--', 'g^', 'ro'])
            
    def drawTable(self, method_vals, pipe_label="", reservoir_1_label="", reservoir_2_label=""):
        table = PrettyTable()
        sn_title = "sn"
        time_title = "curr_time"
        if self.A1 and self.A2: # if we have two reservoirs
            table.field_names = [sn_title, time_title, pipe_label, 
                                 reservoir_1_label, reservoir_2_label]
            table.align = 'l' # align all values to the left
            for sn in (method_vals[self.pipe_label]):
                pipe = method_vals[self.pipe_label]
                reservoir_1 = method_vals[self.reservoir_1_label]
                reservoir_2 = method_vals[self.reservoir_2_label]
                
                t, z_pipe = pipe[sn]
                t, z_reservoir_1 = reservoir_1[sn]
                t, z_reservoir_2 = reservoir_2[sn]
                
                row = [sn, t, z_pipe, z_reservoir_1, z_reservoir_2]
                table.add_row(row)
                
        else:
            table.field_names = [sn_title, time_title, pipe_label]
            table.align = 'l' # align all values to the left
            for sn in (method_vals[self.pipe_label]):
                pipe = method_vals[self.pipe_label]
                t, z_pipe = pipe[sn]
                row = [sn, t, z_pipe]
                table.add_row(row)
        print(table)

if __name__ == '__main__':
    pass


