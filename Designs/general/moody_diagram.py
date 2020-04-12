# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 00:32:13 2019

@author: John Oyegbite
"""

from math import sqrt, log
from prettytable import PrettyTable

class FrictionFactor(object):
    
    def __init__(self, e_D):
        """Assumes e_D is a Real number (conventionaly between 0 and 1)
        
        """
        self.e_D = e_D
        self.laminar_label = "Laminar"
        self.critical_label = "Critical"
        self.turbulent_label = "Turbulent"
    
    def compute_f_ColeBrook(self, R):
        """Assumes R is a Real number
        
        Computes f using the COLEBROOK's equation
        # See Google for COLEBROOK's equation for friction factor
        It's implicit equation that has the friction factor on both the LHS 
        and RHS
        """
        # assume a starting correct value for the "f" on the right hand side (RHS)
        # uses friction factor from Barr's equation as the starting value.
        f_initial = self.compute_f_BARR(R)
        
        relative_roughness = self.e_D
        a = relative_roughness / 3.71
        b = 2.51 / (R * sqrt(f_initial))
    
        # Compute the f on the LHS ------ (1) 
        f_final = 1 / (-2 * log((a + b), 10))**2
        
        # Make sure friction factor is correct to at least 6 decimal place.
        tolerance = 0.0000001
        
        # if the f on the LHS is not within tolerance limit,
        # replace it on the RHS and recompute the f on the LHS till it's within
        # tolerance limit.
        while abs(f_final - f_initial) >= tolerance:
            f_initial = f_final
            b = 2.51 / (R * sqrt(f_initial))
            f_final = 1 / (-2 * log((a + b), 10))**2
            
        return f_final
            
    def compute_f_BARR(self, R):
        """Assumes R is a Real number
        
        Computes f using the BARR's friction equation
        See Google for BARR's equation for friction factor
        """
        upper_LHS_RHS = 5.02 * log(R / (4.518 * log((R / 7), 10)), 10)
        lower_LHS_RHS = R * (1 + (R**0.52) / (29 * (1/self.e_D)**0.7))

        upper_RHS_RHS = 1
        lower_RHS_RHS = 3.7*(1/self.e_D)
        bracket = (upper_LHS_RHS / lower_LHS_RHS) + (upper_RHS_RHS / lower_RHS_RHS)
        
        RHS = (-2)*log(bracket, 10)
        f = (1/RHS)**2

        return f
      
    def interpolate(self, boundary_1, boundary_2,  param_2_middle):
        """Assumes boundary_1 is a list of two Real numbers
                   boundary_2 is a list of two Real numbers
                   param_2_middle is a Real numbers.
        Returns: a Real number that correspond to param_2_middle, 
                which is the interpolation of boundary_1 and boundary_2
        """
        lower, upper = 0, 1
        lower_1, upper_1 = boundary_1[lower], boundary_1[upper]
        lower_2, upper_2 = boundary_2[lower], boundary_2[upper]
        param_1_middle = (((param_2_middle - lower_2)*(upper_1 - lower_1))/(upper_2 - lower_2)) + lower_1
        
        return param_1_middle
        
    def get_f_from_all_laminar_Reynolds_number(self):
        """Returns all the friction factor that correspond to the
        Reynolds number from 0 to 2000.
        
        """
        # laminar Reynolds number span from 0 to 2000 (2000 not inclusive)
        R_fr = 1
        R_to = 2000 + 1
        # note that f = 64/R from Darcy's equation.
        laminar_friction_factor_and_Reynolds = \
            [(64/R, R) if R != 0 else (0, R) for R in range(R_fr, R_to)]
            
        return laminar_friction_factor_and_Reynolds

    def get_f_from_all_critical_Reynolds_number(self):
        """Returns all the friction factor that correspond to the
        Reynolds number from 2000 to 4000.
        
        """
        # critical Reynolds number span from 2000 to 4000 inclusive
        R_from_2000 = 2000
        R_to_4000 = 4000 
        
        # Computes from laminar
        f_laminar_2000 = 64 / R_from_2000
        # Computes from turbulent
        f_turbulent_4000 = self.compute_f_ColeBrook(R_to_4000)
        
        critical_friction_factor_and_Reynolds = []
        for R in range(R_from_2000+1, R_to_4000 + 1):
            f_critical = self.interpolate([f_laminar_2000, f_turbulent_4000],
                                              [R_from_2000, R_to_4000],
                                              R)
            critical_friction_factor_and_Reynolds.append((f_critical, R))
            
        return critical_friction_factor_and_Reynolds
    
    def get_f_from_all_turbulent_Reynolds_number(self):
        """Returns all the friction factor that correspond to the
        Reynolds number from 4000 to 10^6.

        """
        R_fr_4000 = 4000 + 1
        R_to_10_6 = 1e6
        
        turbulent_friction_factor_and_Reynolds = []
        for R in range(R_fr_4000, int(R_to_10_6)+1):
            f_turbulent = self.compute_f_ColeBrook(R)
            turbulent_friction_factor_and_Reynolds.append((f_turbulent, R))
            
        return turbulent_friction_factor_and_Reynolds
    
    def drawMoodyTable(self, f_R):
        """Assumes f_R is a list of tuples
        where the first element is the friction factor and
              the second element is the Reynolds number
        
        """
        print(f"e/D: {self.e_D} \n")
        table = PrettyTable()
        
        sn_title = "sn"
        friction_factor_title = "friction factor"
        reynolds_no_title = "Reynolds number"
        
        # set the header of the table
        table.field_names = [sn_title, friction_factor_title, reynolds_no_title]
        
        # set text-align to left of the columns
        table.align[sn_title] = 'l'
        table.align[friction_factor_title] = 'l'
        table.align[reynolds_no_title] = 'l'
        
        for index in range(len(f_R)):
            f_index, R_index = 0, 1
            table.add_row([index, f_R[index][f_index], f_R[index][R_index]])
            
        print(table)

    def plot_single_graph(self, fig="", _title="", R_label="", f_label="", 
                          R_vals=[], f_vals=[], _label="",
                          color_plus_lineType=['r--', 'b--'], line_thickness=1.0,
                          R_lim = (0, 1e6),
                          f_lim = (0.008, 0.1)):
        '''R_vals: a list of numbers to be plotted
            f_vals:  a list of numbers to be plotted
        '''
        import matplotlib.pyplot as plt
        plt.figure(fig)
        plt.clf()
        plt.xlim(R_lim)
        plt.ylim(f_lim)
        plt.xlabel(R_label)
        plt.ylabel(f_label)
        plt.plot(R_vals, f_vals, color_plus_lineType[0], 
                 label = _label, linewidth = line_thickness)
        plt.legend(loc = 'upper right')
        plt.title(_title)
        plt.show()
            
    def plotMoodyChart(self, f_R, Rlim=(0, 4000), flim=(0.008, 0.1), single=""):
        """Assumes f_R is a list of tuples
        where the first element is the friction factor and
              the second element is the Reynolds number
        
        """
        import numpy as np
        f_R_np = np.array(f_R)

        fvals = f_R_np[:,0] # get all the friction factor values from first index
        Rvals = f_R_np[:,1] # get all the Reynolds number values from second index
        
        non_laminar_text = f" for values of e/D = {self.e_D}"
        
        title = f"Friction factor f Vs {int(Rlim[0])} < R < {int(Rlim[1])}"
        
        if single != self.laminar_label: # include the ending label text 
                                        # for critical and turbulent region
            title += non_laminar_text
            
        title += " || Matric No: 140402032"
        
        self.plot_single_graph(fig=f"{single}", 
                               _title=title,
                               R_label="Reynolds number, R",
                               f_label="friction factor, f", 
                               R_vals=Rvals, 
                               f_vals=fvals,
                              _label=f"{single}", 
                              color_plus_lineType=['b--', 'r-'],
                              line_thickness=1.0,
                              R_lim = Rlim,
                              f_lim = flim)
        
    def moodyDiagram(self):
        f_R_laminar = self.get_f_from_all_laminar_Reynolds_number()
#        R_lim = (1,  2000)
##        f_lim = (0, 1)
#        self.plotMoodyChart(f_R_laminar, R_lim, f_lim, 
#                            self.laminar_label)
        self.drawMoodyTable(f_R_laminar)

#        f_R_critical = self.get_f_from_all_critical_Reynolds_number()
#        R_lim = (2000, 4000)
#        f_lim = (0.008, 0.1)
#        self.plotMoodyChart(f_R_critical, R_lim, f_lim, 
#                            self.critical_label + f" (e/D: {self.e_D})")
#        self.drawMoodyTable(f_R_critical) 

#        f_R_turbulent = self.get_f_from_all_turbulent_Reynolds_number()
#        R_lim = (4000, 1e6)
#        f_lim = (0.008, 0.1)
#        self.plotMoodyChart(f_R_turbulent, R_lim, f_lim, 
#                            self.turbulent_label + f" (e/D: {self.e_D})")
#        self.drawMoodyTable(f_R_turbulent)
      
if __name__ == "__main__":
    e_D = 0.00001   
    F = FrictionFactor(e_D)
    F.moodyDiagram()
   
        
        
        
        

