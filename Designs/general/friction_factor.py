# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 00:32:13 2019

@author: John Oyegbite
"""
from math import sqrt, log

class FrictionFactor(object):
    
    def __init__(self, dzdt, v, e, D):
        """Assumes dzdt, v, e and D are real numbers
        
        where dzdt is the velocity; unit in m/s or ft/s.
                 v is the kinematic viscosity; unit in m^2/s or ft^2/s.
                 e is the pipe roughness height; unit in m or ft.
                 D is the diameter of pipe; unit in m or ft.
        """
        self.V = dzdt
        self.v = v
        self.e = e
        self.D = D
        self.laminar_label = "Laminar"
        self.critical_label = "Critical"
        self.turbulent_label = "Turbulent"

    # Get the Reynolds number  from the formula for Reynolds number.
    def getR(self):
        """Computes Reynolds number from the formula"""
        # Reynolds number uses the absolute value of the velocity
        V = abs(self.V)
        return (V * self.D) / self.v # formula for Reynolds number
    
    # Compute friction factor from COLEBROOK's equation for friction factor.
    def compute_f_ColeBrook(self, R, e, D):
        """Assumes R, e and D are Real numbers
        
        Computes f using the COLEBROOK's equation
        # See Google for COLEBROOK's equation for friction factor
        It's implicit equation that has the friction factor on both the LHS 
        and RHS
        """
        # assume a starting correct value for the "f" on the right hand side (RHS)
        # uses friction factor from Barr's equation as the starting value.
        f_initial = self.compute_f_BARR(R, self.e, self.D)

        relative_roughness = e / D
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
    
    # Compute friction factor, f from BARR's equation for friction factor.
    def compute_f_BARR(self, R, e, D):
        """Assumes R, e and D are Real numbers
        
        Computes f using the BARR's friction equation
        See Google for BARR's equation for friction factor
        """
        upper_LHS_RHS = 5.02 * log(R / (4.518 * log((R / 7), 10)), 10)
        lower_LHS_RHS = R * (1 + (R**0.52) / (29 * (D/e)**0.7))

        upper_RHS_RHS = 1
        lower_RHS_RHS = 3.7*(D/e)
        bracket = (upper_LHS_RHS / lower_LHS_RHS) + (upper_RHS_RHS / lower_RHS_RHS)
        
        RHS = (-2)*log(bracket, 10)
        f = (1/RHS)**2
        
        return f
    
    # Compute friction factor, f from CHEN's equation for friction factor.
    def compute_f_CHEN(self, R, e, D):
        """Assumes R, e and D are Real numbers
        
        Computes f using the CHEN'S friction equation
        See Google for CHEN'S equation for friction factor
        """
        upper_RHS_LHS = e/D
        lower_RHS_LHS = 3.7065
        
        log_square = log((((upper_RHS_LHS)**1.1096)/2.8257) + (7.149/R)**0.8961, 10)
        LHS_RHS_RHS = (5.0452/R) * log_square
        
        curly_bracket = (upper_RHS_LHS / lower_RHS_LHS) - LHS_RHS_RHS
        
        RHS = (-4)*log(curly_bracket, 10)
        f = (1/RHS)**2

        return f
      
    def getF_and_R(self):
        """Computes friction factor based on a computed Reynolds number.
        
        Returns a dictionary of Reynolds number and friction factor.
        
        """
        R = self.getR()
        # Treat the laminar case: f = 64/R
        if R < 2000: # R < 2000 implies flow is laminar
            if R > 0:
                f = 64/R # Darcy's equation for friction factor
                return {"R": R, "f": f}
            elif R == 0:
                return {"R": R, "f": 0}

        # Treat the Turbulent case: use Colebrook's equation.
        elif R > 4000: # R > 4000 implies flow is turbulent.
            # store the final value of "f" after few iterations.
            f_brooks = self.compute_f_ColeBrook(R, self.e, self.D)
            f_chen = self.compute_f_CHEN(R, self.e, self.D)
            f_barr = self.compute_f_BARR(R, self.e, self.D)
            return {"R": R, "f": f_brooks, "f_barr": f_barr, "f_chen": f_chen}
        # Treat the Critical case: interpolate between laminar and turbulent
        # boundaries.
        else: # R is between 2000 and 4000 implies flow is critical.
            # interpolate between the two reynolds number boundaries
            # that is, the laminar and turbulent.
            R_upper_limit, R_middle, R_lower_limit = 4000, R, 2000
            f = 0.025
            f_upper_limit = self.compute_f_ColeBrook(R_upper_limit, self.e, self.D)
            f_lower_limit = 64/R_lower_limit
            f_middle = self.interpolate([f_lower_limit, f_upper_limit],
                                        [R_lower_limit, R_upper_limit],  R_middle)
            
            return {"R": R, "f": f_middle}
      
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
        param_1_middle = (((param_2_middle - lower_2)*(upper_1 - lower_1))/
                          (upper_2 - lower_2)) + lower_1
        return param_1_middle

      
if __name__ == "__main__":
#    dzdt = 2
#    v = 1e-6
#    e = 0.0001524
#    D = 0.8
    
    dzdt = 12
    v = 1e-4
    e = 0.07
    D = 1.398
    
    F = FrictionFactor(dzdt, v, e, D)
    print(F.getF_and_R())

        
        
        
        