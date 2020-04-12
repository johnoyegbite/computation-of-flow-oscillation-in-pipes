# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 09:40:03 2019

@author: USER
"""
import time
import os
import webbrowser
from collections import OrderedDict
import matplotlib.pyplot as plt
from prettytable import PrettyTable

from Designs.general import numerical_method, html_output
from Designs.numerical_methods.multi_step_methods.explicit import adams_bashforth
from Designs.numerical_methods.multi_step_methods.implicit import adams_moulton
from Designs.numerical_methods.single_step_methods.explicit import runge_kutta_fourth_order
from Designs.numerical_methods.single_step_methods.implicit import trapezoidal_predictor_corrector



class FlowOscillation(numerical_method.NumericalMethod):
    
    def __init__(self, D, L, v, e, g=9.81):
        numerical_method.NumericalMethod.__init__(self, D, L, v, e, g)
        
        self.working_directory = os.getcwd()
        
        self.graph_img_loc = os.path.join(self.working_directory, "graphs")
        self.create_folder(self.graph_img_loc)
        
        self.image_format_png_ext = "png"
        self.image_format_pdf_ext = "pdf"
        
        self.html_path = os.path.join(self.working_directory, "results")
        self.create_folder(self.html_path)
        
        self.initial_h_label = "initial_h"
        self.final_h_label = "final_h"
        self.modified_t_z_label = "modified_t_z"
        self.modified_t_z__final_h_label = "modified_t_z__final"
        self.initial_t_z_label = "initial_t_z_label"
        self.total_time = "total_time"
        
        self.modified_t_z_A1_label = "modified_t_z_A1"
        self.modified_t_z_A2_label = "modified_t_z_A2"
    
    def create_folder(self, path):
        try:
            if not os.path.exists(path):
                os.makedirs(path)
        except OSError:
            print("Creation of the directory {path} failed!")

    def get_name_of_graph_folder(self):
        graph_folder_name = ""
        graph_folder_name += f"(D_{self.D})"
        graph_folder_name += f"-(L_{self.L})"
        graph_folder_name += f"-(errTol_{self.err_tolerance})"
        
        if self.A1 and self.A2:
            graph_folder_name += f"-(D1_{self.D1})"
            graph_folder_name += f"-(D2_{self.D2})"
        return graph_folder_name
        
    def set_reservoir_param(self, D1, D2, length_minor_losses, f=None):
        numerical_method.NumericalMethod.set_reservoir_param(self, D1, D2,
                                                             length_minor_losses,
                                                             f)
        self.D1 = D1
        self.D2 = D2
    
    def set_initial_conditions(self, t0, z0, dzdt, h, t, err_tolerance=0.001,
                               refined_tolerance=None):
        numerical_method.NumericalMethod.set_initial_conditions(self, t0, z0, 
                                                                dzdt, h, t)
        # Error tolerance for all method to meet
        self.err_tolerance = err_tolerance
        
        # Error tolerance for exact values which is the refined values 
        # obtained from Runge-Kutta.
        
        # Case 1: Use a constant value 
#        self.err_tolerance_for_refined = 0.000000001
        
        # Case 2: Or reduce the defined method's error tolerance by 10^4
        #         if the user does not input the tolerance that would be used 
        #         for refinining.
        #         If the user does input the refined tolerance and it's less 
        #         than the error tolerance, use it.
        if refined_tolerance and refined_tolerance <= err_tolerance:
            self.err_tolerance_for_refined = refined_tolerance
        else:
            self.err_tolerance_for_refined = self.err_tolerance * 0.0001
        
        self.graph_folder_name = self.get_name_of_graph_folder()
        self.img_path = os.path.join(self.graph_img_loc, self.graph_folder_name)
        self.create_folder(self.img_path)
        
        
    def get_t_z_from_step(self, h, method=""):
        # Depending on the method,
        method_dict = {}
        if method == self.trapezoidal_label:
            tpc = trapezoidal_predictor_corrector.TrapezoidalRule(
                self.D, self.L, self.v, self.e, self.g)
            tpc.set_initial_conditions(self.t0, self.z0, self.V0, h, self.t)
            if self.A1 and self.A2:
                tpc.set_reservoir_param(self.D1, self.D2, self.length_minor_losses,
                                        self.f)
            method_dict = tpc.trapezoidal_solution()[self.pipe_label]
        elif method == self.runge_kutta_label:
            rk = runge_kutta_fourth_order.RungeKutta(self.D, self.L, self.v, self.e, self.g)
            rk.set_initial_conditions(self.t0, self.z0, self.V0, h, self.t)
            if self.A1 and self.A2:
                rk.set_reservoir_param(self.D1, self.D2, self.length_minor_losses,
                                        self.f)
            method_dict = rk.runge_kutta_solution()[self.pipe_label]
        elif method == self.bashforth_label:
            ab = adams_bashforth.AdamsBashForth(self.D, self.L, self.v, self.e, self.g)
            ab.set_initial_conditions(self.t0, self.z0, self.V0, h, self.t)
            if self.A1 and self.A2:
                ab.set_reservoir_param(self.D1, self.D2, self.length_minor_losses,
                                        self.f)
            method_dict = ab.bash_forth_solution()[self.pipe_label]
        elif (method == self.moulton_label) or method == "":
            am = adams_moulton.AdamsMoulton(self.D, self.L, self.v, self.e, self.g)
            am.set_initial_conditions(self.t0, self.z0, self.V0, h, self.t)
            if self.A1 and self.A2:
                am.set_reservoir_param(self.D1, self.D2, self.length_minor_losses,
                                        self.f)
            method_dict = am.moulton_solution()[self.pipe_label]

        # Create an ordered dictionary.
        ordered_list = []
        for index in sorted(method_dict.keys()):
            ordered_list.append(method_dict[index])
        return OrderedDict(ordered_list)
    
    def get_exact_from_(self, method="Runge-Kutta"):
        err_tolerance_for_refined = self.err_tolerance_for_refined
        
        first_step = self.h
        refined_step = (1/2.0) * first_step
        
        # Get time and head from a particular h in a dictionary format
        # where key: time, value: head.
        
        # first_t_z is assumed to be the orginal time-head value.
        # refined_t_z is the refined time-head value by reducing the step size, h
        first_t_z = self.get_t_z_from_step(first_step, method)
        refined_t_z = self.get_t_z_from_step(refined_step, method)
        
        # Store the step you started with and also the time-head dictionary
        # you started with.
        initial_h = first_step
        initial_t_z = first_t_z.copy()
        
        # Keep refining the head (z) till the absolute difference the
        # refined_t_z and the first_t_z is less than the err_tolerance_for_refined.
        not_refined = True
        while not_refined:
            count = len(first_t_z)
            # I use first_t_z as the control for the loop because all the keys (time)
            # in first_t_z is a subset of refined_t_z 
            # i.e. refined_t_z would have keys (time)
            for _time in first_t_z:
                if _time in refined_t_z:
                    error = abs(first_t_z[_time] - refined_t_z[_time])
                    if error > err_tolerance_for_refined:
                        first_step = refined_step
                        refined_step = (1/2.0) * first_step
                        break
                    # Stop if we have reached the last key of the first_t_z
                    # and time-head dictionary is refined.
                    if (error < err_tolerance_for_refined) and count == 1:
                        not_refined = False
                        break
                count -= 1
            first_t_z = self.get_t_z_from_step(first_step)
            refined_t_z = self.get_t_z_from_step(refined_step)

        refined_h = refined_step
        refined_at_original_time = {}
        if self.A1 and self.A2: # store the refined {time: head} for both reservoirs
            refined_A1_at_original_time = {}
            refined_A2_at_original_time = {}
        for _time in initial_t_z:
            if (_time in refined_t_z):
                refined_at_original_time[_time] = refined_t_z[_time]
                if self.A1 and self.A2:
                    refined_A1_at_original_time[_time] = (refined_t_z[_time] * self.A)/self.A1
                    refined_A2_at_original_time[_time] = (refined_t_z[_time] * self.A)/self.A2
        
        if self.A1 and self.A2:
            return {
                    self.initial_h_label: initial_h, 
                    self.final_h_label: refined_h, 
                    self.modified_t_z_label: refined_at_original_time,
                    self.modified_t_z__final_h_label: refined_t_z,
                    self.initial_t_z_label: initial_t_z,
                    self.modified_t_z_A1_label: refined_A1_at_original_time,
                    self.modified_t_z_A2_label: refined_A2_at_original_time
                    }
       
        return {
                self.initial_h_label: initial_h, 
                self.final_h_label: refined_h, 
                self.modified_t_z_label: refined_at_original_time,
                self.modified_t_z__final_h_label: refined_t_z,
                self.initial_t_z_label: initial_t_z
                }


    def is_step_tolerable(self, mid_step, exact_details, method="", 
                          err_tolerance=0.001):
        exact = exact_details[self.modified_t_z_label]
        
#        print(f"\n\t mid step: {mid_step} --- error tol: {err_tolerance}")
        step_dict = self.get_t_z_from_step(mid_step, method)
#        print()
#        print(f"{method}....")
#        print(step_dict)
        
        table = PrettyTable()
        table.field_names = ["sn", "time", "exact", "z_at_{mid_step}", "error"]
        sn = 0
        head_precision = str(len(str(self.err_tolerance)))+"f"
        
        is_tolerable = True
        
        for exact_time, exact_value in exact.items():
            time_list = sorted(step_dict.keys())
#            print(f"exact time: {exact_time}")
            for index, curr_time_step in enumerate(time_list):
                if curr_time_step == exact_time:
                    step_value = step_dict[curr_time_step]
#                    print(".........")
#                    print((exact_time, exact_value, step_value))
                    error = abs(exact_value - step_value)
                    table.add_row([f"{sn}", f"{exact_time}", f"{exact_value:.{head_precision}}", 
                                   f"{step_value:.{head_precision}}", f"{error:.{head_precision}}"])
                    sn += 1
                    if error > err_tolerance:
#                        print(f"{error}, {False}")
                        is_tolerable = False
                    break
                elif curr_time_step > exact_time or (index == len(time_list)-1
                                                     and curr_time_step != exact_time):
                    if index > 0: # to avoid IndexError on the next line.
                        prev_time_step = time_list[index - 1]
                        prev_z_value = step_dict[prev_time_step]
                        curr_z_value = step_dict[curr_time_step]
                        boundary_1 = [prev_z_value, curr_z_value]
                        boundary_2 = [prev_time_step, curr_time_step]
                        param_2 = exact_time
                        correct_z_value = self.interpolate(boundary_1, 
                                                           boundary_2,  param_2)
#                        print(".........")
#                        print((exact_time, exact_value, correct_z_value))
                        error = abs(exact_value - correct_z_value)
                        table.add_row([f"{sn}", f"{exact_time}", f"{exact_value:.{head_precision}}", 
                                   f"{correct_z_value:.{head_precision}}", f"{error:.{head_precision}}"])
                        sn += 1
                        if error > err_tolerance:
#                            print(f"{error}, {False}")
                            is_tolerable = False
                        break

#                print((exact_time, exact_value, step_value))
#        print(f"{mid_step}, {True}")
#        print(table)
#        print(step_dict)
        return is_tolerable

    def compare_exact_with(self, exact_details, method="", err_tolerance=0.001):
        """Compare a numerical method with the exact values gotten from
        Runge-kutta and refine it till it matches the exact values with the
        requried error tolerance.
        
        PS: the step size is controlled by reducing the step size by half (1/2)
            each time till the error tolerance is met then we track back to
            get the actual step size since one might have jumped above the 
            actual step size. 
            This is done by bisection method.
        """
        exact = exact_details[self.modified_t_z_label]
        exact_2 = exact_details[self.modified_t_z__final_h_label]
#        print("\n\nExact...")
#        print(exact)
        
        method_dict = self.get_t_z_from_step(self.h, method)
#        print(f"{method}.....")
#        print(method_dict)
#        print()
        
        step_tolerance = 0.00000001
        low_step = self.h
        high_step = low_step
        
        initial_h = self.h
        not_done = True
        while not_done:
            count = len(exact)
            for exact_time, exact_value in exact.items():
                if exact_time in method_dict:
                    error = abs(exact_value - method_dict[exact_time])
                    if error > err_tolerance:
                        # Also helps to keep track of the step before and after/on
                        # ascertaining
                        low_step = high_step
                        high_step = (1/2.0) * low_step
                        break
                    if (error < err_tolerance) and count == 1:
                        not_done = False
                        # Use bisection method to get the actual step size.
#                        mid_step = (high_step + low_step) / 2.0
                        mid_step = high_step + 0.05
#                        print(f"\n\n{method}-----------------{low_step}, {high_step}")
                        not_seen_actual_step_size = True
##                        while abs(high_step - low_step) > step_tolerance:
#                        while not_seen_actual_step_size:
#                            is_step_tolerable = self.is_step_tolerable(mid_step,
#                                                                       exact_details, 
#                                                                       method, 
#                                                                       err_tolerance)
#                            if is_step_tolerable:
#                                high_step = mid_step
##                                mid_step = (high_step + low_step) / 2.0
#                                mid_step = high_step + 0.05
#                            else:
#                                not_seen_actual_step_size = False
#                                low_step = mid_step
##                                mid_step = (high_step + low_step) / 2.0
#                                mid_step = high_step
##                            print(f"\n-----------------{low_step}, {high_step}")
                        break
                count -= 1
            method_dict = self.get_t_z_from_step(high_step, method)
            
        final_h = high_step # change this later
        modified_at_original_time = {}
        
        if self.A1 and self.A2: # store the refined {time: head} for both reservoirs
            modified_A1_at_original_time = {}
            modified_A2_at_original_time = {}
            
        for _time in exact:
            if (_time in method_dict):
                modified_at_original_time[_time] = method_dict[_time]
                if self.A1 and self.A2:
                    modified_A1_at_original_time[_time] = (method_dict[_time] * self.A)/self.A1
                    modified_A2_at_original_time[_time] = (method_dict[_time] * self.A)/self.A2
                    
        if self.A1 and self.A2:
            return {
                    self.initial_h_label: initial_h, 
                    self.final_h_label: final_h, 
                    self.modified_t_z_label: modified_at_original_time, 
                    self.modified_t_z_A1_label: modified_A1_at_original_time,
                    self.modified_t_z_A2_label: modified_A2_at_original_time
                    }  
        return {
                self.initial_h_label: initial_h, 
                self.final_h_label: final_h, 
                self.modified_t_z_label: modified_at_original_time
                }

    def compare_all_methods(self, tolerance=0.001, to_print=False):
        all_methods_info_dict = {}
        # Get refined runge-kutta as the exact values and store the time
        # required to do that.
#        time_before_refined = time.time()
        refined_details = self.get_exact_from_(self.runge_kutta_label)
#        time_after_refined = time.time()
#        total_time_refined = time_after_refined - time_before_refined
#        refined_details[self.total_time] = total_time_refined
        all_methods_info_dict[self.exact_label] = refined_details
        
        # Compare the trapezoidal method with the exact values gotten from
        # Runge-kutta and refine it till it matches the exact values with the
        # requried error tolerance and also store the time required to do that.
        time_before_trapezoidal = time.time()
        modified_trapezoidal_details = self.compare_exact_with(refined_details, 
                                                       method=self.trapezoidal_label, 
                                                       err_tolerance=tolerance)
        time_after_trapezoidal = time.time()
        total_time_trapezoidal = time_after_trapezoidal - time_before_trapezoidal
        modified_trapezoidal_details[self.total_time] = total_time_trapezoidal
        all_methods_info_dict[self.trapezoidal_label] = modified_trapezoidal_details
        
        # Compare the Adams Bashforth method with the exact values gotten from
        # Runge-kutta and refine it till it matches the exact values with the
        # requried error tolerance and also store the time required to do that.
        time_before_bashforth = time.time()
        modified_bashforth_details = self.compare_exact_with(refined_details, 
                                                       method=self.bashforth_label, 
                                                       err_tolerance=tolerance)
        time_after_bashforth = time.time()
        total_time_bashforth = time_after_bashforth - time_before_bashforth
        modified_bashforth_details[self.total_time] = total_time_bashforth
        all_methods_info_dict[self.bashforth_label] = modified_bashforth_details
#        print(modified_bashforth_details)
        
        # Compare the Runge-Kutta method with the exact values gotten from
        # Runge-kutta and refine it till it matches the exact values with the
        # requried error tolerance and also store the time required to do that.
        time_before_runge_kutta = time.time()
        modified_runge_kutta_details = self.compare_exact_with(refined_details, 
                                                       method=self.runge_kutta_label, 
                                                       err_tolerance=tolerance)
        time_after_runge_kutta = time.time()
        total_time_runge_kutta = time_after_runge_kutta - time_before_runge_kutta
        modified_runge_kutta_details[self.total_time] = total_time_runge_kutta
        all_methods_info_dict[self.runge_kutta_label] = modified_runge_kutta_details
        
        # Compare the Adams Moulton method with the exact values gotten from
        # Runge-kutta and refine it till it matches the exact values with the
        # requried error tolerance and also store the time required to do that.
        time_before_moulton = time.time()
        modified_moulton_details = self.compare_exact_with(refined_details, 
                                                       method=self.moulton_label, 
                                                       err_tolerance=tolerance)
        time_after_moulton = time.time()
        total_time_moulton = time_after_moulton - time_before_moulton
        modified_moulton_details[self.total_time] = total_time_moulton
        all_methods_info_dict[self.moulton_label] = modified_moulton_details
        
  
        if to_print:
            print("Exact...")
            print(refined_details)
            print("--------------------------\n")
            
            print("Trapezoidal...")
            print(modified_trapezoidal_details)
            print(total_time_trapezoidal)
            print("--------------------------\n")

            print("Runge-Kutta...")
            print(modified_runge_kutta_details)
            print(total_time_runge_kutta)
            print("--------------------------\n")
                    
            print("Adams Bashforth...")
            print(modified_bashforth_details)
            print(total_time_bashforth)
            print("--------------------------\n")
                    
            print("Adams Moulton...")
            print(modified_moulton_details)
            print(total_time_moulton)
            print("--------------------------\n")
            
        return all_methods_info_dict
    
    def plot_multiple_graph(self, fig="", 
                            _title="", 
                            t_label="", 
                            z_label="", 
                          numerical_details=[], 
                          line_thickness=1.0):
        numerical_labels = ['r--', 'b--', 'g--', 'k--', 'b^']
        plt.figure(fig)
        plt.clf()
        plt.xlabel(t_label)
        plt.ylabel(z_label)
        
        time_vals = 0
        head_vals = 1
        _label = 2
        for numerical_method_index in range(len(numerical_details)):
            numerical_method = numerical_details[numerical_method_index]
            plt.plot(numerical_method[time_vals],
                     numerical_method[head_vals], 
                     numerical_labels[numerical_method_index], 
                     label = numerical_method[_label], 
                     linewidth = line_thickness)
        plt.legend(loc='upper right')
        plt.title(_title)
        # To save the plot as format in the desired location.
        # bbox_inches: removes undesirable, whitespace around the image.
        img_loc = self.img_path
        img_name = fig
        
        img_png_ext = self.image_format_png_ext
        img_pdf_ext = self.image_format_pdf_ext
        
        png_path = os.path.join(img_loc, f"{img_name}.{img_png_ext}")
        pdf_path = os.path.join(img_loc, f"{img_name}.{img_pdf_ext}")
        
        plt.savefig(png_path, bbox_inches='tight')
        plt.savefig(pdf_path, bbox_inches='tight')
        
        plt.close()
    
    def plot_all_numerical_graph(self, all_numerical_method_dict, 
                                 modified_label, fluid_type="Water", _type=""):
        # tp_t_l implies trapezoidal time list
        # tp_z_l implies trapezoidal head list
        # rk == Runge-Kutta, ab == Adams Bashforth, am == Adams Moulton
        tp_t_l, tp_z_l = [], []
        rk_t_l, rk_z_l = [], []
        ab_t_l, ab_z_l = [], []
        am_t_l, am_z_l = [], []
        
        trapezoidal_dict = all_numerical_method_dict[self.trapezoidal_label][modified_label]
        runge_kutta_dict = all_numerical_method_dict[self.runge_kutta_label][modified_label]
        bashforth_dict = all_numerical_method_dict[self.bashforth_label][modified_label]
        moulton_dict = all_numerical_method_dict[self.moulton_label][modified_label]

        for key in all_numerical_method_dict[self.runge_kutta_label][modified_label]:
            if key in trapezoidal_dict:
                tp_t_l.append(key)
                tp_z_l.append(trapezoidal_dict[key])
            # It is expected that key would be in Runge-Kutta since it is used
            # as the control.
            rk_t_l.append(key)
            rk_z_l.append(runge_kutta_dict[key])
            if key in bashforth_dict:
                ab_t_l.append(key)
                ab_z_l.append(bashforth_dict[key])
            if key in moulton_dict:
                am_t_l.append(key)
                am_z_l.append(moulton_dict[key])
        
        # Store the x-axis, y-axis and label values for plotting
        trapezoidal_plot = [tp_t_l, tp_z_l, self.trapezoidal_label]
        runge_kutta_plot = [rk_t_l, rk_z_l, self.runge_kutta_label]
        bashforth_plot = [ab_t_l, ab_z_l, self.bashforth_label]
        moulton_plot = [am_t_l, am_z_l, self.moulton_label]
        
        numerical_details = [trapezoidal_plot, runge_kutta_plot, 
                             bashforth_plot, moulton_plot]
        
        # Function that plot all numerical methods
        self.plot_multiple_graph(fig=f"{_type}-Numerical methods", 
                               _title=f"Graph of ({fluid_type}) HEAD against TIME",
                               t_label="Time, t", 
                               z_label=f"{_type} - Head, z", 
                               numerical_details=numerical_details
                               )
        
    def draw_table(self, all_numerical_method_dict, modified_label):
        table = PrettyTable()
        table.field_names = ["Sn", "Time, t", f"{self.exact_label} (z)", 
                             self.trapezoidal_label, self.runge_kutta_label,
                             self.bashforth_label, self.moulton_label]

        exact_dict = all_numerical_method_dict[self.exact_label]
        trapezoidal_dict = all_numerical_method_dict[self.trapezoidal_label]
        runge_kutta_dict = all_numerical_method_dict[self.runge_kutta_label]
        bashforth_dict = all_numerical_method_dict[self.bashforth_label]
        moulton_dict = all_numerical_method_dict[self.moulton_label]
        
        i = 0
        head_precision = str(len(str(self.err_tolerance)))+"f"
        time_precision = "1f"
        ex_z, tp_z, ab_z, am_z = -1, -1, -1, -1
        for _time in runge_kutta_dict[modified_label]:
            if _time in exact_dict[modified_label]:
                ex_z = exact_dict[modified_label][_time]
            if _time in trapezoidal_dict[modified_label]:
                tp_z = trapezoidal_dict[modified_label][_time]
            
            rk_t = (_time)
            rk_z = (runge_kutta_dict[modified_label][_time])
            
            if _time in bashforth_dict[modified_label]:
                ab_z = (bashforth_dict[modified_label][_time])
            if _time in moulton_dict[modified_label]:
                am_z = (moulton_dict[modified_label][_time])
            
            sn = i if len(str(i)) > 1 else str(i).zfill(2)
            i += 1
            rk_t = f"{rk_t:.{time_precision}}"
            ex_z = f"{ex_z:.{head_precision}}"
            tp_z = f"{tp_z:.{head_precision}}"
            rk_z = f"{rk_z:.{head_precision}}"
            ab_z = f"{ab_z:.{head_precision}}"
            am_z = f"{am_z:.{head_precision}}"
            
            table.add_row([sn, rk_t, ex_z, tp_z, rk_z, ab_z, am_z])
        seperate = "-"*5
        table.add_row([seperate, seperate, seperate, seperate,
                       seperate, seperate, seperate])
        table.add_row(["-", f"Initial h ({self.dimen_unit}): ", 
                       f"{seperate}",
                       f"{trapezoidal_dict[self.initial_h_label]:.{head_precision}}",
                       f"{runge_kutta_dict[self.initial_h_label]:.{head_precision}}",
                       f"{bashforth_dict[self.initial_h_label]:.{head_precision}}",
                       f"{moulton_dict[self.initial_h_label]:.{head_precision}}"
                       ])
        table.add_row(["-", f"Final h ({self.dimen_unit}): ", 
                       f"{seperate}",
                       f"{trapezoidal_dict[self.final_h_label]:.{head_precision}}",
                       f"{runge_kutta_dict[self.final_h_label]:.{head_precision}}",
                       f"{bashforth_dict[self.final_h_label]:.{head_precision}}",
                       f"{moulton_dict[self.final_h_label]:.{head_precision}}"
                       ])
        table.add_row(["-", f"Time t ({self.time_unit}): ",
                       f"{seperate}",
                       f"{trapezoidal_dict[self.total_time]:.{head_precision}}",
                       f"{runge_kutta_dict[self.total_time]:.{head_precision}}",
                       f"{bashforth_dict[self.total_time]:.{head_precision}}",
                       f"{moulton_dict[self.total_time]:.{head_precision}}"
                       ])

        return table.get_html_string()            
    
    
    def get_all_numerical_values(self, to_print=False, dimen_unit="m", 
                                 fluid_type="Water", time_unit="s"):
        self.dimen_unit = dimen_unit
        self.time_unit = time_unit        
        all_numerical_method_dict = self.compare_all_methods(self.err_tolerance, 
                                                         to_print)
        if to_print:
            print("Pipe...")
        pipe_html_table = self.draw_table(all_numerical_method_dict, 
                                          self.modified_t_z_label)
        self.plot_all_numerical_graph(all_numerical_method_dict,
                                      self.modified_t_z_label,
                                      fluid_type,
                                      _type=self.pipe_label)
        # Plot graphs and generate tables for the reservoir if available.
        reservoir_1_html_table = None
        reservoir_2_html_table = None
        if self.A1 and self.A2:
            if to_print:
                print()
                print("Reservoir 1...")
            reservoir_1_html_table = self.draw_table(all_numerical_method_dict, 
                                                self.modified_t_z_A1_label)
            if to_print:
                print()
                print("Reservoir 2...")
            reservoir_2_html_table = self.draw_table(all_numerical_method_dict, 
                                                self.modified_t_z_A2_label)
            self.plot_all_numerical_graph(all_numerical_method_dict,
                                      self.modified_t_z_A1_label,
                                      fluid_type,
                                      _type=self.reservoir_1_label)
            self.plot_all_numerical_graph(all_numerical_method_dict,
                                      self.modified_t_z_A2_label,
                                      fluid_type,
                                      _type=self.reservoir_2_label)
        
        date_time = time.strftime("%A %B, %Y @ %I:%M %p")

        pipe_img_src = os.path.join(self.img_path, f"{self.pipe_label}-Numerical methods.{self.image_format_png_ext}")
        reservoir1_img_src = os.path.join(self.img_path, f"{self.reservoir_1_label}-Numerical methods.{self.image_format_png_ext}")
        reservoir2_img_src = os.path.join(self.img_path, f"{self.reservoir_2_label}-Numerical methods.{self.image_format_png_ext}")
        
        html_class = html_output.HTML(date_time, self.D, self.L, self.v, 
                                         self.e, self.g, self.t0, self.z0, 
                                         self.V0, self.h, self.t, self.err_tolerance, 
                                         self.err_tolerance_for_refined, self.D1, self.D2, 
                                         self.length_minor_losses, self.f,
                                         pipe_html_table, reservoir_1_html_table,
                                         reservoir_2_html_table, pipe_img_src, 
                                         reservoir1_img_src, reservoir2_img_src,
                                         self.pipe_label, self.reservoir_1_label, 
                                         self.reservoir_2_label,
                                         dimen_unit, time_unit, fluid_type)
        html = html_class.get_html()
        date_time_for_folder = str(time.strftime("%a_%b_%Y-%I_%M_%p")) + "-" + str(time.time())
        html_path = os.path.join(self.html_path, date_time_for_folder)
        self.write_to_html(html, html_path)
        self.open_in_webbrowser(html_path)
        
    def write_to_html(self, html, html_path):
        if not os.path.exists(html_path):
            os.makedirs(html_path)
            file_path = os.path.join(html_path, "flow.html")
            with open(file_path, "w") as file:
                file.write(html)
                return file_path
                
    def open_in_webbrowser(self, html_path):
        new = 2 # open in a new tab, if possible.
        full_path = os.path.abspath(html_path) # get the absolute path of the file.
        file_path = os.path.abspath(os.path.join(html_path, "flow.html"))
        webbrowser.open_new_tab(file_path) # open html file on default browser.
        webbrowser.open(full_path, new=new) # open html file on my own (Windows) computer.
        
        
        
if __name__ == "__main__":     
###### Formulated Test Case
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
#    err_tolerance = 0.001

#    flow = FlowOscillation(D, L, v, e, g)
#    flow.set_initial_conditions(t0, z0, dzdt, h, t, err_tolerance)
         
####### Streeter reservoir example (Ex. 10.15)
#        # pipe properties
#    D = 3 # diameter of pipe; ft
#    L = 2000 # length of pipe; ft
#    v = 1e-4 # kinematic viscosity; f^2/s
#    e = 0.007093 # Roughness height; ft
#        # initial conditions
#    t0 = 0 # time; secs
#    z0 = 1132 # head; ft
#    dzdt = 0 # velocity; ft/s
#    t = 500 # final time; secs
#    h = 5 # time step; secs
#    g = 32.2 # acceleration due to gravity; ft/secs^2
#    err_tolerance = 0.001
##    refined_tolerance = 0.000001
#
#    D1 = 15.958 # diameter of first reservoir; ft
#    D2 = 19.544 # diameter of second reservoir; ft
#    length_minor_losses  = 438 # minor losses; ft
#    
##     MAKE SURE you set the Reservoir parmeters before the initial conditons
#    flow = FlowOscillation(D, L, v, e, g)
#    flow.set_reservoir_param(D1, D2, length_minor_losses)
#    flow.set_initial_conditions(t0, z0, dzdt, h, t, err_tolerance)
##        
###### CEG 848, UNILAG, M.Sc(Civil Engr) Exam, Second Semester 2017/2018 Session, 
#    # Reservoir Question 1b.
#        # pipe properties
    D = 0.8 # diameter of pipe; m
    L = 2000 # length of pipe; m
    v = 1e-4 # kinematic viscosity; m^2/s
    e = 0.007093 # Roughness height; m
        # initial conditions
    t0 = 0 # time; secs
    z0 = 7.8125 # head; m
    dzdt = 0 # velocity; m/s
    t = 1000 # final time; secs
    h = 5 # time step; secs
    g = 9.81 # acceleration due to gravity; m/secs^2
    err_tolerance = 0.001
#      # reservoir properties
    D1 = 10 # diameter of first reservoir; m
    D2 = 8 # diameter of second reservoir; m
    length_minor_losses  = 160 # minor losses == 8% of pipe length; m
#    
    # MAKE SURE you set the Reservoir parmeters before the initial conditons
    flow = FlowOscillation(D, L, v, e, g)
    flow.set_reservoir_param(D1, D2, length_minor_losses)
    flow.set_initial_conditions(t0, z0, dzdt, h, t, err_tolerance)
#
###### Formulated Test Case
#        #<b> pipe properties </b> <br>
#    D = 3 # diameter of pipe; ft <br>
#    L = 2000 # length of pipe; ft <br>
#    v = 0.0001 # kinematic viscosity; f^2/s <br>
#    e = 0.006 # Roughness height; ft <br>
#        #<b> initial conditions </b>  <br>
#    t0 = 0 # time; secs <br>
#    z0 = 0 # head; ft <br>
#    dzdt = 3 # velocity; ft/s <br>
#    t = 50 # final time; secs <br>
#    h = 0.5 # time step; secs <br>
#    g = 32.2 # acceleration due to gravity; ft/secs^2
#    err_tolerance = 0.001
    
#    flow = FlowOscillation(D, L, v, e, g)
#    flow.set_initial_conditions(t0, z0, dzdt, h, t, err_tolerance)
         
###### Streeter Reservoir example (Problem 10.56)
#      i.e hyperbolic theoretical solution. 
#        # pipe properties
#    D = 8 # diameter of pipe; ft
#    L = 3000 # length of pipe; ft
#    v = 0.16146 # kinematic viscosity; ft^2/s
#    e = 0.0005 # Roughness height; ft
#        # initial conditions
#    t0 = 0 # time; secs
#    z0 = 93.75 # head; ft
#    dzdt = 0 # velocity; ft/s
#    t = 500# final time; secs
#    h = 1 # time step; secs
#    g = 32.2 # acceleration due to gravity; ft/secs^2
#    err_tolerance = 0.001
##
#    D1 = 20 # diameter of first reservoir; ft
#    D2 = 20 # diameter of second reservoir; ft
#    length_minor_losses  = 1800 # minor losses == 8% of pipe length; ft
#    
    # MAKE SURE you set the Reservoir parmeters before the initial conditons
#    flow = FlowOscillation(D, L, v, e, g)
#    flow.set_reservoir_param(D1, D2, length_minor_losses)
#    flow.set_initial_conditions(t0, z0, dzdt, h, t, err_tolerance)
#
###### Formulated Test Case
    
    to_print = False
    dimen_unit = "m"
    fluid_type = "Water"
    print(flow.get_all_numerical_values(to_print, dimen_unit, fluid_type))
    
    
    
        
    
