# -*- coding: utf-8 -*-

def validate_compulsory_modules():
    all_installed = True
    
    libraries_installation = {
        "matplotlib": "pip install matplotlib",
        "wx": "pip install -U wxPython",
        "prettytable": "pip install prettytable",
        "webbrowser": "pip install webbrowser",
        "string": "pip install string"
        }
    
    error = "\nNOTE: You cannot use this program if you don't have the following"
    error += " libraries installed.\n\n"
    for library in libraries_installation:
        try:
            __import__(library)
        except ModuleNotFoundError:
            all_installed = False
            error += f"{library} is not installed\n"
            error += f"Please install {library}\n"
            error += f"To install {library}, "
            error += f"just open your command prompt and type:\n"
            error += f"   {libraries_installation[library]}   \n"
            error += f"Please make sure you have an Internet Connection\n\n"
    return all_installed, error



# make sure all important libraries are installed
all_libraries_installed, library_error = validate_compulsory_modules()

if not all_libraries_installed:
    print(library_error)
else: 
    ###############################################################################
    """GUI using wx not still in development
    """
    # import wx
    
    # class FlowPanel(wx.Panel):
    #     def __init__(self, parent):
    #         wx.Panel.__init__(self, parent)
            
    #         main_sizer = wx.BoxSizer(wx.VERTICAL)
            
    #         self.row_obj_dict = {}
    
    #         self.list_ctrl = wx.ListCtrl(
    #             self, size=(-1, 100), 
    #             style=wx.LC_REPORT | wx.BORDER_SUNKEN
    #         )
            
    #         self.list_ctrl.InsertColumn(0, 'Pipe Properties', width=100)
    #         self.list_ctrl.InsertColumn(1, 'Inital Conditions', width=100)
    #         self.list_ctrl.InsertColumn(2, 'Reservoir Parmeters', width=100)
    #         # sizer.Add(window, proportion, flag, margin)
    #         main_sizer.Add(
    #                 # window
    #                 self.list_ctrl, 
    #                 # proportion
    #                 0,         # don't make vertically stretchable
    #                 # flag
    #                 wx.ALL |   # and make border all around
    #                 wx.EXPAND, # make horizontally stretchable
    #                 # border
    #                 0         # set margin or border width to 5
    #                 )
            
    #         design_button = wx.Button(self, label='Design')
    #         design_button.Bind(wx.EVT_BUTTON, self.on_design)
    #         main_sizer.Add(design_button, 0, wx.ALL | wx.CENTER, 10) 
            
    #         self.SetSizer(main_sizer)

    #     def on_design(self, event):
    #         print("DESIGNED!")

    # class FlowFrame(wx.Frame):
    #     def __init__(self):
    #         wx.Frame.__init__(self, parent=None, title='FLOW OSCILLATIONS IN PIPES')       
    #         self.panel = FlowPanel(self)
  
    #         self.Show()
###############################################################################
    """Temporary user interaction
    """
    
    from prettytable import PrettyTable
    from Designs import flow_oscillation
    
    class ValidateInputs(object):
        def __init__(self):
            self.data = {}
            self.errors = []
            
            self.D_name = "pipe's diameter"
            self.L_name = "length of pipe"
            self.v_name = "kinematic viscosity"
            self.e_name = "roughness height"
            
            self.t0_name = "initial time"
            self.z0_name = "initial head"
            self.V0_name = "initial velocity"
            self.h_name = "time step"
            self.t_name = "final time"
            self.g_name = "acceleration due to gravity"
            self.err_tolerance_name = "error tolerance"
            
            self.D1_name = "first reservoir diameter"
            self.D2_name = "second reservoir diameter"
            self.length_minor_losses_name = "length of minor losses"
            self.f_name = "friction factor \n No friction factor? Press N!"
            
            self.dimen_unit_name = "unit of dimension (ft, m, inch, e.t.c)"
            self.fluid_type = "fluid type (water, oil, e.t.c)"
            
            self.data_names = [self.D_name, self.L_name, self.v_name, 
                              self.e_name,
                              self.t0_name, self.z0_name, self.V0_name,
                              self.h_name, self.t_name, self.g_name,
                              self.err_tolerance_name
                              ]
            
            self.reservoir_name = [self.D1_name, self.D2_name, 
                                   self.length_minor_losses_name,
                                   self.f_name]
            self.fluid_unit_name = [self.dimen_unit_name, self.fluid_type]
            
            self.data_names_with_reservoir = self.data_names[:]
            self.data_names_with_reservoir.extend(self.reservoir_name)
            self.data_names_with_reservoir.extend(self.fluid_unit_name)
            self.has_reservoir = False # default
            
            self.data_names.extend(self.fluid_unit_name)
        
        def is_real_number(self, char):
            try:
                float(char)
            except ValueError:
                return False
            return True
        
        def validate_data(self, char, prop_name):
            if not self.is_real_number(char):
                return f"\nError: Please enter a valid '{prop_name}'!"
            char = float(char)
            if(prop_name == self.t0_name or prop_name == self.z0_name or prop_name == self.V0_name):
                if char < 0:
                    return f"\nError: {prop_name} must be positive!"
            else:
                if char <= 0:
                    return f"\nError: {prop_name} must be greater than zero!"
                   
            return None
        
        def confirm_reservoir(self):
            while True:
                is_reservoir = input("\nIs your pipe connected to a reservoir? or Enter 'e' to stop\n (Y/N)?->: ")
                is_reservoir = is_reservoir.lower()
                if is_reservoir == 'e':
                    stop_design = self.ask_to_stop_design()
                    if stop_design:
                        return None
                if is_reservoir == 'y' or is_reservoir == 'yes' or is_reservoir == '1':
                    self.has_reservoir = True
                    data_names = self.data_names_with_reservoir
                    return data_names
                elif is_reservoir == 'n' or is_reservoir == 'no' or is_reservoir == '0':
                    data_names = self.data_names
                    return data_names
                else:
                    if is_reservoir != "e":
                        print("\nError: Wrong input!")
                    
        def print_header(self):
            print("----------------------------------------------------")
            print("Collecting informations...")
            print("----------------------------------------------------")
            
        def validate_friction_factor(self, prop_name, prop):
            if prop == 'n' or prop == 'no' or prop == '0':
                self.data[prop_name] = None
                return False
            else:
                error = (self.validate_data(prop, prop_name))
                if error == None:
                    self.data[prop_name] = float(prop)
                    return False
                else:
                    print(error)
            return True
        
        def ask_to_stop_design(self):
            while True:
                ask_to_stop = input("\nDo you want to stop Design?\n (Y/N)?->: ")
                ask_to_stop = ask_to_stop.lower()
                if ask_to_stop == 'y' or ask_to_stop == 'yes' or ask_to_stop == '1':
                    return True
                elif ask_to_stop == 'n' or ask_to_stop == 'no' or ask_to_stop == '0':
                    return False
                else:
                    print("\nError: Wrong input!")
                    
        def collect_data(self):
            self.print_header()
            data_names = self.confirm_reservoir()

            if data_names == None:
                return data_names
            
            for prop_name in data_names:
                print(f"\n{prop_name.upper()}----------------")
                stay_on_current = True
                while stay_on_current:
                    prop = input(f"Enter {prop_name} or 'e' to stop ->: ")
                    prop = prop.lower()
                    
                    stop_design = False # default value
                    # if the user do not wish to continue with the design at
                    # any point.
                    if prop == 'e':
                        stop_design = self.ask_to_stop_design()
                        if stop_design:
                            break
                    # check if property is friction factor; which can be empty
                    if prop_name == self.f_name:
                        stay_on_current = self.validate_friction_factor(prop_name, prop)
                    
                    # just store if property doesn't have to do with numbers
                    elif prop_name == self.dimen_unit_name or prop_name == self.fluid_type:
                        # if character is too long
                        max_char_len = 30
                        if len(prop) > max_char_len:
                            print("\nError: Too long!")
                        else:
                            self.data[prop_name] = prop
                            stay_on_current = False
                    # validate property with numbers
                    else:
                        error = (self.validate_data(prop, prop_name))
                        if error == None:
                            prop = float(prop)
                            check = True # for all properties
                            
                            # make sure pipe's diameter is less than the
                            # reservoir's diameter.
                            if prop_name == self.D1_name or prop_name == self.D2_name:
                                if prop <= self.data[self.D_name]:
                                    print(f"\nError: {prop_name} must be greater than {self.D_name}!")
                                    print(f"\n{self.D_name} is {self.data[self.D_name]}")
                                    check = False
                            # make sure initial time is less than the
                            # final time.
                            if prop_name == self.t_name:
                                if prop <= self.data[self.t0_name]:
                                    print(f"\nError: {prop_name} must be greater than {self.t0_name}!")
                                    check = False
                            # make sure error tolerance is less than 1.
                            if prop_name == self.err_tolerance_name:
                                if prop >= 1:
                                    print(f"\nError: {prop_name} must be less than one (1)!")
                                    check = False
                                    
                            if check:
                                self.data[prop_name] = prop
                                stay_on_current = False
                        else:
                            print(error)
                if stop_design:
                    return None # no data should be used
                print(f"\n----------END OF {prop_name.upper()}\n")
#            print(self.table_data(self.data))
            return self.data
        
        def table_data(self, data):
            # use PrettyTable() to generate a table
            table = PrettyTable()
            design_prop_title = "Design Properties"
            design_values_title = "Values"
            # Create the header
            table.field_names = [design_prop_title, design_values_title]
            
            # align header: "design_prop_title" values to the left
            table.align[design_prop_title] = 'l'
            # align header: "design_values_title" values to the left
            table.align[design_values_title] = 'l'
            
            for prop_name, value in data.items():
                table.add_row([prop_name, value])
                
            return table
        
        def design_flow(self):
            data = self.collect_data()
            
            if data != None: # if we have data
                table = self.table_data(data)
                print(table)
                print("-----------------------------------------------------------")
                print("\nKindly wait while the computer process your result...")
                print("\nResult would be displayed on your BROWSER soon!")
                
                D = data[self.D_name]
                L = data[self.L_name]
                v = data[self.v_name]
                e = data[self.e_name]
                
                t0 = data[self.t0_name]
                z0 = data[self.z0_name]
                dzdt = data[self.V0_name]
                h = data[self.h_name]
                t = data[self.t_name]
                g = data[self.g_name]
                err_tolerance = data[self.err_tolerance_name]
                
                dimen_unit = data[self.dimen_unit_name]
                fluid_type = data[self.fluid_type]
                
                flow = flow_oscillation.FlowOscillation(D, L, v, e, g)
                if self.has_reservoir:
                    D1 = data[self.D1_name]
                    D2 = data[self.D2_name]
                    length_minor_losses = data[self.length_minor_losses_name]
                    f = data[self.f_name]
                    flow.set_reservoir_param(D1, D2, length_minor_losses, f)
                flow.set_initial_conditions(t0, z0, dzdt, h, t, err_tolerance)
                
                to_print = False # for debugging purpose
                print("\nComputing result...")
                print("\nThis would take a moment...")
                print("\nPlease wait!")
                flow.get_all_numerical_values(to_print, dimen_unit, 
                                              fluid_type)
                print("\n-------------------------------------------------")
                print("DONE!")
                print("\nYour result should be displayed in the BROWSER now!")
                print("Or check the ouput folder and open it on a browser!")
                print("\nThank you for using this program.")
                print("-------------------------------------------------\n")
                
            else:
                print("\n-------------------------------------")
                print("Thank you for using this program.")
                print("\nYou can 'Run' the program again to start!!!")
                print("-------------------------------------\n")
            
            self.print_developer_info()
            
        def print_developer_info(self):
            print("\n-------------------------------------------------")
            print("DEVELOPER INFO------------------------------------")
            print("\nTo contact John Oyegbite for further improvement on this program:")
            print("\nEmail: johnoyegbite@gmail.com")
            print("\nLinkedin: https://www.linkedin.com/in/john-oyegbite-67bb9913a")
            print("\nMobile No.: +2348090565698")
            print("-------------------------------------------------\n")

            
    if __name__ == '__main__':
        # GUI yet to be created.
#        app = wx.App()
#        frame = FlowFrame()
#        app.MainLoop()
        
        # Created a temporary console USER interaction
        inputs = ValidateInputs()
        inputs.design_flow()
