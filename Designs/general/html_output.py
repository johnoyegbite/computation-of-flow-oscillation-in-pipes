# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 00:32:13 2019

@author: John Oyegbite
"""

"""Creates an HTML with values given"""
class HTML(object):
    
    def __init__(self, date_time, D, L, v, e, g ,
                 t0, z0, V0, h, t, err_tolerance, refined_tolerance,
                 D1, D2, length_minor_losses, f,
                 pipe_table, reservoir1_table, reservoir2_table,
                 pipe_img_src, reservoir1_img_src, reservoir2_img_src,
                 pipe_label, reservoir1_label, reservoir2_label,
                 dimen_unit="m", time_unit="s", fluid_type="Water"):
        
        self.date_time = date_time
        
        self.fluid_type = fluid_type
        self.D = D
        self.L = L
        self.v = v
        self.e = e
        self.g = g
        
        self.t0 = t0
        self.z0 = z0
        self.V0 = V0
        self.h = h
        self.t = t
        
        self.err_tolerance = err_tolerance
        self.refined_tolerance = refined_tolerance
        self.D1 = D1
        self.D2 = D2
        self.length_minor_losses = length_minor_losses
        self.f = f
        
        self.pipe_label = pipe_label
        self.reservoir1_label = reservoir1_label
        self.reservoir2_label = reservoir2_label
        
        self.pipe_table = pipe_table
        self.reservoir1_table = reservoir1_table
        self.reservoir2_table = reservoir2_table
        
        self.pipe_img_src = pipe_img_src
        self.reservoir1_img_src = reservoir1_img_src
        self.reservoir2_img_src = reservoir2_img_src
        
        self.dimen = dimen_unit
        self.time_unit = time_unit
    
    def get_start_html(self):
        html = "<!DOCTYPE html>\n"
        html += "<html lang='en'>\n"

        html += "<head>\n"
        
        html += "<meta charset='UTF-8'>\n"
        html += "<meta name='viewport' content='width=device-width, initial-scale=1.0'>\n"
        html += "<meta http-equiv='X-UA-Compatible' content='ie=edge'>\n"
        html += "<title>Flow Oscillations in Pipes Design</title>\n"
        html += "<link rel='stylesheet' href='../../Designs/style/flow.css'>\n"
        html += "<link rel='stylesheet' href='../../style/flow.css'>\n"
        html += "<script src='https://code.jquery.com/jquery-3.3.1.min.js'></script>\n"

        html += "</head>\n"
        html += "<body>\n"
        return html
        
    def get_page_header(self):
        html = "<header class='header-cont'>\n"
        html += "<div class='header center-text'>\n"
        html += "Flow Oscillations in Pipes Design\n"
        html += "</div>\n"
        html += f"<div class='date center-text'>TABLES AND GRAPHS</div>\n"
        html += f"<div class='date center-text'>John Oyegbite || {self.date_time}</div>\n"
        html += "</header>\n"
        return html
        
    def get_aside_section(self):
        # for the navigation drawer
        html = "<aside class='navigation-drawer'>\n"
        html += "<div class='flow-parameters-header'>Flow Parameters</div>\n"
        html += "<div class='navigation-head'>Pipe properties:</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'> D = {self.D} </div>\n"
        html += f"<div class='unit'>Diameter of pipe ({self.dimen})</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'> L = {self.L}</div>\n"
        html += f"<div class='unit'>Length of pipe ({self.dimen})</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>v = {self.v}</div>\n"
        html += f"<div class='unit'>Kinematic viscosity "
        html += f"({self.dimen}<sup>2</sup>/{self.time_unit})</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>e = {self.e}</div>\n"
        html += f"<div class='unit'>Roughness height ({self.dimen})</div>\n"
        html += "</div>\n"
        html += "<hr>\n"
        html += "<div class='navigation-head'>Initial conditions:</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>t0 = {self.t0}</div>\n"
        html += f"<div class='unit'>Initial time ({self.time_unit})</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>z0 = {self.z0}</div>\n"
        html += f"<div class='unit'>Initial head ({self.dimen})</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>V0 = {self.V0} </div>\n"
        html += f"<div class='unit'>Initial velocity ({self.dimen}/{self.time_unit})</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>t = {self.t} </div>\n"
        html += f"<div class='unit'>Final time ({self.time_unit})</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>h = {self.h} </div>\n"
        html += f"<div class='unit'>Initial time step ({self.time_unit})</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>g = {self.g}</div>\n"
        html += f"<div class='unit'>Acceleration due to gravity "
        html += f"({self.dimen}/{self.time_unit}<sup>2</sup>)</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>Error tolerance = {self.err_tolerance}</div>\n"
        html += "<div class='unit'>Error tolerance</div>\n"
        html += "</div>\n"
        html += "<div class='value-unit'>\n"
        html += f"<div class='value'>Refined tolerance = {self.refined_tolerance}</div>\n"
        html += "<div class='unit'>Refined tolerance</div>\n"
        html += "</div>\n"
        if self.D1 and self.D2:
            html += "<hr>\n"
            html += "<div class='navigation-head'>Reservoir parameter:</div>\n"
            html += "<div class='value-unit'>\n"
            html += f"<div class='value'>D1 = {self.D1} </div>\n"
            html += f"<div class='unit'>Diameter of first reservoir ({self.dimen})</div>\n"
            html += "</div>\n"
            html += "<div class='value-unit'>\n"
            html += f"<div class='value'>D2 = {self.D2}</div>\n"
            html += f"<div class='unit'>Diameter of second reservoir ({self.dimen})</div>\n"
            html += "</div>\n"
            if self.f:
                html += "<div class='value-unit'>\n"
                html += f"<div class='value'>f = {self.f} </div>\n"
                html += "<div class='unit'>Constant friction factor</div>\n"
                html += "</div>\n"
            html += "<div class='value-unit'>\n"
            html += f"<div class='value'>Length of minor losses = {self.length_minor_losses} </div>\n"
            html += f"<div class='unit'>Length of minor losses ({self.dimen})</div>\n"
            html += "</div>\n"
        html += "<div class='clear'></div>\n"
        html += "<div class='clear'></div>\n"
        html += "</aside>\n"
        return html
        
    def get_main_section(self, sn, header, table, img_src):
        html = "<div class='flow-container'>\n"
        html += f"<div class='calculation-header'>{header} calculations</div>\n"
        html += "<div class='table-header-cont'>\n"
        html += "<div class='table-header'>\n"
        html += f"<span>{header} Table</span>\n"
        html += f"<p>(error tolerance: {self.err_tolerance}, "
        html += f"refined error tolerance: {self.refined_tolerance})</p>\n"
        html += "</div>\n"
        html += "<div class='table-body'>\n"
        html += f"{table}\n"
        html += "</div>\n"
        html += "</div>\n"
        html += "<div class='graph-header-cont'>\n"
        html += f"<div class='graph-header'>{header} Graph</div>\n"
        html += "<div class='graph-body'>\n"
        html += "<img "
        html += f"src='{img_src}' " 
        html += f"alt='Graph of {header}-Head against Time'>\n"
        html += "<p class='center-text caption'>"
        html += f"Graph of ({self.fluid_type}) Head against Time in {header}"
        html += "</p>\n"
        html += "</div>\n"
        html += "</div>\n"
        html += "</div>\n"
        html += "<div class='clear'></div>\n"
        return html
    
    def get_end_html(self):
        html = "</body>\n"
        html += "</html>\n"
        return html
    
    def get_html(self):
        html = self.get_start_html()
        html += self.get_page_header()
        html += self.get_aside_section()
        html += "<div class='container'>\n"
        
        sections = {}
        if self.D1 and self.D2:
            sections = {
                    self.pipe_label: [self.pipe_table, self.pipe_img_src],
                    self.reservoir1_label: [self.reservoir1_table, self.reservoir1_img_src],
                    self.reservoir2_label: [self.reservoir2_table, self.reservoir2_img_src]
                    }
        else:
            sections = {self.pipe_label: [self.pipe_table, self.pipe_img_src]}
        
        sn = 1
        for section_header in sections:
            table_index = 0
            img_src_index = 1
            table = sections[section_header][table_index]
            img_src = sections[section_header][img_src_index]
            html += self.get_main_section(sn, section_header, table, img_src)
            sn += 1
        html += "</div>"
        html += self.get_end_html()
        return html
 

    

