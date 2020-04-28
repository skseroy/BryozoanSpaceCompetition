# Sasha Seroy and Danny Grunbaum, 2018
# School of Oceanography, Friday Harbor Labs
# University of Washington
# Seattle, WA 98105


# ==========================================================================================
# A SPATIALLY-EXPLICIT SIMULATION TO MODEL BRYOZAON (MEMBRANIPORA MEMBRANACEA) COLONY GROWTH
# Originally coded for Python 2.7

# -------------------------------------------
# Instructions:
# -------------------------------------------

# 1. Designate flags for plotting
#       - Saves parameters.txt and exports data as col_data.csv (area and energy content by colony)
#           and EMD_data.csv (energy, mass, developmental state by zooid within colony)
# 2. Designate zooid geometry parameters to describe the desired zooid size
# 3. Set habitat parameters including domain (kelp blade) size and number of initial settling colonies
# 4. Set desired time period for the duration with start day, end day and time step interval
# 5. Set settlement parameters if you want multiple cohorts to settle after inital settlement
#       - Note: all cohorts currently require the same number of settlers
# 6. Set energy parameters - broken down by which parameters are included in each term (i.e qf, qin etc.)
# 7. Colonies can be given two different growth rates coeeficients, 'cont' and 'pred' (e.g. undefended or defended)
# 8. Colony growth rates are not direct inputs but coefficients (S_c) to modify qf. Some experimenting required to find desired S_cs



# ---------------------------------------------------------------
# Current parameters:
# ---------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.path as mplPath
import time
import random

from math import *
from random import *
from matplotlib import animation
from scipy.integrate import odeint
from sympy.geometry import Polygon, Point
from functools import wraps
from datetime import datetime


# ===========================================================================
# Flags for plotting
# ===========================================================================

t = datetime.now()
time_created = str(t.month)+'-'+str(t.day)+'-'+str(t.hour)+'-'+str(t.minute)+'-'+str(t.second)

interactive_plot = True             # True turns on interactive plotting
plot_nums = False                   # True turns on zooid numbering
save_final_fig = False              # True turns on saving final time step figure
save_movie_frames = False           # True turns on saving all time step frames
save_data = False                   # True turns on saving parameters.txt, col_data and pop_data files
data_name = 'B1_plate_rep_'+time_created

# ===========================================================================
# Zooid geometry parameters
# ===========================================================================

# biologically accurate p (summer 2017 data, no pred pH 7.9) to get average zooid area in mm2 (0.2972 mm2)
p = 0.54041
theta = pi / 6
x1 = p * cos(theta / 2)
y1 = p * sin(theta / 2)

# ===========================================================================
# Habitat parameters
# ===========================================================================

dom = [100, 50]         # habitat dimensions (mm2 x mm2)

# ===========================================================================
# Time step conditions
# ===========================================================================

# can be edited to reflect actual time in the field
start_day=0.            # start time of simulation

end_day=81.             # end time of simulation

time_steps_per_day = 2. # number of time steps per day (i.e 2. = each time step is 1/2 day)

num_intervals = time_steps_per_day*(end_day + 1.)  # number of time steps during simulation

# ===========================================================================
# Settlement parameters - designate settlement time and position
# ===========================================================================

# List settlement times and positions in the same order

# Settlement times of colonies in chronological order
settlement_days = [2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 9, 11, 14, 23]

settlement_times = []
for i in range(len(settlement_days)):
    a = settlement_days[i]*time_steps_per_day + 1
    settlement_times.append(a)

# Enter settlement positions in order to correspond with settlement days
settlement_positions = [[68.8, 47.0],
			            [73.3, 13.4],
			            [78.7, 37.0],
			            [65.9, 39.6],
			            [86.4, 37.3],
			            [50.5, 15.6],
			            [35.9, 34.5],
			            [39.6, 3.0],
			            [3.0, 9.2],
			            [82.4, 47.0],
			            [17.8, 17.2],
			            [8.2, 12.7],
			            [14.6, 47.0],
			            [35.0, 3.0],
			            [50.3, 47.0],
			            [50.3, 38.5]]
'''
Note, this version of the code allows for colonies to plotted on top of each other
if the settlement positions are too close to each other. Make sure to input positions
such that colonies will not overlap.
'''
# ===========================================================================
# Energy input from feeding (qf) parameters
# ===========================================================================

alg = 50000.            # algal density for feeding (cells/ml)

cl = 22.08*0.7          # zooid clearance rate (ml/day)

F = 1.15e-6             # food quality of algae (Joules / cell)

start_dev = 0.1         # zooids start with this initial dev state (no units)

dev_rate = 0.3          # rate of development of any given zooid (1/day)

Dev_max = 1.0           # maximum developmental state - adult zooid (no units)

# all colonies recieve the same growth rate
gro = 1.0222441     # growth rate coefficient - value from delta q graph

# make a list of growth rate coefficients for all simulation colonies
S_cs = np.repeat(gro, len(settlement_days))

# ===========================================================================
# Energy output to growth (qg) parameters
# ===========================================================================

start_mass = 1.0e-8     # starting mass of a zooid (gC)

M_max = 1.96e-7         # maximum mass of a zooid (gC)

Q_m = 4.184e4           # mass-energy conversion (Joules / gC)

g = 0.5                 # growth rate of tissue addition (1/day)


# ===========================================================================
# Energy output to metabolic rate (qm) parameters
# ===========================================================================

q_0 = 0.05              # basal metabolic rate (Joules/day)


# ===========================================================================
# Energy translocation (qin/qout) parameters
# ===========================================================================

# Note: o_a + o_l must = 0.5

o_a = 0.25 * 2.         # axial energy translocation coeffiecient (none)
o_a_rate = 1.           # axial energy translocation rate (1/day)

o_l = 0.25 * 6.         # lateral translocation coefficient - (none)
o_l_rate = 1.           # lateral energy translocation rate (1/day)

start_E = 0.001         # zooid starting energy (Joules)

E_r = 7.                # threshold energy below which zooid cannot send energy downstream

E_r_smoothing_param=0.1


# ===========================================================================
# Store parameters in an associated .txt file
# ===========================================================================

if save_data == True:
    f = open('parameters_' + data_name + '.txt', 'w') #Creates a txt file to store parameter values for the run

    print >>f, "\ninit_set", init_set, "\ndom", dom  # send param info to txt file
    print >>f, "\nstart_day", start_day, '\nend_day', end_day, '\nnum_intervals', num_intervals

    if len(settlement_times) != 0:
        print >>f, "\ncohort_settlers", cohort_settlers, "\nsettlement_times", settlement_times

    print >>f, "\nalg", alg, "\ncl", cl, "\nF", F, "\nstart_dev", start_dev, "\ndev_rate", dev_rate, "\nDev_max", Dev_max

    if len(S_cs) != 0:
        print >>f, "\nS_cs", S_cs
    else:
        print >>f, '\nS_c', S_cs

    if len(S_cs2) != 0:
        print >>f, "\nS_cs2", S_cs2
    else:
        print >>f, '\nS_c', S_cs

    print >>f, "\nstart_mass", start_mass, "\nM_max", M_max, "\nQ_m", Q_m, "\ng", g, "\nstart_E", start_E
    print >>f, "\nq_0", q_0
    print >>f, "\no_a", o_a, "\no_a_rate", o_a_rate, "\no_l", o_l, "\no_l_rate", o_l_rate, "\nE_r", E_r, "\nE_r_smoothing_param", E_r_smoothing_param
    f.close()

# ===========================================================================
# Additional parameters
# ===========================================================================

rnd_digits = 5            # round values to 5 decimals

# ===================================================================================================================
# Define a zooid class
# ===================================================================================================================

class Zooid():

    # helper functions to calculate position of zooids
    x1 = p * cos(theta / 2)
    y1 = p * sin(theta / 2)

    # define a zooid with the following characteristics
    def __init__(self, p=p, theta=theta, ID=0, orientation=0, midPoint=[0., 0.],
                 vertices=[(-x1 - p / 2, 0), (-p / 2, y1), (p / 2, y1), (x1 + p / 2, 0), (p / 2, -y1), (-p / 2, -y1)],
                 orig_verts = [(-x1 - p / 2, 0), (-p / 2, y1), (p / 2, y1), (x1 + p / 2, 0), (p / 2, -y1), (-p / 2, -y1)],
                 #shape = Polygon((-x1 - p / 2, 0), (-p / 2, y1), (p / 2, y1), (x1 + p / 2, 0), (p / 2, -y1), (-p / 2, -y1)),
                 ax_mom1 = None, ax_mom2 = None, ax_dtr1 = None, ax_dtr2 = None, lat_neighb1 = None, lat_neighb2 = None, potential_neighbs = [],
                 ztype=None, E_tot=0.0, deltaE_tot = 0.0, area=(4 * (0.5 * x1 * y1) + (p * 2 * y1)),
                 spines=0, mass=start_mass, dev_state=start_dev, q_in=0.0, q_out=0.0, q_f=0.0, q_g=0.0, q_m=q_0, dEdt=0.0, neighb_area = 0.0,
                 open_edge_dtr1 = True, open_edge_dtr2 = True, axial_dtr1_z = None, axial_dtr2_z = None, col_ID = None, col_cond = 'c'):
        self.p = p
        self.theta = theta
        self.ID = ID
        self.orientation = orientation
        self.midPoint = midPoint
        self.vertices = vertices
        self.ax_mom1 = ax_mom1
        self.ax_mom2 = ax_mom2
        self.ax_dtr1 = ax_dtr1
        self.ax_dtr2 = ax_dtr2
        self.lat_neighb1 = lat_neighb1
        self.lat_neighb2 = lat_neighb2
        self.potential_neighbs = potential_neighbs
        self.ztype = ztype
        self.E_tot = E_tot
        self.deltaE_tot = deltaE_tot
        self.area = area
        self.spines = spines
        self.mass = mass
        self.dev_state = dev_state
        self.q_in = q_in
        self.q_f = q_f
        self.q_out = q_out
        self.q_g = q_g
        self.q_m = q_m
        self.dEdt = dEdt
        self.neighb_area = neighb_area
        self.axial_dtr1_z = axial_dtr1_z # actual daughter zooid instances stored
        self.axial_dtr2_z = axial_dtr2_z # actual daughter zooid instances stored
        self.open_edge_dtr1 = open_edge_dtr1
        self.open_edge_dtr2 = open_edge_dtr2
        self.round_mp = [round(self.midPoint[0], rnd_digits), round(self.midPoint[1], rnd_digits)]
        self.col_ID = col_ID
        self.orig_verts = orig_verts
        self.col_cond = col_cond

    # Plot zooids based of colony growth rate and energy available
    def Plot(self, ID_z=''):
        if self.col_cond == 'c':
            if self.E_tot > E_r:
                color = 'cornflowerblue'
            else:
                color = 'lightblue'
        else:
            if self.E_tot > E_r:
                color = 'mediumvioletred'
            else:
                color = 'plum'

        plt_verts=[[self.vertices[0,ipv],self.vertices[1,ipv]] for ipv in range(6)]
        zooid = patches.Polygon(plt_verts, color=color, alpha=0.50)
        ax.add_patch(zooid)

    # locates and creates on-axis zooids, adds them to zooid attribute on-axis daughter list

    def Get_axial_dtr1_midPoint(self, ID=''):
        x1 = p * cos(theta / 2)
        y1 = p * sin(theta / 2)
        theta_0 = self.orientation
        
        # make on-axis daughter 1
        wrk1=(sqrt((x1 + p) ** 2 + y1 ** 2))
        wrk2=theta_0 + atan(y1 / (x1 + p))
        xD1 = round(self.midPoint[0] + wrk1 * cos(wrk2), rnd_digits)
        yD1 = round(self.midPoint[1] + wrk1 * sin(wrk2), rnd_digits)
        return [xD1, yD1]

    def Get_axial_dtr2_midPoint(self, ID=''):
        x1 = p * cos(theta / 2)
        y1 = p * sin(theta / 2)
        theta_0 = self.orientation

        # make on-axis daughter 2
        wrk1=sqrt((x1 + p) ** 2 + y1 ** 2)
        wrk2=theta_0 - atan(y1 / (x1 + p))
        xD2 = round(self.midPoint[0] + wrk1 * cos(wrk2), rnd_digits)
        yD2 = round(self.midPoint[1] + wrk1 * sin(wrk2), rnd_digits)
        return [xD2, yD2]

    def Get_axial_dtr1(self, ID=''):
        x1 = p * cos(theta / 2)
        y1 = p * sin(theta / 2)
        theta_0 = self.orientation
        
        # make on-axis daughter 1
        wrk1=(sqrt((x1 + p) ** 2 + y1 ** 2))
        wrk2=theta_0 + atan(y1 / (x1 + p))
        xD1 = self.midPoint[0] + wrk1 * cos(wrk2)
        yD1 = self.midPoint[1] + wrk1 * sin(wrk2)

        # make rotation matrix for daughter 2
        c, s = np.cos(theta_0), np.sin(theta_0)
        rotate = np.matrix([[c, -s], [s, c]])
        mp1 = np.matrix([[xD1],[yD1]])
        orig_verts = [(-x1 - p / 2, 0), (-p / 2, y1), (p / 2, y1), (x1 + p / 2, 0), (p / 2, -y1), (-p / 2, -y1)]
        vertices = np.matrix([[v[0] for v in orig_verts],[v[1] for v in orig_verts]])
        
        # create daughter 1
        daughter1 = Zooid(midPoint=[xD1, yD1], orientation=self.orientation, ztype='dtr', E_tot = start_E, vertices = (rotate*vertices) + mp1)
        self.axial_dtr1_z = daughter1


    def Get_axial_dtr2(self, ID=''):
        x1 = p * cos(theta / 2)
        y1 = p * sin(theta / 2)
        theta_0 = self.orientation

        # make on-axis daughter 2
        wrk1=sqrt((x1 + p) ** 2 + y1 ** 2)
        wrk2=theta_0 - atan(y1 / (x1 + p))
        xD2 = self.midPoint[0] + wrk1 * cos(wrk2)
        yD2 = self.midPoint[1] + wrk1 * sin(wrk2)
        
        # make rotation matrix for daughter 2
        c, s = np.cos(theta_0), np.sin(theta_0)
        rotate = np.matrix([[c, -s], [s, c]])
        mp2 = np.matrix([[xD2],[yD2]])
        orig_verts = [(-x1 - p / 2, 0), (-p / 2, y1), (p / 2, y1), (x1 + p / 2, 0), (p / 2, -y1), (-p / 2, -y1)]
        vertices = np.matrix([[v[0] for v in orig_verts],[v[1] for v in orig_verts]])
        
        # create daughter 2
        daughter2 = Zooid(midPoint=[xD2, yD2], orientation=self.orientation, ztype='dtr', E_tot = start_E, vertices = (rotate*vertices) + mp2)
        self.axial_dtr2_z = daughter2

    def oldGet_axial_dtr1(self, ID=''):
        x1 = p * cos(theta / 2)
        y1 = p * sin(theta / 2)
        theta_0 = self.orientation

        # make on-axis daughter 1
        xD1 = self.midPoint[0] + (sqrt((x1 + p) ** 2 + y1 ** 2)) * cos(theta_0 + atan(y1 / (x1 + p)))
        yD1 = self.midPoint[1] + (sqrt((x1 + p) ** 2 + y1 ** 2)) * sin(theta_0 + atan(y1 / (x1 + p)))
        mp1 = [xD1, yD1]
        daughter1 = Zooid(midPoint=[xD1, yD1], orientation=self.orientation, ztype='dtr', E_tot = start_E)
        rt = daughter1.shape.rotate(daughter1.orientation).translate(daughter1.midPoint[0], daughter1.midPoint[1])
        daughter1.shape = rt

        self.axial_dtr1_z = daughter1

    def Get_neighbor_space(self, ID=''):
        area = mplPath.Path.circle(center=(self.midPoint[0], self.midPoint[1]), radius=3 * p)
        self.neighb_area = area


    def Point_within_radius(self,P, ID='',radius=3 * p):
        d2=(self.midPoint[0]-P[0])**2+(self.midPoint[1]-P[1])**2
        if d2<=radius**2:
            return True
        else:
            return False
        

# ===================================================================================================================
# Define a colony class
# ===================================================================================================================

class Colony():

    # define a colony class with the following characteristics
    def __init__(self, ancestrula=int(2 * pi / theta), theta=theta, ID=0, midPoint=[0., 0.],
                 zooid_list=[],new_zooids_to_plot= [], perimeter_zooids = [], colony_area=0., colony_energy=0., Es=[], Devs = [], Masses =[],
                 dEs_dt=[], dDs_dt = [], dMs_dt = [], dAs_dt=[],
                 potential_neighb_cols = [], spine_coeff = 1.0, spine_cond = 'c'):
        self.ancestrula = int(2 * pi / theta)
        self.theta = theta
        self.ID = ID
        self.midPoint = midPoint
        self.round_mp = [round(self.midPoint[0], rnd_digits), round(self.midPoint[1], rnd_digits)]
        self.zooid_list = []
        self.perimeter_zooids = []
        self.colony_area = colony_area
        self.colony_energy = colony_energy
        self.Es = Es
        self.Devs = Devs
        self.Masses = Masses
        self.dEs_dt = dEs_dt
        self.dDs_dt = dDs_dt
        self.dMs_dt = dMs_dt
        self.dAs_dt = dAs_dt
        self.new_zooids_to_plot= []
        self.potential_neighb_cols = []
        self.potential_neighb_cols_IDs = []
        self.x_min = midPoint[0]
        self.y_min = midPoint[1]
        self.x_max = midPoint[0]
        self.y_max = midPoint[1]
        self.spine_coeff = spine_coeff
        self.spine_cond = spine_cond
        


    # locates and creates the zooids for the first ancestrula ring (zooids 1- 12...if theta = pi/6)
    def Get_ancestrula(self, ID=''):
        for i in range(self.ancestrula):
            x1 = p * cos(self.theta / 2)
            y1 = p * sin(self.theta / 2)
            x_in = self.midPoint[0] + ((x1 + p / 2) * cos(self.theta * (i)))
            y_in = self.midPoint[1] + ((x1 + p / 2) * sin(self.theta * (i)))
            mp1 = np.matrix([[x_in],[y_in]])
            c, s = np.cos(self.theta*i), np.sin(self.theta*i)

            rotate = np.matrix([[c, -s], [s, c]])
            orig_verts = [(-x1 - p / 2, 0), (-p / 2, y1), (p / 2, y1), (x1 + p / 2, 0), (p / 2, -y1), (-p / 2, -y1)]
            vertices = np.matrix([[v[0] for v in orig_verts],[v[1] for v in orig_verts]])
                    
            new_zooid1 = Zooid(ID=i, midPoint= [x_in,y_in], col_ID = self.ID, orientation=theta * i, ztype='anc', E_tot=5.0,
                               vertices = (rotate*vertices) + mp1, col_cond = self.spine_cond)
            self.zooid_list.append(new_zooid1)


    # locates and creates the zooids for the second ancestrula ring (zooids 13-24...if theta = pi/6)

        for i in range(self.ancestrula):
            x1 = p * cos(self.theta / 2)
            y1 = p * sin(self.theta / 2)
            x_out = self.midPoint[0] + ((x1 + 3 * p / 2) * cos(self.theta * (i) + self.theta / 2))
            y_out = self.midPoint[1] + ((x1 + 3 * p / 2) * sin(self.theta * (i) + self.theta / 2))
            mp2 = np.matrix([[x_out],[y_out]])
            c, s = np.cos(self.theta*(i+0.5)), np.sin(self.theta*(i+0.5))

            rotate = np.matrix([[c, -s], [s, c]])
            orig_verts = [(-x1 - p / 2, 0), (-p / 2, y1), (p / 2, y1), (x1 + p / 2, 0), (p / 2, -y1), (-p / 2, -y1)]
            vertices = np.matrix([[v[0] for v in orig_verts],[v[1] for v in orig_verts]])
            
            new_zooid2 = Zooid(ID=i + self.ancestrula, midPoint=[x_out, y_out], col_ID = self.ID, orientation=theta * (i) + theta / 2, ztype='anc', E_tot=5.0,
                              vertices = (rotate*vertices) + mp2, col_cond = self.spine_cond)
            self.zooid_list.append(new_zooid2)
            

        # Update the max/min geometrical characteristics for the colony
        # assign colony max/min values with maximum buffer to account for using midpoints instead of vertices
        for iz in range(len(self.zooid_list)):
            self.x_max = max(self.x_max,self.zooid_list[iz].midPoint[0]+p*.5)
            self.x_min = min(self.x_min,self.zooid_list[iz].midPoint[0]-p*.5)
            self.y_max = max(self.y_max,self.zooid_list[iz].midPoint[1]+p*.5)
            self.y_min = min(self.y_min,self.zooid_list[iz].midPoint[1]-p*.5)


    # plots zooids of type ancestrula
    def Plot_ancestrula(self, ID=''):
        for zooid in self.zooid_list:
            if zooid.ztype == 'anc':
                zooid.Plot()
                if plot_nums == True:
                    ax.text(zooid.midPoint[0], zooid.midPoint[1], zooid.ID)

    def Find_perimeter_zooids(self, ID=''): #new version, uses if dtr neighbor is designated or not
        for zooid in self.zooid_list:
            if zooid.ax_dtr1 == None or zooid.ax_dtr2 == None:
                self.perimeter_zooids.append(zooid)

    def Clear_perimeter_zooids(self, ID=''):
        del self.perimeter_zooids[:]

    def Clear_new_zooids_to_plot(self, ID=''):
        del self.new_zooids_to_plot[:]

    def Get_colony_area(self, ID=''):
        self.colony_area = 0.
        for zooid in self.zooid_list:
            self.colony_area = self.colony_area + zooid.area

    def Get_colony_energy(self, ID=''):
        self.colony_energy = 0.
        for zooid in self.zooid_list:
            self.colony_energy = self.colony_energy + zooid.E_tot


    def Find_neighbor_space(self, ID=''):
        for z in self.zooid_list:
            area = mplPath.Path.circle(center=(z.midPoint[0], z.midPoint[1]), radius=3 * p)
            z.neighb_area = area

    # locate potential neighbors within a colony
    def Locate_potential_neighbors(self, ID=''):
        for i in range(len(self.zooid_list)):
            potential_neighbs = []
            for j in range(len(self.zooid_list)):
                if self.zooid_list[i].ID != self.zooid_list[j].ID:
                    if self.zooid_list[i].Point_within_radius(self.zooid_list[j].midPoint):
                        potential_neighbs.append(self.zooid_list[j])
            self.zooid_list[i].potential_neighbs = potential_neighbs


    def Assign_neighbors(self, ID=''):
        for i in range(len(self.zooid_list)):
            i_ro_verts = [[round(self.zooid_list[i].vertices[0,ii],rnd_digits),round(self.zooid_list[i].vertices[1,ii],rnd_digits)] for ii in range(6) ]
            for p_neighbor in self.zooid_list[i].potential_neighbs:
                j_ro_verts = [[round(p_neighbor.vertices[0,ii],rnd_digits),round(p_neighbor.vertices[1,ii],rnd_digits)] for ii in range(6) ]

                # test potential neighbors for AXIAL MOM 1
                if self.zooid_list[i].ztype == 'dtr' and i_ro_verts[0] in j_ro_verts and i_ro_verts[1] in j_ro_verts:
                    self.zooid_list[i].ax_mom1 = p_neighbor.ID

                # test potential neighbors for AXIAL MOM 2
                elif self.zooid_list[i].ztype == 'dtr' and i_ro_verts[0] in j_ro_verts and i_ro_verts[5] in j_ro_verts:  # ans. zooids do not have mom
                    self.zooid_list[i].ax_mom2 = p_neighbor.ID

                # test potential neighbors for LATERAL NEIGHBOR 1
                elif i_ro_verts[1] in j_ro_verts and i_ro_verts[2] in j_ro_verts:
                    self.zooid_list[i].lat_neighb1 = p_neighbor.ID

                # test potential neighbors for LATERAL NEIGHBOR 2
                elif i_ro_verts[4] in j_ro_verts and i_ro_verts[5] in j_ro_verts:
                    self.zooid_list[i].lat_neighb2 = p_neighbor.ID

                # test potential neighbors for AXIAL DAUGHTER 1
                elif i_ro_verts[2] in j_ro_verts and i_ro_verts[3] in j_ro_verts:
                    self.zooid_list[i].ax_dtr1 = p_neighbor.ID
                    self.zooid_list[i].open_edge_dtr1=False

                # test potential neighbors for AXIAL DAUGHTER 2
                elif i_ro_verts[3] in j_ro_verts and i_ro_verts[4] in j_ro_verts:
                    self.zooid_list[i].ax_dtr2 = p_neighbor.ID
                    self.zooid_list[i].open_edge_dtr2=False



# ===================================================================================================================
# Define a habitat class
# ===================================================================================================================

class Habitat():

    # define a habitat class to contain the following characteristics
    def __init__(self, ID=1, domain=[50, 100],
        initial_settlement=10, colony_mps = [], colony_list=[], pop_area = 0.0, settlement_attempts = 100):
        self.ID = ID
        self.domain = domain
        self.shape = Polygon((0,0), (self.domain[0], 0), (self.domain[0], self.domain[1]), (0, self.domain[1]))
        self.colony_list = []
        self.nudi_list = []
        self.pop_area = pop_area
        self.colony_mps = []
        self.settlement_attempts= settlement_attempts

    def Settle(self, x, y):
        new_colony = Colony(ID= len(self.colony_list), midPoint=[x, y])
        new_colony.Get_ancestrula()

        if interactive_plot == True:
            new_colony.Plot_ancestrula()
        new_colony.Es = np.zeros([len(new_colony.zooid_list)])
        new_colony.Devs = np.zeros([len(new_colony.zooid_list)])
        new_colony.Masses = np.zeros([len(new_colony.zooid_list)])

        for z in new_colony.zooid_list:
            idz = z.ID
            new_colony.Es[idz] = z.E_tot
            new_colony.Devs[idz] = z.dev_state
            new_colony.Masses[idz] = z.mass
        self.colony_list.append(new_colony)
        self.colony_mps.append([x, y])


    def Plot_settlement(self, ID=''):
        for colony in self.colony_list:
            colony.Plot_ancestrula()

    def Plot_all(self, ID=''):
        for col in self.colony_list:
            for z in col.zooid_list:
                z.Plot()
                if plot_nums == True:
                    ax.text(z2.midPoint[0], z2.midPoint[1], z2.ID)

    def Get_pop_area(self, ID=''):
        self.pop_area = 0
        for colony in self.colony_list:
            colony.Get_colony_area()
            self.pop_area = self.pop_area + colony.colony_area


# ===================================================================================================================
# Energy function to be integrated:
# ===================================================================================================================


def Energy_func(A,t,j):

    n = len(A)
    dAdts = np.zeros([n]) # derivatives of all values in A

    n1 = n/3
    dEdts = np.zeros([n1])
    dMdts = np.zeros([n1])
    dDdts = np.zeros([n1])

    tmpEs=A[0 : len(A)/3] # new energy values are the first
    tmpMasses = A[len(A)/3 : 2*len(A)/3]
    tmpDevs = A[2*len(A)/3 : len(A)]
    qins = np.zeros([n1])
    qouts = np.zeros([n1])
    qfs = np.zeros([n1])
    qgs = np.zeros([n1])

    col = kelp.colony_list[j]

    for iz in range(n1):
        if tmpEs[iz] - E_r > 0.:
            q = (tmpEs[iz] - E_r) / (tmpEs[iz] - E_r + E_r_smoothing_param)

            # =================================================
            # Calculate FLUX OUT (q out) through neighbor pores
            # zooids own daughter pores, any pores with flux out
            # check to make sure not negative

            # axial daughter pore 1:
            if col.zooid_list[iz].ax_dtr1 is not None:
                axd1 = col.zooid_list[iz].ax_dtr1  # ID no. of axial mom 1
                aflux_out1 = o_a * q * max(0., tmpEs[iz] - E_r)
                qouts[iz] += aflux_out1  # add to qout of iz
                qins[axd1] += aflux_out1  # add to qin of dtr

            # axial daughter pore 2:
            if col.zooid_list[iz].ax_dtr2 is not None:
                axd2 = col.zooid_list[iz].ax_dtr2  # ID no. of axial mom 2
                aflux_out2 = o_a * q * max(0., tmpEs[iz] - E_r)
                qouts[iz] += aflux_out2  # add to qout of iz
                qins[axd2] += aflux_out2  # add to qin of dtr


            # lateral neighbor pore 1:
            if col.zooid_list[iz].lat_neighb1 is not None:
                ltnb1 = col.zooid_list[iz].lat_neighb1
                lflux_out1 = o_l * q * max(0., tmpEs[iz] - E_r)
                qouts[iz] += lflux_out1  # add to qout of iz
                qins[ltnb1] += lflux_out1  # add to qin of neighbor


            # lateral neighbor pore 2:
            if col.zooid_list[iz].lat_neighb2 is not None:
                ltnb2 = col.zooid_list[iz].lat_neighb2  # ID no. of lat neighbor 2
                lflux_out2 = o_l * q * max(0., tmpEs[iz] - E_r)
                qouts[iz] += lflux_out2  # add to qout of iz
                qins[ltnb2] += lflux_out2  # add to qin of neighbor


    for iz in range(n1):

        # =================================================
        # Calculate energy consumed by growth
        qgs[iz] = g * Q_m * (M_max - tmpMasses[iz])

        # Calculate energy obtained by feeding
        qfs[iz] = col.spine_coeff * alg * cl * F * tmpDevs[iz]  # col.spine_coeff is a constant
        qm = q_0

        if tmpEs[iz] > 0.:  # if there is some energy available proceed
            dEdts[iz] = qins[iz] + qfs[iz] - qouts[iz] - qgs[iz] - qm
        else:  # if there is no energy available only feed and recieve E from neighbors
            dEdts[iz] = qins[iz] + qfs[iz]

    # calculate ODE for updating mass
    for iz in range(n1):
        if tmpMasses[iz] < M_max and tmpEs[iz] > 0.:  # mass only advances if the zooid is below the adult zooid mass
            dMdts[iz] = qgs[iz] / Q_m
        else:
            dMdts[iz] = 0.

    for iz in range(n1):
        if qfs[iz] > 0.:  # If the colony is feeding, advance dev_state
            if tmpDevs[iz] < Dev_max and tmpEs[iz] > 0.:  # dev_state only advances if the zooid is still developing
                dDdts[iz] = dev_rate
            else:
                dDdts[iz] = 0.  # if dev_state is already at the maximum,  dev_state no longer changes

    #update ind. attribute with new values
    for iz in range(n1):
        #energy flow values
        col.zooid_list[iz].q_in = qins[iz]
        col.zooid_list[iz].q_out = qouts[iz]
        col.zooid_list[iz].q_f = qfs[iz]
        col.zooid_list[iz].q_g = qgs[iz]

        col.zooid_list[iz].dEdt = dEdts[iz]
        col.zooid_list[iz].dMdt = dMdts[iz]
        col.zooid_list[iz].dDdt = dDdts[iz]

        col.zooid_list[iz].E_tot = tmpEs[iz]
        col.zooid_list[iz].mass = tmpMasses[iz]
        col.zooid_list[iz].dev_state = tmpDevs[iz]

    dAdts = np.concatenate((dEdts, dMdts, dDdts))

    return dAdts


# ===================================================================================================================
# Integrate in small time steps
# ===================================================================================================================

def Habitat_integrate(t0, t1):
    print('beginning time step', t1)
    for j in range(len(kelp.colony_list)):
        i = kelp.colony_list[j]
        # integrate change of state variables, Etot
        atol = 5e-3
        rtol = 1e-3
        hmax = 0.25
        hmin = 0.00005

        master_list = np.concatenate((i.Es, i.Masses, i.Devs))
        all = odeint(Energy_func, master_list, [t0, t1], args=(j,), atol=atol, rtol=rtol, hmax=hmax, hmin=hmin)

        new_array = all[-1]
        i.Es = new_array[0: len(new_array) / 3]  # new energy values are the first
        i.Masses = new_array[len(new_array) / 3: 2 * len(new_array) / 3]
        i.Devs = new_array[2 * len(new_array) / 3: len(new_array)]

        for j2 in i.zooid_list:
            j2.E_tot = i.Es[j2.ID]
            j2.dev_state = i.Devs[j2.ID]
            j2.mass = i.Masses[j2.ID]

    # establish new perimeter zooid list for each colony
    for i2 in kelp.colony_list:
        i2.Clear_perimeter_zooids()
        i2.Find_perimeter_zooids()

    # find new potential neighbor colonies
    for c1 in kelp.colony_list:
        for c2 in kelp.colony_list:
            if c1.ID != c2.ID:  # if it is not the same colony
                if c2.ID not in c1.potential_neighb_cols_IDs:  # if it is not already listed in neighb col IDs
                    # find potential overlapping colonies to search
                    if c1.ID != c2.ID and c1.x_min - 3 * p < c2.x_max and \
                                            c1.x_max + 3 * p > c2.x_min and c1.y_min - 3 * p < c2.y_max and c1.y_max + 3 * p > c2.y_min:
                        c1.potential_neighb_cols_IDs.append(c2.ID)

    # create new daughter for each colony
    col_index = range(len(kelp.colony_list))
    shuffle(col_index)  # make a new list with the shuffled colonies
    for i_shuff in col_index:
        i3 = kelp.colony_list[i_shuff]
        midPoints = []  # appears to be a list of midpoints of new axial daughter zooids
        for z in i3.perimeter_zooids:
            # does it have enough energy to make one new zooid?
            if z.open_edge_dtr1 == False or i3.Es[z.ID] < E_r + start_E:  # if it does, try to make new axial daughters
                continue
            axial_dtr1_midPoint = z.Get_axial_dtr1_midPoint()
            # Test whether too close to edge of habitat
            if axial_dtr1_midPoint[0] > kelp.domain[0] - 3. * p or axial_dtr1_midPoint[0] < 3. * p or \
                            axial_dtr1_midPoint[1] > kelp.domain[1] - 3. * p or axial_dtr1_midPoint[1] < 3. * p:
                z.open_edge_dtr1 = False  # this spot will be forever occluded
            elif axial_dtr1_midPoint in midPoints:  # check dtr1: is it a duplicate of

                z.open_edge_dtr1 = False  # this spot will be forever occluded
            else:
                for i4 in i3.potential_neighb_cols_IDs:  # only check potential neighbor colonies
                    a = kelp.colony_list[i4].perimeter_zooids + kelp.colony_list[i4].new_zooids_to_plot
                    for p1 in a:
                        if (axial_dtr1_midPoint[0] - p1.midPoint[0]) ** 2 + (
                            axial_dtr1_midPoint[1] - p1.midPoint[1]) ** 2 < 9. * p ** 2:
                            z.open_edge_dtr1 = False  # this spot will be forever occulded an not in perimeter zooids
                            break  # this breaks from the p1 loop
                    if z.open_edge_dtr1 == False:
                        break  # this breaks from the i4 loop
            if z.open_edge_dtr1 == False:
                continue  # skip to the next zooid in the z loop

            # print 'Getting axial dtr1'
            z.Get_axial_dtr1()
            dtr1 = z.axial_dtr1_z
            dtr1.Get_neighbor_space()

            i3.Es[z.ID] = i3.Es[z.ID] - start_E

            z.open_edge_dtr1 = False
            dtr1.ID = max(zooid.ID for zooid in i3.zooid_list) + 1
            dtr1.col_ID = i3.ID
            dtr1.col_cond = i3.spine_cond
            # Add mother/daughter connection
            dtr1.ax_mom2 = z.ID
            z.ax_dtr1 = dtr1.ID
            z.axial_dtr1_z = dtr1
            i3.new_zooids_to_plot.append(dtr1)
            i3.zooid_list.append(dtr1)
            i3.Es = np.concatenate([i3.Es, [dtr1.E_tot]])  # add its energy to E array
            i3.Masses = np.concatenate([i3.Masses, [dtr1.mass]])
            i3.Devs = np.concatenate([i3.Devs, [dtr1.dev_state]])
            midPoints.append(dtr1.round_mp)  #

            i3.x_max = max(i3.x_max, dtr1.midPoint[0] + p * .5)
            i3.x_min = min(i3.x_min, dtr1.midPoint[0] - p * .5)
            i3.y_max = max(i3.y_max, dtr1.midPoint[1] + p * .5)
            i3.y_min = min(i3.y_min, dtr1.midPoint[1] - p * .5)

        for z in i3.perimeter_zooids:
            # does it have enough energy to make another new zooid?
            if z.open_edge_dtr2 == False or i3.Es[z.ID] < E_r + start_E:  # if it does, try to make new axial daughters
                continue  # skip to the next zooid in the z loop
            axial_dtr2_midPoint = z.Get_axial_dtr2_midPoint()

            # Test whether too close to edge of habitat
            if axial_dtr2_midPoint[0] > kelp.domain[0] - 3. * p or axial_dtr2_midPoint[0] < 3. * p or \
                            axial_dtr2_midPoint[1] > kelp.domain[1] - 3. * p or axial_dtr2_midPoint[1] < 3. * p:
                z.open_edge_dtr2 = False  # this spot will be forever occluded
            elif axial_dtr2_midPoint in midPoints:  # check dtr2: is it a duplicate of
                z.open_edge_dtr2 = False  # this spot will be forever occluded
            else:
                # print 'checking for overlaps with neighboring colonies, zooid ',z.ID
                for i4 in i3.potential_neighb_cols_IDs:  # only check potential neighbor colonies
                    a = kelp.colony_list[i4].perimeter_zooids + kelp.colony_list[i4].new_zooids_to_plot
                    for p1 in a:
                        if (axial_dtr2_midPoint[0] - p1.midPoint[0]) ** 2 + (
                            axial_dtr2_midPoint[1] - p1.midPoint[1]) ** 2 < 9. * p ** 2:
                            z.open_edge_dtr2 = False  # this spot will be forever occulded an not in perimeter zooids
                            break  # this breaks from the p1 loop
                    if z.open_edge_dtr2 == False:
                        break  # this breaks from the i4 loop
            if z.open_edge_dtr2 == False:
                continue

            z.Get_axial_dtr2()
            dtr2 = z.axial_dtr2_z
            dtr2.Get_neighbor_space()

            # add the zooid, update mom's E
            i3.Es[z.ID] = i3.Es[z.ID] - start_E
            z.open_edge_dtr2 = False
            dtr2.ID = max(zooid.ID for zooid in i3.zooid_list) + 1
            dtr2.col_ID = i3.ID
            dtr2.col_cond = i3.spine_cond
            # Add mother/daughter connection
            dtr2.ax_mom1 = z.ID
            z.ax_dtr2 = dtr2.ID
            z.axial_dtr2_z = dtr2

            i3.new_zooids_to_plot.append(dtr2)
            i3.zooid_list.append(dtr2)
            i3.Es = np.concatenate([i3.Es, [dtr2.E_tot]])  # add its energy to E array
            i3.Devs = np.concatenate([i3.Devs, [dtr2.dev_state]])
            i3.Masses = np.concatenate([i3.Masses, [dtr2.mass]])

            midPoints.append(dtr2.round_mp)

            # Update the max/min geometrical characteristics for the colony
            # assign colony max/min values with maximum buffer to account for using midpoints instead of vertices
            i3.x_max = max(i3.x_max, dtr2.midPoint[0] + p * .5)
            i3.x_min = min(i3.x_min, dtr2.midPoint[0] - p * .5)
            i3.y_max = max(i3.y_max, dtr2.midPoint[1] + p * .5)
            i3.y_min = min(i3.y_min, dtr2.midPoint[1] - p * .5)

    # Designate the potential neighbor space, find neighbors and store them for all zooids in the colony
    for i8 in kelp.colony_list:
        i8.Find_neighbor_space()
        i8.Locate_potential_neighbors()
        i8.Assign_neighbors()


def Get_EMD_data(EMD_data, time_pt):
    # calculate colony areas at various time points
    for icol in kelp.colony_list:
        for jz in range(len(icol.zooid_list)):
            rowc = [time_pt, icol.ID, icol.zooid_list[jz].ID, icol.Es[jz], icol.Masses[jz], icol.Devs[jz]]
            EMD_data.append(rowc)
    return EMD_data



# collect data on colony size and after every time step - adds fine within loop but when printing master list, hasn't added anything
def Get_col_data(col_data, time_pt):
    # calculate colony areas at various time points
    for i in kelp.colony_list:
        i.Get_colony_area()
        i.Get_colony_energy()
        rowc = [time_pt, i.ID, i.spine_coeff, len(i.zooid_list), i.colony_area, i.colony_energy]
        col_data.append(rowc)
    return col_data


def Get_pop_data(pop_data, time_pt):
    kelp.Get_pop_area()
    rowp = [time_pt, len(kelp.colony_list), kelp.pop_area]
    pop_data.append(rowp)
    return pop_data


# ============================================================================================================
# Plot simulation and collect data
# ============================================================================================================


# set up interactive plot if designated true
if interactive_plot == True:
    plt.ion()

# set up the plot
fig = plt.figure()
ax = fig.add_subplot(111)

# create the habitat
kelp = Habitat(domain=dom)

# get neighbor space and locate neighbors for ancestrula zooids
for colny in kelp.colony_list:
    colny.Get_ancestrula()

for colny in kelp.colony_list:
    colny.Find_neighbor_space()
    colny.Locate_potential_neighbors()
    colny.Assign_neighbors()

# Plot ancestrula
if interactive_plot == True:
    kelp.Plot_settlement()

# Set up the plot
axes = plt.gca()
axes.set_xlim([0, kelp.domain[0]])
axes.set_ylim([0, kelp.domain[1]])
plt.gca().set_aspect('equal')


# Infrastructure to collect data:
kelp.Get_pop_area()
start_pop_area = kelp.pop_area

# set up an array to collect colony size data (time, col ID, # of zooids, col area)
# start with t0 data
col_data = []
for ic in kelp.colony_list:
    col_data.append([0.0, ic.ID, ic.spine_coeff, len(ic.zooid_list), ic.colony_area, ic.colony_energy])
    
# set up an array to collect total population data (time, # of colonies, tot area of all colonies)
# start with t0 data
pop_data = [[0.0, 0, start_pop_area]]

# establish the time duration and time steps of the simulation
ts = np.linspace(start_day,end_day,num_intervals)

# make energy array for initial ancestrula zooids - input for Habitat_integrate2()
for coln in kelp.colony_list:
    coln.Es = np.zeros([len(coln.zooid_list)])
    for z in coln.zooid_list:
        idz = z.ID
        coln.Es[idz] = z.E_tot

# make dev_state array for initial ancestrula zooids
for coln2 in kelp.colony_list:
    coln2.Devs = np.zeros([len(coln2.zooid_list)])
    for z2 in coln2.zooid_list:
        idz2 = z2.ID
        coln2.Devs[idz2] = z2.dev_state

# make mass array for initial ancestrula zooids
for coln3 in kelp.colony_list:
    coln3.Masses = np.zeros([len(coln3.zooid_list)])
    for z3 in coln3.zooid_list:
        idz3 = z3.ID
        coln3.Masses[idz3] = z3.mass

# set up an array to collect colony energy, mass and devs data (time, col ID, zooid ID, E, Mass, Dev)
# start with t0 data
EMD_data = []
for icoln in kelp.colony_list:
    id = icoln.ID
    for jz in range(len(icoln.zooid_list)):
        EMD_data.append([0.0, id, icoln.zooid_list[jz].ID, icoln.Es[jz], icoln.Masses[jz], icoln.Devs[jz]])

# Grow colonies in response to how much energy zooids have
for it in range(len(ts)-1):
    for s in range(len(settlement_times)):
        if it == settlement_times[s]:
            sp = settlement_positions[s]
            kelp.Settle(x = sp[0], y = sp[1]) # settle and plot new colonies
            max_colid = max(colony.ID for colony in kelp.colony_list)
            kelp.colony_list[max_colid].spine_coeff = S_cs[s]
    Habitat_integrate(ts[it],ts[it+1])
    Get_col_data(col_data, ts[it + 1])  # collect data after each time step
    Get_pop_data(pop_data, ts[it + 1])  # collect data after each time step
    Get_EMD_data(EMD_data, ts[it + 1])  # collect data after each time step

    if interactive_plot == True: # need to clear the plot here

        ax.cla()
        axes.set_xlim([0, kelp.domain[0]])
        axes.set_ylim([0, kelp.domain[1]])
        plt.gca().set_aspect('equal')
        plt.xlabel('Kelp length (mm)')
        plt.ylabel('Kelp width (mm)')

        ax.text(kelp.domain[0] - 16, kelp.domain[1] + 1, 't = ')
        ax.text(kelp.domain[0] - 10, kelp.domain[1] + 1, round(ts[it + 1], 1))
        ax.text(kelp.domain[0] - 2, kelp.domain[1] + 1, 'd')

        axis_text1 = str(str(len(settlement_positions)) + ' colonies')
        ax.text(0, kelp.domain[1] + 1, axis_text1)

        for i9 in kelp.colony_list:
            for z2 in i9.zooid_list:
                z2.Plot()
                if plot_nums == True:
                    ax.text(z2.midPoint[0], z2.midPoint[1], z2.ID)
            i9.Clear_new_zooids_to_plot()
        
        fig.canvas.draw()
        fig.show()

    if save_movie_frames ==  True:
        frame_num = "%it.png" % it
        fig.savefig(frame_num, format='png')

# save simulation data as csv files
if save_data == True:
    np.savetxt("col_data_" + data_name + ".csv", col_data, delimiter=",")
    np.savetxt("pop_data_" + data_name + ".csv", pop_data, delimiter=",")
    np.savetxt("EMD_data_" + data_name + ".csv", EMD_data, delimiter=",")

# plot all at the end if you do not want interactive plotting -- good for replicate runs
if interactive_plot == False:
    kelp.Plot_all() # plots all zooids of all colonies

# save final figure as png file
if save_final_fig == True:
    fig.savefig("final_fig_" + data_name + ".png", format='png', dpi=300)

print('Done', data_name)
