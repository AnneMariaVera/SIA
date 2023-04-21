"""
@author: Anne Hermann

Overview: 
    Set up variables for the SIA model.
    Create the geometry on an artifical ice sheet and discretize it. 
    Plot initial velocity field and heights over time.
    
    Assumption: 1. velocity at base is equal to zero.
                2. the bed is constant along the x axis.
                3. height at left and right boundary is zero
    
"""

import SIA as sia
import matplotlib.pyplot as plt
import numpy as np

# --------------------------  initialize variables  -------------------------- 
n = 3           # Glens flow law
rho = 910       # kg m^-3
A = 10E-16      # a^-1 Pa^-3 
g = 9.81        # m s^-2

# -------------------------      Discretization      ------------------------- 
# grid in x dir.
delta_x = 5
x_lim = 100
bed = np.arange(-x_lim,x_lim+1,delta_x)

# time discretization
t_0=0
N=50
delta_t=100

# ----------------------------     Elevation      ----------------------------
#TODO: look for a better initial ice sheet
surface = -0.001*bed**2+10  

# --------------------------   SIA Solution Plot   --------------------------- 
# calculate solution
#TODO: change boundary cond. in SIA file!
h,v = sia.solution(surface, delta_x, delta_t, t_0, N, n, rho, A, g)

# Create Plots
v_min = np.min(v[0])
v_max = np.max(v[0])
# create plots
for i in range(0,len(h),2):
    fig , ax = plt.subplots(2)
    # create plot
    ax[0].plot(bed,surface,"k-")
    ax[0].plot(bed,h[i],"c.",markersize=3)
    # add labels and title
    ax[0].set(title = f"time = {t_0+i*delta_t} [a]",ylabel=r"$h\,[m]$")
    plt.rc('axes',titlesize=12)
    # switch off borders in plot
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[1].set_ylim([v_min, v_max])
    ax[1].plot(bed[1:len(bed)],v[i])
    ax[1].set(ylabel=r"$v_x\,[ma^{-1}]$",xlabel=r"$x\,[m]$")
