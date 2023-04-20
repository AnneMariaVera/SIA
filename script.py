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
# grid size in x dir.
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

# --------------------  SIA velocity cal. of initial state ------------------- 
# calculate the velocity 
v_x = sia.velocity(surface,delta_x,n,rho,A,g)

# -------------------------   Plot initial velocity    ----------------------- 
fig , ax = plt.subplots(2)
# increase spacing between subplots
fig.tight_layout(pad=3.0)
# create plot
ax[0].plot(bed[1:len(bed)],surface[1:len(surface)])
ax[1].plot(bed[1:len(bed)],v_x,"r.",markersize=1)
# label axes
ax[0].set( ylabel = r"$h\, [m]$")
ax[1].set(xlabel = r"$x\, [m]$", ylabel = r"$v_x\,[ma^{-1}]$")
# switch off borders in plot
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)



# --------------------------   SIA Solution Plot   --------------------------- 
# calculate solution
#TODO: change boundary cond. in SIA file!
sol = sia.solution(surface, delta_x, delta_t, t_0, N, n, rho, A, g)

# create plots
for i in range(0,len(sol),2):
    fig , ax = plt.subplots()
    # create plot
    ax.plot(bed,surface,"k-")
    ax.plot(bed,sol[i],"c.",markersize=3)
    # add labels and title
    ax.set(title = f"time = {t_0+i*delta_t} [a]",xlabel=r"$x\,[m]$"
           ,ylabel=r"$h\,[m]$")
    plt.rc('axes',titlesize=12)
    # switch off borders in plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
