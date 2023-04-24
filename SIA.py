"""
@author: Anne Hermann

Overview: 
    Calculate the velocities in x direction at the surface of an ice sheet 
    along a flow line. 
    
    Gov. equ.: see Greve and Blater eq. (3.84), p. 79. 
    Assumption: 1. velocity at base is equal to zero.
                2. the bed is constant along the x axis.
        
Input variables: 
    base        Location of the base of the ice sheet
    surface     Location of the surface of the ice sheet 
    boundary    contains information about the boundary condition
    a           Accumultaion rate
    delta_x     Size of the grid in x dir.
    t_0         start time
    N           Number of time steps
    delta_t     step size time
    n           Glen flow law exponent (3 or 4)
    rho         ice density
    A           flow rate factor (Temp. dependent)
    g           gravity 
"""
# TODO: Add boundary condition and base    
# TODO: Add accumulation
# TODO: Adjust for Flux

import numpy as np

def velocity(base,surface,boundary,delta_x,n,rho,A,g): 
    # Calculate velocity at the top
    H = surface-base
    # velocity at the base is zero
    v_x = -2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i]-surface[i-1])/delta_x)**n*H[i]**(n+1) 
           for i in range(1,len(surface))])
    return v_x

def flux(base,surface,boundary,delta_x,delta_t,n,rho,A,g):
    H = surface-base
    Q = 2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i+1]-surface[i])/delta_x)**n*((H[i+1]+H[i])/2)**(n+2) 
         for i in range(0,len(surface)-1)])  
    Q=np.append(boundary[0],Q)
    Q=np.append(Q,boundary[1])
    return Q
    

def height(base,surface,boundary,a,delta_x,delta_t,n,rho,A,g):
    Q = flux(base,surface,boundary,delta_x,delta_t,n,rho,A,g)
    h_new = np.array([surface[i]+delta_t/delta_x*(Q[i+1]-Q[i]) 
                    for i in range(0,len(Q)-1)])+delta_t*a
    h_new = [0 if h_new[i]-base[i]<0 else h_new[i] for i in range(0,len(base))] 
    return h_new

def solution(base,surface,boundary,a,delta_x,delta_t,t_0,N,n,rho,A,g):
    # save initial condition in solution list sol
    sol_h = [surface]
    sol_v = [velocity(base,surface,boundary,delta_x,n,rho,A,g)]
    Q=[flux(base, surface, boundary, delta_x, delta_t, n, rho, A, g)]
    time = t_0 
    # Iterate through time and append solution to sol list
    while time<=t_0+(N-1)*delta_t:
        time = time + delta_t
        sol_h.append(height(base,sol_h[len(sol_h)-1],boundary,a,delta_x,delta_t,n,rho,A,g))
        sol_v.append(velocity(base,sol_h[len(sol_h)-1],boundary,delta_x,n,rho,A,g)) 
        Q.append(flux(base, surface, boundary, delta_x, delta_t, n, rho, A, g))

    return [sol_h,sol_v,Q]