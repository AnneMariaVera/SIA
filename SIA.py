"""
@author: Anne Hermann

Overview: 
    Calculate the velocities in x direction at the surface of an ice sheet 
    along a flow line. 
    
    Gov. equ.: see Greve and Blater eq. (3.84), p. 79. 
    Assumption: 1. velocity at base is equal to zero.
                2. the bed is constant along the x axis.
        
Input variables: 
    h           Heights of the ice sheet 
    delta_x     Size of the grid in x dir.
    t_0         start time
    N           Number of time steps
    delta_t     step size time
    n           Glen flow law exponent (3 or 4)
    rho         ice density
    A           flow rate factor (Temp. dependent)
    g           gravity 
"""
      

import numpy as np

def velocity(h,delta_x,n,rho,A,g): 
    # Calculate velocity at the top
    v_x = np.array([2*(rho*g)**3*((h[i-1]-h[i])/delta_x)**3*A*1/4*
                    (h[i])**4 for i in range(1,np.size(h))])
    return v_x

def height(h,delta_x,delta_t,n,rho,A,g):
    # Diffusion coef. 
    D = [2*A/(n+2)*(rho*g)**n*((h[i+1]-h[i])/delta_x)**(n-1)*
         ((h[i]+h[i+1])/2)**(n+2) for i in range(0,np.size(h)-1)]
    
    # height 
    h_new = [h[i]+delta_t/(delta_x)**2*(D[i]*(h[i+1]-h[i])-
                        D[i-1]*(h[i]-h[i-1])) for i in range(1,np.size(h)-1)]
    
    return np.array(h_new)

def solution(h,delta_x,delta_t,t_0,N,n,rho,A,g):
    # save initial condition in solution list sol
    sol_h = [h]
    sol_v = [velocity(h,delta_x,n,rho,A,g)]
    # First time step
    surface = height(h, delta_x, delta_t, n, rho, A, g)
    time = t_0 + delta_t

    #TODO: boundary condition
    # h=0 at the left and right boundary right now
    surface = np.append(np.append([0],surface),[0])
    sol_h.append(surface)
    sol_v.append(velocity(surface, delta_x, n, rho, A, g))
    
    # Iterate through time and append solution to sol list
    while time<=t_0+(N-1)*delta_t:
        time = time + delta_t
        surface = height(surface, delta_x, delta_t, n, rho, A, g)
        surface = np.append(np.append([0],surface),[0])
        sol_h.append(surface)
        sol_v.append(velocity(surface, delta_x, n, rho, A, g))
        
    return [sol_h,sol_v]