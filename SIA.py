"""
@author: Anne Hermann

Overview: 
    Calculate the velocities in x direction at the surface of an ice sheet 
    along a flow line. 
    
    The boundary condition is a prescribed flux condition. 
    
    Gov. equ.: see Greve and Blater eq. (3.84), p. 79. 
    Assumption: 1. velocity at base is equal to zero.
                2. the bed is constant along the x axis.
        
Input variables: 
    bed         Location of the bed of the ice sheet
    surface     Location of the surface of the ice sheet 
    boundary    contains information about the flux at the boundary
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

import numpy as np

def velocity_00(bed, surface, boundary, delta_x, n, rho, A, g):
    # Calculate the mean velocity
    # the velocity at the bed is zero
    H = surface-bed    # ice thickness
    v_x = -2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i+1]-surface[i])/delta_x)**n*H[i]**(n+1)
            for i in range(0, len(surface)-1)])
    return v_x

def velocity_11(bed, surface, boundary, delta_x, n, rho, A, g):
    # Calculate the mean velocity
    # the velocity at the bed is zero
    H = surface-bed    # ice thickness
    v_x = -2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i+1]-surface[i])/delta_x)**n*H[i]**(n+1)
            for i in range(0, len(surface)-1)])
    if H[len(H)-1]!=0:
        v_x = np.append(v_x,boundary[1,0]/H[len(H)-1])
    return v_x

def velocity_01(bed, surface, boundary, delta_x, n, rho, A, g):
    # Calculate the mean velocity
    # the velocity at the bed is zero
    H = surface-bed    # ice thickness
    v_x = -2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i+1]-surface[i])/delta_x)**n*H[i]**(n+1)
            for i in range(0, len(surface)-1)])
    if H[len(H)-1]!=0:
        v_x = np.append(v_x,boundary[1,0]/H[len(H)-1])
    return v_x

def velocity_10(bed, surface, boundary, delta_x, n, rho, A, g):
    # Calculate the mean velocity
    # the velocity at the bed is zero
    H = surface-bed    # ice thickness
    v_x = -2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i]-surface[i-1])/delta_x)**n*H[i]**(n+1)
            for i in range(1, len(surface))])
    if H[0]!=0:
        v_x = np.append(boundary[0,0]/H[0],v_x)
    return v_x

def flux_11(bed, surface, boundary, delta_x, delta_t, n, rho, A, g):
    H = surface-bed    # ice thickness
    Q = 2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i+1]-surface[i])/delta_x)**n*((H[i+1]+H[i])/2)**(n+2)
         for i in range(0, len(surface)-1)])
    # Append the flux at the boundaries
    Q = np.append(boundary[0,1], Q)
    Q = np.append(Q, boundary[1,1])
    return Q

def flux_10(bed, surface, boundary, delta_x, delta_t, n, rho, A, g):
    H = surface-bed    # ice thickness
    Q = 2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i+1]-surface[i])/delta_x)**n*((H[i+1]+H[i])/2)**(n+2)
         for i in range(0, len(surface)-1)])
    # Append the flux at the boundaries
    Q = np.append(boundary[0,1], Q)
    return Q

def flux_01(bed, surface, boundary, delta_x, delta_t, n, rho, A, g):
    H = surface-bed    # ice thickness
    Q = 2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i+1]-surface[i])/delta_x)**n*((H[i+1]+H[i])/2)**(n+2)
         for i in range(0, len(surface)-1)])
    # Append the flux at the boundaries
    Q = np.append(Q, boundary[1,1])
    return Q

def flux_00(bed, surface, boundary, delta_x, delta_t, n, rho, A, g):
    H = surface-bed    # ice thickness
    Q = 2*A/(n+2)*(rho*g)**n*np.array(
        [((surface[i+1]-surface[i])/delta_x)**n*((H[i+1]+H[i])/2)**(n+2)
         for i in range(0, len(surface)-1)])
    return Q

# def velocity_inner(bed, surface, delta_x,nz, n, rho, A, g):
#     # TODO: Think about grid in ice sheet (bottom to top)
#     # Pattyn
#     # Calculate the inner velocity
#     # the velocity at the base is zero
#     H = surface-bed    # ice thickness
#     v_inner = np.array([np.linspace(bed[i],surface[i],nz)])
#     v_x = -2*A/(n+2)*(rho*g)**n*np.array(
#         [((surface[i+1]-surface[i])/delta_x)**n*H[i]**(n+1)
#             for i in range(0, len(surface)-1)])
#     return v_x

def height_11(bed, surface, boundary, a, delta_x, delta_t, n, rho, A, g):
    Q = flux_11(bed, surface, boundary, delta_x, delta_t, n, rho, A, g)
    h_new = np.array([surface[i]+delta_t/delta_x*(Q[i+1]-Q[i])
                      +a*delta_t for i in range(0, len(Q)-1)])
    h_new = [bed[i] if h_new[i]-bed[i] < 0 else h_new[i]
             for i in range(0, len(bed))]
    return h_new

def height_10(bed, surface, boundary, a, delta_x, delta_t, n, rho, A, g):
    Q = flux_10(bed, surface, boundary, delta_x, delta_t, n, rho, A, g)
    h_new = np.array([surface[i]+delta_t/delta_x*(Q[i+1]-Q[i])
                      +a*delta_t for i in range(0, len(surface)-1)])
    h_new = np.append(h_new,boundary[1,1])
    h_new = [bed[i] if h_new[i]-bed[i] < 0 else h_new[i]
             for i in range(0, len(bed))]
    return h_new

def height_01(bed, surface, boundary, a, delta_x, delta_t, n, rho, A, g):
    Q = flux_01(bed, surface, boundary, delta_x, delta_t, n, rho, A, g)
    h_new = np.array([surface[i]+delta_t/delta_x*(Q[i]-Q[i-1])
                      +a*delta_t for i in range(1, len(surface)-1)])
    h_new = np.append(boundary[0,1],h_new)
    h_new = [bed[i] if h_new[i]-bed[i] < 0 else h_new[i]
             for i in range(0, len(bed))]
    return h_new

def height_00(bed, surface, boundary, a, delta_x, delta_t, n, rho, A, g):
    Q = flux_00(bed, surface, boundary, delta_x, delta_t, n, rho, A, g)
    h_new = np.array([surface[i]+delta_t/delta_x*(Q[i]-Q[i-1])
                      +a*delta_t for i in range(1, len(surface)-1)])
    h_new = np.append(boundary[0,1],h_new)
    h_new = np.append(h_new,boundary[1,1])
    h_new = [bed[i] if h_new[i]-bed[i] < 0 else h_new[i]
             for i in range(0, len(bed))]
    return h_new


def solution(bed, surface, boundary, a, delta_x, delta_t, t_0, N, n, rho, A, g):
    if (boundary[0,0],boundary[1,0])==(0,0):
        flux = flux_00
        height = height_00
        velocity=velocity_00
    elif (boundary[0,0],boundary[1,0])==(0,1):
        flux = flux_01
        height = height_01
        velocity = velocity_01
    elif (boundary[0,0],boundary[1,0])==(1,0):
        flux = flux_10
        height = height_10
        velocity = velocity_10
    elif (boundary[0,0],boundary[1,0])==(1,1):
        flux = flux_11
        height = height_11
        velocity = velocity_11


    # save initial condition in solution list sol_h, sol_v and Q
    sol_h = [surface]
    sol_v = [velocity(bed, surface, boundary, delta_x, n, rho, A, g)]
    Q = [flux(bed, surface, boundary, delta_x, delta_t, n, rho, A, g)]
    timesteps = 1
    sol_h.append(height(bed, sol_h[len(sol_h)-1], boundary, a, delta_x, delta_t, n, rho, A, g))
    sol_v.append(velocity(bed, sol_h[len(sol_h)-1], boundary ,delta_x, n, rho, A, g))
    Q.append(flux(bed, sol_h[len(sol_h)-1],boundary, delta_x, delta_t, n, rho, A, g))
    timesteps = timesteps+1
    # Iterate through time and append solution to sol list
    while timesteps <= N and max(abs(np.array(sol_h[len(sol_h)-1])-np.array(sol_h[len(sol_h)-2])))>=0.00001:
        sol_h.append(
            height(bed, sol_h[len(sol_h)-1], boundary, a, delta_x, delta_t, n, rho, A, g))
        sol_v.append(
            velocity(bed, sol_h[len(sol_h)-1], boundary ,delta_x, n, rho, A, g))
        Q.append(flux(bed, sol_h[len(sol_h)-1],
                 boundary, delta_x, delta_t, n, rho, A, g))
        timesteps = timesteps+1

    return [sol_h, sol_v, Q,timesteps]
