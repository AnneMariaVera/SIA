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

def sigma_coords(bed, surface,nz):
    """ Returns a mesh between the bed and surface with nz number of grids in the z direction. 
    The corresponding sigma mesh is also evalutated"""
    H = surface-bed    # ice thickness
    Z = np.linspace(bed,surface,nz)
    sigma = np.transpose(np.array([(surface[i]-Z[:,i])/H[i] for i in range(len(surface))]))    
    return Z,sigma

def velocity_x_inner(bed, surface, delta_x,nz, n, rho, A, g, end):
    # TODO: Think about grid in ice sheet (bottom to top)
    # Pattyn
    # Calculate the inner velocity
    # the velocity at the base is zero
    H = surface-bed    # ice thickness
    Z,sigma=sigma_coords(bed,surface,nz)
    v_x = np.array([2*A/(n+1)*(rho*g)**n*H[:end-1]**(n+1)*((surface[1:]-surface[:end-1])/delta_x)**(n-1)*
                    (surface[1:]-surface[:end-1])/delta_x*(sigma[i,:end-1]**(n+1)-1) for i in range(nz)])
    return v_x,Z,sigma

def velocity_z_inner(v_x,surface,bed,sigma,delta_x,end):
    endz = len(sigma[:,0])
    endx = end
    H = surface-bed
    v_z_inner = np.array([H[0:endx-2]*(sigma[i,:endx-2]-1)*(v_x[i,1:]-v_x[i,:endx-2])/delta_x for i in range(endz)])
    return v_z_inner

def velocity_inner(bed, surface,A,g,rho, delta_x,nz, n, end):
    v_x,Z,sigma = velocity_x_inner(bed, surface, delta_x,nz, n, rho, A, g, end)
    v_z =  velocity_z_inner(v_x,surface,bed,sigma,delta_x,end)
    return v_x,v_z,Z,sigma

def flux_11(bed, surface, boundary, delta_x, n,end, coef):
    H = surface-bed    # ice thickness
    Q = coef*np.array(((surface[1:]-surface[:end-1])/delta_x)**n*((H[1:]+H[:end-1])/2)**(n+2))
    # Append the flux at the boundaries
    Q = np.append(boundary[0,1], Q)
    Q = np.append(Q, boundary[1,1])
    return Q

def flux_10(bed, surface, boundary, delta_x, n, end, coef):
    H = surface-bed    # ice thickness
    Q = coef*np.array(((surface[1:]-surface[:end-1])/delta_x)**n*((H[1:]+H[:end-1])/2)**(n+2))
    # Append the flux at the boundaries
    Q = np.append(boundary[0,1], Q)
    return Q

def flux_01(bed, surface, boundary, delta_x, n, end, coef):
    H = surface-bed    # ice thickness
    Q = coef*np.array(((surface[1:]-surface[:end-1])/delta_x)**n*((H[1:]+H[:end-1])/2)**(n+2))
    # Append the flux at the boundaries
    Q = np.append(Q, boundary[1,1])
    return Q

def flux_00(bed, surface, boundary, delta_x, n, end, coef):
    H = surface-bed    # ice thickness
    Q = coef*((surface[1:]-surface[:end-1])/delta_x)**n*((H[1:]+H[:end-1])/2)**(n+2)
    return Q

def height_11(bed, surface, boundary, a, delta_x, delta_t, n, end, coef):
    Q = flux_11(bed, surface, boundary, delta_x, n,end, coef)
    h_new = np.array(surface[:]+delta_t/delta_x*(Q[1:end+1]-Q[:end])
                      +a*delta_t)
    h_new = np.array([bed[i] if h_new[i]-bed[i] < 0 else h_new[i]
             for i in range(len(bed))])
    return h_new,Q

def height_10(bed, surface, boundary, a, delta_x, delta_t, n, end, coef):
    Q = flux_10(bed, surface, boundary, delta_x, n, end, coef)
    h_new = np.array(surface[:end-1]+delta_t/delta_x*(Q[1:end]-Q[:end-1])
                      +a*delta_t)
    h_new = np.append(h_new,boundary[1,1])
    h_new = np.array([bed[i] if h_new[i]-bed[i] < 0 else h_new[i]
             for i in range(len(bed))])
    return h_new,Q

def height_01(bed, surface, boundary, a, delta_x, delta_t, n, end, coef):
    Q = flux_01(bed, surface, boundary, delta_x, n, end, coef)
    h_new = np.array(surface[1:end]+delta_t/delta_x*(Q[1:end]-Q[0:end-1])
                      +a*delta_t)
    h_new = np.append(boundary[0,1],h_new)
    h_new = np.array([bed[i] if h_new[i]-bed[i] < 0 else h_new[i]
             for i in range(len(bed))])
    return h_new,Q

def height_00(bed, surface, boundary, a, delta_x, delta_t, n, end, coef):
    Q = flux_00(bed, surface, boundary, delta_x, n, end, coef)
    h_new = surface[1:end-1]+delta_t/delta_x*(Q[1:end-1]-Q[:end-2])+a*delta_t
    h_new = np.append(boundary[0,1],h_new)
    h_new = np.append(h_new,boundary[1,1])
    # check if the ice sheet melted to the bottom, h_new is then equal to the bed.
    h_new = np.array([bed[i] if h_new[i]-bed[i] < 0 else h_new[i]
             for i in range(len(bed))])
    return h_new,Q

def solution(bed, surface, boundary, a, delta_x, nz, delta_t, t_0, N, n, rho, A, g):
    end = len(surface)
    coef = 2*A/(n+2)*(rho*g)**n
    if (boundary[0,0],boundary[1,0])==(0,0):
        flux = flux_00
        height = height_00
    elif (boundary[0,0],boundary[1,0])==(0,1):
        flux = flux_01
        height = height_01
    elif (boundary[0,0],boundary[1,0])==(1,0):
        flux = flux_10
        height = height_10
    elif (boundary[0,0],boundary[1,0])==(1,1):
        flux = flux_11
        height = height_11

    # save initial condition in solution list sol_h, sol_v and Q
    sol_h = [surface]
    v_x,v_z,Z,sigma = velocity_inner(bed, surface,A,g,rho, delta_x,nz, n, end)
    sol_vx=[v_x]
    sol_vz=[v_z]
    Zp=[Z]
    sigmap=[sigma]
    timesteps = 0
    h_new,Q_new = height(bed, sol_h[len(sol_h)-1], boundary, a, delta_x, delta_t, n, end, coef)
    sol_h.append(h_new)
    Q=[Q_new]
    timesteps = timesteps+1

    # Iterate through time and append solution to sol list
    while timesteps <= N and max(abs(np.array(sol_h[len(sol_h)-1])-np.array(sol_h[len(sol_h)-2])))>=0.00001:
        h_new,Q_new = height(bed, sol_h[len(sol_h)-1], boundary, a, delta_x, delta_t, n, end, coef)
        sol_h.append(h_new)
        Q.append(Q_new)
        v_x,v_z,Z,sigma = velocity_inner(bed, sol_h[len(sol_h)-1],A,g,rho, delta_x,nz, n, end)
        sol_vx.append(v_x)
        sol_vz.append(v_z)
        Zp.append(Z)
        sigmap.append(sigma)
        timesteps = timesteps+1

    Q.append(flux(bed,sol_h[len(sol_h)-1],boundary, delta_x, n, end, coef))
    return [sol_h, sol_vx,sol_vz, Q, Zp, sigmap,timesteps]