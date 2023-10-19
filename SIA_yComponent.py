"""
@author: Anne Hermann

Overview: 
    Calculate the velocity field and the geometry evolution of an Ice sheet 
    along a cross-section acounting for transversal flow.
    
    The boundary condition is a prescribed flux condition. 
    
    Gov. equ.: see Greve and Blater eq. (3.84), p. 79. 
    Assumption: velocity at the base is equal to zero (no basal sliding).
        
Input variables: 
    bed         [m]         3D bed of the ice sheet 
    surface     [m]         3D surface of the ice sheet 
    boundary                contains information about the flux at the boundary
    a           [m/a]       2D Accumultaion rate
    dx          [m]         Size of the grid in x dir.
    nz                      number of gridpoints in z direction 
    t_0         [a]         start time
    N                       Number of time steps
    dt          [a]         step size time
    n                       Glen flow law exponent (3 or 4)
    rho         [kg/m^3]    ice density
    A           [Pa/s]      flow rate factor (Temp. dependent)
    g           [m/s^2]     gravity 
"""

import numpy as np

def velocity_field(bed,surface,a,dx,dy,nz,n,coef):
    H=surface-bed
    Z = np.linspace(bed[1,:],surface[1,:],nz)
    vx =coef*(n+2)/(n+1)*(((surface[1,1:]-surface[1,:-1])/dx)**2+((surface[2,:-1]-surface[0,:-1])/(2*dy))**2)**((n-1)/2)*(surface[1,1:]-surface[1,:-1])/dx*np.array([((surface[1,:-1]-Z[i][:-1])**(n+1)-H[1,:-1]**(n+1)) for i in range(nz)]) 
    vy =coef*(n+2)/(n+1)*(((surface[1,1:]-surface[1,:-1])/dx)**2+((surface[2,:-1]-surface[0,:-1])/(2*dy))**2)**((n-1)/2)*(surface[2,:-1]-surface[0,:-1])/(2*dx)*np.array([((surface[1,:-1]-Z[i][:-1])**(n+1)-H[1,:-1]**(n+1)) for i in range(nz)]) 
    return vx,vy,Z

def flux(bed,surface,a,dx,dy,n,coef):
    H = surface-bed    # ice thickness  

    # Boundary condition 
    xdiv = np.argmax(surface[1,:])
    rB = sum(a[1,xdiv:]*dx)
    lB = -sum(a[1,:xdiv]*dx)

    xdiv0 = np.argmax(surface[0,:])
    rB0 = sum(a[0,xdiv0:]*dx)
    lB0 = -sum(a[0,:xdiv0]*dx)
    
    xdiv2 = np.argmax(surface[2,:])
    rB2 = sum(a[2,xdiv2:]*dx)
    lB2 = -sum(a[2,:xdiv2]*dx)

    # Calculate fluxes
    Qx = np.zeros((3,np.size(H[0,:])+1))
    Qx[0,1:-1] = -coef*(((surface[0,1:]-surface[0,:-1])/dx)**2+((surface[1,:-1]-surface[0,:-1])/(dy))**2)**((n-1)/2)*(surface[0,1:]-surface[0,:-1])/dx*((H[0,1:]+H[0,:-1])/2)**(n+2)
    Qx[1,1:-1] = -coef*(((surface[1,1:]-surface[1,:-1])/dx)**2+((surface[2,:-1]-surface[0,:-1])/(2*dy))**2)**((n-1)/2)*(surface[1,1:]-surface[1,:-1])/dx*((H[1,1:]+H[1,:-1])/2)**(n+2)
    Qx[2,1:-1] = -coef*(((surface[2,1:]-surface[2,:-1])/dx)**2+((surface[2,:-1]-surface[1,:-1])/(dy))**2)**((n-1)/2)*(surface[2,1:]-surface[2,:-1])/dx*((H[2,1:]+H[2,:-1])/2)**(n+2)
    Qx[0,0] = lB0
    Qx[0,-1] = rB0
    Qx[1,0] = lB
    Qx[1,-1] = rB
    Qx[2,0] = lB2
    Qx[2,-1] = rB2

    Qy = np.zeros((4,np.size(H[0,:])))
    Qy[1,:-1] = -coef*(((surface[1,1:]-surface[1,:-1])/dx)**2+((surface[1,:-1]-surface[0,:-1])/(dy))**2)**((n-1)/2)*(surface[1,:-1]-surface[0,:-1])/dy*(H[0,:-1]+H[1,:-1])/2
    Qy[2,:-1] = -coef*(((surface[1,1:]-surface[1,:-1])/dx)**2+((surface[2,:-1]-surface[1,:-1])/(dy))**2)**((n-1)/2)*(surface[2,:-1]-surface[1,:-1])/dy*(H[2,:-1]+H[1,:-1])/2

    return Qx,Qy

def height(bed,surface,dx,dy,dt,a,n,coef):
    Qx,Qy = flux(bed,surface,a,dx,dy,n,coef)
    h_new = np.zeros((3,np.size(bed[0,:])))
    h_new[0,:] = surface[0,:]-dt/dx*(Qx[0,1:]-Qx[0,:-1])-dt/dy*(Qy[1,:]-Qy[0,:])+a[0,:]*dt
    h_new[1,:] = surface[1,:]-dt/dx*(Qx[1,1:]-Qx[1,:-1])-dt/dy*(Qy[2,:]-Qy[1,:])+a[1,:]*dt
    h_new[2,:] = surface[2,:]-dt/dx*(Qx[2,1:]-Qx[2,:-1])-dt/dy*(Qy[3,:]-Qy[2,:])+a[2,:]*dt
    h_new[0,:] = np.array([bed[0,i] if h_new[0,i]-bed[0,i] < 0 else h_new[0,i]
             for i in range(len(bed[0,:]))])
    h_new[1,:] = np.array([bed[1,i] if h_new[1,i]-bed[1,i] < 0 else h_new[1,i]
             for i in range(len(bed[1,:]))])
    h_new[2,:] = np.array([bed[2,i] if h_new[2,i]-bed[2,i] < 0 else h_new[2,i]
             for i in range(len(bed[2,:]))])
    return h_new

def solution(bed, surface, a, dx, dy, nz, dt, N, n, rho, A, g):
    # end = len(surface)
    coef = 2*A/(n+2)*(rho*g)**n

    # save initial condition in solution list sol_h, sol_v and Q
    sol_h = [surface]
    vx,vy,Z = velocity_field(bed,surface,a,dx,dy,nz,n,coef)
    sol_vx=[vx]
    sol_vy=[vy]
    Zgrid=[Z]
    timesteps=0
    # calculate first step

    h_new = height(bed, sol_h[-1], dx, dy, dt, a, n, coef)
    sol_h.append(h_new)
    vx,vy,Z = velocity_field(bed,sol_h[-1],a,dx,dy,nz,n,coef)
    Zgrid.append(Z)
    sol_vx.append(vx)
    sol_vy.append(vy)

    timesteps = 1

    # Iterate through time and append solution to sol list
    while timesteps < N and max(abs(sol_h[-1][1,:]-sol_h[-2][1,:]))>=0.00001:
        h_new = height(bed, sol_h[-1], dx, dy, dt, a, n, coef)
        sol_h.append(h_new)
        vx,vy,Z = velocity_field(bed,sol_h[-1],a,dx,dy,nz,n,coef)
        sol_vx.append(vx)
        sol_vy.append(vy)
        Zgrid.append(Z)
        timesteps = timesteps+1

    return [sol_h, sol_vx, sol_vy, Zgrid,timesteps]
