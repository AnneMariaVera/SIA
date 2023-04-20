import numpy as np

def velocity(h,delta_x,n,rho,A,g):
    """
    @author: Anne Hermann

    Overview: 
        Calculate the velocities in x direction at the surface of an ice sheet 
        along a flow line. 
    
        Gov. equ.: see Greve and Blater eq. (3.84), p. 79. 
        Assumption: 1. velocity at base is equal to zero.
                    2. the bed is constant along the x axis.
        
    Input: 
        h           Heights of the ice sheet 
        delta_x     Size of the grid in x dir.   
        n           Glen flow law exponent (3 or 4)
        rho         ice density
        A           flow rate factor (Temp. dependent)
        g           gravity 
        
    Output: 
        velocity in x direction 
    """
       
    # Calculate velocity at the top
    v_x = np.array([2*(rho*g)**3*((h[i-1]-h[i])/delta_x)**3*A*1/4*
                    (h[i])**4 for i in range(1,np.size(h))])
    return v_x

def height(h,delta_x,delta_t,n,rho,A,g):
    """
    @author: Anne Hermann

    Overview: 
        Calculate the new height of the ice sheet along a flow line according 
        to the velocity.
    
        Gov. equ.: see Pattyn eq. (3.43), p.17. 
        Discretization: see Pattyn eq. (3.45) and (3.46)
        Assumption: 1. velocity at base is equal to zero.
                    2. the bed is constant along the x axis.
        
    Input: 
        h           Heights of the ice sheet 
        delta_x     Size of the grid in x dir. 
        delta_t     time step size 
        n           Glen flow law exponent (3 or 4)
        rho         ice density
        A           flow rate factor (Temp. dependent)
        g           gravity 
        
    Output: 
        new height  
    """ 
    
    # Diffusion coef. 
    D = [2*A/(n+2)*(rho*g)**n*((h[i+1]-h[i])/delta_x)**(n-1)*
         ((h[i]+h[i+1])/2)**(n+2) for i in range(0,np.size(h)-1)]
    
    # height 
    h_new = [h[i]+delta_t/(delta_x)**2*(D[i]*(h[i+1]-h[i])-
                        D[i-1]*(h[i]-h[i-1])) for i in range(1,np.size(h)-1)]
    
    return np.array(h_new)

def solution(h,delta_x,delta_t,t_0,N,n,rho,A,g):
    """
    @author: Anne Hermann

    Overview: 
        Calculate the solutions for different time steps.  
    
        Assumption: 1. velocity at base is equal to zero.
                    2. the bed is constant along the x axis.
        
    Input: 
        h           Initial heights of the ice sheet 
        delta_x     Size of the grid in x dir.   
        delta_t     time steps
        t_0         Start time 
        N           Number of time steps
        n           Glen flow law exponent (3 or 4)
        v_bx        ice velocity at the base
        rho         ice density
        A           flow rate factor (Temp. dependent)
        g           gravity 
        
    Output: 
        velocity in x direction 
    """
    # save initial condition in solution list sol
    sol = [h]
    
    # First time step
    surface = height(h, delta_x, delta_t, n, rho, A, g)
    time = t_0 + delta_t

    #TODO: boundary condition
    # h=0 at the left and right boundary right now
    surface = np.append(np.append([0],surface),[0])
    sol.append(surface)
    
    # Iterate through time and append solution to sol list
    while time<=t_0+(N-1)*delta_t:
        time = time + delta_t
        surface = height(surface, delta_x, delta_t, n, rho, A, g)
        surface = np.append(np.append([0],surface),[0])
        sol.append(surface)
        
    return sol