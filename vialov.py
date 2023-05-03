"""
@author: Anne Hermann

Overview: Calculate the Vialov Profile from Greve (eq. 5.117 and 5.119).

Input:
    L = extent of ice sheet to the left and right from the peak at x = 0
    a = accumulation rate 
    delta_x = grid size in x dir.
    n = 3 or 4, Glen's flow law
    A = flow rate factor
    rho = density of the ice 
    g = gravity
    
"""

import numpy as np
def vialov(L,a,delta_x,n,A,rho,g):
    A_0=2*A*(rho*g)**n/(n+2)
    h_0=2**(n/(n+2))*(a/A_0)**(1/(2*n+2))*L**(1/2)
    x = np.arange(-L,L+delta_x,delta_x)
    h = h_0*(1-(abs(x)/L)**((n+1)/n))**(n/(2*n+2))
    return h