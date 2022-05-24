#!/usr/bin/env python3
"""
Created on Fri May 2 12:42:32 2022

@author: Dean Thomas

Test Biot-Savart integration in ionosphere_integrals
"""

import logging
from magnetopost import ionosphere_integrals as ii
import numpy as np
import matplotlib.pyplot as plt

def biot_savart_disk(sigma, omega, z, R, num_pts, axis=2):
    """ Compare analytic solution for a spinning, charged disk to the Biot-Savart
    numeric integration in ionosphere_integrals
    
    Analytic solution for the magnetic field along z-axis of a
    spinning, charged disk in x-y plane, centered at the origin
    
    Bz = (mu0 sigma omega)/2 * ((2 z^2 + R^2)/(z^2 + R^2)^1/2 - 2 abs(z))
    
    Arguments:
        
    R = current loop radius
    z = height above loop where B is measured
    sigma = surface charge density of disk
    omega = angular velocity of disk 
    mu0 = 1 in algorithm being tested
    num_pts = disk is centered in grid of num_pts^2 points
    axis = omega along 0 (x), 1 (y), or 2 (z) axis
    
    Returns the numeric integration of Biot-Savart and analytic result
    """
        
    # Analytic solution from Biot-Savart
    Banalytic = (sigma * omega)/2 * ((2*z**2 + R**2)/np.sqrt(z**2 + R**2) \
        - 2*np.abs(z))

    # Do numeric integration

    # To ensure that we include the origin, make sure num_pts is odd
    n_pts = num_pts
    if( np.mod(n_pts,2) == 0 ):
        n_pts = n_pts + 1

    # Create a grid n_pts by n_pts

    # Initialize x-y-z arrays
    x_arr = np.zeros(n_pts**2)
    y_arr = np.zeros(n_pts**2)
    z_arr = np.zeros(n_pts**2)
    
    # Initialize Jx, Jy, and Jz arrays
    Jx_arr = np.zeros(n_pts**2)
    Jy_arr = np.zeros(n_pts**2)
    Jz_arr = np.zeros(n_pts**2)
   
    # measure = dx * dy * dz 
    dx = np.float64( 2*R/(n_pts - 1) )
    dy = dx
    dz = dx  
    measure = np.full(n_pts**2, dx*dy*dz)
    
    # Determine which grid points are on disk
    # If not on disk, J = 0, otherwise J is non-zero
    i = 0;
    while i < n_pts:
        j = 0;
        while j < n_pts:
            x_arr[i*n_pts + j] = -R + j*dx
            y_arr[i*n_pts + j] = -R + i*dy
            r_xy = np.sqrt(x_arr[i*n_pts + j]**2 + y_arr[i*n_pts + j]**2)
            
            # See if we're on the disk 
            if( R >= r_xy ):
                # Clockwise current loop in x-y plane when viewed from positive z
                # Magnitude of J = rho r omega = sigma/dz r omega 
                # J = (sigma/dz) (omega z-hat) x (x i-hat + y j-hat)
                Jx_arr[i*n_pts + j] = - sigma * omega * y_arr[i*n_pts + j] / dz
                Jy_arr[i*n_pts + j] =   sigma * omega * x_arr[i*n_pts + j] / dz 
            
            # Increment counters
            j += 1
        i += 1            
    
    # Observation point, the z in the analytic solution 
    # always along axis of rotation (aka axis)
    x0 = np.zeros(3)
    x0[axis] = z
    
    if( axis == 0 ):
        Bnum = ii._integral_bs(z_arr, x_arr, y_arr, x0, Jz_arr, Jx_arr, Jy_arr, measure)
    elif ( axis == 1 ):
        Bnum = ii._integral_bs(y_arr, z_arr, x_arr, x0, Jy_arr, Jz_arr, Jx_arr, measure)
    elif ( axis == 2 ):
        Bnum = ii._integral_bs(x_arr, y_arr, z_arr, x0, Jx_arr, Jy_arr, Jz_arr, measure)
    else:
        Bnum = 0
        
    logging.info('Test ionosphere Biot-Savart spinning charged disk')
    logging.info('Num points 1-D: {:}'.format(num_pts))
    logging.info('Axis: {:}'.format(axis))
    logging.info('Numeric estimate: {:}'.format(Bnum))
    logging.info('Analytic estimate: {:}'.format(Banalytic))
    
    # Need to select which B component is seen at observation point (aka axis) 
    return Bnum[axis], Banalytic

def biot_savart_line(I, z, L, num_pts, axis = 0):
    """Compare analytic solution for a current in a straight wire to the 
    Biot-Savart numeric integration in ionosphere_integrals
    
    Analytic solution for the magnetic field along z-axis for current in a 
    finite, straight wire along the x-axis, centered at the origin

    Bz = mu0 I/(4 pi z) * L/sqrt(L^2/4 + z^2)
    
    Arguments:
        
    z = height above wire where B is measured
    I = current in wire
    L = total length of wire, symmetric about origin along x-axis
    mu0 = 1 in algorithm being tested
    num_pts = L is divided into num_pts segments
    axis = wire along x (0), y (1), or z (2) axis
    
    Returns the numeric integration of Biot-Savart and analytic result
    """
      
    # Analytic solution based on Biot-Savart
    Banalytic = - I/4/z/np.pi * L/np.sqrt(L**2/4+z**2)
    
    # To ensure that we include origin, make sure num_pts is odd
    n_pts = num_pts
    if( np.mod(n_pts,2) == 0 ):
        n_pts = n_pts + 1

    # x-y-z pts
    cnt_arr = np.arange(0,n_pts)
    x_arr = -L/2 + L/n_pts * cnt_arr[:]
    y_arr = np.zeros(n_pts)
    z_arr = np.zeros(n_pts)
    
    # dx, dy, and dz
    dx_arr = np.full(n_pts, L/(n_pts-1))
    dy_arr = dx_arr
    dz_arr = dx_arr

    # measure = dx * dy * dz 
    measure = dx_arr * dy_arr * dz_arr
    
    # Current density vector at each point
    # Magnitude is I / perpendicular area
    Jx_arr = I / (dy_arr * dz_arr)
    Jy_arr = np.zeros(n_pts)
    Jz_arr = np.zeros(n_pts)
    
    # Observation point, the z in the analytic solution
    # Need to rotate observation point as axis changes
    # x-axis -> z obs. pt, y-axis -> x obs. pt, z-axis -> y obs. pt
    x0 = np.zeros(3)
    x0[np.mod(axis+2,3)] = z
    
    if (axis == 0):
        Bnum = ii._integral_bs(x_arr, y_arr, z_arr, x0, Jx_arr, Jy_arr, Jz_arr, measure)
    elif (axis == 1):
        Bnum = ii._integral_bs(z_arr, x_arr, y_arr, x0, Jz_arr, Jx_arr, Jy_arr, measure)
    elif (axis == 2):
        Bnum = ii._integral_bs(y_arr, z_arr, x_arr, x0, Jy_arr, Jz_arr, Jx_arr, measure)
    else:
        Bnum = 0.
        
    logging.info('Test ionosphere Biot-Savart for line current')
    logging.info('Num points 1-D: {:}'.format(num_pts))
    logging.info('Axis: {:}'.format(axis))
    logging.info('Numeric estimate: {:}'.format(Bnum))
    logging.info('Analytic estimate: {:}'.format(Banalytic))
    
    # Need to select which B component is seen at observation point  
    # x-axis -> y component, y-axis -> z component, z-axis -> x component
    return Bnum[np.mod(axis+1,3)], Banalytic

def loop_thru_tests( test_type, axis, num = 6 ):
    """Loop through test cases looking at a grid size varying from 10^1 to 10^num
    
    Arguments:
        
    test_type = type of Biot-Savart test to execute, i.e., line of current (1)
                or spinning charged disk (2)
    axis = vary orientation, i.e., change axis for line current or spin axis 
            of disk
    num = number of orders of magnitude to vary grid size. 
    
    Note: 1-D current line grows as 10^num while 2-D spinning disk grows as 10^(2*num)
    """
    
    Bnum = np.zeros(num)
    Banalytic = np.zeros(num)
    cnt = np.arange(1,num+1)
    
    logging.info('Test ionosphere loop...')

    i = 0
    while i < num:
        if( test_type == 1):
            Bnum[i], Banalytic[i] = biot_savart_line(1, 10, 256, 10**(i+1), axis)
        else: 
            Bnum[i], Banalytic[i] = biot_savart_disk(1, 1, 100, 1, 10**(i+1), axis)
        i += 1
        
    fig, ax = plt.subplots()
    ax.plot(cnt, abs(Bnum[:] - Banalytic[:]))
    if (test_type == 1):
        plt.title('Line current, Axis {:}'.format(axis))
    else:
        plt.title('Spinning charged disk, Axis {:}'.format(axis))
    plt.yscale('log')
    plt.ylabel('log(Bnumeric - Banalytic)')
    plt.xlabel('log(Num Pts in 1-D)')
    plt.show()

def test_loops():
    """Test ionosphere integrals with line current and spinning charged disk
    with current and disk aligned along all three axes x (0), y (1), and z (2)
    """
    
    logging.info('Test multiple ionosphere cases...')

    # line current
    loop_thru_tests( 1, 0, 7 )
    loop_thru_tests( 1, 1, 7 )
    loop_thru_tests( 1, 2, 7 )
    # spinning charged disk
    loop_thru_tests( 2, 0, 3 )    
    loop_thru_tests( 2, 1, 3 )    
    loop_thru_tests( 2, 2, 3 )    

if __name__ == "__main__":
    test_loops()
