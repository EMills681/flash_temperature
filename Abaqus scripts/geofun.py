""" geofun - a module of geometry functions"""
#this will need to be saved in the Abaqus working directory for input_flash_temp to work
# This code is used for the geometry sketch when creating a part in Abaqus
import math
from os import X_OK
import numpy as np

####################################################################

### Variables
## | Variable   | Description                                   | ##
## |------------|-----------------------------------------------| ##
## |     m      | module (mm)                                   | ##
## |     z      | number of teeth (mm)                          | ##
## |    psi     | pressure angle   (radians)                   | ##
## |     r      | radius (mm)                                  | ##
## |     d      | pitch circle diameter d = m*z   (mm)         | ##
## |     db     | base circle diameter db = d*cos(psi) (mm)    | ##  
## |     da     | addendum circle diameter da = d +2m (mm)     | ##  
## |     rb     | base circle radius rb = db/2 (mm)            | ## 
 
## |     x      | x-coordinate                                  | ##  
## |     y      | y-coordinate                                  | ##

####################################################################

def pitch_circle(m, z):
    """
    Returns pitch circle diameter 
    """   
    d = m * z  
    return d

def base_circle(m, z, psi):
    """
    Returns base circle diameter 
    """     
    db = pitch_circle(m, z) * math.cos(psi)
    return db

def ad_circle(m, z):
    """
    Returns addendum circle diameter 
    """     
    da = pitch_circle(m, z) + 2 * m
    return da

def dedendum_rad(m, z):
    """
    Returns dedendum circle radius 
    """ 
    r = 0.5 * pitch_circle(m, z) - 1.25 * m
    return r

def clearance(m):
    """
    Returns clearance 
    """ 
    c = 0.25 * m
    return c


def inv_a(m, z, psi, r):
    """
    Returns involute angle  
    """ 
    alpha = math.acos(base_circle(m, z, psi) / 2 / r) 
    inv = math.tan(alpha) - alpha
    return inv

def spline1(m, z, psi):
    """
    Returns spline 1 coordinates as an array 
    """ 
    btpstep = (pitch_circle(m, z) - base_circle(m, z, psi)) / 20
    ptastep = (ad_circle(m, z) - pitch_circle(m, z)) / 20

    if 0.5*base_circle(m, z, psi) - dedendum_rad(m, z) > clearance(m):
        r = np.zeros(22)

        for i in range(1,12):       
            r[i] = (base_circle(m, z, psi) / 2) + (i-1) * btpstep
        for i in range(12,22):
            r[i] = (pitch_circle(m, z) / 2) + (i - 11) * ptastep

        spline_1 = np.zeros((22,2))
        for i in range(1,22):   
            inv = inv_a(m, z, psi, r[i])
            spline_1[i,0] = r[i] * math.cos(inv)
            spline_1[i,1] = r[i] * math.sin(inv)       
    
        r[0] = dedendum_rad(m, z) + clearance(m)
        inv = -inv_a(m, z, psi, base_circle(m, z, psi)-r[0])
        spline_1[0,0] = r[0] * math.cos(inv)
        spline_1[0,1] = r[0] * math.sin(inv)         
    else:
        r = np.zeros(21)
        for i in range(11):       
            r[i] = (base_circle(m, z, psi) / 2) + i * btpstep
        for i in range(11,21):
            r[i] = (pitch_circle(m, z) / 2) + (i - 10) * ptastep

        spline_1 = np.zeros((21,2))
        for i in range(21):   
            inv = inv_a(m, z, psi, r[i])
            spline_1[i,0] = r[i] * math.cos(inv)
            spline_1[i,1] = r[i] * math.sin(inv)
    return spline_1


def reflection_ang(m, z, psi):
    """
    Returns reflection angle from x axis
    """ 
    t = math.pi * m / 2             #tooth thickness
    ang = t / pitch_circle(m, z)       #angle (s = r theta - radians)
    r = pitch_circle(m, z) / 2
    r_angle = inv_a(m, z, psi, r) + ang    
    return r_angle

def reflect_xcoords(m, z, psi, x, y):
    """
    Returns reflected x coordinates (other side of tooth) 
    """ 
    grad = math.tan(reflection_ang(m, z, psi))
    xr = ((1 - grad**2) * x + 2 * grad * y) / (1 + grad**2)
    return xr

def reflect_ycoords(m, z, psi, x, y):
    """
    Returns reflected y coordinates (other side of tooth) 
    """    
    grad = math.tan(reflection_ang(m, z, psi))
    yr = ((grad**2 - 1) * y + 2 * grad * x) / (1 + grad**2)
    return yr

def spline2(m, z, psi):
    """
    Returns spline 2 coordinates as an array 
    """ 
    spline_1 = spline1(m, z, psi)
    if 0.5*base_circle(m, z, psi) - dedendum_rad(m, z) > clearance(m):
        n = 22
    else:
        n = 21
    spline_2 = np.zeros((n,2))
    for i in range(n):   
        spline_2[i,0] = reflect_xcoords(m, z, psi, spline_1[i,0], spline_1[i,1])
        spline_2[i,1] = reflect_ycoords(m, z, psi, spline_1[i,0], spline_1[i,1])
    return spline_2

def array_to_tuple(array):
    array_of_tuples = map(tuple, array)
    a_t = tuple(array_of_tuples)
    return  a_t

def dedendum_coords(m,z, psi):
    """
    Returns coordinates to plot dedendum arcs
    [0,0],[0,1]: nearside, arc by spline 1 (x&y)
    [1,0],[1,1]: farside, arc by spline 1 (x&y)
    [2,0],[2,1]: nearside, arc by spline 2 (x&y)
    [3,0],[3,1]: farside, arc by spline 2 (x&y)
    """  
    r = dedendum_rad(m, z)     #dedendum radius
    refl = reflection_ang(m, z, psi)
    coords = np.zeros((4,2))

    coords[0,0] = r
    coords[0,1] = 0.0
    coords[1,0] = r * math.cos(-refl)
    coords[1,1] = r * math.sin(-refl)
    coords[2,0] = reflect_xcoords(m, z, psi, coords[0,0], coords[0,1])
    coords[2,1] = reflect_ycoords(m, z, psi, coords[0,0], coords[0,1])
    coords[3,0] = reflect_xcoords(m, z, psi, coords[1,0], coords[1,1])
    coords[3,1] = reflect_ycoords(m, z, psi, coords[1,0], coords[1,1])
 
    return coords