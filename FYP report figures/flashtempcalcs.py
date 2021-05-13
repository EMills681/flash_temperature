""" Flash temperature functions - Functions used to predict Flash Temperature"""
import math
import numpy as np

#######################################################################

### Gear geometry
## |  Variable  | Description                                      | ##
## |------------|--------------------------------------------------| ##
## |     mg     | gear ratio                                       | ##
## |     m      | module (mm)                                      | ##
## |     z      | number of teeth                                  | ##
## |     z1     | number of teeth (pinion)                         | ##
## |    psi     | pressure angle (radians)                         | ##
## |     b      | face width (mm)                                  | ##
## |     d      | pitch circle diameter d = m*z            (mm)    | ##
## |     db     | base circle diameter db = d*cos(psi)     (mm)    | ##  
## |     da     | addendum circle diameter da = d +2m      (mm)    | ##  
## |     rb     | base circle radius rb = db/2             (mm)    | ##
## |     r1     | pitch point radius (pinion) (mm)                 | ##
## |     r      | radius at arbitary point (mm)                    | ##   

## |   step    | Number of steps (odd to include pitch point)      | ##  

### Material properites
## |  Variable  | Description                                      | ##
## |------------|--------------------------------------------------| ##
## |     E     | Young's modulus (MPa, N/mm^2)                     | ##  
## |     v     | Poisson's ratio                                   | ## 
## |    hc     | Heat conductivity (N/s/C)                         | ##  
## |   rhoM    | Density (kg/m3)                                   | ##  
## |    cM     | Specific heat per unit mass (J/kg/C)              | ##   

### Operating conditions
## |  Variable  | Description                                      | ##
## |------------|--------------------------------------------------| ##
## |     N1    | Rotational speed of pinion (rpm)                  | ## 
## |     T1    | Toruqe on pinion (Nm)                             | ##

### Air properites
## |  Variable  | Description                                      | ##
## |------------|--------------------------------------------------| ##
## |  rhoair    | Density of air (kg/m3) @ 20 oC                   | ##  
## |   vair    | Dynamic viscocity of air (kg/m/s) @ 20 oC         | ## 
## |  Cpair    | Specific heat capacity of air (J/kg/C) @ 20 oC    | ##  
## |   kair    | Thermal conductivity of air (W/m/C) @ 20 oC       | ##  

#######################################################################

##Standard Gear Geometry Functions##
def pitch_circle_d(m, z):
    """
    Returns pitch circle diameter (mm)
    """   
    d = m * z  
    return d

def base_circle_d(m, z, psi):
    """
    Returns base circle diameter (mm)
    """     
    db = pitch_circle_d(m, z) * math.cos(psi)
    return db

def ad_circle_d(m, z):
    """
    Returns addendum circle diameter (mm)
    """     
    da = pitch_circle_d(m, z) + 2 * m
    return da

def centre_distance(z1, mg, m):
    """
    Returns centre distance (mm)
    """     
    a = z1 * (1 + 1 / mg) * m / 2
    return a

def load(T1, r):
    """
    Returns Ft: nominal tangential load (N)
    """         
    F = 1000 * T1 / r
    return F

def unit_load(T1, b, r1):
    """
    Returns wt: unit load (N/mm)
    """ 
    ft = load(T1, r1)
    wt = 1000 * ft / b
    return wt

def velocity(d, N1):
    """
    Returns vt: velocity (m/s) (pitch line velocity)
    """ 
    v = d / 2000 * 2 * math.pi * N1 / 60
    return v

##Path of Contact Functions##
def rotation_angles(m, z1, psi, mg):
    """
    rotation_angles[0]: angle of approach (rad)
    rotation_angles[1]: angle of recess (rad)
    rotation_angles[2]: total angle of rotation (rad)
    """ 
    z2 = z1*mg

    ra1 = ad_circle_d(m, z1) / 2
    rb1 = base_circle_d(m, z1, psi) / 2
    r1 = pitch_circle_d(m, z1) / 2
    ra2 = ad_circle_d(m, z2) / 2
    rb2 = base_circle_d(m, z2, psi) / 2
    r2 = pitch_circle_d(m, z2) / 2
    theta_a = psi - math.atan((r1*math.sin(psi)-(math.sqrt(ra2**2 - rb2**2)-r2*math.sin(psi)))/rb1)
    theta_r = math.acos(rb1/ra1) - psi
    theta_tot = theta_a + theta_r
    return theta_a, theta_r, theta_tot

def line_of_action(m, z1, psi, step):
    """
    Assuming gear ratio = 1
    line of action: straight line between -1 and 1 across roll angle, it will be 0 at pitch point.
    line_of_action[0,0]: parameter at start of line of action
    line_of_action[0,1]: angle of rotation at start of line of action (angle of approach) (rad)
    line_of_action[-1,0]: parameter at end of line of action
    line_of_action[-1,1]: angle of rotation at end of line of action (angle of recess) (rad)
    """ 
    ra1 = ad_circle_d(m, z1) / 2

    alpha_E = math.acos(0.5* base_circle_d(m, z1, psi) / ra1) 
    te = math.tan(alpha_E)/math.tan(psi) - 1
    ta = - te

    increment = (-ta + te)/(step-1)
    L = np.zeros((step,2))
    for i in range(step):
        L[i,0] = ta + increment*i
        L[i,1] = math.atan((L[i,0]+1)*math.tan(psi))
    return L

def radius_of_curvature(m, z1, psi, step, mg):
    """
    Calculated using line of action values and BS ISO 6336:20
    radius_of_curvature[i]: respective local relative radius of curvature (mm)
    """ 
    # rho12[i,0]: rhoy1, rho12[i,1]: rhoy2 
    rho12 = np.zeros((step,2))
    rho = np.zeros((step,1))
    # line_of_action[:,0]: linear parameter for line of action
    loa = line_of_action(m, z1, psi,step)[:,0]
    #centre distance
    a = centre_distance(z1, mg, m)
    for i in range(step):
        rho12[i,0] = (1 + loa[i])/(1 + mg) * a * math.sin(psi)
        rho12[i, 1] = (mg - loa[i])/(1 + mg) * a * math.sin(psi)  
        rho[i] = (rho12[i,0] * rho12[i, 1]) / (rho12[i,0] + rho12[i, 1])
    return rho

## Load sharing factor ##
def load_sharing_f(m, z1, psi, step):
    """
    Returns load sharing factor values along line of action
    """ 
    L = line_of_action(m, z1, psi, step)[:,0]
    x = np.zeros(step)
    for i in range(step):       
        x[i] = -0.8 * L[i]**2 / (L[-1]*1.15)**2 + 0.8
    return x 

## Normal unit load ##
def normal_unit_load(psi, T1, b, r1):
    """
    Returns normal unit load, wn: N/mm
    """ 
    wt = unit_load(T1, b, r1)
    wn = wt / math.cos(psi)
    return wn

## Semi-width of Hertzian contact band ##
def reduced_modulus(E,v):
    """
    Returns reduced modulus of Elasticity (Pa)
    """ 
    Er = E / (1-v)**2
    return Er

def Hertz_stress(m, z1, psi, b, step, mg, E, v, T1):
    Er = reduced_modulus(E,v)
    r = pitch_circle_d(m, z1)/2
    Ft = load(T1, r)
    X = load_sharing_f(m, z1, psi, step)
    rel = radius_of_curvature(m, z1, psi, step, mg)
    pnom = np.zeros(step)
    for i in range(step):    
        pnom[i] = math.sqrt(Er/2*math.pi) * math.sqrt(Ft * X[i]/(b * rel[i] * math.cos(psi)))
    return pnom

def Hertz_contact(m, z1, psi, b, step, mg, E, v, T1):
    """
    Returns Hertzian contact band width bh across meshing path (mm)
    """ 
    rel = radius_of_curvature(m, z1, psi, step, mg)
    p_y = Hertz_stress(m, z1, psi, b, step, mg, E, v, T1)            #could add factors iso 6336
    Er = reduced_modulus(E,v)
    bh = np.zeros(step)
    for i in range(step):    
        bh[i] = 4 * rel[i] * p_y[i] / Er
    return bh

## Tangential velocities ## 
def tang_velocities(m, z1, psi, mg, N1, step):
    """
    Returns v1 and v2: velocity of pinion and wheel at each contact position (m/s)
    v12[:,0] = tang_velocities(m, z1, psi, mg, step)[:,0] : velocity of pinion at point y
    v12[:,1] = tang_velocities(m, z1, psi, mg, step)[:,1] : velocity of wheel at point y
    """ 
    z2 = z1*mg

    d1 = pitch_circle_d(m, z1)
    db1 = base_circle_d(m, z1, psi)
    d2 = pitch_circle_d(m, z2)
    db2 = base_circle_d(m, z2, psi)

    a = centre_distance(z1, mg, m)    
    L = line_of_action(m, z1, psi, step)[:,0]

    v12 = np.zeros((step,2))

    for i in range(step):
        y1b1 = 2*(1 + L[i])/(1 + mg) * a * math.sin(psi)
        v12[i, 0] = 2 * math.pi * N1 / 60 * d1 / 2000 * math.sin(psi)* y1b1 / (math.sqrt(d1**2 - db1**2))
        y2b2 = 2*(mg - L[i])/(1 + mg) * a * math.sin(psi)
        v12[i, 1] = 2 * math.pi * N1 / 60 / mg * d2 / 2000 * math.sin(psi)* y2b2 / (math.sqrt(d2**2 - db2**2))
    return v12

## Thermal contact coefficient
def thermal_coefficent(hc, rhoM, cM):
    """
    Returns thermal coefficient from material properties(BS ISO 6336:20)
    """ 
    Bm = (0.001*hc*rhoM*cM)**0.5
    return Bm

## FLASH TEMPERATURE ##
def flash_temp(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu):
    """
    Returns flash temperature, T (oC)
    Array of size 'step' for temperature at each point along line of action
    """ 
    T = np.zeros(step)
    r1 = pitch_circle_d(m, z1) / 2

    x = load_sharing_f(m, z1, psi, step)
    wn = normal_unit_load(psi, T1, b, r1)
    bh = Hertz_contact(m, z1, psi, b, step, mg, E, v, T1)
    v12 = tang_velocities(m, z1, psi, mg, N1, step)
    Bm = thermal_coefficent(hc, rhoM, cM)

    for i in range(step):  
        T[i] = 1.110 * (mu * x[i] * wn / math.sqrt(2*bh[i])) * (abs(v12[i,0] - v12[i,1])/ Bm) / (math.sqrt(v12[i,0]) + math.sqrt(v12[i,1]))
    return T

#### Misc ######

def datum_points(m, z1, psi, step):
    """
    Find datum points to plot between each step
    """ 
    Y = line_of_action(m, z1, psi, step)[:,1]
    rb = base_circle_d(m, z1, psi)/2
    ry =  np.zeros(step)
    datums = np.zeros((step,2))
    for i in range(step):
        ry[i] = rb / math.cos(Y[i])
        inv = math.tan(Y[i]) - Y[i]
        datums[i,0] = ry[i] * math.cos(inv)
        datums[i,1] = ry[i] * math.sin(inv)
    return datums

def temps_and_time(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu):
    """
    Find temp and time between each step (seconds) 
    temps_and_time(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM)[:,0] :temperature between steps
    temps_and_time(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM)[:,1] : time between steps
    """ 
    t =  np.zeros((step-1, 2))
    T = flash_temp(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu)
    Y = line_of_action(m, z1, psi, step)[:,1]
    
    for i in range(step-1):
        t[i,0] = (T[i+1] + T[i]) / 2
        t[i,1] =  30 * (Y[i+1] - Y[i]) /(math.pi * N1)   
    return t
    
def tmax(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu):
    tmax = np.amax(flash_temp(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu))
    return tmax

def rotation_time(m, z1, psi, step, N1):
    """
    remainding rotation time before meshing again
    """ 
    total_t = 60 / N1
    tooth_rotation = (line_of_action(m, z1, psi, step)[-1,1] - line_of_action(m, z1, psi, step)[0,1])
    tooth_t = 30 * tooth_rotation / (math.pi * N1)
    rotation_t = total_t - tooth_t
    return rotation_t

def heat_trans_co(m, z1, N1, rhoair, vair, Cpair, kair):
    """
    heat transfer coefficients at gear tooth surfaces, h mW/mm2/C)
    [0]: top surface
    [1]: gear sides
    [2]: non-meshing gear surfaces
    """     
    Pr = vair * Cpair / kair
    ht = 0.0197 * kair * 3**0.2 * Pr**0.6 * (rhoair * N1 / vair)**0.8 * (0.5 * ad_circle_d(m, z1) * 10**-3)**0.6 * 10**-3
    hs = 0.0197 * kair * 3**0.2 * Pr**0.6 * (rhoair * N1 / vair)**0.8 * (0.5 * pitch_circle_d(m, z1) * 10**-3)**0.6 * 10**-3
    hr = hs / 3
    return ht, hs, hr

def meshing_position(m, z1, psi, step):
    """
    Find datum points to plot between each step
    """ 
    Y = line_of_action(m, z1, psi, step)[:,1]
    rb = base_circle_d(m, z1, psi)/2
    ry =  np.zeros(step)
    datums = np.zeros((step,2))
    for i in range(step):
        ry[i] = rb / math.cos(Y[i])
        inv = math.tan(Y[i]) - Y[i]
        datums[i,0] = ry[i] * math.cos(inv)
        datums[i,1] = datums[i,0] - datums[0,0]
    return datums