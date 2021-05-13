""" Gear Functions - Functions used to Calculate Flash Temperature"""
#Figures in FYP report

import math
import flashtempcalcs as ftc
import matplotlib.pyplot as plt

####################################
#rpm variations
t1 = ftc.flash_temp(m = 2, z1 = 30, T1 = 10, b = 17.4, N1 = 500, psi =  math.radians(20), step = 61, mg = 1, E = 2600, v = 0.3, hc = 0.230, rhoM = 1410, cM = 1470, mu = 0.29)
t2 = ftc.flash_temp(m = 2, z1 = 30, T1 = 10, b = 17.4, N1 = 1000, psi =  math.radians(20), step = 61, mg = 1, E = 2600, v = 0.3, hc = 0.230, rhoM = 1410, cM = 1470, mu = 0.32)
t3 = ftc.flash_temp(m = 2, z1 = 30, T1 = 10, b = 17.4, N1 = 1500, psi =  math.radians(20), step = 61, mg = 1, E = 2600, v = 0.3, hc = 0.230, rhoM = 1410, cM = 1470, mu = 0.35)
mp = ftc.meshing_position(m=2, z1=30, psi=math.radians(20), step=61)[:,1]
plt.plot(mp,t1,'--', label="500 rpm")
plt.plot(mp,t2,'-.', label="1000 rpm")
plt.plot(mp,t3,'-', label="1500 rpm")
plt.xlim(left = 0)
plt.ylim(bottom = 0)
plt.legend(loc="lower right")
plt.xlabel('Meshing position (mm)')
plt.ylabel('Flash Temperatures ($^\circ$C)')
plt.grid()
plt.rc('font', size = 20)
plt.show()

# #torque variations
t0 = ftc.flash_temp(m = 2, z1 = 30, T1 = 3, b = 17.4, N1 = 1000, psi =  math.radians(20), step = 61, mg = 1, E = 2600, v = 0.3, hc = 0.230, rhoM = 1410, cM = 1470, mu = 0.55)
t1 = ftc.flash_temp(m = 2, z1 = 30, T1 = 5, b = 17.4, N1 = 1000, psi =  math.radians(20), step = 61, mg = 1, E = 2600, v = 0.3, hc = 0.230, rhoM = 1410, cM = 1470, mu = 0.42)
t2 = ftc.flash_temp(m = 2, z1 = 30, T1 = 7, b = 17.4, N1 = 1000, psi =  math.radians(20), step = 61, mg = 1, E = 2600, v = 0.3, hc = 0.230, rhoM = 1410, cM = 1470, mu = 0.35)
t3 = ftc.flash_temp(m = 2, z1 = 30, T1 = 10, b = 17.4, N1 = 1000, psi =  math.radians(20), step = 61, mg = 1, E = 2600, v = 0.3, hc = 0.230, rhoM = 1410, cM = 1470, mu = 0.32)
mp = ftc.meshing_position(m=2, z1=30, psi=math.radians(20), step=61)[:,1]
plt.plot(mp,t0,'--', label="3 N m")
plt.plot(mp,t1,'-.', label="5 N m")
plt.plot(mp,t2,'-', label="7 N m")
plt.plot(mp,t3,':', label="10 N m")
plt.xlim(left = 0)
plt.ylim(bottom = 0)
plt.legend(loc="lower right")
plt.xlabel('Meshing position (mm)')
plt.ylabel('Flash Temperatures ($^\circ$C)')
plt.grid()
plt.rc('font', size = 20)
plt.show()

# #load sharing factor
plt.plot(ftc.meshing_position(m=2, z1=30, psi=math.radians(20), step=61)[:,1], ftc.load_sharing_f(m=2, z1=30, psi=math.radians(2), step = 61))
plt.xlabel('Meshing position (mm)')
plt.ylabel('Load sharing factor')
plt.grid()
plt.xlim(left = 0)
plt.ylim(0,1)
plt.rc('font', size = 20)
plt.show()

#partitions vs continuous temp
surface_temps = ftc.flash_temp(m = 2, z1 = 30, T1 = 10, b = 17.4, N1 = 1000, psi =  math.radians(20), step = 61, mg = 1, E = 3590, v = 0.34, hc = 0.230, rhoM = 1420, cM = 1285, mu=0.32)+20
continu_mesh = a_p = ftc.meshing_position(m=2, z1=30, psi=math.radians(20), step=61)[:,1]
a_t = ftc.temps_and_time(m=2, z1=30, T1=10, b=17.4, N1=1000, psi=math.radians(20), step=11, mg=1, E = 3590, v = 0.34, hc = 0.230, rhoM = 1420, cM = 1285, mu=0.32)[:,0]+20
a_p = ftc.meshing_position(m=2, z1=30, psi=math.radians(20), step=11)[:,1]
plot_temps = a_t[0],a_t[0],a_t[1],a_t[1],a_t[2],a_t[2],a_t[3],a_t[3],a_t[4],a_t[4],a_t[5],a_t[5],a_t[6],a_t[6],a_t[7],a_t[7],a_t[8],a_t[8],a_t[9],a_t[9]
plot_positions= a_p[0],a_p[1],a_p[1],a_p[2],a_p[2],a_p[3],a_p[3],a_p[4],a_p[4],a_p[5],a_p[5],a_p[6],a_p[6],a_p[7],a_p[7],a_p[8],a_p[8],a_p[9],a_p[9],a_p[10]
plt.plot(plot_positions,plot_temps, continu_mesh, surface_temps)
plt.xlabel('Meshing position (mm)')
plt.ylabel('Surface Temperatures ($^\circ$C)')
plt.grid()
plt.xlim(left = 0)
plt.ylim(0, 100)
plt.rc('font', size = 20)
plt.show()

#Mao(2007) comparison
plt.plot(ftc.meshing_position(m=2, z1=30, psi=math.radians(20), step=61)[:,1], ftc.flash_temp(m = 2, z1 = 30, T1 = 10, b = 17.4, N1 = 1000, psi =  math.radians(20), step = 61, mg = 1, E = 2600, v = 0.3, hc = 0.230, rhoM = 1410, cM = 1470, mu = 0.32))
plt.xlabel('Meshing position (mm)')
plt.ylabel('Flash Temperatures ($^\circ$C)')
plt.grid()
plt.xlim(left = 0)
plt.ylim(0, 70)
plt.rc('font', size = 20)
plt.show()

#Hertz contact
plt.plot(ftc.meshing_position(m=2, z1=30, psi=math.radians(20), step=61)[:,1],ftc.Hertz_contact(m=2, z1=30, psi=math.radians(20), b=17.4, step=61, mg=1, E=3590 , v =0.34, T1=10))
plt.xlabel('Meshing position (mm)')
plt.ylabel('Hertzian half width (mm)')
plt.grid()
plt.xlim(left = 0)
plt.ylim(bottom = 0.2)
plt.rc('font', size = 20)
plt.show()


############################################################
# # INPUTS:
# #   gear geometry
# mg = 1
# m = 2
# z1 = 30
# pressure_angle = 20
# psi = math.radians(pressure_angle)
# b = 17.4

# step = 11           # must be an odd number to include pitch point

# #   material properties
# E = 3590                ## Modulus of elasticity (N/mm2 = MPa)
# v = 0.34                  ## Poisson ratio
# hc = 0.230                  ##heat conductivity
# rhoM = 1420                     ## density
# cM = 1285          ##specific heat per unit mass

# #   Operating conditions
# T1 = 8      #(Nm)
# N1 = 1000   #(rpm)
#############################################################
