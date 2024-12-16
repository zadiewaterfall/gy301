#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:50:04 2024

@author: libbywaterfall
"""

import numpy as np
import matplotlib.pyplot as plt

# SETTING UP THE PARAMETERS
dx = 20 * 1609.34  # 20 miles in meters 
L = 300 * 1609.34  # 300 miles in meters (total LENGTH of grid)
D = 1.12e-8 * 60 * 60 * 24 * 365  # diffusivity rate of crude oil in water in m^2/year
u = 0.8 * 60 * 60 * 24 * 365  # advection velocity in m/year (0.8 is the average velocity of the Loop current in m/s)



# DEFINING VARIABLES AND INITIAL CONCENTRATIONS 
BPD = 60000  # BARRELS of oil released per DAY
VPB = 159  # VOLUME of oil in LITERS per BARREL
DAYS = 87  # DURATION of oil release in DAYS 
LITERS_TO_CUBIC_METERS = 1e-3  # conversion factor from LITERS to CUBIC METERS

# MY CALCULATIONS 
total_oil_liters = BPD * VPB * DAYS  # total VOLUME of oil in LITERS
total_oil_cubic_meters = total_oil_liters * LITERS_TO_CUBIC_METERS  # total VOLUME of oil in CUBIC METERS



#SPACIAL STEPS 
dx = 32186.8  # spatial step in METERS
x = np.arange(0, L, dx)  # grid points
nodes = len(x)


dt = 1.1e-15 # time step in YEARS (further reduced for stability)
totaltime = 14  # total simulation time in YEARS 


# INITIAL CONDITIONS
oil = np.zeros(nodes)  # oil concentration in m^3/m
source = 0  # source at the wellhead (x = 0)
oil[source] = total_oil_liters / (87 * 365 / (dt * 365))  # distribute oil release over time



# CHECK STABILITY
courant = dt * u / dx
s = dt * D / dx**2

import sys #stop code from running if it is not stable
if courant**2 > 2 * s:
    print('Unstable #1') 
    sys.exit() #code will stop running if unstable 



if s + courant / 4 > 0.5:
    print('Unstable #2')   
    sys.exit() #code will stop running if unstable 

# BUILD THE A MATRIX
A = np.zeros((nodes, nodes))
for i in range(1, nodes - 1):
    A[i, i + 1] = s - courant / 2
    A[i, i] = 1 - 2 * s
    A[i, i - 1] = s + courant / 2

# BOUNDARY CONDITIONS
A[0, 0] = 1
A[-1, -1] = 1

# TIME LOOP
fig, ax = plt.subplots(1, 1)
ax.plot(x / 1609.34, oil, '--k', label='Initial')

time = 0
while time < time:
    new_oil = np.dot(A, oil)  # apply diffusion and advection
    oil[:] = new_oil
    oil[source] += total_oil_liters / (87 * 365 / (dt * 365)) * dt  # continue source emission 
    time += dt

# PLOT RESULTS
ax.plot(x / 1609.34, oil, label='Final')
ax.set_xlabel('Distance from Wellhead (miles)')
ax.set_ylabel('Oil Concentration (m^3/m)')
ax.set_title('Diffusion and Advection of Oil from Deepwater Horizon')
ax.legend()
plt.show()
