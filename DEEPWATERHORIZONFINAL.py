#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:35:19 2024

@author: libbywaterfall
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


## INITIAL CONDITIONS
##--------------------
dx = 12614400 # meters (calculated for stability) 
x = np.arange(1, 402336, dx/1000)*1000 # grid from 1 m to 402336 m (*1000 because stabilize )
nodes = len(x) # calculating number of grid nodes 
dt = 0.25/1.8 # time step divided by MAX velocity (for stability)



## AVERAGE VELOCITY OF THE LOOP CURRENT IN THE GULF OF MEXICO IS 0.8 METERS PER SECOND
##------------------------------------------------------------------------------------
loop_current_velocityMS = 0.8  # AVERAGE velocity in METERS/SECOND
loop_current_velocityMY = loop_current_velocityMS * 60 * 60 * 24 * 365 #velocity of the Loop Current in the Gulf of Mexico in METERS/YEAR 



#SCENARIO 1: SPATIALLY VARYING VELOCITIES 
##---------------------------------------
increased_velocityMY = np.ones(nodes) * loop_current_velocityMY # initialize w/ 0.8 METERS/SECOND converted to METERS/YEAR


increased_velocityMS = 1.6 # is some areas of the Loop Current the velocity increased (METERS per SECOND)
increased_velocityMY[(x >= 1e8) & (x <= 1.5e8)] = increased_velocityMS * 60 * 60 * 24 * 365 #1.6 METERS/SECOND converted to METERS/YEAR

increased_velocityMS = 1.8 # velocity at 1.8 METERS/SECOND
increased_velocityMY[(x >= 2e8) & (x <= 2.5e8)] = increased_velocityMS * 60 * 60 * 24 * 365  #1.8 METERS/SECOND converted to METERS/YEAR

increased_velocityMS = 1 # velocity at 1 METER/SECOND
increased_velocityMY[(x >= 3e8) & (x <= 3.5e8)] = increased_velocityMS * 60 * 60 * 24 * 365 #1 METER/SECOND converted to METERS/YEAR


"""
Since the supply in the case of the Deepwater Horizon oil spill was added only
at one location, the location where the connection of the oil pipeline to the 
rig was severed, the initial concentration everywhere else on the grid should 
be zero. The supply of oil is added and then moved along the grid by varying
velocities.
"""


## INITIAL CONDITION FOR CONCENTRATION OF CRUDE OIL  
##-------------------------------------------------
C = np.zeros(nodes) #initialize concentration everywhere as 0 
C[x <= 1] = 60000  # when x is less than or equal to 1 METER, the concentration is 60000 BARRELS per DAY
C[x > 1] = 0  # when x is greater than 1 meter set concentration equal to 0 



## OIL SUPPLY
##------------
R = 60000  # supply rate of crude oil in barrels per day from the severed pipeline
S = np.zeros(nodes)
S[1] = R  # supply rate being added ONLY at the first location 



## STABILITY CHECK 
##----------------
courant = dt * increased_velocityMY / dx  #stability based on timestep(dt), velocity(m/y), and spatial step(dx) 
if np.max(courant) > 1: # if the max courant number is greater than 1 it is unstable 
    print('Unstable #1')  # Print 'Unstable #1' if the simulation is unstable
print(dt) # current time step for reference
print(courant) # courant number for spacial elements for stability 


    
## CREATE THE A MATRIX 
##-----------------------
A = np.zeros((nodes, nodes)) # initialize empty matrix and then index

for i in range (1,nodes): # loop through nodes 
    A[i,i] = 1 - courant[i] # diagonal 
    A[i,i-1] = courant[i] # left diagonal 

A[0,0] = 1 # boundary condition that first node stays the same 

file = np.genfromtxt('Data_for_Deepwater_Horizon_Oil_Spil.txt', skip_header = 1) # load txt file - ignore row 1/the headings for data 
# we only want the numerical data - once it is here it can be explained in the code

concentration = interp1d(file[:,0], file[:,1]) # interp function for txt file data 


totaltime = 8 # years 
time = 1 # years 


while time < totaltime: # run code until totaltime is reached 
    interpolated_concentration = concentration(time) # at current time step call for concentration data 
    new_C = np.dot(A, C) # new concentration(C) w/ regard to flux(A)   
    C[:] = new_C # replaces values in concentration vector C w/ new values
    C += S * dt # incorperates source term(S), scaled by timestep (dt), into concentration vector C
    time += dt # increases the time by the duration of the timestep(dt) 
    C[C > loop_current_velocityMY] =  loop_current_velocityMY #limiting concentration based on velocity 
    

## SCENARIO 2: CONSTANT VELOCITY
##-------------------------------
"""
It is unlikely that in an environment heavily affected by ocean currents there
would be a scenario where the velocity of such current is constant however, 
plotting this against statislly varying velocities as a theorical scenario will 
help us understand how varied velocities impact the spread of crude oil. This 
scenario runs with the velocity as 0.8 m/s which is the average speed of 
The Loop Current in the Gulf of Mexico. (0.8 m/s = 12,614,400 m/y).
"""
constant_velocityMY = np.ones(nodes) * loop_current_velocityMY  # set all velocity values to 0.8 METERS/SECOND, converted to METERS/YEAR


## STABILITY CHECK FOR A CASE WITH A CONSTANT VELOCITY 
##-----------------------------------------------------
courant_constant = dt * constant_velocityMY / dx # find courant number for the constant velocity scenario
A_constant = np.zeros((nodes, nodes)) # create empty matrix A_constant of size (nodes, nodes)

for i in range(1, nodes): # loop through nodes from second node (i=1)
    A_constant[i, i] = 1 - courant_constant[i] # diagonal 
    A_constant[i, i-1] = courant_constant[i]  # left diagonal

A_constant[0, 0] = 1  # boundary condition (no change)

C_constant = np.zeros(nodes) # create concentration vector C_constant w/ zeros
C_constant[x <= 1] = 60000  # setting initial concentration at first node


time = 1 # years 
while time < totaltime: # loop til totaltime
    new_C_constant = np.dot(A_constant, C_constant)  # advection step for constant velocity
    C_constant[:] = new_C_constant # update concentration vector C_constant w/ new values
    C_constant += S * dt # add source term(S) multiplied by time step(dt) to concentration vector C_constant
    time += dt # increases the time by the duration of the timestep(dt) 
    C_constant[C_constant > loop_current_velocityMY] = loop_current_velocityMY  # limiting concentration based on velocity
    





## PLOTTING 
##-----------------------------------------------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 1, figsize=(8, 6), constrained_layout=True)  # create a figure and axis - adjusting the figure

ax.plot(x, C_constant, color='red', linewidth = 2, linestyle='--', label='Oil Movement - Constant Velocity (0.8 m/s)') # plot oil concentration over distance for constant velocity 
ax.plot(x, C, color = 'black', linewidth = 2, label='Oil Movement - Varied Velocities') # plot oil concentration over distance for varied velocities 


## HIGHLIGHT REGIONS OF INCREASED VELOCITIES 
##-------------------------------------------
ax.axvspan(1e8, 1.5e8, color='blue', alpha=0.2, label='Velocity = 1.6 m/s') #shown on the graph in m/s to minimize the amount of visual distractions on the graph
ax.axvspan(2e8, 2.5e8, color='blue', alpha=0.4, label='Velocity = 1.8 m/s') #shown in m/s on graph 
ax.axvspan(3e8, 3.5e8, color='blue', alpha=0.6, label='Velocity = 1 m/s') #shown in m/s on graph 


## SET LABELS, TITLE, AND LEGEND 
##-------------------------------
ax.set_xlabel('Distance from Wellhead (meters)')
ax.set_ylabel('Oil Concentration (barrels/meter)')
ax.set_title('Advection of Oil from the Deepwater Horizon Explosion by the Loop Current')
ax.legend()

## GRAPH!
##--------
plt.show()















