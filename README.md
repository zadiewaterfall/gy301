#Deepwater Horizon Oil Spill Simulation:

"""This model looks at an Earth Science problem related to the advection of ocean currents in the Gulf of Mexico 
that transport crude oil from the Deepwater Horizon oil spill. The Deepwater Horizon oil spill occurred on 
April 10th, in 2010. The explosion of the rig led to the release of 60,000 barrels of oil per day into the Gulf of Mexico. 
The main ocean current that runs through the Gulf of Mexico is called the Loop Current and it moves at approximately 
0.8 meters per second, though it is not constant. The velocity of the Loop Current varies between 0.6 and 1.8 meters per second. 
The simulation models the movement of oil under the influence of the Loop Current using both spatially varying and constant velocity scenarios. 
The severed pipeline connection to the Deepwater Horizon rig was finally plugged after 87 days of 60,000 barrels per day being emptied into the Gulf.
"""

#Quick Aspects:
- Models oil advection with spatially varying velocities.
- Includes a theoretical scenario with constant velocity for comparison.
- Visualizes oil concentration over distance from the wellhead.
- Highlights areas of increased velocity for spatial analysis.


##Python Packages Needed:
- Python 3.x
- NumPy
- Matplotlib
- SciPy

##Usage:
- Place the provided data file Data_for_Deepwater_Horizon_Oil_Spil.txt in the same directory as the script.
- Run the Python script: DEEPWATERHORIZONFINAL.py
- The script generates a plot comparing oil movement under spatially varying and constant velocity scenarios.

##Assumptions in this Model:
- Initial oil concentration is zero everywhere except at the wellhead.
- Oil is supplied at a constant rate of 60,000 barrels per day.
- The velocity of the Loop Current is modeled both as spatially varying and constant for comparison.

##Plot Features:
- Oil Movement (Constant Velocity): Represented by a red dashed line.
- Oil Movement (Varied Velocities): Represented by a black solid line.
- Highlighted Regions: Indicate areas of increased velocity with varying transparency.

#Output:
- Plot showing distance from the wellhead (in meters) on the x-axis.
- Oil concentration (in barrels per meter) on the y-axis.
- A legend indicating velocity conditions and highlighted regions.

#Author:
Zadie Waterfall 














