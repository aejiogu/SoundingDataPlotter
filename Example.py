# -*- coding: utf-8 -*-
"""
Created on Tues Nov 09 18:31:54 2021

@author: Amanze Ejiogu

Name: Example.py
Description: Example script meant to showcase how to use SoundingFunctions
"""
import SoundingFunctions as sf
from datetime import datetime
import matplotlib.pyplot as plt

#Initializing constants
date = datetime(2021,9,1,12) #Sounding from 9/01/2021 at 12z
location = 'IAD' #IAD= NWS Sterling, VA

#Retrieve sounding data
data = sf.pullData(date,location)
print(data)

"""
Here, data is a Pandas dataframe that contains altitude(Alt), potential temperature (theta_k),
relative humidity (RH), temperature (temp), specific humidity (q), atmopsheric 
refractivity (N), and moist adiabatic lapse rate (malr) as columns.
"""

#Plot data with PBL heights and return detected values.
figure, detection_values = sf.detectPBL(date,location)
plt.show(figure)
print(detection_values)

"""
Here, the figure, detection altitudes and gradients of each profile, 
except for Richardson index, are returned. See 'Example.jpg' for the figure.
"""
