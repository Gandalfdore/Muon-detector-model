# -*- coding: utf-8 -*-
"""
COMPLETE MODEL FOR MUON FLUX DETECTION

@author: MoonPenguin
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as integrate



plt.style.use('coolplot.mplstyle')

df = pd.read_csv('muon_counter_data.csv')  #,header=1,nrows=4)

datapoints = 10 #how many points to generate for the model (# of the x axis values) 

beta = np.linspace(0,89,datapoints) # defining the points between 0 and 89 degrees (we don't use 90deg because of numerical overflow problems) 



################### Dimension of the scintillators ##################### 

F0=0.498 # the constant of the flux (this is a parameter to be sweeped)

height = 1.5  # all units in "cm"  of the detector, height is the distance between the scintillators
width = 3   
length = 21 


############################## Aperture angles definition ##################

phi = 2*np.arctan(width/height) # gives the aperture angles of the detector in radians    
print ('The angle of the width aperture is:',phi*180/np.pi,"degrees")
alpha = 2*np.arctan(length/height)                                                      
print ('The angle of the length aperture is:',alpha*180/np.pi,"degrees")



################### Correction of the aperture (in case of need) ##############
angle_correction = -50. # IN DEGREES
angle_correction = angle_correction*np.pi/180

alpha = alpha + angle_correction  #the arc of vision of the detector in degrees ALPHA - WIDE SIDE, PHI - NARROW SIDE of the detectors
phi = phi + angle_correction

print (f'Corrected angle of the width aperture is:{phi*180/np.pi} degrees')                                                   
print (f'Corrected angle of the length aperture is:{alpha*180/np.pi} degrees')

############################################



def cos_func_squared (y,x,zenith_angle):
    
    a = 10 # in km  atmosphere height
    r = 6400 # in km  radius of the planet earth
    
    
    # new_cos = abs(np.cos(y/2) * abs(np.cos(zenith_angle + x )))
    new_cos = abs(np.cos(np.arcsin(np.sin(zenith_angle + x ) * np.cos(y/2))))   ###test
    # output = (new_cos**2)   ###test
    # output = 1
    output = (a/(-r*new_cos + np.sqrt(a**2 + 2*a*r + (r*new_cos)**2)))**2
    return output


##########################################



def integral_of_the_flux (zenith_angle, alpha, phi, height, width, length,F0):
    """Here we define the integrand and we integrate to get the flux through the detectors for a given zenith angle."""
    
    alim = -alpha/2
    blim = alpha/2
    
    alim_2 = -phi/2
    blim_2 = phi/2
        
    f = lambda y,x: np.sin(x + np.pi/2) * np.cos(y/2) * (width - (height*abs(np.tan(y)))) * (length - (height*abs(np.tan(x)))) * cos_func_squared (y,x,zenith_angle)
    
    result = F0 * integrate.dblquad(f, alim, blim, alim_2, blim_2)[0]
    
    return result#, Angular_integral[0], Area_integral[0]

###############################




"""This part of the code contains the parameters to calculate the muon rate for a thin 1x1 cm chip"""
  

height2 = 0.0001  # all units in "cm"  of the detector, height is the distance between the scintillators
width2 = 1
length2 = 1 

phi2 = 2*np.arctan(width2/height2) # gives the aperture angles of the detector in radians    #126.9
print ('The angle of the width aperture for a chip is:',phi2*180/np.pi,"degrees")
alpha2 = 2*np.arctan(length2/height2)                                                        #171.8
print ('The angle of the length aperture for a chip is:',alpha2*180/np.pi,"degrees")

##############################



######################### Plotting code ##############

integral_table1 = np.linspace(0,1,datapoints)
integral_table2 = np.linspace(0,1,datapoints)

beta1=(np.pi/180)*beta # convert the beta degrees to radians

for i in range(len(beta)):
    integral_table1[i] = integral_of_the_flux (beta1[i], alpha, phi, height, width, length,F0)
    integral_table2[i] = integral_of_the_flux (beta1[i], alpha2, phi2, height2, width2, length2,F0)
    print("Progress:", int(100*(i+1)/len(beta)),"%")
    
scaling_factor =  width*length  # normalization factor for the detector's surface area (basically its own area seen directly from above)

fig1, (ax1) = plt.subplots(1, 1, figsize = (5, 5), dpi = 400)
ax1.errorbar(df['angle'], df['average_cpm_coincidental']/scaling_factor,yerr=df['average_cpm_coincidental_error']/scaling_factor,xerr=0.3, fmt=".",linewidth=4, label='Experimental data')
ax1.plot(beta,integral_table1/scaling_factor,linewidth=3,label = "fit to detector data" )
ax1.plot(beta,integral_table2,linewidth=3,label = "chip " )


ax1.set_xlabel('$Zenith$ $angle$ ($^{\circ}$)', fontsize=12)
ax1.set_ylabel('$Absolute$ $count$', fontsize=12)
ax1.set_title('$Muon$ $impacts$ $per$ $minute$ $for$ $1x1 cm$ $chip$', fontdict=None, loc='center', pad=None)

# ax1.set_ylim([0,2.1])
ax1.legend()

plt.plot()
