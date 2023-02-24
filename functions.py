# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 17:34:56 2023

@author: QCTlab
"""
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt


def fudger (height, width, length, fudge_distance1, fudge_distance2):
    """This function accounts for part of the volume of the scintillators next to the surface which isn't active
    
    INPUT: height, width, length - of the scinitillators
    fudge_distance1 - the dead zone distance in scintillator 1
    fudge_distance2 - the dead zone distance in scintillator 2
    
    OUTPUTS: modified (height, width, length)
        aperture angle range 1
        aperture angle range 2
    """
    
    l = fudge_distance1 
    m = fudge_distance2
    
    height +=  l + m
    width -= 2 * ( l + m ) / 2
    length -= 2 * ( l + m ) / 2
    
    phi = 2 * np.arctan ( width / height ) # gives the aperture angles of the detector in radians    
    print ('The angle of the width aperture is:', phi*180/np.pi, "degrees")
    alpha = 2 * np.arctan ( length / height )                                                      
    print ('The angle of the length aperture is:', alpha*180/np.pi, "degrees")
    
    return height, width, length, phi, alpha


##############################################


def generator_func_arr (argument_array, F0, fudge_distance, height, width, length):
    """ This hunction geneerates prediction values. For Arrays!
    INPUT: F0 - the intensity of the flux
        fudge distance - the amount of dead zone arounf the detector
        height, width, length - geometrical dimemsions of the detector
        argument_array - an array with the feature (angles in this case)
    OUTPUT: an array of the predicted values"""
    
    height, width, length, phi, alpha = fudger(height, width, length, fudge_distance) # transforms the geometrical parameters and gives derived geometrical values
    
    integral_table = np.zeros(len(argument_array))  # table to store the computed vales
    argument_array=(np.pi/180)*argument_array # the argument values in radians
    
    for i in range(len(argument_array)):
        integral_table[i] = integral_of_the_flux (argument_array[i], alpha, phi, height, width, length, F0) # compute the value, based on the given parameters
    
    return integral_table


##############################################


def chassis_tilt (max_shift_ratio, array):
    '''
    This function compensates for the tilt of the chassis of the detector during the measurments. 
    The discrepency in the real and the measured data was shown to be between [0% and 6% ] for the range [0 to 90deg] in our experiment.
    
    INPUT: 
        max_shift_ratio - the ration between the not maximum not tilted and tilted values
        array - the array to be manipulated
        
    OUTPUT: 
        manipulated array
    '''
    
    black_adder = np.linspace(1, max_shift_ratio, array.size)
    array *= black_adder
    
    return array


##############################################

def cos_func_squared (y, x, zenith_angle, clear_sky = True):
    """Function to calculate the modified cosinus law of muons. Uses a transformed angle, check the documetation pdf for more info.
    
    INPUT:
        swwep angle 1 - spehere zenith angle
        sweep angle 2 - spehere equatorial angle
        global zenith angle
        clear_sky - whether we assume buildings around the area we are measuring or not
        
    OUTUT: 
        the modified cosie law
        """
    
    a = 10 # in km  atmosphere height
    r = 6400 # in km  radius of the planet earth
    discrim_angle = 25 # in degrees
    
    if clear_sky == True:
        new_cos = abs (np.cos ( np.arcsin ( np.sin ( zenith_angle + x ) * np.cos( y/2 ) ) ) )
    else:
        if abs((x + zenith_angle)*180/np.pi) <= (90 - discrim_angle) or abs((x + zenith_angle)*180/np.pi) >= (90 + discrim_angle) :  # giving a tolerance of 10deg over the horizon, due to the buildings blocking the detector's "sight" and blocking part of the muons, thus rendering the data incomplete
            new_cos = abs (np.cos ( np.arcsin ( np.sin ( zenith_angle + x ) * np.cos( y/2 ) ) ) )
        else:
            new_cos = 0.00001
            
    
    output = ((a/(-r*new_cos + np.sqrt(a**2 + 2*a*r + (r*new_cos)**2)))**2)  # (1/np.exp(1 + zenith_angle))*np.exp(1 + 0.1*zenith_angle)
    
    return output

##############################################

def integral_of_the_flux (zenith_angle, alpha, phi, height, width, length, F0, clear_sky = True):
    """Here we define the integrand and we integrate to get the flux through the detectors for a given zenith angle.
    
    INPUT: 
        global zenith angle
        aperture angle range 1
        aperture angle range 2
        height, width, length of the scinitllators
        F0 - is the intensity of the flux (check the documentation)
        clear_sky - whether we assume buildings around the area we are measuring or not
        
    OUTPUT:
        the resulting integral (check the documentation)
    """
    
    alim = -alpha/2
    blim = alpha/2
    
    alim_2 = -phi/2
    blim_2 = phi/2
    
    # f = lambda y,x: np.sin(x + np.pi/2) * np.cos(y/2) * (width - (height*abs(np.tan(y)))) * (length - (height*abs(np.tan(x)))) * cos_func_squared (y,x,zenith_angle)
    f = lambda y,x: np.sin(x + np.pi/2) * (width - (height*abs(np.tan(y)))) * (length - (height*abs(np.tan(x)))) * cos_func_squared (y,x,zenith_angle)
    
    if clear_sky == True:
        f = lambda y,x: np.sin(x + np.pi/2) * (width - (height*abs(np.tan(y)))) * (length - (height*abs(np.tan(x)))) * cos_func_squared (y, x, zenith_angle, clear_sky = True)
    else:
        f = lambda y,x: np.sin(x + np.pi/2) * (width - (height*abs(np.tan(y)))) * (length - (height*abs(np.tan(x)))) * cos_func_squared (y, x, zenith_angle, clear_sky = False)
    
    result = F0 * integrate.dblquad(f, alim, blim, alim_2, blim_2)[0]
    
    return result
    

################################################

def namestr(obj, namespace = globals()):
    namespace
    """This function returns the name of an array"""
    return [name for name in namespace if namespace[name] is obj]

##################################################

def plotter (x_labels, y_labels, x_predictions, y_predictions):
    """This funtion plots, two set of data sets, on one axis"""
    
    fig1, (ax1) = plt.subplots(1, 1, figsize = (5, 5), dpi = 400)
    ax1.plot(x_predictions, y_predictions, linewidth=3, label = "fit to detector data", zorder=1)
    
    ax1.scatter(x_labels, y_labels, c='r', linewidths=5, label='Experimental data', zorder=2)
    
    ax1.set_xlabel('$Zenith$ $angle$ ($^{\circ}$)', fontsize=12)
    ax1.set_ylabel('$Absolute$ $count$', fontsize=12)
    # ax1.set_title('$Muon$ $impacts$ $per$ $minute$ $for$ $1x1 cm$ $chip$', fontdict=None, loc='center', pad=None)
    
    ax1.legend()
    
    plt.plot()
    
    return None

#################################

def min_value_parser (threshold, array, sweep_array1, sweep_array2):      #
    """This function takes a 2D array and gives one the values and the indeces of where its values are smaller then "threshold" """
    
    bb = np.flatnonzero(array < threshold)
    
    if bb.size == 0:
        print ("===No values below this threshold===")
    
    for i in bb:
        print (i)
        index_2 = i % array.shape[1]
        index_1 = int((i - index_2) / array.shape[1])
        print ("index_1:",index_1,"index_2:",index_2)
        print ('\n',"indexes:", index_1, index_2, "| RMSE value:", array[index_1, index_2])
        print ("length value:", sweep_array1[index_1], "| Width value:", sweep_array2[index_2])
    
    return None

#########################

def least_sqares_method (y_pred, y_label):
    """takes values of two arrays (of equal size) and computes the least squares metric"""
    sum = 0
    
    for i in range(len(y_pred)):
        dif = (y_pred[i] - y_label[i])**2  
        sum = sum + dif
    
    sum = sum/len(y_pred)
    sum = np.sqrt (sum)
    
    return sum      
#######################

def fitter_func(argument_array, F0, fudge_distance):
    """Just a special case of ===generator_func_arr==="""
    
    height = 1.5
    width = 3   
    length = 21 
    
    return generator_func_arr (F0, fudge_distance, height, width, length, argument_array) #generator_func_single (argument, F0, fudge_distance, height, width, length)

#################