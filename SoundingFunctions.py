# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 16:34:08 2021

@author: Amanze Ejiogu, Rahim Karmara

Name: SoundingFunctions.py
Description: Fully enclosed set of functions that can be used to process UWYO 
soundings using the methodology set in "Estimation of atmospheric mixing layer 
height from radiosonde data. Atmospheric Measurement Techniques, Wang, X. Y., 
& Wang, K. C. (2014)."

Notes: To import, just type "import SoundingFunctions". Then, you can use the functions
however you please.
"""
#All the functions that will be used 
#from datetime import datetime, timedelta #To log time and date stuff
import pandas as pd #Used for dataframes, which hold data
import matplotlib.pyplot as plt #Plotting
import numpy as np #Used for calculations and arrays
import metpy.calc as mpcalc #Used to calulate variables missing from NWS
from metpy.units import units #Adds Units
from siphon.simplewebservice.wyoming import WyomingUpperAir #Pulls data from UWYO site
import scipy.signal #Used to process the profile gradients
import matplotlib.patches as mpatches

#Ignore this, just to ensure plots look good
plt.rcParams['axes.spines.top'] = True
plt.rcParams['axes.spines.right'] = True


##Constants
dalr = 9.8 #degrees C/km
min_height = 600 #Minimum altitude of allowed PBL detection
max_height = 3000 #Maximum altitude of allowed PBL detection
sort_num = 10  #Number of max/min gradients to keep track of
height_delta = 50 #Maximum altitude difference between PBL heights
#default_filt = [0.25,0.5,0.25] #The filter used for the smoothing
ep = 15 #Using an epsilon here to determine the number of meters to average for stability
colnames = ['Theta_Height', 'RH_Height', 'q_Height'
            , 'N_Height', 'dTheta', 'dRH', 'dq', 'dN']
stab_vals = ['U','S','CU', 'N/A']
stab_colors = ['red', 'green', 'yellow', 'gray']
    
#Setting up the legend
unstable = mpatches.Patch(color='red', label='Unstable')
stable = mpatches.Patch(color='green', label='Stable')
con_stab = mpatches.Patch(color='yellow', label='Conditionally Unstable')
non = mpatches.Patch(color='gray', label='N/A')

def pullData(date, station):
    """ 
    Pulls data from the UWYO servers and returns what the other functions
    need

    Parameters
    ----------
    date : datetime
        Date of sounding (with time) in datetime format. Hours can only be 
        0z or 12z.
    station : string
        3-5 digit ground station identifier.

    Returns
    -------
    final_data : Pandas dataframe
        Processed data that contains altitude(Alt), potential temperature (theta_k),
        relative humidity (RH), temperature (temp), specific humidity (q), 
        atmopsheric refractivity (N), and moist adiabatic lapse rate (malr) as columns.

    """
    try: #Try to get data from server
        raw_data = WyomingUpperAir.request_data(date, station)
    
    except:
        pass #If there's an error, pass it onto
    
    #These parts pull the data that we want from the sounding
    h = raw_data['height'].values * units(raw_data.units['height']) 
    P = raw_data['pressure'].values *units(raw_data.units['pressure']) 
    T_degC = raw_data['temperature'].values * units(raw_data.units['temperature']) 
    Td_degC = raw_data['dewpoint'].values * units(raw_data.units['dewpoint']) 
    
    #We'll return to this part soon (pulls lat/long and wind direction data)
    u = raw_data['u_wind'].values * units(raw_data.units['u_wind'])
    v = raw_data['v_wind'].values * units(raw_data.units['v_wind']) 
    lat = raw_data['latitude'].values * units(raw_data.units['latitude'])
    lon = raw_data['longitude'].values * units(raw_data.units['longitude']) 
    
    #Converting the temeratures to Kelvin
    T_kelvin = T_degC.to('kelvin') 
    Td_kelvin = Td_degC.to('kelvin')
    
    #Calculating the other variables that we need
    Theta = mpcalc.potential_temperature(P, T_kelvin)
    RH = mpcalc.relative_humidity_from_dewpoint(T_kelvin, Td_kelvin)*100 
    q = mpcalc.specific_humidity_from_dewpoint(P, Td_kelvin)
    mr = mpcalc.mixing_ratio_from_specific_humidity(q)
    e = mpcalc.vapor_pressure(P, mr)
    malr = mpcalc.moist_lapse(P,T_degC[0])
    Ri = mpcalc.gradient_richardson_number(h, Theta, u, v)

    
    N = (77.6*(P/T_kelvin)).magnitude + ((3.73E5)*(e/(T_kelvin**2))).magnitude
    
    data = {'Alt': h,
            'Theta_K': Theta,
            'RH' : RH, 
            'Temp' : T_kelvin,
            'q' : q,
            'N' : N,
            'Ri':Ri,
            'malr':malr}
    #Returns data as a Pandas dataframe
    final_data = pd.DataFrame(data)
    
    return final_data


def inf2nan(array):
    """
    For internal use. For converting the infinite values in an array to nan.

    Parameters
    ----------
    array : NumPy array
        A NumPy array of data.

    Returns
    -------
    fixed_array : NumPy array
        Array with all np.inf replaced with nans.

    """
    array[abs(array) == np.inf] = np.nan
    fixed_array = array

    
    return fixed_array

def sortlist(order,array, N):    
    """
    For internal use. Returns the N number of max/min values in an array w/o nans.
    Parameters
    ----------
    order : string
        "min" or "max". The ordering of the sorted array (decending or ascending).
    array : NumPy array
        A 1-D NumPy array of data.
    N : int
        The number of sorted values to be returned.

    Returns
    -------
    sorted_array: NumPy array
        A 1-D NumPy array of sorted data.

    """
    if order == "max":
        sorted_indexes = (-array).argsort()
        
    elif order == "min":
        sorted_indexes = array.argsort()

    return sorted_indexes[:N] 

#Trying to calculate the filter to be used for sigfilter
def binomcoeffs(n): 
    """
    For internal use. Calculate a binomial signal filter

    Parameters
    ----------
    n : int
        Order of binomial filter.

    Returns
    -------
    filt: list
        Signal filter of desired order.

    """
    filt = (np.poly1d([0.5, 0.5])**n).coeffs
    return filt

default_filt = binomcoeffs(2).tolist() #The filter used for the smoothing


#1D signal filter using filtfit. 
def sigfilter(filt , data):
    """
    For internal use. Filters the data using a filter.

    Parameters
    ----------
    filt : list
        Filter to use on the data.
    data : NumPy array
        1-D NumPy array of filter.

    Returns
    -------
    filtered_signal : NumPy array
        Filtered 1-D signal.

    """
    filtered_signal = scipy.signal.filtfilt(filt,1,data)
    
    return filtered_signal

def detectPBL(date, loc, plot= True, return_data= True):
    """
    Detects the PBL height using the methodology adapted from Wang et. al, 2014.
    Cutoff is set to 10,000 meters by default. See "Constants" section to adjust.
    Parameters
    ----------
    date: datetime
        Datetime of the sounding. Must be 0 or 12 UTC.
        
    loc: string
        3-5 digit ground station identifier.
    
    plot: bool
        Default: True. Whether or not the data should be plotted and returned.
    
    return_data: bool
        Default: True. Whether or not the function should return detected values.
    

    Returns
    -------
    comp_fig : Matplotlib figure object
        Only returns if plot = True. Plotted data figure.
    vals: Dataframe
        Only return_data = True. Detection altitudes and values for all data.

    """
    if date.hour != 0 and date.hour != 12:
       raise ValueError('Error: Hour value must be 0z or 12z.')
       
    time = date.strftime('%m/%d/%Y, %H:%M')
    
    data = pullData(date, loc)
    alt = data["Alt"]
    theta_k = data["Theta_K"]
    rel_h = data["RH"]
    temp = data["Temp"]
    q = data["q"]
    N = data["N"]
    malr = data['malr']
    
    #Missing data requires that we replace missing data
    Ri = inf2nan(data["Ri"])

    ##Detecting PBL

    min_index = 0
    max_index = 0
    
    for i in range(0, len(alt)):
        current_alt = alt[i]
        
        if current_alt >= min_height and min_index == 0:
            min_index = i 
            continue
        
        if current_alt >= max_height:
            max_index = i
            break
        
    #Calculate gradients wrt altitudes
    theta_height = inf2nan(np.diff(theta_k)/np.diff(alt/1000))
    rel_h_height = inf2nan(np.diff(rel_h)/np.diff(alt/1000))
    q_height = inf2nan(np.diff(q)/np.diff(alt/1000))
    N_height = inf2nan(np.diff(N)/np.diff(alt/1000))
    
    
    
    #Ri_height = inf2nan(np.diff(Ri)/np.diff(alt/1000))

    
    #Now we put this 1D signal through a filter
    theta_height =  sigfilter(default_filt,theta_height)
    rel_h_height =  sigfilter(default_filt,rel_h_height)
    q_height =  sigfilter(default_filt,q_height)
    N_height =  sigfilter(default_filt,N_height)
    
    #Ri_height = sigfilter(default_filt, Ri_height)
    
    
    """
    Changing this up to better reflect the methodology used in the paper
    
    """
    
    theta_PBL_indexes = sortlist("max", theta_height[min_index:max_index], sort_num) + min_index
    rh_PBL_indexes = sortlist("min", rel_h_height[min_index:max_index], sort_num) + min_index
    
    q_PBL_indexes = sortlist("min", q_height[min_index:max_index], sort_num) + min_index
    N_PBL_indexes = sortlist("min", N_height[min_index:max_index], sort_num) + min_index
    
    #Ri_PBL_indexes = sortlist()
    
    
    
    #Here, we run through the theta_height array to check if any of the values are within
    #the height difference set above
    
    pot_heights = pd.DataFrame({'rh': alt[rh_PBL_indexes].tolist(),
                                'q': alt[q_PBL_indexes].tolist(),
                                'N': alt[N_PBL_indexes].tolist()})
    temp_indexes = []
    
    for i in range(0,len(theta_PBL_indexes)):
        #print(loc + time)
        current_theta_height = alt[theta_PBL_indexes[i]]
        
        #Basically, we subtract the current height and find the minimum value
        #of each column. If 2 or more of them are within range, then the indexes
        #are collected
        height_diff = abs(pot_heights - current_theta_height)
        min_diff = np.asarray(height_diff.min())
        
        if np.count_nonzero(min_diff <= height_delta) >= 2 :
            
            #NOTE: This is the index relative to the PBL indexes of before 
            #Use these values to get the raw indexes from the previous arrays
            temp_indexes = height_diff.idxmin()
            temp_indexes['theta'] = i
            
            break
        else:
            #Here we test if the three categories on their own are within range
            current_rh_height = pot_heights.loc[i,"rh"]
            height_diff = abs(pot_heights - current_rh_height)
            rh_diff = np.asarray(height_diff.min())
            
            if np.count_nonzero(rh_diff <= height_delta) == 3:
                temp_indexes = height_diff.idxmin()
                temp_indexes['theta'] = i
                #break
        
    #Convert the temp indexes back into the regular indexes
    PBL_lines = [theta_PBL_indexes[temp_indexes['theta']],
                 rh_PBL_indexes[temp_indexes['rh']],
                 q_PBL_indexes[temp_indexes['q']],
                 N_PBL_indexes[temp_indexes['N']]]
    

    ##Determining stability
    lapse_rate = -inf2nan(np.diff(temp)/np.diff(alt/1000))
    stability_vals = []
    for i in range(0, len(alt)):
        
        #For each altitude, establish the index range that is at least epsilon 
        #wide
        plus = 1
        
        if i == 0: #Different rules for the first index; only consider values above
            minus = 0
            while abs(alt[i + plus] - alt[i]) < ep:
                plus += 1
                
        elif i == len(alt) - 1: #Different rules for the last index; only consider values below
            plus = 0
            while abs(alt[i] - alt[i - minus]) < ep:
                minus += 1
                
        else:
            minus = 1
            while abs(alt[i + plus] - alt[i - minus]) < (ep*2):
                plus += 1
                minus += 1
                
                if (i + plus) >= len(alt) - 1 or (i - minus) <= 0:
                    break
                
        #Take the average of all the values within this range
        avg_malr = np.mean(inf2nan(malr[i-minus:i+plus]))
        avg_lapse_rate = np.mean(lapse_rate[i-minus:i+plus])
        
        #If the average lapse rate is greater than the average malr...
        if avg_lapse_rate >= avg_malr:
            
            #...but is lower than the dalr, then this altitude is conditionally unstable
            if avg_lapse_rate <= dalr:
                rating = "CU"

                
            #...but is greater than the dalr, then this altitude is stable
            else:
                rating = "S"                
                
        #If the average lapse rate is less than the average malr
        else:
            rating = "U"
        
        #Trying to figure out the bounds of the highlighted sections
        
        if i == 0:
            lower = 0
            upper = np.nanmean(alt[i:i+1])

            
        elif i == len(alt):
            upper = alt[i]
            lower = np.nanmean(alt[i-1:i])

        else:
            lower = np.nanmean(alt[i-1:i])
            upper = np.nanmean(alt[i:i+1])

        
        stability_vals.append([alt[i], rating, lower, upper])
        
        if np.nan in stability_vals:
            stability_vals[2] = "N/A"

        
    sections = pd.DataFrame(stability_vals, columns= ['center','rating', 'lower', 'upper'] )

    vals = ([ alt[PBL_lines[0]], alt[PBL_lines[1]], alt[PBL_lines[2]],
                 alt[PBL_lines[3]], theta_height[PBL_lines[0]], rel_h_height[PBL_lines[1]],
                 q_height[PBL_lines[2]], N_height[PBL_lines[3]]])    
    
    vals_df = pd.DataFrame([vals], columns= colnames)
    
    
    if plot == True:
        ##Here the main chart is plotted
        plt.figure(figsize= [10,7.5])
        comp_fig, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(1,6)
        comp_fig.suptitle(loc + ", " + time + " UTC ")
        
        #Plotting Richardson Number
        ax1.set_ylim([0,max_height])
        ax1.plot(Ri,alt, color='deeppink')
        ax1.set_xlabel('Ri', fontsize = 10)
        ax1.set_title("Est. PBL Height:", fontsize = "small")
        ax1.set_ylabel("Altitude (m)")
        ax1.set_ylim([0,10000])

        #Adding Potential Temperature Values
        ax2.plot(theta_k,alt,color="blue")
        ax2.set_xlabel("θ (K)", fontsize = 10)#, color="blue")
        ax2.get_yaxis().set_visible(False)
        ax2.axhline(y = alt[PBL_lines[0]], xmin = 0, xmax = 1, color = "blue", 
                    linestyle = "--")
        ax2.set_title(str(alt[PBL_lines[0]]) + " m", fontsize = "small", color = "blue")
        ax2.set_ylim([0,10000])
        
        #Adding Relative Humidity
        ax3.plot(rel_h,alt, color="green")
        ax3.set_xlabel("RH (%)", fontsize = 10)
        ax3.get_yaxis().set_visible(False)
        ax3.axhline(y = alt[PBL_lines[1]], xmin = 0, xmax = 1, color = "green", 
                    linestyle = "--")
        ax3.set_title(str(alt[PBL_lines[1]]) + " m", fontsize = "small", color = "green")
        ax3.set_ylim([0,10000])
        
        #Adding Specific Humidity
        ax4.plot(q,alt, color="purple")
        ax4.set_xlabel("SH (g/kg)", fontsize = 10)
        ax4.get_yaxis().set_visible(False)
        ax4.axhline(y = alt[PBL_lines[2]], xmin = 0, xmax = 1, color = "purple", 
                    linestyle = "--")
        ax4.set_title(str(alt[PBL_lines[2]]) + " m", fontsize = "small", color = "purple")
        ax4.set_ylim([0,10000])
        
        #Adding Atmospheric Refractivity
        ax5.plot(N,alt, color="goldenrod")
        ax5.set_xlabel("N", fontsize = 10)
        ax5.get_yaxis().set_visible(False)
        ax5.axhline(y = alt[PBL_lines[3]], xmin = 0, xmax = 1, color = "goldenrod", 
                    linestyle = "--")
        ax5.set_title(str(alt[PBL_lines[3]]) + " m", fontsize = "small", color = "goldenrod")
        ax5.set_ylim([0,10000])
        
        #Plotting stability values
        ax6.set_xticks([])
        ax6.get_yaxis().set_visible(False)
        ax6.set_xlabel("Stability", fontsize = 10)
        ax6.set_ylim([0,10000])
        ax6.set_xlim([0,1])
        ax6.legend(handles= [stable, con_stab, unstable, non],  prop={'size': 6}, 
                   loc = 'upper left',  bbox_to_anchor=(1, 1))

    
        for i in range(0,len(sections)):
            #Figure out the color
            ind = stab_vals.index(sections['rating'][i])
            color = stab_colors[ind]        
            if sections['lower'][i] == np.nan or sections['upper'][i] == np.nan:
                continue
                
            ax6.axhspan(sections['lower'][i], sections['upper'][i],
                        color= color)
            
        return comp_fig, vals
    