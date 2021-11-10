# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 21:41:15 2021

@author: Amanze Ejiogu

Name: OzoneDataPlotter.py
Description:This portion of the program will plot datasets for comparison. It
will attempt to derive the planetary boundary layer height from those days
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import scipy.signal
import metpy.calc as mpcalc
from metpy.units import units
#from scipy import signal

#Enter the location of the processed data filename list
filename_loc = "./Filenames.csv"
file_list = pd.read_csv(filename_loc, header = None)

plt.rcParams['axes.spines.top'] = True
plt.rcParams['axes.spines.right'] = True


#Other constants
dalr = 9.8 #degrees C/km
min_height = 600 #Minimum altitude of allowed PBL detection
max_height = 3000 #Maximum altitude of allowed PBL detection
sort_num = 10  #Number of max/min gradients to keep track of
height_delta = 50 #Maximum altitude difference between PBL heights
dz = 15 #The height, in meters, above and below which should be used to sample the atmopsphere
#default_filt = [0.25,0.5,0.25] #The filter used for the smoothing
ep = 15 #Using an epsilon here to determine the number of meters to average for stability
stab_vals = ['U','S','CU', 'N/A']
stab_colors = ['red', 'green', 'yellow', 'gray']
    
#Setting up the legend
unstable = mpatches.Patch(color='red', label='Unstable')
stable = mpatches.Patch(color='green', label='Stable')
con_stab = mpatches.Patch(color='yellow', label='Conditionally Unstable')
non = mpatches.Patch(color='gray', label='N/A')

#For converting the infinite values to nan
def inf2nan(array):
    array[abs(array) == np.inf] = np.nan
    
    return(array)

#A bit hacky, but returns the N number of max/min values in an array w/o nans
def sortlist(order,array, N):        
    if order == "max":
        sorted_indexes = (-array).argsort()
        
    elif order == "min":
        sorted_indexes = array.argsort()

    return(sorted_indexes[:N])

#Trying to calculate the filter to be used for sigfilter
def binomcoeffs(n): 
    return (np.poly1d([0.5, 0.5])**n).coeffs

default_filt = binomcoeffs(2).tolist() #The filter used for the smoothing


#1D signal filter using filtfit. Defaults to a 1-2-1 filter
def sigfilter(filt , data):
    filtered_signal = scipy.signal.filtfilt(filt,1,data)
    
    return filtered_signal

stability_vals = []

#For now, we skip files with insufficient data
skip_list = [11]

colnames = ['Date_Time_UTC', 'Location', 'Theta_Height', 'RH_Height', 'q_Height'
            , 'N_Height', 'dTheta', 'dRH', 'dq', 'dN']

vals = []

for j in range(0, len(file_list)):
    #j = 40
    if j in skip_list:
        continue
    
    
    current_datanum = j
    data_loc = "./"+ file_list[0][current_datanum]  
    time = file_list[1][current_datanum]
    loc = file_list[2][current_datanum]
    
    #def OzoneDataPlotter(data_loc):
    data = pd.read_csv(data_loc)
    ozone = data['O3']
    alt = data['Alt']
    theta_k = data['Theta_K']
    rel_h = data['RH']
    temp = data['Temp']
    malr = data['malr']
    
    q = data["q"]
    N = data["N"]
#    u = data['u'].values * units(data.units['u'])
#    v =  units.Quantity(data['v'], 'm/s')
    #Calculate Ri
#    Ri = mpcalc.gradient_richardson_number(units.Quantity(alt, "meters")
#                                           , units.Quantity(theta_k, "kelvin")
#                                           , u, v)
    
    
    
    #Missing data requires that we replace missing data
    #Ri = inf2nan(data["Ri"])

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
                
    #If there is no consensus, default to the regular method
    """
    #For the other measures, the minimum slope is the PBL
    rh_PBL_index = np.nanargmin(rel_h_height[min_index:max_index])

    #Trying *any* large perterbation
    q_PBL_index = np.nanargmax(abs(q_height[min_index:max_index]))
    N_PBL_index = np.nanargmax(abs(N_height[min_index:max_index]))
    """
        
    #Convert the temp indexes back into the regular indexes
    PBL_lines = [theta_PBL_indexes[temp_indexes['theta']],
                 rh_PBL_indexes[temp_indexes['rh']],
                 q_PBL_indexes[temp_indexes['q']],
                 N_PBL_indexes[temp_indexes['N']]]
    
    
    """
    #For potential temperature, the maximum slope is the PBL
    theta_PBL_index = np.nanargmax(theta_height[min_index:max_index])
    
    #For the other measures, the minimum slope is the PBL
    rh_PBL_index = np.nanargmin(rel_h_height[min_index:max_index])
    
    #Trying *any* large perterbation
    q_PBL_index = np.nanargmax(abs(q_height[min_index:max_index]))
    N_PBL_index = np.nanargmax(abs(N_height[min_index:max_index]))

    Remember, all of these index values are relative to the *subscript* so
    add the index to the min_index to get the actual value.
    
    
    #Determine where the marks should be placed
    PBL_lines = np.asarray([theta_PBL_index,rh_PBL_index,q_PBL_index,N_PBL_index]) 
                           + min_index
            
    """
    
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
    
    """
    #Here, we compress these values so that the final display step is quicker
    quick
    for i in range(0,len(sections)):
        
        current_rating = sections['rating'][i]
    """
    
    ##Here the main chart is plotted
    plt.figure(figsize= [10,7.5])
    #plt.figure(figsize= [16,12])
    comp_fig, (ax1,ax2,ax3,ax4,ax5, ax6) = plt.subplots(1,6)
    comp_fig.suptitle(loc + ", " + time + " UTC")
    
    
    #Plotting O3
    ax1.plot(theta_k,alt,color="red")
    ax1.set_xlabel("O3 (ppbv)", fontsize = 10)#, color="red")
    ax1.set_ylabel("Altitude (m)")
    ax1.set_title("Est. PBL Height:", fontsize = "small")
    ax1.set_ylim([0,10000])
    #ax1.set_title("Theta K & Mixing Ratio vs Altitude ")
    
    #Adding Potential Temperature Values
    ax2.plot(theta_k,alt,color="blue")
    ax2.set_xlabel("θ (K)", fontsize = 10)#, color="blue")
    ax2.get_yaxis().set_visible(False)
    ax2.axhline(y = alt[PBL_lines[0]], xmin = 0, xmax = 1, color = "blue", 
                linestyle = "--")
    ax2.set_title(str(alt[PBL_lines[0]]) + " m", fontsize = "small", color = "blue")
    ax2.set_ylim([0,10000])

    #ax2.set_ylabel("Altitude (m)")
    #ax1.set_xlim([250,400])
    #ax2.set_title("Theta K & Mixing Ratio vs Altitude ")
    
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
    #break
    #ax6.legend()
    
    #Adding Richardson Number
    """
    ax6.plot(Ri,alt, color = 'deeppink')
    ax6.set_xlabel('Ri', fontsize = 10)
    ax6.get_yaxis().set_visible(False)
    ax6.set_ylim([0,10000])
    """
    #plt.tight_layout()    
    plt.savefig("./PBL Detection/" + loc + time.replace(":","")
                +  ".jpg", dpi = 400,bbox_inches='tight') 
    plt.close()
    
    sections.to_csv('./Stability/' + loc + time.replace(':','') + '.csv' )
    
    vals.append([time, loc, alt[PBL_lines[0]], alt[PBL_lines[1]], alt[PBL_lines[2]],
                 alt[PBL_lines[3]], theta_height[PBL_lines[0]], rel_h_height[PBL_lines[1]],
                 q_height[PBL_lines[2]], N_height[PBL_lines[3]]])

vals_df = pd.DataFrame(vals, columns = colnames)
vals_df.to_csv("./PBL Detection/PBL Height Detection Values.csv")
#help_me = np.asarray(stability_vals).reshape(46,2)
#plt.show()
    