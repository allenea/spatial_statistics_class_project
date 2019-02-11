#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 17:37:14 2019

@author: allenea
Author: Eric Allen, University of Delaware, allenea@udel.edu

LAST MODIFIED 2/7/2019 at 2PM

Completed for GEOL659 (spatial statistics, F18) project comparing model evaluation results of traditional
methods with the radical approach of using interpolated data for verification to give a bigger N in the analysis.


This program combines several programs developed in my GEOL_Project directory and an intermidate step in the
interpolation-for-verification program which involved:

    1- Gathering Verification Data and Metadata (NJMET, DelDOT, DEOS, NCEI LCD, NDBC, (optional Mesowest))
    
    2- Process, Quality Control, Adjust winds from anemometer heights to 10m (need to add in Monin-Obukhov Stability Parameter)
    
    3- Extract case study data by case and by observing station
    
    4- Averaging routine (by station) - Forward avg 30 mins so that at analysis time only data leading up to that time is considered.
    
    5- Trim hour average to the case study window
    
    6- Map (cartopy) locations of observations being used
    
    7- Plot data in time-series and do analysis (line,dot, and box)... Find outliers or bad stations and remove and/or document.
    
    8- Wrote a program for statistical functions for analysis
    
    9- Created a OO class to parse namelist wps files to extract the information to be able to recreate the grids for each domain
        and create a cartopy map that can be used to plot data
        
    10- Created a wrf model validation program to perform statistical analysis on the performance of the model in two ways
        first: by case study, second: by independent variable (in this case data assimilation type). Creates line graphs to visualize how the models compare.
            output files store information that can be extracted and put in a table or in the future use that data for more analaysis.
            **** I want to add more statistics like Willmott's Index of Agreement or the Nash-Sutcliffe statistic.... But I need to have a mean of observations
                for each timestep and that's hard to do with how I have it coded now and data stored by station.I'd need to iterate through twice.
                Wonder if there is a better way. 
              
    11- Created this Kriging Analysis program that creates an interpolated grid using the verification data and extracts random locations within the domain
        (the domain which happens to be a subset of the planned WRF domain but for this project using an old one it's actually bigger...since I didn't have the
        time or server spave to rerun the models before the spring semester although it should just be ready to run.). Creates plots for interpolated data and
         variance plus csv file with data.
        Output data in trim_hr_avg file format except it's all combined into one file.
        
    12- Created a program to tell about the spatial distribution of data including the validation data set, and on a sample of 10,000 simulations.
        Plus this program gets its own analysis on the N sample size

Combined with research that involved:
    
    1 - Gathering DEOS data (20 stations) and CMLF data (water and atmoshere)
    
    2 - Process, Quality Control* (different than above), Adjust winds from anemometer heights to 10m (need to add in Monin-Obukhov Stability Parameter)
    
    3 - All data to little_r format at original sensor heights (added to elevation for WRF). CMLF sensor height missing. Units converted as needed. 
        There is also a data file not used with 10m heights. Both of these datasets can be used in the above parts 3 - 7. Not used for validation.
    
    4 - Program to help in the identifying of case studies when we have CMLF data 
        
    5 - Configuring WRF runs and getting initial runs (last year - all of which need to be redone). Made adjustments based on that to improve the model
        and have it work for the coupled ocean/atmosphere COAWST model still being configured.
        
    6 - Set up and connected the data assimilation routines and set parameterizations in OBSGRID and WRF
    
    7 - Created a program to store data from little_r to OBS files for OBSGRID.
    
    8 - NCL scripts to plot data. 
    
    9 - Script to analyze all WRF/OBSGRID text files for errors and warnings and get statistics (numbers) on assimilation which can be easily stored in a spreadsheet.
    
    10 - Program to plot radar data of sea breeze
    
    11 - Moved model outputs to local computer from FARBER to run above #10. These runs are useless and bad data was assimilated. 
        Time on all DEOS data was off because of an innapropraite time conversion based off information given about the data. It was not in fact in local time. Those runs should perform worse
        in analysis until model can be re-run on all 3 cases.
     
TO-DO
    - Re-run WRF model on the 3 cases with the new configuration
    - Expand to 5-10 cases (including 1-2 null sea breeze cases)
    - Redo analysis
    - Compare statistics of RU-WRF and other operational models to evaluate the performance of WRF.
    - Re-run and compare analysis to the other runs.
    - PAPER - END WORK ON THIS - Interpolation will not be used unless it yeilds interesting results.
    
    - Select whether using marine based observations is beneficial. If they are add NDBC buoy data to the list of available observations plus
        some NJ Mesonet stations to fill the other side of the bay for COAWST.
        
    - Implement what works best in the coupled ocean/atmosphere COAWST model and run with that assimilation setup and without any assimilation.
    - Statistical Analysis on COAWST runs compared to all above statistic results.
    
    - Acknowledge the other research that shows how it can be improve (coastline, different sensitivity tests, starting the runs at different time of the day).
        Someone else can try to improve that stuff in the model.
    
    - Can this be the masters. Otherwise I want to stop doing the modeling stuff ASAP and look at SB Convection and Radar Analysis/Prediction
    - Modeling SB Convection Events identified by Dan Moore.
    
ASSETS:
    - Dan Moore Automatic SB detection algorithm

    
?Problems? With this Program:
    
    Coded to use Lambert Projections when taking i,j to lat/lon or lat/lon to i,j, fix for other projections.
    
    Is PyKrige-UK actually taking latitude and longitude into account the right way for trend surface analysis? I think so....
    
    PyKrige doesn't have geographic coordinates for UK - Euclidean should be "okay" for the purposes of GEOL659.
    
    I suspect some bugs in PyKrige that still need ironing out based on the documentation and what I see in their github code.
    
    I believe that UK.execute("grid",long_grid,lat_grid) should be the full grid and not just a slice but it crashes my computer.......
           PyKrige has poor documentation and support. Their example uses linspace to create the list of longitudes and latitudes.
           grid_lon = np.linspace(0.0, 360.0, 1); grid_lat = np.linspace(-90.0, 90.0, 1)

    Choose Better Colormaps - https://matplotlib.org/users/colormaps.html
    
    Map the random lat/lon locations eventually.
    Intermediate script to store in format for new_wrf_model_analysis or plots
    
"""

#SPECIAL IMPORTS
import pykrige as pk
import cartopy
import wrf

# BASIC IMPORTS
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import datetime
import os
import sys
import calendar
import random



# FUNCTIONS
# =============================================================================
#
from Randomness_Analysis_Part1 import data_randomness
import namelist_plot as nplt
#
# =============================================================================
# Make Interpolation Grid
def make_interp_grid(namelist_wps_filename,plot_domain_number):
    """
    Use namelist file to recreate the grid and store it. 
    returns ny by nx array for Latitude and Longitude at every gridpoint. 
    returns the domain information in wps and the plt_idx or specific grid being used.
    
    """
    wps = nplt.wps_info(namelist_wps_filename)
    plt_idx = wps.plot_domain_number(plot_domain_number)
    wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()
    
    corner_x, corner_y = nplt.wps_info.reproject_corners(corner_lon_full[plt_idx,:], corner_lat_full[plt_idx,:], wpsproj, latlonproj)
    
    nplt.wps_info.print_info()

    nx = wps.e_we[plt_idx]   # x direction #nx = 256 #% Number of grid points in the x direction
    ny = wps.e_sn[plt_idx]   # y direction #ny = 256 #%*2; % Number of grid points in the y direction
    print()
    n = int(ny * nx)
    grid_y = np.empty(n)
    grid_x = np.empty(n)
    
    print("nx: ",nx," by ny: ",ny)
    count1 = 0
    for i in range(ny):
        for j in range(nx): # x
            [grid_y[count1],grid_x[count1]] = wrf.xy_to_ll_proj(j, i, meta=False,  map_proj=1, truelat1=wps.truelat1, truelat2=wps.truelat1, stand_lon=wps.stand_lon,\
                    ref_lat=wps.ref_lat, ref_lon=wps.ref_lon, pole_lat=wps.pole_lat, pole_lon=wps.pole_lon, known_x=wps.ref_x, known_y=wps.ref_y, dx=wps.dx[plt_idx], dy=wps.dy[plt_idx]) #Lat/Lon
            count1 +=1
            
    # Reshape the lat/lon
    newLon = np.reshape(grid_x,[ny,nx]); #ROW BY COLUMN... ROW = Constant Lat 
    newLat = np.reshape(grid_y,[ny,nx]); # Row By column.... Row = Changing Lon.

    return newLon,newLat,wps,plt_idx
# =============================================================================


# =============================================================================
# Universal Kriging 
def Ukriging(xcord,ycord,zvalue,lonGrid,latGrid,interpolation_type,plot_variogram=False):
    """
    Includes Trend Surface Analysis, but is currently lacking the feature to handle geographic coordinates
    instead using euclidean. Haven't explored if this is the best variogram model or if it is the best drift term.
    
    Uses the top row and column to recreate the grid.
    
    RETURNS: kriged field, variance field, and the kriging object (for statistics)
    
    class pykrige.uk.UniversalKriging(x, y, z, variogram_model=’linear’, variogram_parameters=None, variogram_function=None,
                        nlags=6, weight=False, anisotropy_scaling=1.0,
                        anisotropy_angle=0.0, drift_terms=None,
                        point_drift=None, external_drift=None, external_drift_x=None, external_drift_y=None,
                        specified_drift=None, functional_drift=None, verbose=False, enable_plotting=False)
    
    
    VARIOGRAM MODEL - VARIOGRAM PARAMETERS
    linear - {'slope': slope, 'nugget': nugget}
    power - {'scale': scale, 'exponent': exponent, 'nugget': nugget}
    gaussian - {'sill': s, 'range': r, 'nugget': n}      OR      {'psill': p, 'range': r, 'nugget':n}
    spherical - {'sill': s, 'range': r, 'nugget': n}     OR       {'psill': p, 'range': r, 'nugget':n}
    exponential - {'sill': s, 'range': r, 'nugget': n}   OR       {'psill': p, 'range': r, 'nugget':n}
    hole-effect - {'sill': s, 'range': r, 'nugget': n}
    """
    
    # Create the ordinary kriging object. Required inputs are the X-coordinates of
    # the data points, the Y-coordinates of the data points, and the Z-values of the
    # data points. Variogram is handled as in the ordinary kriging case.
    # drift_terms is a list of the drift terms to include; currently supported terms
    # are 'regional_linear', 'point_log', and 'external_Z'. Refer to
    # UniversalKriging.__doc__ for more information.
    UK = pk.UniversalKriging(xcord, ycord, zvalue, variogram_model=interpolation_type,
                drift_terms=['regional_linear'],verbose=False, enable_plotting=plot_variogram)
    # Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
    # grid of points, on a masked rectangular grid of points, or with arbitrary points.
    # (See UniversalKriging.__doc__ for more information.)
    z_krige, ss1 = UK.execute('grid', lonGrid[0,:], latGrid[:,0])

    return z_krige, ss1, UK
# =============================================================================


# =============================================================================
# Original Kriging 
def Okriging(xcord,ycord,zvalue,lonGrid,latGrid,interpolation_type,plot_variogram=False):
    """
    NO Trend Surface Analysis, but it is capable of handling geographic coordinates
    instead using euclidean. Haven't explored if this is the best variogram model or if it is the best drift term.
    
    Uses the top row and column to recreate the grid.
    
    (list - float) xcord, ycord, zvalue
    (array - ny by nx) lonGrid, latGrid
    (str) interpolation_type
    (bool) plot_variogram
    
    RETURNS: kriged field, variance field, and the kriging object (for statistics)
    """
    # Create ordinary kriging object:
    OK = pk.OrdinaryKriging(xcord,ycord , zvalue, variogram_model=interpolation_type, verbose=False,
                         enable_plotting=plot_variogram, coordinates_type='geographic')
    
    # Execute on grid:
    z_krige, ss1 = OK.execute('grid', lonGrid[0,:], latGrid[:,0])
    
    return z_krige, ss1, OK
# =============================================================================


# =============================================================================
# Get analysis window
def get_analysis_window(case_time, analysis_start_hour_UTC, analysis_length_hrs):
    """
    Calculates the analysis window. Doesn't go multiple months or years.
    
    (str) Case_Time: i.e. '2014-06-04_06:00'
    (int) analysis_start_hour_UTC:  12
    (int) analysis_length_hrs: i.e. 18
    
    RETURNS datetime tuple for start and end of analysis
    """
    
    dtuple = datetime.datetime.strptime(case_time, "%Y-%m-%d_%H:%M")
    yrCH = dtuple.year
    monCH = dtuple.month
    dayCH = dtuple.day
    
    
    analysis_start = datetime.datetime(dtuple.year, dtuple.month, dtuple.day,analysis_start_hour_UTC)
    # CHANGES DAY
    if analysis_start_hour_UTC + analysis_length_hrs > 23:
        # END HOUR
        analysis_end_hour_UTC = (analysis_start_hour_UTC + analysis_length_hrs)%24
        
        ## NUMBER OF ADDITIONAL DAYS
        analysis_end_day_UTC = (analysis_start_hour_UTC + analysis_length_hrs)//24  

        if dtuple.month == 1 or dtuple.month == 3 or  dtuple.month == 5 or  dtuple.month == 7 or  dtuple.month == 8 or  dtuple.month == 10:
            
            if dtuple.day + analysis_end_day_UTC > 32 and (dtuple.day + analysis_end_day_UTC) <= 60:
                monCH = dtuple.month +1
                dayCH =  (dtuple.day + analysis_end_day_UTC) - 32
                
            elif dtuple.day + analysis_end_day_UTC <= 32:
                dayCH =  dtuple.day + analysis_end_day_UTC 
            else:
                print("ANALYSIS WINDOW PROBLEM")
                sys.exit(0)
                
        elif dtuple.month == 2 and calendar.isleap(dtuple.year) == False:
            #and dtuple.day == 28
            if dtuple.day + analysis_end_day_UTC > 28 and (dtuple.day + analysis_end_day_UTC) <= 60:
                monCH = dtuple.month +1
                dayCH =  (dtuple.day + analysis_end_day_UTC) - 28
                
            elif dtuple.day + analysis_end_day_UTC <= 28:
                dayCH =  dtuple.day + analysis_end_day_UTC 
            else:
                print("ANALYSIS WINDOW PROBLEM")
                sys.exit(0)

        elif dtuple.month == 2 and calendar.isleap(dtuple.year) == True:
            # and dtuple.day == 29
            if dtuple.day + analysis_end_day_UTC > 29 and (dtuple.day + analysis_end_day_UTC) <= 61:
                monCH = dtuple.month +1
                dayCH =  (dtuple.day + analysis_end_day_UTC) - 29
                
            elif dtuple.day + analysis_end_day_UTC <= 29:
                dayCH =  dtuple.day + analysis_end_day_UTC 
            else:
                print("ANALYSIS WINDOW PROBLEM")
                sys.exit(0)
                
        elif dtuple.month == 4 or dtuple.month == 6 or dtuple.month == 9 or dtuple.month ==11:
            if dtuple.day + analysis_end_day_UTC > 30 and (dtuple.day + analysis_end_day_UTC) <= 60:
                monCH = dtuple.month +1
                dayCH =  (dtuple.day + analysis_end_day_UTC) - 30
                
            elif dtuple.day + analysis_end_day_UTC <= 30:
                dayCH =  dtuple.day + analysis_end_day_UTC 
            else:
                print("ANALYSIS WINDOW PROBLEM")
                sys.exit(0)
                
        elif dtuple.month == 12:
            if dtuple.day + analysis_end_day_UTC > 30 and (dtuple.day + analysis_end_day_UTC) <= 62:
                yrCH = dtuple.year + 1
                monCH = 1
                dayCH =  (dtuple.day + analysis_end_day_UTC) - 30
            elif dtuple.day + analysis_end_day_UTC <= 30:
                dayCH =  dtuple.day + analysis_end_day_UTC 
            else:
                print("ANALYSIS WINDOW PROBLEM")
                sys.exit(0)
                
                
        analysis_end = datetime.datetime(yrCH, monCH, dayCH, analysis_end_hour_UTC)

    else:
        analysis_end_hour_UTC = analysis_start_hour_UTC + analysis_length_hrs
        analysis_end = datetime.datetime(dtuple.year, dtuple.month, dtuple.day, analysis_end_hour_UTC)  

    return analysis_start, analysis_end  

# =============================================================================

def get_value_interp(lon,lat,kriged,wps,plt_idx):
    """
    Takes a lat/lon combo from the list of random locations and extracts the data from that point.
    
    (float) Latitude
    (float) Longitude
    (namelist_plot.wps_info) wps
    (int) plt_idx
    
    RETURNS (float) value
    """
    [lat_x,lon_y]= wrf.ll_to_xy_proj(lat, lon, meta=False,  map_proj=1, truelat1=wps.truelat1, truelat2=wps.truelat1, stand_lon=wps.stand_lon,\
                               ref_lat=wps.ref_lat, ref_lon=wps.ref_lon, pole_lat=wps.pole_lat, pole_lon=wps.pole_lon, known_x=wps.ref_x, known_y=wps.ref_y, dx=wps.dx[plt_idx], dy=wps.dy[plt_idx])
    # kriged = NY, NX
    value = kriged[lon_y,lat_x]
    return value

# =============================================================================

def randomx(number_validate,min_lon,max_long):
    """  
    Generate a list of random longitudes.
    
    (int) Number of random validation points
    (float) Longitude Minimum
    (float) Longitude Maximum
    
    RETURNS (list) Longitude
    """
    
    lstx = []
    for i in range(number_validate):
        x = random.uniform(min_lon,max_long)
        lstx.append(x)
    return lstx

# =============================================================================

def randomy(number_validate,min_lat,max_lat):
    """  
    Generate a list of random latitudes
    
    (int) Number of random validation points
    (float) Latitude Minimum
    (float) Latitude Maximum
    RETURNS (list) Latitude
    """
    
    lsty = []
    for i in range(number_validate):
        x = random.uniform(min_lat,max_lat)
        lsty.append(x)
    return lsty

# =============================================================================
# INTERPOLATION OPTIONS - PyKrige
def get_interpolation_option(option):
    """  
    Get the interpolation type that will be used in this program
    
    (int) Option Number
    RETURNS (str) option
    """
    if option == 1: return'linear'
    elif option == 2: return'power'
    elif option == 3: return'gaussian'
    elif option == 4: return'spherical'
    elif option == 5: return'exponential'
    elif option == 6: return'hole-effect'
    else:
        print(option," - Invalid Option: Try Again \n\nValid Options: Option 1: Linear, "+\
              "Option 2: Power, Option 3: Gaussian, Option 4: Spherical, "+\
              "Option 5: Exponential, Option 6: Hole-Effect")

        option = input("Select Another Interpolation Type (INTEGER 1 - 6):")
        try:
            option = int(option)
            return option
        except:
            try:
                option = int(input("Select Another Interpolation Type (INTEGER 1 - 6):"))
                return option
            except:
                print("Second Invalid Interpolation Type Attempt.... Exiting progam.")
                sys.exit(0)
        return get_interpolation_option(option)
# =============================================================================
#####                           END SET FUNCTIONS                          ####
# =============================================================================

# =============================================================================
# =============================================================================
# =============================================================================
#####                         MODIFY AS NEEDED                            #####
# =============================================================================
# =============================================================================
# =============================================================================

def fmt_data_path(case, subdirlistIN,subdirlistOUT):
    """
    SETS FILE FORMAT FOR CASES AND THEIR VARIATIONS
    Files should be arranged and named so you can simply iterate through.
    
    
    # DIRECTORY MAKER: Build Your Own Naming Mechanism
    
    ROOTDIR =  Macintosh HD⁩ ▸ ⁨Users⁩ ▸ NAME ▸ ⁨Documents⁩ ▸ MAIN_DIR ▸
                /Users/myusername/Documents/MAIN_DIR
    PATH = Sub_DIR ▸ Sub_2_Dir ▸ Sub_3_Dir  ....
                /Verification_Data/interpolation_pre/10m
    case[:10] = specific case
                /2014-06-04
    filename  = filename with file type (.csv)
                /all_hr_avg_trim.csv
    
    For Example: 
    # INPUT DATA FILE EXAMPLE
    /Users/myusername/Documents/MAIN_DIR/Verification_Data/interpolation_pre/10m/2014-06-04/all_hr_avg_trim.csv
    /Users/myusername/Documents/MAIN_DIR/Verification_Data/interpolation_pre/10m/2014-06-08/all_hr_avg_trim.csv
    /Users/myusername/Documents/MAIN_DIR/Verification_Data/interpolation_pre/10m/2015-08-14/all_hr_avg_trim.csv
    
    # OUTFILE FILE EXAMPLE
    /Users/myusername/Documents/MAIN_DIR/Verification_Data/interpolation_post/10m/2014-06-04/interpolated_verification_data.csv
    /Users/myusername/Documents/MAIN_DIR/Verification_Data/interpolation_post/10m/2014-06-08/interpolated_verification_data.csv
    /Users/myusername/Documents/MAIN_DIR/Verification_Data/interpolation_post/10m/2015-08-14/interpolated_verification_data.csv
    
    #FIGURE DIRECTORY EXAMPLE
    /Users/myusername/Documents/MAIN_DIR/Python_Programs/GEOL_Spatial_Stats/figs/
    /Users/myusername/Documents/MAIN_DIR/Python_Programs/GEOL_Spatial_Stats/figs/
    /Users/myusername/Documents/MAIN_DIR/Python_Programs/GEOL_Spatial_Stats/figs/
    """
    
    #ROOT DIRS
    data_dir = os.path.abspath('../..')
    fig_dir = os.getcwd()
    
    
    filenameIN = "/all_hr_avg_trim.csv"
    filenameOUT = "/interpolated_verification_data.csv"
    fig_dir_name = "figs"
    
    # SET INPUT DATA PATH
    pathIN = "/"
    for subdir_in in subdirlistIN:
        pathIN = pathIN +subdir_in + "/"
    instring = data_dir+pathIN+case[:10]+filenameIN
    
    # CHECK TO MAKE SURE THE PATH EXISTS
    if not os.path.exists(instring): print("INVALID DATA PATH"); sys.exit(0)

    ## SET AND CREATE OUTPUT DATA FILE AND PATH
    pathOUT = "/"
    for subdir_out in subdirlistOUT:
        pathOUT = pathOUT +subdir_out + "/"
    outdir = data_dir+pathOUT+case[:10]
    outstring = outdir+filenameOUT
    
    #Checks to see if path exits, if not create it
    if not os.path.exists(outdir): os.makedirs(outdir)

    # CREATE DIRECTORY FOR FIGURES CREATED IN THIS PROGRAM
    fig_path =fig_dir+"/"+fig_dir_name+"/"
    if not os.path.exists(fig_path): os.makedirs(fig_path)

    return instring, outstring, fig_path
# =============================================================================
# =============================================================================
####                        USER DEFINED VARIABLES                         ####
# =============================================================================
# =============================================================================

# THE USER SHOULD MODIFY THIS SECTION TO MEET THEIR REQUIREMENTS

# ADMINISTRATIVE - os
sub1dirIN = "Verification_Data"; sub2dirIN = "interpolation_pre"; sub3dirIN = "10m";
subdirlistIN = [sub1dirIN,sub2dirIN,sub3dirIN]
sub1dirOUT = "Verification_Data"; sub2dirOUT = "interpolation_post"; sub3dirOUT = "10m";
subdirlistOUT = [sub1dirOUT,sub2dirOUT,sub3dirOUT]

# CASE STUDY/ANALYSIS
casestudy_time = ['2014-06-04_06:00','2014-06-08_06:00','2015-08-14_06:00']

analysis_start_hour_UTC = 10 
analysis_length_hrs = 18

timesteps_per_hour = 2
run_hours = 24
preset_time_len = (timesteps_per_hour * run_hours) + 1


# INTERPOLATION OPTIONS
interpolation_type = get_interpolation_option(3)    #"gaussian"  -- listed above
my_coordinates_type='geographic'                    # 'euclidean' ## GEO IS INVALID IN UNIVERSAL KRIGING - Just leave as geographic I guess
plot_variogram = False                              #Toggle on/Off based on if you want a variogram figure

# Number of interpolated observations being extracted and saved
number_validate = 250

## GRID TO INTERPOLATE
interp_namelist_str= "namelist.wps"
interp_domain_number = 1

#PLOTTING GRID - If different
plot_namelist_str='namelist.wps_MASTER' #UNUSED
real_plot = 3 #UNUSED
  
# =============================================================================
####                               CONSTANTS                               ####
# =============================================================================

# VARIABLES - Be careful - You will have to make changes in the main program if you change this
variables= ['Wind_Speed (m/s)','Wind_Direction (deg)','Air_Temperature (K)','Dewpoint_Temperature (K)','Relative Humidity (%)','Pressure (Pa)']  

sec = 60 # 60 seconds in a minute

HEADER_INTP = ["ID_String","DATE"]+variables+["Latitude","Longitude"]

# =============================================================================
# =============================================================================
####                            MAIN PROGRAM                               ####
# =============================================================================
# =============================================================================          
#%% Initial Step - Create grid and get random locations(lat/lon) - takes a few seconds
lon_grid, lat_grid,wpsINTP,plt_idxINTP = make_interp_grid(interp_namelist_str,interp_domain_number)
valid_lon = randomx(number_validate,np.min(lon_grid[0,:]),np.max(lon_grid[0,:]))
valid_lat = randomy(number_validate,np.min(lat_grid[:,0]),np.max(lat_grid[:,0]))

#%% START OF INTERPOLATION
# Calculate randomness of the random data in the domain.
OGMean, OGMeanNearest,OGStandardDeviation, OGSTDNearest = data_randomness(valid_lat,valid_lon,geo_fmt="degrees")
 
print("\n","Spatial Interpolation: ", interpolation_type)

for case_time in casestudy_time:
    firstPass = True # keeps from repeatedly filling the first 3 columns in the interpolated data csv
    
    # Admin - Set Input file, output file and figure directory paths
    data_file,outstring,fig_dir = fmt_data_path(case_time,subdirlistIN,subdirlistOUT)

    analysis_start, analysis_end  = get_analysis_window(case_time, analysis_start_hour_UTC, analysis_length_hrs)
    start_dte = datetime.datetime(int(case_time[0:4]),int(case_time[5:7]),int(case_time[8:10]),int(case_time[11:13]), int(case_time[14:16]))
    end_dte = datetime.datetime(int(case_time[0:4]),int(case_time[5:7]),int(case_time[8:10])+1,int(case_time[11:13]), int(case_time[14:16]))
        
    storeMatrix = np.zeros(((len(valid_lon)*int(preset_time_len)),len(HEADER_INTP))).astype(object)
    
    for variable in variables:
        i_loc = 0
        print("\n",case_time,"\t",variable)
        timeZeus = []
        # Store Statistics
        Q1 = []; Q2 = []; cR = []

        # Read In Data   
        data = pd.read_csv(data_file, low_memory=False)
        data.columns.tolist()
        data = data.mask(data == " ",other = np.nan)
        data = data.mask(data == "",other = np.nan)
        
       # Load In Data
        ID_String = np.array(data['ID_String'])
        WSpd= np.array(data['Wind_Speed (m/s)'])
        WDir = np.array(data['Wind_Direction (deg)'])
        temp = np.array(data['Air_Temperature (K)']) 
        dewpt = np.array(data['Dewpoint_Temperature (K)']) 
        rh = np.array(data['Relative Humidity (%)']) 
        press = np.array(data['Pressure (Pa)']) 
        lat = np.array(data['Latitude']) 
        long = np.array(data['Longitude']) 
        
        #elev = np.array(data['Elevation (m)']) 
        #name = np.array(data['Name_string']) 
        #fm = np.array(data['FM_string']) 
        #source = np.array(data['Source_string']) 
        year =np.array(data['YEAR']) 
        month =np.array(data['MONTH']) 
        day = np.array(data['DAY']) 
        hour =np.array(data['HOUR']) 
        minute = np.array(data['MINUTE']) 
        
        Date = [0] * len(year)
        for i in range(len(year)):
            Date[i] = datetime.datetime(int(year[i]),int(month[i]),int(day[i]),int(hour[i]),int(minute[i]))
    
        del minute, hour,day,month,year   
                
        timesteps = list(set(Date))
        timesteps.sort()
        

        for i in range(len(timesteps)):
            
            #store data for interpolation
            pass_lat = []; pass_lon = []; pass_value = []
            
            if timesteps[i] >= analysis_start and timesteps[i] <= analysis_end:
                #print("Time: ", timesteps[i])
                step_Date = list(np.where(Date == np.datetime64(timesteps[i]))[0])
                for j in step_Date:
                    
                    # PER VARIABLE - EXTRACT FOR ANALYSIS  (Verification Data)
                    if variable == 'Wind_Speed (m/s)':
                        if np.isnan(WSpd[j]) == False:
                            pass_lat.append(lat[j])
                            pass_lon.append(long[j])
                            pass_value.append(WSpd[j]) 
                            o_str = "WSPD"
                            cmap1 = matplotlib.cm.viridis
                            maxvalue= 10; minvalue = 0
                    elif variable == 'Wind_Direction (deg)': 
                        if np.isnan(WDir[j]) == False:
                            pass_lat.append(lat[j])
                            pass_lon.append(long[j])
                            pass_value.append(WDir[j])   
                            o_str = "WDIR"
                            cmap1 = matplotlib.cm.cubehelix # DONT EVER USE INTERPLOTED WIND DIRECTION!!
                            maxvalue = 360; minvalue = 0
                    elif variable == 'Air_Temperature (K)':
                        if np.isnan(temp[j]) == False:
                            pass_lat.append(lat[j])
                            pass_lon.append(long[j])
                            pass_value.append(temp[j]) 
                            o_str = "2T"
                            cmap1 = matplotlib.cm.YlOrRd     
                            maxvalue = 320; minvalue = 285
                    elif variable == 'Dewpoint_Temperature (K)':
                        if np.isnan(dewpt[j]) == False:
                            pass_lat.append(lat[j])
                            pass_lon.append(long[j])
                            pass_value.append(dewpt[j])  
                            o_str = "2Td"
                            cmap1 = matplotlib.cm.YlOrRd
                            maxvalue = 320; minvalue = 285
                    elif variable == "Relative Humidity (%)":
                        if np.isnan(rh[j]) == False:
                            pass_lat.append(lat[j])
                            pass_lon.append(long[j])
                            pass_value.append(rh[j]) 
                            o_str = "RH"
                            cmap1 = matplotlib.cm.cool
                            maxvalue = 100; minvalue = 0
                    elif variable == "Pressure (Pa)":
                        if np.isnan(press[j]) == False:
                            pass_lat.append(lat[j])
                            pass_lon.append(long[j])
                            pass_value.append(press[j]) 
                            o_str = "PRESS"
                            cmap1 = matplotlib.cm.rainbow_r
                            maxvalue = 95000; minvalue = 110000

                    else: 
                        print ("INVALID VARIABLE OPTION")
                        sys.exit(0)
                    
            if len(pass_lat) != len(pass_value):
                print ("DATA NOT THE SAME LENGTH")
                sys.exit(0)
                
            if len(pass_lon) != 0:
                #z_krige,ss1,Kr = Okriging(pass_lon,pass_lat,pass_value,lon_grid,lat_grid,interpolation_type,plot_variogram)
                z_krige,ss1,Kr = Ukriging(pass_lon,pass_lat,pass_value,lon_grid,lat_grid,interpolation_type,plot_variogram)
    
                if Kr.get_statistics()[0] > 1.0:
                    print("----- VALUE ERROR ",timesteps[i].strftime('%-m-%-d-%Y %H:%M')," -----")
                    Kr.print_statistics()
                    print("--------------------------------------------\n")

                    # THIS IS NEEDED TO CONTINUE TO FILL THE SPREADSHEET EVEN WHEN WE DON'T WANT TO USE THE DATA FOR THIS TIMESTEP
                    z_krige[:,:] = np.nan; ss1[:,:] = np.nan
                    
                else:
                    #Store Interpolation Statistics
                    timeZeus.append(timesteps[i].strftime('%-m-%-d-%Y %H:%M'))
                    Q1.append(Kr.get_statistics()[0])
                    Q2.append(Kr.get_statistics()[1])
                    cR.append(Kr.get_statistics()[2])
                    
                    ## GET/STORE WPS Data, Get Plot Domain, Calculate Domain Info For Plotting
                    #wps2 = nplt.wps_info(file=plot_namelist_str)
                    #plt_idx2 = wps2.plot_domain_number(real_plot)
                    # define the coordinate system that the grid lons and grid lats are on
                    #lambert = cartopy.crs.LambertConformal(central_longitude=wps2.ref_lon, central_latitude=wps2.ref_lat,\
                    #                        false_easting=wps2.dx[plt_idx2], false_northing=wps2.dy[plt_idx2], secant_latitudes=None, standard_parallels=None,\
                    #                        globe=None)
                    
                    ## SET UP PLOT 1 - Interpolated Data Values
                    fig1, ax = plt.subplots(1, 1, figsize=(10,10), subplot_kw={'projection': cartopy.crs.PlateCarree()})
                    cbar_ax = fig1.add_axes([0, 0, 0.1, 0.1])
         
                    ## Map Features
                    ax.coastlines('10m','black')
                    ax.add_feature(cartopy.feature.STATES.with_scale('10m'))
                    title = str(interpolation_type+' Interpolation on '+timesteps[i].strftime('%-m-%-d-%Y  %H')+"Z").upper()
                    ax.set_title(title,fontsize=16)
                    ax = nplt.wps_info.set_lambert_ticks(ax,1.0,1.0,14,14)

                    fig1.canvas.draw()

                    #Plot Data
                    arr2 = np.array(pass_lon)
                    area = (arr2)/(arr2) * 125
                    pc1 =ax.pcolormesh(lon_grid,lat_grid,z_krige,cmap=cmap1)####, vmin=minvalue, vmax=maxvalue)
                    ax.scatter(pass_lon,pass_lat,c=pass_value,s=area, cmap=cmap1,edgecolor='black')####, vmin=minvalue, vmax=maxvalue)

                    # Adding the colorbar
                    plt.draw()
                    posn = ax.get_position()
                    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0, 0.04, posn.height])
                    plt.colorbar(pc1, cax=cbar_ax)
                    #Save and close figure
                    plt.savefig(fig_dir+"Interpolation_"+o_str+"_"+interpolation_type+"_"+timesteps[i].strftime('%m%d%Y')+"_hr"+str(timesteps[i].hour)+".png")
                    plt.close()   
        
        
                    # SET UP PLOT 2 - Variance
                    # define the coordinate system that the grid lons and grid lats are on
                    #lambert = cartopy.crs.LambertConformal(central_longitude=wps2.ref_lon, central_latitude=wps2.ref_lat,\
                    #                        false_easting=wps2.dx[plt_idx2], false_northing=wps2.dy[plt_idx2], secant_latitudes=None, standard_parallels=None,\
                    #                        globe=None)
                    
                    
                    fig2, ax2 = plt.subplots(1, 1, figsize=(10,8), subplot_kw={'projection': cartopy.crs.PlateCarree()})
                    cbar_ax2 = fig2.add_axes([0, 0, 0.1, 0.1])
                    
                    ## Map Features
                    title = str(interpolation_type+' Interpolation-Variance on '+timesteps[i].strftime('%-m-%-d-%Y  %H')+"Z").upper()
                    ax2.set_title(title,fontsize=16)
                    ax2.coastlines('10m','black')
                    ax2.add_feature(cartopy.feature.STATES.with_scale('10m'))
                    ax2 = nplt.wps_info.set_lambert_ticks(ax2,xskip=1.,yskip=1.,x_thickness=14,y_thickness=14)
                    fig2.canvas.draw()

                    # Plot Data
                    pc2 =ax2.pcolormesh(lon_grid,lat_grid,ss1,cmap=cmap1)#####, vmin=minvalue, vmax=maxvalue)

                    # Adding the colorbar
                    plt.draw()
                    posn2 = ax2.get_position()
                    cbar_ax2.set_position([posn2.x0 + posn2.width + 0.01, posn2.y0, 0.04, posn2.height])
                    plt.colorbar(pc2, cax=cbar_ax2)
                    
                    #Save and close figure
                    plt.savefig(fig_dir+"Variance_"+o_str+"_"+interpolation_type+"_"+timesteps[i].strftime('%m%d%Y')+"_hr"+str(timesteps[i].hour)+".png")
                    plt.close()
                    continue
                ## STORING DATA NEEDS TO HAPPEN HERE
                storeCount = 1 # Used in naming locations
                for i2 in range(len(valid_lon)):
                    #Only store on the first pass through (first variable for each case)
                    if firstPass == True:
                        storeMatrix[i_loc,0] = "LOC"+str(storeCount)
                        storeMatrix[i_loc,1] = timesteps[i].strftime('%Y%m%d%H%M%S')
                        storeMatrix[i_loc,8]   = valid_lat[i2]
                        storeMatrix[i_loc,9]   = valid_lon[i2]
                    if variable == 'Wind_Speed (m/s)':
                        storeMatrix[i_loc,2]   = get_value_interp(valid_lon[i2],valid_lat[i2],z_krige,wpsINTP,plt_idxINTP)
                    elif variable == 'Wind_Direction (deg)': 
                        storeMatrix[i_loc,3]   = get_value_interp(valid_lon[i2],valid_lat[i2],z_krige,wpsINTP,plt_idxINTP)
                    elif variable == 'Air_Temperature (K)':
                        storeMatrix[i_loc,4]   = get_value_interp(valid_lon[i2],valid_lat[i2],z_krige,wpsINTP,plt_idxINTP)
                    elif variable == 'Dewpoint_Temperature (K)':
                        storeMatrix[i_loc,5]  = get_value_interp(valid_lon[i2],valid_lat[i2],z_krige,wpsINTP,plt_idxINTP)
                    elif variable == "Relative Humidity (%)":
                        storeMatrix[i_loc,6]  = get_value_interp(valid_lon[i2],valid_lat[i2],z_krige,wpsINTP,plt_idxINTP)
                    elif variable == "Pressure (Pa)":
                        storeMatrix[i_loc,7]  = get_value_interp(valid_lon[i2],valid_lat[i2],z_krige,wpsINTP,plt_idxINTP)
                    else: 
                        print ("INVALID VARIABLE OPTION")
                        sys.exit(0)
                    storeCount +=1
                    i_loc +=1
            #SAVE DATA FOR SINGLE TIMESTEP
            
            else:
                # If not in analysis window then store as np.nan values to keep file formats the same and easy for model analysis.
                storeCount = 1
                for i2 in range(len(valid_lon)):
                    if firstPass == True:
                        storeMatrix[i_loc,0] = "LOC"+str(storeCount)
                        storeMatrix[i_loc,1] = timesteps[i].strftime('%Y%m%d%H%M%S')
                        storeMatrix[i_loc,8]   = valid_lat[i2]
                        storeMatrix[i_loc,9]   = valid_lon[i2]
                    if variable == 'Wind_Speed (m/s)':
                        storeMatrix[i_loc,2]   = np.nan
                    elif variable == 'Wind_Direction (deg)': 
                        storeMatrix[i_loc,3]   = np.nan
                    elif variable == 'Air_Temperature (K)':
                        storeMatrix[i_loc,4]   = np.nan
                    elif variable == 'Dewpoint_Temperature (K)':
                        storeMatrix[i_loc,5]  = np.nan
                    elif variable == "Relative Humidity (%)":
                        storeMatrix[i_loc,6]  = np.nan
                    elif variable == "Pressure (Pa)":
                        storeMatrix[i_loc,7]  = np.nan
                    else: 
                        print ("INVALID VARIABLE OPTION")
                        sys.exit(0)
                    storeCount +=1
                    i_loc +=1
            
            #SAVE DATA FOR A SINGLE TIMESTEP
        # Calculate and print statistics for each variable
        #if variable == 'Wind_Speed (m/s)' or variable == 'Air_Temperature (K)' or variable == 'Dewpoint_Temperature (K)': 
        try:
            Q1_MEAN = np.nanmean(Q1)
            Q2_MEAN = np.nanmean(Q2)
            cRMean = np.nanmean(cR)
            print(case_time[:10]+" - " + variable+'\n\t%-20s'%'Date/Time'+'\t\tQ1\t\tQ2\t\tcR')
            for rats in range(len(cR)):
                print('\t%-20s'%timeZeus[rats],'\t\t%3.3f'%Q1[rats],'\t\t%3.3f'%Q2[rats],'\t\t%3.3f'%cR[rats])
                
            print('\t%-20s'%"AVERAGE VALUES",'\t\t%3.3f'%Q1_MEAN,'\t\t%3.3f'%Q2_MEAN,'\t\t%3.3f'%cRMean)
        
        except:
            pass
        
        #SAVE DATA FOR A SINGLE VARIABLE - no longer fill the first 3 columns (name,lat,lon)
        firstPass = False
        
    #END OF TIME STEPS FOR ONE CASE - ALL VARIABLES
    station_data = pd.DataFrame(storeMatrix, index=None, columns = HEADER_INTP)
    
    # DO NOT SAVE ONCE WE HAVE THE DATA WE WANT
    #station_data.to_csv(outstring,index=False)
    
    
    
    
    #END OF VARIABLES
#END OF CASES   
        
