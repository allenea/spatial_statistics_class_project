#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 19:08:15 2019

@author: allenea
"""

# =============================================================================
# =============================================================================
# =============================================================================
#####                         MODIFY AS NEEDED                            #####
# =============================================================================
# =============================================================================
# =============================================================================
import os
import sys
import numpy as np
import pandas as pd

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
    /Users/myusername/Documents/MAIN_DIR/Verification_Data/interpolation_post/10m/2014-06-04/interpolated_verification_data.csv
    /Users/myusername/Documents/MAIN_DIR/Verification_Data/interpolation_post/10m/2014-06-08/interpolated_verification_data.csv
    /Users/myusername/Documents/MAIN_DIR/Verification_Data/interpolation_post/10m/2015-08-14/interpolated_verification_data.csv
    
    """
    
    #ROOT DIRS
    data_dir = os.path.abspath('../')
    
    filenameIN =  "/interpolated_verification_data.csv"
    
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
    outdir = data_dir+pathOUT+case[:10]+"/"
    outstring = outdir
    
    #Checks to see if path exits, if not create it
    if not os.path.exists(outdir): os.makedirs(outdir)

    return instring, outstring

HEADER_OUT  = ["ID_String", 'DATE','Wind_Speed (m/s)','Wind_Direction (deg)','Air_Temperature (K)',"Dewpoint_Temperature (K)",\
              'Relative Humidity (%)','Pressure (Pa)','Latitude','Longitude']

# ADMINISTRATIVE - os

sub1dirIN = "Verification_Data"; sub2dirIN = "interpolation_post"; sub3dirIN = "10m";
subdirlistIN =  [sub1dirIN,sub2dirIN,sub3dirIN] 
sub1dirOUT = "Verification_Data"; sub2dirOUT = "interpolation_split"; sub3dirOUT = "10m";
subdirlistOUT = [sub1dirOUT,sub2dirOUT,sub3dirOUT]
# CASE STUDY/ANALYSIS
casestudy_time = ['2014-06-04_06:00','2014-06-08_06:00','2015-08-14_06:00']
num_random_locs = 250

for case_time in casestudy_time:
     # Admin - Set Input file, output file and figure directory paths
    data_file,outstring = fmt_data_path(case_time,subdirlistIN,subdirlistOUT)
    
    # Read In Data   
    data = pd.read_csv(data_file, low_memory=False)
    data.columns.tolist()
    data = data.mask(data == " ",other = np.nan)
    data = data.mask(data == "",other = np.nan)
    dataf = np.array(data)
    for num in range(1,num_random_locs + 1):
        ID_String_IDX = "LOC"+str(num)
        output = outstring+ID_String_IDX +"_hr_avg_trim.csv"
        stationData = []
        for row in dataf:
            if ID_String_IDX == row[0]:
                stationData.append(row)
        df= pd.DataFrame(stationData,columns=HEADER_OUT)
        df.to_csv(output, index=False)
