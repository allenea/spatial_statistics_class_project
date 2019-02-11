#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 10:30:06 2018

@author: allenea


"""

#%%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import os
import pandas as pd
import namelist_plot as nplt
import cartopy.crs as ccrs
import cartopy
import cartopy.io.img_tiles as cimgt

#%%

# =============================================================================
# =============================================================================
# =============================================================================
#####                         MODIFY AS NEEDED                            #####
# =============================================================================
# =============================================================================
# =============================================================================
import sys

def fmt_data_path(case, subdirlistIN):
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
    data_dir = os.path.abspath('../..')
    
    filenameIN =  "/all_hr_avg_trim.csv"
    
    # SET INPUT DATA PATH
    pathIN = "/"
    for subdir_in in subdirlistIN:
        pathIN = pathIN +subdir_in + "/"
    instring = data_dir+pathIN+case[:10]+filenameIN
    # CHECK TO MAKE SURE THE PATH EXISTS
    if not os.path.exists(instring): print("INVALID DATA PATH"); sys.exit(0)


    return instring

HEADER_OUT  = ["ID_String", 'DATE','Wind_Speed (m/s)','Wind_Direction (deg)','Air_Temperature (K)',"Dewpoint_Temperature (K)",\
              'Relative Humidity (%)','Pressure (Pa)','Latitude','Longitude']

# ADMINISTRATIVE - os

sub1dirIN = "Verification_Data"; sub2dirIN = "interpolation_pre"; sub3dirIN = "10m";
subdirlistIN =  [sub1dirIN,sub2dirIN,sub3dirIN] 

# CASE STUDY/ANALYSIS
casestudy_time = ['2014-06-04_06:00','2014-06-08_06:00','2015-08-14_06:00']

for case_time in casestudy_time:
     # Admin - Set Input file, output file and figure directory paths
    data_file = fmt_data_path(case_time,subdirlistIN)
    
    # Read In Data   
    data = pd.read_csv(data_file, low_memory=False)
    data.columns.tolist()
    short_data = data[:len(list(set(data.ID_String)))]
    ID_String = short_data.ID_String
    Latitude = short_data.Latitude
    Longitude = short_data.Longitude
      
    ## SET INPUT DATA
    wpsfile = 'namelist.wps_MASTER'
    plot_domains= 3
    
    ## GET/STORE WPS Data, Print Data, Get Plot Domain, Calculate Domain Info
    wps = nplt.wps_info(wpsfile)
    #wps.print_info()
    plt_idx = wps.plot_domain_number(plot_domains)
    wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()
    
    ## SET UP PLOT
    fig2 = plt.figure(figsize=(15,20))
    ax1 = plt.axes(projection=wpsproj)
    ax1 = nplt.wps_info.set_lambert_ticks(ax1,1.0,1.0,14,14)

    ## ADDING FEATURES
    ax1.coastlines('10m','black')
    ax1.add_feature(cartopy.feature.STATES.with_scale('10m'))
    
    ## REPORJECT
    corner_x, corner_y = nplt.wps_info.reproject_corners(corner_lon_full[plt_idx,:], corner_lat_full[plt_idx,:], wpsproj, latlonproj)
    
    ### ZOOM FUNCTION
    ns_zoom = 0.0
    we_zoom = 0.0
    corner_x, corner_y = nplt.wps_info.plot_zoom_out(corner_x, corner_y,ns_zoom,we_zoom)
    ########
    
    ## SET DOMAIN LIMITS TO ZOOM 
    ax1.set_xlim([corner_x[0], corner_x[3]])
    ax1.set_ylim([corner_y[0], corner_y[3]])
    
    ax1.set_title("Locations of Observation Data for "+case_time[:10]+" Simulation" ,va='bottom', fontsize=20)
    
    
    request = cimgt.OSM()
    #####request = cimgt.GoogleTiles(style="satellite")
    ax1.add_image(request, 9)
    fig2.canvas.draw()

    #####[longitude start, longitude end, latitude start, latitude end]
    #ax1.background_img(resolution='full')
    #ax1.add_wmts(arcmapurl,"ESRI_Imagery_World_2D/MapServer/tile/0/0/0")
    
    
    ax1.plot(Longitude, Latitude, 'o', color='r',markersize=10, transform=ccrs.PlateCarree(),label="Observations")
        
       
    
    legend = ax1.legend(loc='best', frameon = True, shadow=True, fontsize=16,title = "Data Sources:")
    legend.get_title().set_fontsize(20)

    plt.savefig(os.getcwd()+'/Final_OBS_MAP'+case_time[:10].replace("-","_")+'.png')
    plt.show()
        
                