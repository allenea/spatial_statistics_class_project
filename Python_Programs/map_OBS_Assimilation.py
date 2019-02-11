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

#%%
# SAMPLE: 2014-06-04_06
casestudy_time = ['2014-06-04_06:00','2014-06-08_06:00','2015-08-14_06:00']
data_type = ["data_assimilation_original"]######,"data_assimilation_10m","verification_original", "verification_10m"]
runTime_D = 1
mapdir = os.getcwd()+"/Maps/"

float64 = np.float64

dtypesDict={"ID_String": object, "DATE":float64, "Wind_Speed (m/s)":float64, "Wind_Direction (deg)":float64,\
            "Air_Temperature (K)"  : float64, "Dewpoint_Temperature (K)" :float64, "Relative_Humidity (%)" :float64,\
            "Pressure (Pa)":float64, "Latitude"   :float64, "Longitude" :float64, "Elevation_SensorHeight (m)":float64,\
            "Name_string":object, "FM_string":object, "FM_Code":object, "Source_string":object}

HEADER_OUT = ["ID_String", 'DATE','Wind_Speed (m/s)','Wind_Direction (deg)','Air_Temperature (K)',"Dewpoint_Temperature (K)",\
              'Relative Humidity (%)','Pressure (Pa)','Latitude','Longitude','Elevation (m)','Name_string','FM_string','FM_Code','Source_string']

isLabelList = [True,False]           
for dtype in data_type:              
    #OPEN/READ DATA FILE - GET FARBER LOCATION
    data_dir = os.path.abspath('../')
    
    if dtype == "data_assimilation_original":
        ### STANDARD ASSIMILATION DATA (CMLF, DEOS)
        data_file = data_dir+'/Ferry_Data_Farber/data_CMLF2011_2016_D2010_2017.txt'
        outdf = '_'

    elif dtype ==  "data_assimilation_10m":
        ### STANDARD ASSIMILATION DATA INTERPOLATED TO 10m WINDS (CMLF, DEOS)
        data_file = data_dir+'/Ferry_Data_Farber/Assimilation_Data_10m_DEOS_CMLF.txt'
        outdf = '_10mda_'

    elif dtype == "verification_original":
        ### Verification Data Archive, for NJ and DE, (DEOS, NDBC, NCEI, DELDOT, NJMesonet,Meoswest**) - can be used for assimilation
        data_file = data_dir+'/Verification_Data/all_verification_OBSdata_BADARC.txt'
        outdf = '_Verify_'

    elif dtype == "verification_10m":
        ### Verification Data Archive, for NJ and DE, INTERPOLATED TO 10m WINDS (DEOS, NDBC, NCEI, DELDOT, NJMesonet,Meoswest**)
        data_file = data_dir+'/Verification_Data/all_verification_10mdata_BADARC.txt'
        outdf = '_10mVerify_'
    
    
    infile = pd.read_table(data_file, delim_whitespace=True,dtype = dtypesDict,header=None, names = HEADER_OUT)
    infile["FM_string"] = infile["FM_string"] +" "+ infile["FM_Code"]
    readData = infile.drop(['FM_Code'],axis=1)
    
    #print readData
    data = np.array(readData)
    del infile, readData


    for case_time in casestudy_time:
            print (case_time)   # FOR FARBER
                
            year = case_time[0:4]
                
            month = case_time[5:7]
            intMonth = int(month)
                
            day = case_time[8:10] 
            intDay = int(day)
            
            endDay = intDay + runTime_D
            
            if (endDay < 10):
                endDay = '0'+ str(endDay)    
            else:
                endDay = str(endDay)
                
                
            hour = case_time[11:13]
            intHour = int(hour)
            
            minute = case_time[14:16]
            intMinute = int(minute)
            ## for each file iterate hour, hour +1 , hour +2 and create output filess
                    
            print ("Year:",year, ". Month:",month,". Day:",day,". Hour:",hour,".Minute:", minute)
            
            # TIME STEP-WHAT DO I WANT TO USE?
            Start_Time = year+month+day+hour + minute +'00'
            fmtStart = month+'-'+day+'-'+year +" "+hour+":"+minute+ "Z"
            End_Time = year+month+endDay+hour + minute +'00'
            fmtEnd = month+'-'+endDay+'-'+year +" "+hour+":"+minute + "Z"
            
            print (Start_Time, "   ", End_Time)
            
            
            #%% Iterate backwards through the data.
            name = []; lat = []; lon = []; CMLF_LAT = []; CMLF_LON =[];
            for row in data:
                
                time = str(row[1])
                if (time[0:4] != year): continue
                elif (intMonth ==12 or intMonth == 1): print ("*************YEAR TROUBLE*********************")
                if int(time[4:6]) < intMonth: continue
                elif (((intMonth==2 and intDay == 28) or (intMonth==2 and intDay == 29)) or intDay==30 or intDay ==31): print ("*************END MONTH TROUBLE*********************")
                if (int(time[6:8]) < intDay): continue
                if (int(time[6:8])) > int(endDay):continue
                if(time >= Start_Time and time<= End_Time):
                    name.append(row[0])
                    lat.append(row[8])
                    lon.append(row[9])
                    if row[0] == 'CMLF':
                        CMLF_LAT.append(row[8])
                        CMLF_LON.append(row[9])
                elif(time >= End_Time):
                    break
            #%%    
               # TRIM
            stations = np.array(name)
            stationName, indicies = np.unique(stations,return_index=True)
            DuniqueLat = []; DuniqueLon =[];DuniqueName=[]
            FuniqueLat = []; FuniqueLon =[];FuniqueName=[]
            
            for index in indicies:
                if stations[index] == "CMLF":
                    for index2 in range(len(CMLF_LAT)):
                        if CMLF_LAT[index2] == '-888888.00000' or CMLF_LON[index2]=='-888888.00000':
                            continue
                        else:
                            FuniqueName.append("CMLF")
                            FuniqueLat.append(float(CMLF_LAT[index2]))
                            FuniqueLon.append(float(CMLF_LON[index2]))
                else:
                    DuniqueName.append(stations[index])
                    DuniqueLat.append(float(lat[index]))
                    DuniqueLon.append(float(lon[index])) 
            
            
            #%%
            for isLabel in isLabelList:
      
                ## SET INPUT DATA
                wpsfile = 'namelist.wps'
                plot_domains= 3
                
                ## GET/STORE WPS Data, Print Data, Get Plot Domain, Calculate Domain Info
                wps = nplt.wps_info(wpsfile)
                #wps.print_info()
                plt_idx = wps.plot_domain_number(plot_domains)
                wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()
                
                ## SET UP PLOT
                fig2 = plt.figure(figsize=(15,20))
                ax1 = plt.axes(projection=wpsproj)
                
                ## ADDING FEATURES
                ax1.coastlines('10m','black')
                ax1.add_feature(cartopy.feature.STATES.with_scale('10m'))
                
                ## REPORJECT
                corner_x, corner_y = nplt.wps_info.reproject_corners(corner_lon_full[plt_idx,:], corner_lat_full[plt_idx,:], wpsproj, latlonproj)
                
                ### ZOOM FUNCTION
                ns_zoom =35
                we_zoom = 45
                corner_x, corner_y = nplt.wps_info.plot_zoom_in(corner_x, corner_y,ns_zoom,we_zoom)
                ########
                
                ## SET DOMAIN LIMITS TO ZOOM 
                ax1.set_xlim([corner_x[0], corner_x[3]])
                ax1.set_ylim([corner_y[0], corner_y[3]])
                
                ax1.set_title("Locations of Assimilated Observation Data \n"+fmtStart + " until " + fmtEnd ,va='bottom', fontsize=20)
                
                
                request = cimgt.OSM()
                #####request = cimgt.GoogleTiles(style="satellite")
                ax1.add_image(request, 9)
                
                #####[longitude start, longitude end, latitude start, latitude end]
                #ax1.background_img(resolution='full')
                #ax1.add_wmts(arcmapurl,"ESRI_Imagery_World_2D/MapServer/tile/0/0/0")
                
                
                ax1.plot(DuniqueLon, DuniqueLat, 'o', color='r',markersize=10, transform=ccrs.PlateCarree(),label="Delaware Mesonet (DEOS)")
                ax1.plot(FuniqueLon, FuniqueLat, 'o', color='g',markersize=4, transform=ccrs.PlateCarree(),label="Cape May-Lewes Ferry (CMLF)")
                
                #fig1.canvas.draw()
                
                for label,xpt,ypt in zip(DuniqueName,DuniqueLon,DuniqueLat):
                    if label == "DBBB" and isLabel == True:
                        print (label)
                        ax1.text(xpt, ypt,"  "+label, fontsize=12, fontweight='bold',color = 'r',va='top', zorder=10, transform=ccrs.PlateCarree())
                    elif label == "DBNG" and isLabel == True:
                        print (label)
                        ax1.text(xpt, ypt,"  "+label, fontsize=12, fontweight='bold',color = 'r',va='bottom', zorder=10, transform=ccrs.PlateCarree())
                    elif label == "DDFS"and isLabel == True:
                        ax1.text(xpt, ypt,"  "+label, fontsize=12, fontweight='bold',color = 'r',va='bottom', zorder=10, transform=ccrs.PlateCarree())
                        print(label)
                    elif label == "DHAR"and isLabel == True:
                        ax1.text(xpt, ypt,"  "+label, fontsize=12, fontweight='bold',color = 'r',va='bottom', zorder=10, transform=ccrs.PlateCarree())
                        print(label)       
                    elif isLabel == True:
                        print (label)
                        ax1.text(xpt, ypt, "  "+label, fontsize=12, fontweight='bold',color = 'r',va='top', zorder=10, transform=ccrs.PlateCarree())
                    elif isLabel == False:
                        continue
                
                #ax1.set_title("\n   Cape May-Lewes Ferry  \n Observations",va='top',color='w', fontweight='bold')

                legend = ax1.legend(loc='upper left', frameon = True, shadow=True, fontsize=16,title = "Data Sources:")
                legend.get_title().set_fontsize(20)

    
                #axins = zoomed_inset_axes(ax2, 2, loc=1)
                
                ## INSERT AXIS ZOOM SET

                
                if isLabel == True:
                    plt.savefig(mapdir+'Assimilated_OBS_Map'+outdf+Start_Time[:8]+'_label.png')
                    #plt.show()
                else:
                    plt.savefig(mapdir+'Assimilated_OBS_Map'+outdf+Start_Time[:8]+'.png')
                    #plt.show()
                
                
                
        
        
#%%
"""
    ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
    ax1.coastlines('50m')
    ax1.set_extent([-45, 55, 20, 80], ccrs.PlateCarree())
    ax1.quiver(x, y, u, v, transform=vector_crs)

    ax2 = plt.subplot(2, 1, 2, projection=ccrs.PlateCarree())
    
    
    
fig = plt.figure(figsize=(15,20))
#m = Basemap(projection='lcc',resolution=None, width = 1E6, height = 1E6, lat_0 = 38.8, lon_0 = -75.0)
m = Basemap(projection='lcc', width = 8E6, height = 8E6, llcrnrlon = -76., llcrnrlat=38., urcrnrlon= -74.4, urcrnrlat=40., lon_0= -75.04, lat_0=38.88,epsg=4326)
#m.etopo()   #resolution="f"
m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)

#m.drawlsmask(land_color='#cc9955',ocean_color="aqua")
#m.drawcoastlines()
#m.drawstates()
#m.drawcounties()
# labels = [left,right,top,bottom]
parallels = np.arange(0.,81,1.)
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,2.)
m.drawmeridians(meridians,labels=[True,False,False,True])

x,y = m(DuniqueLat,DuniqueLon)
x2,y2 = m(FuniqueLat,FuniqueLon)

m.plot(x,y,'or',markersize=10)
m.plot(x2,y2,'ow',markersize=4)
for label,xpt,ypt in zip(DuniqueName,x,y):
    if label == "DBBB":
        plt.text(xpt,ypt," "+label, fontsize=12,fontweight='bold',color = 'r',va='bottom')
    elif label == "DBNG":
        plt.text(xpt,ypt,"    "+label, fontsize=12,fontweight='bold',color = 'r',va='top')
    else:
        plt.text(xpt,ypt,"  "+label, fontsize=12,fontweight='bold',color = 'r')
        
plt.text(x2[len(x2)/2],y2[len(y2)/2],"       Cape May-Lewes Ferry", fontsize=12,fontweight='bold',color = 'w')

plt.title("Location of Observations",fontsize=20)
plt.xlabel("Longitude", fontsize=12)
plt.ylabel("Latitude", fontsize=12)
plt.savefig('Observation Map.png')
plt.show()

"""
"""
sys.exit(0)
#m = Basemap(projection='lcc',resolution=None, width = 1E6, height = 1E6, lat_0 = 38.8, lon_0 = -75.0)
m = Basemap(projection='lcc', width = 8E6, height = 8E6, llcrnrlon = -76., llcrnrlat=38.25, urcrnrlon= -74.4, urcrnrlat=40., lon_0= -75.04, lat_0=38.88,epsg=4326)
#m.etopo()   #resolution="f"
m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)

#m.drawlsmask(land_color='#cc9955',ocean_color="aqua")
#m.drawcoastlines()
#m.drawstates()
#m.drawcounties()
# labels = [left,right,top,bottom]
parallels = np.arange(0.,81,1.)
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,2.)
m.drawmeridians(meridians,labels=[True,False,False,True])


x,y = m(DuniqueLat,DuniqueLon)
x2,y2 = m(FuniqueLat,FuniqueLon)

m.plot(x,y,'or',markersize=10)
m.plot(x2,y2,'ow',markersize=3)
for label,xpt,ypt in zip(DuniqueName,x,y):
    if label == "DBBB":
        plt.text(xpt,ypt," "+label, fontsize=12,fontweight='bold',color = 'r',va='bottom')
    elif label == "DBNG":
        plt.text(xpt,ypt,"    "+label, fontsize=12,fontweight='bold',color = 'r',va='top')
    else:
        plt.text(xpt,ypt,"  "+label, fontsize=12,fontweight='bold',color = 'r')
        
#plt.text(x2[len(x2)/2],y2[len(y2)/2],"       Cape May-Lewes Ferry", fontsize=12,fontweight='bold',color = 'w')
plt.title("Location of Observations \n"+fmtStart + " until " + fmtEnd ,va='bottom', fontsize=20)

ax = fig.add_subplot(111)
axins = zoomed_inset_axes(ax, 2, loc=1)
axins.set_xlim(-75.14, -74.94)
axins.set_ylim(38.75, 39.05)

map2 = Basemap(llcrnrlon = -75.14, llcrnrlat=38.75, urcrnrlon= -74.94, urcrnrlat=39.05,  ax=axins)
map2.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)
map2.plot(x2,y2,'ow',markersize=0.5)
plt.title("\n   Cape May-Lewes Ferry  \n Observations",va='top',color='w', fontweight='bold')
mark_inset(ax, axins, loc1=2, loc2=4)


plt.savefig('Observation Map'+Start_Time[:8]+'.png')
plt.show()
"""