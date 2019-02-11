#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 10:30:27 2018

@author: allenea
"""
import glob
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
import sys
import os

casestudy_time = ['2014-06-04_06:00','2014-06-08_06:00','2015-08-14_06:00']
dtype_subusage = ['original','10m']
variables = ['Wind_Speed (m/s)','Wind_Direction (deg)','Air_Temperature (K)','Dewpoint_Temperature (K)','Relative Humidity (%)','Pressure (Pa)']
isRaw = False

mydir = os.path.abspath('../')
data_dir = os.path.abspath('../../')

for case_time in casestudy_time:
    for variable in variables:
        print(variable)
        #for usage in dtype_usage:
        for dtype1 in dtype_subusage:
            if dtype1 == "original":
                if isRaw == True:
                    data_file = data_dir+"/Ferry_Data/case_study_data/"+dtype1+"/"+case_time[0:10]+'/' ### RAW
                else:
                    data_file = data_dir+"/Ferry_Data/hr_avg_trim/"+dtype1+"/"+case_time[0:10]+'/' ### TRIM

    
            elif dtype1 == "10m":
                if variable !='Wind_Speed (m/s)':
                    break
                else:
                    if isRaw == True:
                        data_file = data_dir+"/Ferry_Data/case_study_data/"+dtype1+"/"+case_time[0:10]+'/' #### RAW
                    else:
                        data_file = data_dir+"/Ferry_Data/hr_avg_trim/"+dtype1+"/"+case_time[0:10]+'/' #### TRIM
                        
            if isRaw == True:
                outdir =  mydir+"/plots/Ferry_Only/OBS/"
            elif isRaw == False:
                outdir =  mydir+"/plots/Ferry_Only/AVG/"

    
            if not os.path.exists(outdir):
                os.makedirs(outdir)
    
            plt.figure(figsize=(15,7))
    
    
            for file in glob.glob(data_file+'CMLF*.csv'):
   
                data = pd.read_csv(file, low_memory=False)
                data.columns.tolist()
                data = data.mask(data == " ",other = np.nan)
                data = data.mask(data == "",other = np.nan)
                data = data.mask(data == -888888.0, other = np.nan)
                data['Wind_Speed (m/s)'] = data['Wind_Speed (m/s)'].mask(data['Wind_Speed (m/s)'] < 0, other = np.nan)
    
                #%% Load In Data
                ID_String = np.array(data['ID_String'],dtype=str)
    
                if variable == 'Wind_Speed (m/s)':
                    variable_data= np.array(data['Wind_Speed (m/s)'])
                elif variable == 'Wind_Direction (deg)':
                    variable_data = np.array(data['Wind_Direction (deg)'])
                elif variable == 'Air_Temperature (K)':
                    variable_data = np.array(data['Air_Temperature (K)'])
                elif variable == 'Dewpoint_Temperature (K)':
                    variable_data = np.array(data['Dewpoint_Temperature (K)']) 
                elif variable == "Relative Humidity (%)":
                    variable_data = np.array(data['Relative Humidity (%)'])
                elif variable == "Pressure (Pa)":
                    variable_data = np.array(data['Pressure (Pa)']) 
                else:
                    sys.exit(0)
                
                lat = np.array(data['Latitude']) 
                long = np.array(data['Longitude']) 
                elev = np.array(data['Elevation (m)']) 
                name = np.array(data['Name_string']) 
                fm = np.array(data['FM_string']) 
                source = np.array(data['Source_string']) 
                
                
                time = []
    
                if isRaw == False:
                    year =np.array(data['YEAR'],dtype=int) 
                    month =np.array(data['MONTH'],dtype=int) 
                    day = np.array(data['DAY'],dtype=int) 
                    hour =np.array(data['HOUR'],dtype=int) 
                    minute = np.array(data['MINUTE'],dtype=int) 

                    for idx in range(len(year)):
                        utc_dt= datetime.datetime(year[idx],month[idx],day[idx],hour[idx],minute[idx])
                        time.append(utc_dt.strftime('%d/%H:%M'))
                elif isRaw == True:
                        Date = np.array(data['DATE'],dtype=str)
                        for dte in Date:
                            utc_dt= datetime.datetime(int(dte[0:4]),int(dte[4:6]),int(dte[6:8]),int(dte[8:10]),int(dte[10:12]),int(dte[12:14]))
                            #dt= utc_dt.strftime('%Y-%m-%d %H:%M:%S')
                            dt= utc_dt.strftime('%d/%H:%M')
                            time.append(dt)
    
                
                if ID_String[0] == "CMLF":
                    print (ID_String[0], np.nanmax(variable_data))
                    plt.plot(time,variable_data,'green',linewidth=1,label="CMLF")
                    

                plt.xlabel('* Time (UTC)', fontsize=16,fontweight='bold')
                plt.ylabel(variable, fontsize=16,fontweight='bold',)
                
    
                if isRaw == True:
                    if "Wind Speed" in variable:
                        if dtype1 == "10m":
                            plt.title("Observed "+dtype1+" Wind Speed By The Ferry  "+case_time, fontsize=20)
                        else:
                            plt.title("Observed Wind Speed By The Ferry  "+case_time, fontsize=20)

                    else:
                        plt.title("Observed "+variable+" Observed Near The Ferry  "+case_time, fontsize=20)
                else:
                
                    if "Wind Speed" in variable:
                        plt.title("30 Minute Averaged "+dtype1+" Wind Speed Observed By The Ferry  "+case_time, fontsize=20)
                    else:
                        plt.title("30 Minute Averaged "+variable+" Observed By The Ferry  "+case_time, fontsize=20)
                
                
                #TAB OVER FOR ALL ON 1 plot
                plt.legend(loc="best", ncol=3, fontsize=14, labelspacing = -0.1)
                #N = 6 # we will plot only every 5 date from date_list
                plt.x_labels = time
                #plt.x_labels_major = time[::N]
                if isRaw == True:
                    plt.xticks(time[::75],fontsize=14)
                else:
                    plt.xticks(time[::2],fontsize=14)
    
                plt.gcf().autofmt_xdate()
                #plt.minorticks_on()
                #plt.text(0.5, -4.0, '* Time of Every 75th Observation', fontsize=12,fontweight='bold')
                #plt.title("\n" + name[0] + " ("+ID_String[0]+" ) Wind Speed Analysis on "+ case_time[0:10], fontsize=12,fontweight='bold',)
                #plt.grid(True)
                
                if isRaw == True:
                     plt.savefig(outdir+"CMLF_Time_Series_Obs_"+dtype1+"_"+variable[:-6]+"_"+case_time[0:10]+".png")
                else:
                     plt.savefig(outdir+"CMLF_Time_Series_Avg_"+dtype1+"_"+variable[:-6]+"_"+case_time[0:10]+".png")

                plt.close()