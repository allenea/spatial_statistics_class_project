#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 23:01:19 2018

@author: allenea
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:41:51 2018

@author: allenea
"""

import glob
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import sys


casestudy_time = ['2014-06-04_06:00','2014-06-08_06:00','2015-08-14_06:00']
dtype_subusage = ['10m','original']
dtype_usage = ["Verification", "Assimilation"]
variables = ['Wind_Speed (m/s)','Wind_Direction (deg)','Air_Temperature (K)','Dewpoint_Temperature (K)','Relative Humidity (%)','Pressure (Pa)']

isRaw = False

mydir = os.path.abspath('../')
data_dir = os.path.abspath('../../')

for variable in variables:
    print (variable)
    for case_time in casestudy_time:
        #for usage in dtype_usage:
        for dtype1 in dtype_subusage:
            filelst = []
            for datype in dtype_usage:
                if dtype1 == "original" and datype =="Verification":
                    if isRaw == True:
                        data_file = data_dir+"/Verification_Data/verify_case_study_data/"+dtype1+"/"+case_time[0:10]+'/' ### RAW
                    else:
                        data_file = data_dir+"/Verification_Data/hr_avg_trim/"+dtype1+"/"+case_time[0:10]+'/' ### TRIM

        
                elif dtype1 == "10m" and datype =="Verification":
                    if variable !='Wind_Speed (m/s)':
                        break
                    else:
                        if isRaw == True:
                            data_file = data_dir+"/Verification_Data/verify_case_study_data/"+dtype1+"/"+case_time[0:10]+'/' #### RAW
                        else:
                            data_file = data_dir+"/Verification_Data/hr_avg_trim/"+dtype1+"/"+case_time[0:10]+'/' #### TRIM
                        
                elif dtype1 == "original" and datype =="Assimilation":
                    if isRaw == True:
                        data_file = data_dir+"/Ferry_Data/case_study_data/"+dtype1+"/"+case_time[0:10]+'/' ### RAW
                    else:
                        data_file = data_dir+"/Ferry_Data/hr_avg_trim/"+dtype1+"/"+case_time[0:10]+'/' ### TRIM

        
                elif dtype1 == "10m" and datype =="Assimilation":
                    if variable !='Wind_Speed (m/s)':
                        break
                    else:
                        if isRaw == True:
                            data_file = data_dir+"/Ferry_Data/case_study_data/"+dtype1+"/"+case_time[0:10]+'/' #### RAW
                        else:
                            data_file = data_dir+"/Ferry_Data/hr_avg_trim/"+dtype1+"/"+case_time[0:10]+'/' #### TRIM
                            
                        
                if isRaw == True:
                    outdir =  mydir+"/plots/Near_Ferry/OBS/"
                elif isRaw == False:
                    outdir =  mydir+"/plots/Near_Ferry/AVG/BOX/"
        
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                
                filelst.append(data_file)
        
            if len(filelst) == 0:
                continue
            
            fig, ax = plt.subplots(sharex=True, figsize=(15,7)) 
    
            curMax = 0
            legendNames=[]
            allWinds = []
            landWinds = []
            waterWinds = []
            nearFerry = []

            isDEOS1DA = False; isDEOS2V = False; isDelDOT = False; isASOS = False; isNDBC = False
            isCMLF = False; isNJMET = False;
            print(filelst)
            
            for fileDIR in filelst:
                for file in glob.glob(fileDIR+'*.csv'):

                    data = pd.read_csv(file, low_memory=False)
                    data.columns.tolist()
                    data = data.mask(data == " ",other = np.nan)
                    data = data.mask(data == "",other = np.nan)
                    data = data.mask(data == -888888.0, other = np.nan)
                    data['Wind_Speed (m/s)'] = data['Wind_Speed (m/s)'].mask(data['Wind_Speed (m/s)'] < 0, other = np.nan)
        
                    #%% Load In Data
                    
                    ID_String = np.array(data['ID_String'],dtype=str)
                    if variable == 'Wind_Speed (m/s)':
                        variable_data= np.array(data['Wind_Speed (m/s)'],dtype=float) 
                    elif variable == 'Wind_Direction (deg)':
                        variable_data = np.array(data['Wind_Direction (deg)'],dtype=float) 
                    elif variable == 'Air_Temperature (K)':
                        variable_data = np.array(data['Air_Temperature (K)'],dtype=float) 
                    elif variable == 'Dewpoint_Temperature (K)':
                        variable_data = np.array(data['Dewpoint_Temperature (K)'],dtype=float) 
                    elif variable == "Relative Humidity (%)":
                        variable_data = np.array(data['Relative Humidity (%)'],dtype=float) 
                    elif variable == "Pressure (Pa)":
                        variable_data = np.array(data['Pressure (Pa)'],dtype=float) 
                    else:
                        sys.exit(0)
                    
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
        
                    if ID_String[0] == "CMLF" or ID_String[0] == "BRND1" or ID_String[0] == "CMAN4" or ID_String[0] == "KWWD" or ID_String[0] == "SJSN4" or \
                        ID_String[0] == "44009" or ID_String[0] == "LWSD1" or ID_String[0] == "DCPH" or ID_String[0] == "DLEW" or  ID_String[0] == "DEL03" or\
                         ID_String[0] == "DRHB" or ID_String[0] == "MWCM":
                            print (ID_String[0])
                            if ID_String[0] == "CMLF":
                                ax.plot(time,variable_data,marker= '*',linestyle=None,color ='green',linewidth = 4, label=ID_String[0])
                            else:
                                nearFerry.append(variable_data.tolist())
                                ax.plot(time,variable_data,marker= '*',linestyle="",label=ID_String[0])
        
        
        
            tmpWindsnearFerry = np.array(nearFerry)
            avgWindsNF = np.nanmean(tmpWindsnearFerry,axis=0)
        
            ax.plot(time,avgWindsNF,'black',linewidth=4,label="Average")
            #ax.plot(time,avgWindsLand,'k',linewidth=4,label="Average Land Obs")
            #ax.plot(time,avgWindsWater,'c',linewidth=4,label="Average Marine Obs")
            
            df = pd.DataFrame(nearFerry,columns=time)

            if "Wind_Speed (m/s)" == variable:
                boxplot = df.boxplot(column=time,ax=ax,grid=False,return_type='axes',rot=45, fontsize=14)
                ax.set_yticks(np.arange(0, 13, step=0.5))
                ax.set_yticklabels(np.arange(0, 13, step=0.5), fontsize=14)
                atime = np.arange(len(time))
                ax.set_xticks(atime[::2])
                ax.set_xticklabels(time[::2], fontsize=14)
                ax.set_ylabel(variable, fontsize=16,fontweight='bold')
                ax.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')
                
            elif 'Air_Temperature (K)' == variable:
                boxplot = df.boxplot(column=time,ax=ax,grid=False,return_type='axes',rot=45, fontsize=14)
                ax.set_yticks(np.arange(280, 315, step=5))
                ax.set_yticklabels(np.arange(280, 315, step=5), fontsize=14)
                atime = np.arange(len(time))
                ax.set_xticks(atime[::2])
                ax.set_xticklabels(time[::2], fontsize=14)
                ax.set_ylabel(variable, fontsize=16,fontweight='bold')
                ax.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')
                
            elif 'Dewpoint_Temperature (K)' == variable:
                boxplot = df.boxplot(column=time,ax=ax,grid=False,return_type='axes',rot=45, fontsize=14)
                ax.set_yticks(np.arange(280, 315, step=5))
                ax.set_yticklabels(np.arange(280, 315, step=5), fontsize=14)
                atime = np.arange(len(time))
                ax.set_xticks(atime[::2])
                ax.set_xticklabels(time[::2], fontsize=14)
                ax.set_ylabel(variable, fontsize=16,fontweight='bold')
                ax.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')  
          
            elif 'Relative Humidity (%)' == variable:
                boxplot = df.boxplot(column=time,ax=ax,grid=False,return_type='axes',rot=45, fontsize=14)
                ax.set_yticks(np.arange(0, 101, step=10))
                ax.set_yticklabels(np.arange(0, 101, step=10), fontsize=14)
                atime = np.arange(len(time))
                ax.set_xticks(atime[::2])
                ax.set_xticklabels(time[::2], fontsize=14)
                ax.set_ylabel(variable, fontsize=16,fontweight='bold')
                ax.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')  
                
            elif 'Wind_Direction (deg)' == variable:
                boxplot = df.boxplot(column=time,ax=ax,grid=False,return_type='axes',rot=45, fontsize=14)
                #http://snowfence.cfans.umn.edu/Components/winddirectionanddegreeswithouttable3.htm
                atime = np.arange(len(time))
                awdir = np.arange(-11.25,371.26,22.5)
                text_wdir = ["N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW","N"]
    
                ax.set_xticks(atime[::2])
                ax.set_xticklabels(time[::2], fontsize=14)
                ax.set_yticks(awdir[::2])
                ax.set_yticklabels(text_wdir[::2], fontsize=14)
                ax.set_ylabel(variable, fontsize=16,fontweight='bold')
                ax.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')
                
            else:
                boxplot = df.boxplot(column=time,ax=ax,grid=False,return_type='axes',rot=45, fontsize=14)
                atime = np.arange(len(time))
                ax.set_xticks(atime[::2])
                ax.set_xticklabels(time[::2], fontsize=14)
                ax.set_yticklabels(ax.get_yticks(), fontsize=14)
                ax.set_ylabel(variable, fontsize=16,fontweight='bold')
                ax.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')
            
            plt.legend(loc="best", ncol=4, fontsize=14, labelspacing = -0.1, columnspacing = 0.5)
            plt.xlabel('Time (UTC)', fontsize=16,fontweight='bold')
            plt.ylabel(variable, fontsize=16,fontweight='bold')
  
            plt.x_labels = time
            plt.xticks(time[::2],fontsize=14)
            plt.gcf().autofmt_xdate()
            
            if isRaw == True:
                if "Wind_Speed (m/s)" in variable and dtype1 == "10m":
                    plt.title("Observed "+dtype1+" Wind Speeds (m/s) - "+case_time, fontsize=20)
                    
                elif "Wind_Speed (m/s)" in variable and dtype1 == "original":
                    plt.title("Observed Wind Speeds (m/s) - "+case_time, fontsize=20)
                    
                else:
                    plt.title("Observed "+variable+"  -  "+case_time, fontsize=20)
            else:
            
                if "Wind_Speed (m/s)" in variable and dtype1 == "10m":
                    plt.title("30 Minute Averaged "+dtype1+" Wind Speeds (m/s) Observed Near Ferry On "+case_time, fontsize=20)
                    
                elif "Wind_Speed (m/s)" in variable and dtype1 == "original":
                    plt.title("30 Minute Averaged Wind Speeds (m/s) Observed Near Ferry On "+case_time, fontsize=20)
                
                else:
                    plt.title("30 Minute Averaged "+variable+" Observed Near Ferry On "+case_time, fontsize=20)
                     
                     
            if isRaw == True:
                plt.savefig(outdir+"NEAR_FERRY_Time-Series_plot_"+dtype1+"_"+variable[:-6]+"_"+case_time[0:10]+"BOX.png")
            else:
                plt.savefig(outdir+"NEAR_FERRY_Time-Series_plot_"+dtype1+"_"+variable[:-6]+"_"+case_time[0:10]+"BOX.png")
            plt.close()