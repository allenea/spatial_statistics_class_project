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
import math

 
casestudy_time = ['2014-06-04_06:00','2014-06-08_06:00','2015-08-14_06:00'] #
dtype_subusage = ['original'] #'10m',
dtype_usage = ["Assimilation"]
isRaw = False
variables = ['Wind_Speed (m/s)','Wind_Direction (deg)']
mydir = os.path.abspath('../')
data_dir = os.path.abspath('../../')

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
                if 'Wind_Speed (m/s)' not in variables:
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
                if 'Wind_Speed (m/s)' not in variables:
                    break
                else:
                    if isRaw == True:
                        data_file = data_dir+"/Ferry_Data/case_study_data/"+dtype1+"/"+case_time[0:10]+'/' #### RAW
                    else:
                        data_file = data_dir+"/Ferry_Data/hr_avg_trim/"+dtype1+"/"+case_time[0:10]+'/' #### TRIM
                        
                    
            if isRaw == True:
                outdir =  mydir+"/plots/Delaware_Surface_Winds/WSPD_WDIR/"
            elif isRaw == False:
                outdir =  mydir+"/plots/Delaware_Surface_Winds/WSPD_WDIR/"
    
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            
            filelst.append(data_file)
    
        if len(filelst) == 0:
            continue
        
        fig, (ax1, ax2,ax3) = plt.subplots(3, sharex=True, sharey=False,  figsize=(15,15)) 

        isDEOS1DA = False; isDEOS2V = False; isDelDOT = False; isASOS = False; isNDBC = False
        isCMLF = False; isNJMET = False;
        print(filelst)
        count1 =0; count2=0;count3=0;
        for fileDIR in filelst:
            for file in glob.glob(fileDIR+'*.csv'):

                data = pd.read_csv(file, low_memory=False)
                data.columns.tolist()
                data = data.mask(data == " ",other = np.nan)
                data = data.mask(data == "",other = np.nan)
                data = data.mask(data == -888888.0, other = np.nan)
                data['Wind_Speed (m/s)'] = data['Wind_Speed (m/s)'].mask(data['Wind_Speed (m/s)'] < 0, other = np.nan)
                data['Wind_Direction (deg)'] = data['Wind_Direction (deg)'].mask(data['Wind_Direction (deg)'] < 0, other = np.nan)

                #%% Load In Data
                
                ID_String = np.array(data['ID_String'],dtype=str)
                ws = np.array(data['Wind_Speed (m/s)'],dtype=float) 
                ws = ws * 1.9438444924574
                wdir = np.array(data['Wind_Direction (deg)'],dtype=float) 
                
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

                #%%
                colors = ['b','r', 'g','c','m','y','pink','purple']
                                
                u = [-abs(ws[i])*math.sin((math.pi/180)*wdir[i])for i in range(0,len(ws))]
                u = np.array(u)
                v = [-abs(ws[j])*math.cos((math.pi/180)*wdir[j])for j in range(0,len(ws))]
                v = np.array(v)
                z = [-2.25 for b in range(len(v))]
            
                #print cond_time
                print (ID_String[0])
                if ID_String[0] == "CMLF":
                    continue
                elif ID_String[0] == "DBNG" or ID_String[0]=="DBBB" or ID_String[0]=="DIRL" or ID_String[0]=="DRHB" or ID_String[0]=="DSJR" or ID_String[0]=="DWAR":
                    plotted = True
            
                    ax1.set_title("Coastal Stations", fontsize=12,fontweight='bold')
                    ax1.plot(time, ws, colors[count1],label=ID_String[0])
                    ax1.barbs(time,z,u,v,color=colors[count1],barb_increments=dict(half=1, full=10),rounding=False)
            
                    #ax1.grid(True)
            
                    ax1.set_yticks(np.arange(0, 16, step=1))
                    ax1.set_yticklabels(np.arange(0, 16, step=1), fontsize=14)
                    
                    atime= np.arange(len(time))
                    ax1.set_xticks(atime[::2])
                    ax1.set_xticklabels(time[::2], fontsize=14)
                    #ax1.minorticks_on()
                    
                    ax1.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')
                    ax1.tick_params(axis='both', which='major', length=14)
                    ax1.tick_params(axis='both', which='minor', length=8)
                    ax1.autoscale(enable=True, axis='x', tight=True)
                    count1 +=1
            
                elif ID_String[0]=="DELN" or ID_String[0]=="DGES" or ID_String[0]=="DJCR" or ID_String[0] == "DDFS" or ID_String[0]=="DGUM" or  ID_String[0]=="DDFS" or ID_String[0]=="DSMY":
                    plotted = True

                    ax2.set_title("Central Inland Delaware Stations", fontsize=12,fontweight='bold')
                    ax2.plot(time, ws, colors[count2],label=ID_String[0])
                    ax2.barbs(time,z,u,v,color=colors[count2],barb_increments=dict(half=1, full=10),rounding=False)
            
                    #ax2.grid(True)            
                    ax2.set_yticks(np.arange(0, 16, step=1))
                    ax2.set_yticklabels(np.arange(0, 16, step=1), fontsize=14)
                    
                    atime= np.arange(len(time))
                    ax2.set_xticks(atime[::2])
                    ax2.set_xticklabels(time[::2], fontsize=14)
                    #ax2.minorticks_on()
                    
                    ax2.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')
                    ax2.tick_params(axis='both', which='major', length=14)
                    ax2.tick_params(axis='both', which='minor', length=8)
                    ax2.autoscale(enable=True, axis='x', tight=True)
                    count2 +=1
            
            
                elif ID_String[0]=="DBLK" or ID_String[0]=="DADV" or ID_String[0]=="DSEA" or ID_String[0]=="DLAU" or ID_String[0]=="DBRG" or ID_String[0]=="DSND":
                    plotted = True
                    print ("WEST",count3,ID_String[0])
                    ax3.set_title("Western Delaware Stations", fontsize=12,fontweight='bold')
                    ax3.plot(time, ws, colors[count3],label=ID_String[0])
                    ax3.barbs(time,z,u,v,color=colors[count3],barb_increments=dict(half=1, full=10),rounding=False)
            
                    #ax3.grid(True)
            
                    ax3.set_yticks(np.arange(0, 16, step=1))
                    ax3.set_yticklabels(np.arange(0, 16, step=1), fontsize=14)
                    
                    atime= np.arange(len(time))
                    ax3.set_xticks(atime[::2])
                    ax3.set_xticklabels(time[::2], fontsize=14)
                    #ax3.minorticks_on()
                    ax3.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')
                    ax3.tick_params(axis='both', which='major', length=14)
                    ax3.tick_params(axis='both', which='minor', length=8)
                    ax3.autoscale(enable=True, axis='x', tight=True)
                
                    count3 +=1
            
                else:
                   print ("Not Plotted")
                
                plt.gcf().autofmt_xdate()
                
            ax2.set_ylabel('[Top] Wind Speed (kt) / [Bottom] Wind Direction', fontsize=16,fontweight='bold')

            ax1.axhline(0,linewidth=2,color='k')
            ax2.axhline(0,linewidth=2,color='k')
            ax3.axhline(0,linewidth=2,color='k')
            
            
            ax1.axhline(10,linewidth=2,linestyle=":", color='k', label="full barb-10kt")
            ax2.axhline(10,linewidth=2,linestyle=":", color='k', label="full barb-10kt")
            ax3.axhline(10,linewidth=2,linestyle=":", color='k', label="full barb-10kt")
            
            
            ax1.axhline(1,linewidth=0.5,linestyle="--", color='k', label="half barb-1kt")
            ax2.axhline(1,linewidth=0.5,linestyle="--", color='k', label="half barb-1kt")
            ax3.axhline(1,linewidth=0.5,linestyle="--", color='k', label="half barb-1kt")
            
            ax1.legend(loc="upper left", ncol=4, fontsize=12, labelspacing = -0.1)
            ax2.legend(loc="upper left", ncol=4, fontsize=12, labelspacing = -0.1)
            ax3.legend(loc="upper left", ncol=4, fontsize=12, labelspacing = -0.1)

      
            print (len(time), len(ws), max(ws))
            plt.suptitle('\n 30-Minute Delaware Surface Wind Observations \n(Assimilation Data)\n'+case_time, fontsize=20, fontweight='bold', verticalalignment='baseline', y=0.925,horizontalalignment='center') # or plt.suptitle('Main title')
            
            plt.savefig(outdir + "DEOS_Time_Series_Windplots_w_Direction"+case_time[:10].replace("-","_")+"_NEW.png")
            #plt.show()
            plt.close()