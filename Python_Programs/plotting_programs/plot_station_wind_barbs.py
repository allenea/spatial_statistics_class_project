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
                outdir =  mydir+"/plots/Delaware_Surface_Winds/W_BARBS/"
            elif isRaw == False:
                outdir =  mydir+"/plots/Delaware_Surface_Winds/W_BARBS/"
    
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            
            filelst.append(data_file)
    
        if len(filelst) == 0:
            continue
        
        fig, ax1 = plt.subplots(1,figsize=(18,22))

        isDEOS1DA = False; isDEOS2V = False; isDelDOT = False; isASOS = False; isNDBC = False
        isCMLF = False; isNJMET = False;
        print(filelst)
        count1 =0; count2=0;count3=0;
        names = ["DRHB", "DIRL", "DBBB","DBNG", "DWAR", "DSJR", "DSMY", "DDFS", "DGUM", "DJCR", "DGES", "DELN", "DMIL", "DHAR", "DBRG", "DBLK", "DSND", "DADV", "DSEA", "DLAU"]
        stations =  names[::-1]
        yLabels = [];tTicks=[]
        count1 = 0  ## FOR y axis values
        for fileDIR in filelst:
            for name in stations:
                file = glob.glob(fileDIR+name+'*.csv')[0]

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
                        
    
                u = [-abs(ws[i])*math.sin((math.pi/180)*wdir[i])for i in range(0,len(ws))]
                u = np.array(u)
                v = [-abs(ws[j])*math.cos((math.pi/180)*wdir[j])for j in range(0,len(ws))]
                v = np.array(v)
                
                z = [0+count1 for b in range(len(v))]
                count1 +=1
                
                #print (ID_String[0])
                if ID_String[0] == "CMLF":
                    continue
            
                elif ID_String[0] in stations:
                    ax1.barbs(time,z,u,v,color='k',barb_increments=dict(half=1, full=10),rounding=False)
                    ax1.set_title("DEOS Station Wind-Barb Time-Series Analysis", fontsize=12,fontweight='bold')
                    tTicks.append(z[0])
                    yLabels.append(ID_String[0])
                else:
                   print ("Not Plotted")
                    
            yax = np.arange(0,len(tTicks),1)
            ax1.set_ylabel("Western Delaware (Inland)                              Central Delaware                        Eastern Delaware                          Coastal", fontsize=14)
            ax1.set_yticks(yax)
            ax1.set_yticklabels(list(yLabels), fontsize=14)

            atime= np.arange(len(time))
            ax1.set_xticks(atime[::2])
            ax1.set_xticklabels(time[::2], fontsize=14)
            
            ax1.set_xlabel('Time (UTC)', fontsize=16,fontweight='bold')
            #print (len(time), len(ws), max(ws))
            ax1.tick_params(axis='x', which='major', length=14)
            ax1.tick_params(axis='x', which='minor', length=8)
            ax1.autoscale(enable=True, axis='x', tight=True)
            plt.gcf().autofmt_xdate()

            ax1.axhline(6.5,linewidth=2,color='r')
            ax1.axhline(11.5,linewidth=2,color='r')
            ax1.axhline(15.5,linewidth=2,color='r')
            plt.suptitle('\n Delaware Surface Wind Observations\n(Assimilation Data)\n'+case_time, fontsize=20, fontweight='bold', verticalalignment='baseline', y=0.925,horizontalalignment='center') # or plt.suptitle('Main title')
            
            plt.savefig(outdir+"Station_Wind_Barbs"+case_time[:10]+".png")
            #plt.show()
            plt.close()
# =============================================================================
#                  elif ID_String[0] == "DBNG" or ID_String[0]=="DBBB" or ID_String[0]=="DIRL" or ID_String[0]=="DRHB" or ID_String[0]=="DSJR" or ID_String[0]=="DWAR":
#                elif ID_String[0]=="DELN" or ID_String[0]=="DGES" or ID_String[0]=="DJCR" or ID_String[0] == "DDFS" or ID_String[0]=="DGUM" or  ID_String[0]=="DDFS" or ID_String[0]=="DSMY":
#                elif ID_String[0]=="DBLK" or ID_String[0]=="DADV" or ID_String[0]=="DSEA" or ID_String[0]=="DLAU" or ID_String[0]=="DBRG" or ID_String[0]=="DSND":
          
            
# "DRHB"
# "DIRL"
# "DBBB"
# "DBNG"

# "DWAR"
# "DSJR"
# "DSMY"
# "DDFS"

# "DGUM"
# "DJCR"
# "DGES"
# "DELN"
# "DMIL"


# "DHAR"
# "DBRG"
# "DBLK"
# "DSND"
# "DADV"
# "DSEA"
# "DLAU"
# =============================================================================