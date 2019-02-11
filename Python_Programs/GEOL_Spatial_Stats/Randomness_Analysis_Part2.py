#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 15:13:44 2018

@author: ericallen

Calculates the Spatial Distribution of a dataset.



Todo: take difference between xlim max and min then do a conversion factor and set that as the width???
"""
from math import sin, cos, sqrt, asin,pi,acos
import numpy as np
import random
import matplotlib.pyplot as plt
import namelist_plot as nplt
from Randomness_Analysis_Part1 import data_randomness
import os
import sys
import pandas as pd


def calculate_randomness_distribution(Simulations, data_file, namelist_wps_filename, domain_number, geo_fmt="degrees",figure_name_save="Data_Randomness"):
    if os.path.isfile(data_file) == True:
        if ".txt" in data_file:
            infile = pd.read_table(data_file, low_memory=False)
            infile.columns.tolist()
            infile.columns = [(col.strip()).upper() for col in infile.columns]
            #print("HEADER:  ",infile.columns)
            StationFile = np.array(infile['ID_STRING'])
            LatitudeFile = np.array(infile['LATITUDE'])
            LongitudeFile = np.array(infile['LONGITUDE'])
            
            station_unique, indicies = np.unique(np.array(StationFile),return_index=True)

            Latitude_IN = []; Longitude_IN = [];Station_IN=[]
            for index in indicies:         
                Station_IN.append(StationFile[index])
                Latitude_IN.append(LatitudeFile[index])
                Longitude_IN.append(LongitudeFile[index]) 
                print(StationFile[index],LongitudeFile[index] ,LatitudeFile[index])
                    
        elif ".csv" in data_file:
            infile = pd.read_csv(data_file, low_memory=False)
            infile.columns.tolist()
            infile.columns = [(col.strip()).upper() for col in infile.columns]
            #print("HEADER:  ",infile.columns)
            StationFile = np.array(infile['ID_STRING'])
            LatitudeFile = np.array(infile['LATITUDE'])
            LongitudeFile = np.array(infile['LONGITUDE'])
            
            station_unique, indicies = np.unique(np.array(StationFile),return_index=True)

            Latitude_IN = []; Longitude_IN = [];Station_IN=[]
            for index in indicies:         
                Station_IN.append(StationFile[index])
                Latitude_IN.append(LatitudeFile[index])
                Longitude_IN.append(LongitudeFile[index]) 
        else:
            sys.exit(0)
            #Latitude_IN = np.array(["41.8.44","42.51.0","41.19.0","44.16.58","41.35.27","44.47.48","41.30.51","41.15.48","43.1.40"])
            #Longitude_IN = np.array(["104.48.7","106.19.30","105.35.0","105.30.19","109.13.21","106.57.32","109.27.54","110.57.53","108.23.42"])
    else:
        #Latitude_IN = np.array(["41.8.44","42.51.0","41.19.0","44.16.58","41.35.27","44.47.48","41.30.51","41.15.48","43.1.40"])
        #Longitude_IN = np.array(["104.48.7","106.19.30","105.35.0","105.30.19","109.13.21","106.57.32","109.27.54","110.57.53","108.23.42"])
        sys.exit(0)

    
    n = len(Longitude_IN)
    
    OGMean, OGMeanNearest, OGStandardDeviation, OGSTDNearest =  data_randomness(Latitude_IN,Longitude_IN,geo_fmt=geo_fmt)
    
    
    #FUNCTIONS
    def haversine(lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points 
        on the earth (specified in radians... already converted from DD)
        
        return distance in meters
        """
        # convert decimal degrees to radians 
        #lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
        # haversine formula 
        dlon = lon2 - lon1 
        dlat = lat2 - lat1 
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a)) 
        r = 6371  #* 1000 # Radius of earth in kilometers. multiply by 1000 for meters. Use 3956 for miles
        return c * r



    def random_xy(namelist_wps_filename,domain_number):
        ## Initialize the grid
        wps = nplt.wps_info(namelist_wps_filename)
        plt_idx = wps.plot_domain_number(domain_number)
        wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()
        
        ulat= corner_lat_full[plt_idx,3]
        llat= corner_lat_full[plt_idx,0]
        rlon= corner_lon_full[plt_idx,3]
        llon= corner_lon_full[plt_idx,0]
        
        # =========================TEST====================================================
        #ulat= (45 + 0 / 60.0 + 0 / 3600.0) #*degrad
        #llat= (41 + 0 / 60.0 + 0 / 3600.0)  #*degrad
        #rlon= (-1. * (104 + 0 / 60.0 + 0 / 3600.0) ) #*degrad
        #llon= ( -1. * (111 + 0 / 60.0 + 0 / 3600.0) ) #*degrad
        # =============================================================================
        Latitude = [random.uniform(llat, ulat) for x in range(n)]
        Latitude = [l*degrad for l in Latitude]         # put in radians for law of cosine
        
        Longitude = [random.uniform(llon, rlon) for x in range(n)]
        Longitude = [ll * degrad for ll in Longitude]   # put in radians for law of cosine
        
        return Longitude, Latitude
    
    def dist_calc(rlat,rlat2,rlong,rlong2):
        d = acos(((sin(rlat)*sin(rlat2))+(cos(rlat)*cos(rlat2)*cos(rlong-rlong2))))
        return d
    
    
    #Initialize Arrays
    Mean   = np.zeros(Simulations)
    MeanNearest = np.zeros(Simulations)
    StandardDeviation    = np.zeros(Simulations)
    STDNearest  = np.zeros(Simulations)
    
    #CONSTANT
    degrad = pi / 180.0
    
    for i in range(Simulations):
        dist_list = []
        nearestNeighbor = []
        Longitude, Latitude = random_xy(namelist_wps_filename,domain_number) # As radians
        for j in range(n):
            rlong = Longitude[j]
            rlat = Latitude[j]
            shortList = []
            for k in range(n):
                if j != k:
                    rlong2 = Longitude[k]
                    rlat2 = Latitude[k]
                    #d = dist_calc(rlat,rlat2,rlong,rlong2)
                    d = haversine(rlong, rlat, rlong2, rlat2)
                    dist_list.append(d)
                    shortList.append(d)
                    
            nearestNeighbor.append(min(shortList))
    
        Mean[i] = np.mean(dist_list)
        MeanNearest[i] = np.mean(nearestNeighbor)
        StandardDeviation[i] = np.std(dist_list)
        STDNearest[i] = np.std(nearestNeighbor)
        
    #ALL MEAN
    All_Mean = np.mean(Mean)
    All_Standard_Deviation = np.std(Mean)
    
    #MEAN STD
    All_STD_Mean = np.mean(StandardDeviation)
    All_STD_STD = np.std(StandardDeviation)
    
    #MEAN NEAREST
    All_Near_Mean_Mean = np.mean(MeanNearest)
    All_Near_Mean_STD = np.std(MeanNearest)
    
    #STD STD
    All_Near_STD_Mean_STD = np.mean(STDNearest)
    All_Near_STD_STD = np.std(STDNearest)
    
    print()
    print (Simulations , " Simulations")
    print ("Mean Distance:         %5.5f"%All_Mean)
    print ("STD Distance:          %5.5f"%All_Standard_Deviation)
    
    print ("Mean STD:              %5.5f"%All_STD_Mean)
    print ("STD of STD:            %5.5f"%All_STD_STD)
    
    print ("Mean NN Distance:      %5.5f"%All_Near_Mean_Mean)
    print ("Mean NN STD:           %5.5f"%All_Near_Mean_STD)
    
    print ("Mean Nearest STD:      %5.5f"%All_Near_STD_Mean_STD)
    print ("Mean Nearest STD STD:  %5.5f"%All_Near_STD_STD)
    
    
    #PLOTTING  
    fig, axarr = plt.subplots(2, 2, figsize=(16,16))
    plt.subplots_adjust(hspace=0.25, wspace=0.25)

    axarr[0, 0].hist(Mean, bins=20,color='b')
    axarr[0, 0].set_title('Mean Distances',fontsize=16)
    axarr[0, 0].annotate(str(round(OGMean, 3)),
                    xy=(round(OGMean, 3), 0), xycoords='data',
                    xytext=(axarr[0, 0].get_xlim()[1] * 0.8, axarr[0, 0].get_ylim()[1] * 0.99), textcoords='data',
                    arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",facecolor='red',
                                connectionstyle="angle3,angleA=0,angleB=-90"))
    axarr[0,0].autoscale(enable=True, axis='x', tight=True)
    axarr[0, 0].bar(All_Mean, axarr[0, 0].get_ylim()[1], width=(0.025*axarr[0, 0].get_xlim()[1])/3, color='r', alpha=0.5, label='Average')
    axarr[0, 0].set_xlim(min([axarr[0, 0].get_xlim()[0],OGMean]), axarr[0, 0].get_xlim()[1])
    
    
    axarr[0, 1].hist(MeanNearest, bins=20,color='b')
    axarr[0, 1].set_title('Mean Nearest Neighbor',fontsize=16)
    axarr[0, 1].annotate(str(round(OGMeanNearest, 3)),
                    xy=(round(OGMeanNearest, 3), 0), xycoords='data',
                    xytext=(axarr[0, 1].get_xlim()[1] * 0.8, axarr[0, 1].get_ylim()[1] * 0.99), textcoords='data',
                    arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",facecolor='red',
                                connectionstyle="angle3,angleA=0,angleB=-90"))
    axarr[0, 1].autoscale(enable=True, axis='x', tight=True)
    axarr[0, 1].bar(All_Near_Mean_Mean, axarr[0, 1].get_ylim()[1],width=0.025*axarr[0, 1].get_xlim()[1], color='r', alpha=0.5, label='Average')
    axarr[0, 1].set_xlim(min([axarr[0, 1].get_xlim()[0],OGMeanNearest]), axarr[0, 1].get_xlim()[1])

    
                            
    axarr[1, 0].hist(StandardDeviation, bins=20,color='b')
    axarr[1, 0].set_title('Standard Deviation Mean',fontsize=16)
    axarr[1, 0].annotate(str(round(OGStandardDeviation, 3)),
                    xy=(round(OGStandardDeviation, 3), 0), xycoords='data',
                    xytext=(axarr[1, 0].get_xlim()[1] * 0.8, axarr[1, 0].get_ylim()[1] * 0.99), textcoords='data',
                    arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",facecolor='red',
                                connectionstyle="angle3,angleA=0,angleB=-90"))
    axarr[1, 0].autoscale(enable=True, axis='x', tight=True)
    axarr[1, 0].bar(All_STD_Mean, axarr[1, 0].get_ylim()[1],width=(0.025*axarr[1, 0].get_xlim()[1])/3, color='r', alpha=0.5, label='Average')
    axarr[1, 0].set_xlim(min([axarr[1, 0].get_xlim()[0],OGStandardDeviation]), axarr[1, 0].get_xlim()[1])

    
    
    axarr[1, 1].hist(STDNearest, bins=20,color='b')
    axarr[1, 1].set_title('Standard Deviation Nearest Neighbor',fontsize=16)
    axarr[1, 1].annotate(str(round(OGSTDNearest, 3)),
                    xy=(round(OGSTDNearest, 3), 0), xycoords='data',
                    xytext=(axarr[1, 1].get_xlim()[1] * 0.8, axarr[1, 1].get_ylim()[1] * 0.99), textcoords='data',
                    arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",facecolor='red',
                                connectionstyle="angle3,angleA=0,angleB=-90"))
    axarr[1, 1].autoscale(enable=True, axis='x', tight=True)    
    axarr[1, 1].bar(All_Near_STD_Mean_STD, axarr[1, 1].get_ylim()[1],width=0.025*axarr[1, 1].get_xlim()[1], color='r', alpha=0.5, label='Average')    
    axarr[1, 1].set_xlim(min([axarr[1, 1].get_xlim()[0],OGSTDNearest]), axarr[1, 1].get_xlim()[1])


    
    plt.suptitle("Spatial Distribution Analysis",fontsize=24)

    fig.text(0.5, 0.05,"# Locations: %5.3f"%len(Latitude_IN) + "      Mean Distance: %5.3f"%OGMean +"      Mean Nearest Neighbor: %5.3f"%OGMeanNearest + "      Standard Deviation Distance: %5.3f"%OGStandardDeviation +\
             "      Standard Deviation Nearest Neightbor: %5.3f"%OGSTDNearest,ha='center',fontsize=12)
    
    fig.text(0.5, 0.03,"Units are in Kilometers.  Analysis on 10,000 Simulations.  Distance is calculated using the Haversine formula with Re = 6371km.",ha='center',fontsize=12)

    plt.savefig(figure_name_save+".png")
    plt.close()


Simulations = 10000
namelist_wps_filename = "namelist.wps"
domain_number =1
cases = ["2014-06-04","2014-06-08","2015-08-14"]
data_dir = os.path.abspath('../..')
for i in range(len(cases)):
    data_file = data_dir+"/Verification_Data/interpolation_pre/10m/"+cases[i]+"/all_hr_avg_trim.csv"
    calculate_randomness_distribution(Simulations, data_file, namelist_wps_filename, domain_number, geo_fmt="degrees",figure_name_save="Data_Randomness_"+cases[i].replace("-","_"))