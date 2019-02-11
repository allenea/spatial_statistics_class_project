#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 15:13:44 2018

@author: ericallen
"""
import numpy as np
from math import sin, cos, sqrt, asin,pi, acos

def data_randomness(Latitude,Longitude,geo_fmt="degrees"):
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

    #CONSTANT
    degrad = pi / 180.0 
    n = len(Latitude)
    
    #INITIALIZE
    d_list = []
    nearestNeighbor = []
    
    for k in range(n):   ## range looping is new here
        if geo_fmt == "dms":
            #loc = Station[k]
            ltdeg =float(Latitude[k].split(".")[0])
            ltmin = float(Latitude[k].split(".")[1])
            ltsec = float(Latitude[k].split(".")[2])
            lgdeg = float(Longitude[k].split(".")[0])
            lgmin = float(Longitude[k].split(".")[1])
            lgsec = float(Longitude[k].split(".")[2])
            ddlong = lgdeg + lgmin / 60.0 + lgsec / 3600.0
            ddlat = ltdeg + ltmin / 60.0 + ltsec / 3600.0
        
            rlong = ddlong * degrad
            rlat = ddlat * degrad
        elif geo_fmt == "degrees":
            #loc = Station[k]
            rlong = Longitude[k] * degrad
            rlat = Latitude[k] * degrad
    
    
        stDistance=[]
        for j in range(n):
            if k != j:
                if geo_fmt == "dms":
                    #loc2 = Station[j]
                    ltdeg2 =float(Latitude[j].split(".")[0])
                    ltmin2 = float(Latitude[j].split(".")[1])
                    ltsec2 = float(Latitude[j].split(".")[2])
                    lgdeg2 = float(Longitude[j].split(".")[0])
                    lgmin2 = float(Longitude[j].split(".")[1])
                    lgsec2 = float(Longitude[j].split(".")[2])
                    ddlong2 = lgdeg2 + lgmin2 / 60.0 + lgsec2 / 3600.0
                    ddlat2 = ltdeg2 + ltmin2 / 60.0 + ltsec2 / 3600.0
                    
                    rlong2 = ddlong2 * degrad
                    rlat2 = ddlat2 * degrad
                elif geo_fmt == "degrees":
                    #loc2 = Station[j]
                    rlong2 = Longitude[j] * degrad
                    rlat2 = Latitude[j] * degrad
                    
                #Haversine
                d = haversine(rlong, rlat, rlong2, rlat2)
                
                #Law of Cosine
                #d = acos(((sin(rlat)*sin(rlat2))+(cos(rlat)*cos(rlat2)*cos(rlong-rlong2))))
                
                
                stDistance.append(d)
                d_list.append(d)
            
        nearestNeighbor.append(min(stDistance))
    Mean = np.mean(d_list)
    MeanNearest = np.mean(nearestNeighbor)
    StandardDeviation = np.std(d_list)
    STDNearest = np.std(nearestNeighbor)
    print (len(Latitude)," locations")
    print ("Mean Distance:                         %10.3f"%Mean)
    print ("Mean Nearest Neighbor:                 %10.3f"%MeanNearest)
    print ("Standard Deviation Distance:           %10.3f"%StandardDeviation)
    print ("Standard Deviation Nearest Neightbor:  %10.3f"%STDNearest)
    
    return Mean, MeanNearest,StandardDeviation, STDNearest