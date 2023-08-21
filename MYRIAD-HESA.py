
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 13:22:22 2023

@author: Judith Claassen - Institute for Environmental studies, VU University Amsterdam
"""

#MYRIAD-HESA 
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
#pd.options.mode.chained_assignment = None  # default='warn'
from itertools import compress
import time
import geopandas as gpd
import shapely
#%%
start_time1 = time.time()
acronym_map = {'eq': 'earthquake', 'fl' : 'flood', 'dr' : 'drought', 'ts' : 'tsunami',  'vo' : 'volcano', 'tc' : 'tropical cyclone', 'ls' : 'landslide',  'hw': 'heatwave', 'ew' : 'extreme wind', 'cw' : 'coldwave', 'wf' : 'wildfire'} 
hazards =  ["fl","ls"] #select the hazards to be includes based on acronyms
path = 'C:/.../polygondata/' #path to folder containing the data
timelag = False #if timelag is used timelag = True, if not timelag = False  
lag = 3 #time-lag in days
Area_interest = True #if an area of interest is used (for example a continent or country) Area_interest = True. else Area_interest = False. You can specify the area of interest with the minimum and maximum latitutudes and logitudes of a geographic region, or with a shapefile
#if you specify the area of interest (areai) using latitutes and longtitudes define areai here:
min_lon, min_lat, max_lon, max_lat =165, -20.2, 170, -12.5 #define the minimum and maximum latitude and longitude values of area of interest
areai = shapely.box(min_lon, min_lat, max_lon, max_lat) #create a rectangle geometry from the minimum and maximum latitude and longitude values
#if you want to use your own area of interest based on a shapefile, uncomment the following lines to upload your own shapefile as areai:
#shapefile_path = 'C:/pathtoshapefile/shapefile.shp'#path to the shapefile
#areai = gpd.read_file(shapefile_path).geometry[0] #loading shapefile
#%%Load data single hazard data
mintime = [] #first date registered in single hazard data
maxtime = [] #last date registered in single hazard data
for i in hazards:
    path2 = path + i + '.csv' #creating the path to single hazard data formated as "acronym.csv" for example "ls.csv" for landslides
    locals()[i] = pd.read_csv(path2) #load single hazard data
    locals()[i]['Geometry'] = shapely.from_wkt(locals()[i]['Geometry']) #format Geometry
    locals()[i]['starttime'] = pd.to_datetime(locals()[i]['starttime'], utc=True) #format starttime
    locals()[i]['endtime'] = pd.to_datetime(locals()[i]['endtime'], utc=True) #format endtime
    mintime.append(min(locals()[i].starttime)) #all min starttimes
    maxtime.append(max(locals()[i].endtime)) #all max endtimes
    locals()[i] = locals()[i].loc[:, ~locals()[i].columns.str.contains('^Unnamed')] #remove unnamed columns
    

#%% Selecting only overlapping timeframes & adding hazardtypes
mintime = max(mintime) #select the latest minimum time of all the hazards
maxtime = min(maxtime) #select the earliest minimum time of all the hazards 

for i in hazards:
    locals()[i] =  locals()[i][( locals()[i]['starttime'] >= mintime) & ( locals()[i]['endtime'] <= maxtime)].reset_index(drop=True) #remove the hazards that fall outside of the mintime and maxtime 
    locals()[i]['code'] = i+  locals()[i].index.astype(str) #unique hazard id
    locals()[i]['code2'] = i #hazard type acronym 
#%% Loading dynamic hazards (hazards with individual time-steps)
hazardtime = [] #dynamic hazards
hazardsi = [ 'hw', 'cw', 'ew', 'wf', 'dr']  #dynamic hazards acronyms 
hazardsi = list(set(hazards) & set(hazardsi)) #dynamic hazards loaded

for item in hazardsi:
    hazardtime.append(str(item) + 'active')

for i in hazardtime:  
    path2 = path + i + '.csv' #path to active data
    locals()[i] = pd.read_csv (path2) #load active hazard data formated as "acronymactive.csv" for example "hwactive.csv" for heatwaves
    locals()[i]['Geometry'] = shapely.from_wkt(locals()[i]['Geometry'])  #format Geometry
    locals()[i]['starttime'] = pd.to_datetime(locals()[i]['starttime'], utc=True) #format starttime
    locals()[i]['endtime'] = pd.to_datetime(locals()[i]['endtime'], utc=True) #format endtime
    locals()[i]['endtimelag'] = locals()[i]['endtime'] + timedelta(days = lag) #adding time-lag to an endtime
    unid = np.unique(locals()[i[:2]].Id) #selecting only dynamic hazard information for hazard events includded
    locals()[i] = locals()[i].loc[locals()[i]['Id'].isin(unid)] #selecting only dynamic hazard information for hazard events includded
    idd = np.unique( locals()[i[0:2]].Id) #selecting only dynamic hazard information for hazard events includded
    locals()[i] = locals()[i].loc[locals()[i]['Id'].isin(idd)] #selecting only dynamic hazard information for hazard events includded
    locals()[i] = locals()[i].loc[:, ~locals()[i].columns.str.contains('^Unnamed')] #remove unnamed column
#%% combining all hazards into one dataframe

for i in hazards:
    if i == hazards[0]:
        dfsingle =  locals()[i]
    else:
        dfsingle  = pd.concat([dfsingle, locals()[i]]).reset_index(drop=True)
    
dfsingle.sort_values(by='starttime', inplace=True) #sort by starttime
dfsingle = dfsingle.reset_index(drop=True) #reset index after sorting
dfsingle['index1'] = dfsingle.index #index of the hazards needed for pairs
dfsingle['endtimelag'] =dfsingle['endtime'] + timedelta(days = lag) #create endtime with time-lag
#%% Selecting data that intersects with the area of interest and clipping the data to the area of interest, if applicable

clippin = False #indicate if you want to clip the data to the area of interest. If True, the parts of the hazard polygon that fall outside of the area of interest will be removed.

if Area_interest == True:
   tree = shapely.STRtree(dfsingle.Geometry.values) #make tree too see Geometry overlap
   arr = np.transpose(tree.query(areai,predicate='intersects'))  #find intersecting hazards with the area of interest
   dfsingle = dfsingle.loc[np.sort(arr)].reset_index(drop=True) #remove the hazards that do not intersect with the area of interest
   dfsingle['index1'] = dfsingle.index #index of the hazards needed for pairs
   if clippin == True:
       dfsingle['Geometry'] = shapely.intersection(areai, dfsingle.Geometry) #remove part of the hazard geometry that falls outside of the area of interest

#%% uncomment if you want to check whether the geometries are still valid
#validgeom = []
#for i in dfsingle.Geometry:
#    validgeom.append(shapely.is_valid(i))

#if sum(validgeom) == len(dfsingle):
#    print('All geometries are valid')
#%% function to check overlap of dynamic hazards 

def intersectingdyn(activeii, activej, timelag): 
    try:
        t = False
        tdelta = []
        if timelag == True:
            activejendtime = activej.endtimelag
            activeiiendtime = activeii.endtimelag
        else:
            activejendtime = activej.endtime
            activeiiendtime = activeii.endtime    
        for k in range(0, len(activeii)):            
            activeji = activej[(((activej.starttime >= activeii.starttime[k]) & (activej.starttime <= activeiiendtime[k])) | ((activejendtime >= activeii.starttime[k]) & (activejendtime <= activeiiendtime[k])) | ((activej.starttime <= activeii.starttime[k]) & (activejendtime >= activeii.starttime[k])))].reset_index(drop=True) # select overlapping dates
            if len(activeji) > 0:
                tree2 = shapely.STRtree(activeji.Geometry)
                arr2 = np.transpose(tree2.query(activeii.Geometry[k],predicate='intersects')) 
                if len(arr2) > 0:
                    t = True                
            if t == True:
                tdelta = tdelta + list(((activeii.endtime[k] - activeji.starttime.loc[arr2]).dt.days * -1))
            if tdelta == []:
                t = False
            else:
                t = True
    except:
        t == False
    return t, tdelta   



#%% Making hazard pairs, list of two hazards that overlap in space and time
start_time = time.time() #starting timer
finalpairs = np.array([0,0,0,0]) #starttingpoint of final pairs
for i in range(len(dfsingle)): #loop over all the hazards
    if timelag == True:
        endtime = dfsingle.endtimelag #selecting correct endtime
    else:
        endtime = dfsingle.endtime #selecting correct endtime
    if i % 1000 == 0: #how often to print the time passed
        print('i = {}'.format(i), 'seconds', (time.time() - start_time)) #tracking how fast the pairs are being created
    dftimeoverlap = dfsingle[(((dfsingle.starttime >= dfsingle.starttime[i]) & (dfsingle.starttime <= endtime[i])) | ((endtime >= dfsingle.starttime[i]) & (endtime <= endtime[i])) | ((dfsingle.starttime <= dfsingle.starttime[i]) & (endtime >= dfsingle.starttime[i]))) & (dfsingle.index1 > dfsingle.index1[i])].reset_index(drop=True) # select overlapping dates with hazard event of interest
    tree = shapely.STRtree(dftimeoverlap.Geometry.values) #make tree too see geometry overlap
    arr = np.transpose(tree.query(dfsingle.Geometry[i],predicate='intersects')) #find intersecting hazards
    s = np.delete(arr, np.where(dftimeoverlap.index1[arr] == i), 0) #remove intersection with itself
    arr = np.array(dftimeoverlap.index1[arr])    #obtain the index of overlapping hazards in dfsingle        
    if len(arr) > 0: #check overlaps with dynamic hazards
        timedeltas = [] #timedelta of hazardi to hazardj
        timedeltas2 = [] #timedelta of hazardj to hazardi 
        if (dfsingle.code[i][:2] in hazardsi ): #if hazardi is dynamic            
            activei = locals()[dfsingle.code[i][:2] + 'active'].loc[locals()[dfsingle.code[i][:2] + 'active']['Id'] == dfsingle.Id[i]].reset_index(drop=True) #active polygons of hazardi
            if timelag == True:
                activeiendtime = activei.endtimelag #endtime of hazardi
            else:
                activeiendtime = activei.endtime #endtime of hazardi
        k = 0        
        while k < len(arr):
            j = arr[k]
            if (dfsingle.code[i][:2] in hazardsi ) and (dfsingle.code[j][:2] in hazardsi ): #if both hazardi and hazardj are dynamic
                activej = locals()[dfsingle.code[j][:2] + 'active'].loc[locals()[dfsingle.code[j][:2] + 'active']['Id'] == dfsingle.Id[j]] #active polygons arr
                activeii = activei[(((activei.starttime >= dfsingle.starttime[j]) & (activei.starttime <= endtime[j])) | ((activei.endtime >= dfsingle.starttime[j]) & (activei.endtime <= endtime[j])) | ((activei.starttime <= dfsingle.starttime[j]) & (activei.endtime >= dfsingle.starttime[j])))].reset_index(drop=True) #active polygon i with overlapping time with arr
                iss, timede =  intersectingdyn(activeii, activej, timelag) #see if dynamic hazards intersect in space and time
                if iss == True: #see the time delta if they do not intersect on same day
                    if 0 in timede: #if there is an overlap check the time-lag 
                        timedeltas.append(0)
                        timedeltas2.append(0)
                    elif any(n > 0 for n in timede) & any(n < 0 for n in timede):
                        negtimede = [x for x in timede if x < 0]
                        postimede = [x for x in timede if x > 0]
                        timedeltas.append(min(postimede))
                        timedeltas2.append(max(negtimede)*-1)
                    elif any(n > 0 for n in timede):
                        postimede = [x for x in timede if x > 0]
                        timedeltas.append(min(postimede))
                        timedeltas2.append(999)
                    elif  any(n < 0 for n in timede):
                        negtimede = [x for x in timede if x < 0]
                        
                        timedeltas.append(999)
                        timedeltas2.append(max(negtimede)*-1)                  
                    k+=1 #move to next itteration
                else:
                    arr = np.delete(arr, k) #if they do not intersect, remove pair from arr array
                    k = k #repeat itteration                                                        
            elif (dfsingle.code[i][:2] not in hazardsi ) and (dfsingle.code[j][:2] in hazardsi ): #if hazardi is not dynamic and hazardj is dynamic
                activej = locals()[dfsingle.code[j][:2] + 'active'].loc[locals()[dfsingle.code[j][:2] + 'active']['Id'] == dfsingle.Id[j]] #active polygons of hazardj
                if timelag == True:
                    activejendtime = activej.endtimelag #endtime of hazardj
                else:
                    activejendtime = activej.endtime #endtime of hazardj
                activej = activej[(((activej.starttime >= dfsingle.starttime[i]) & (activej.starttime <= endtime[i])) | ((activejendtime >= dfsingle.starttime[i]) & (activejendtime <= endtime[i])) | ((activej.starttime <= dfsingle.starttime[i]) & (activejendtime >= dfsingle.starttime[i])))].reset_index(drop=True) # select overlapping dates
                tree2 = shapely.STRtree(activej.Geometry)
                arr2 = np.transpose(tree2.query(dfsingle.Geometry[i],predicate='intersects')) #check intersecting polygons
                if len(arr2) > 0:   #if there is an overlap check the time-lag                                     
                    timedeltas.append(min((dfsingle.endtime[i] - activej.starttime.loc[arr2]) * -1).days)
                    if min((dfsingle.endtime[i] - activej.starttime.loc[arr2])* -1).days <=0:                        
                        timedeltas2.append(min((dfsingle.endtime[i] - activej.starttime.loc[arr2]) * -1).days)
                    else:
                        timedeltas2.append(999)
                    k+=1 #move to next itteration
                else:
                    arr = np.delete(arr, k) #if they do not intersect, remove pair from arr array
                    k = k  #repeat itteration                 
            elif (dfsingle.code[i][:2] in hazardsi ) and (dfsingle.code[j][:2] not in hazardsi ): #if hazardi is dynamic but hazardj is not
                activei = activei[(((activei.starttime >= dfsingle.starttime[j]) & (activei.starttime <= endtime[j])) | ((activeiendtime >= dfsingle.starttime[j]) & (activeiendtime <= endtime[j])) | ((activei.starttime <= dfsingle.starttime[j]) & (activeiendtime >= dfsingle.starttime[j])))].reset_index(drop=True) # select overlapping dates
                tree2 = shapely.STRtree(activei.Geometry)
                arr2 = np.transpose(tree2.query(dfsingle.Geometry[j],predicate='intersects')) #check intersecting polygons
                if len(arr2) > 0: 
                    timede = list((activei.endtime.loc[arr2] - dfsingle.starttime[j]).dt.days * -1)
                    ib = list(activei.starttime.loc[arr2] > dfsingle.endtime[j])
                    tfn =[]
                    tfp =[]
                    ttn =[]
                    for x in range(0, len(ib)): #if there is an overlap check the time-lag 
                        tfn.append((timede[x] < 1) & (ib[x] == False) )
                        tfp.append((timede[x] > 0) & (ib[x] == False ))
                        ttn.append((timede[x] < 0 )&( ib[x] == True ))
                    if any(tfn):
                        timedeltas.append(min(timede))
                        timedeltas2.append(min(timede))
                    elif any(tfp) == True & any(ttn) == True:
                        timedeltas2.append(max(list(compress(timede, ttn)))*-1)
                        timedeltas.append(min(list(compress(timede, tfp))))
                    elif any(tfp) == True:
                        timedeltas.append(min(list(compress(timede, tfp))))
                        timedeltas2.append(999)
                    elif any(ttn) == True:
                        timedeltas2.append(max(list(compress(timede, ttn)))*-1)
                        timedeltas.append(999)                    
                    k+=1 #move to next itteration
                else:
                    arr = np.delete(arr, k) #if they do not intersect, remove pair from arr array
                    k = k   #repeat itteration                              
            else: #neither hazardj or hazardi are dynamic
                timedeltas.append((dfsingle.endtime[i] - dfsingle.starttime[j]).days * -1)
                if (dfsingle.endtime[i] - dfsingle.starttime[j]).days * -1 <=0:
                    
                    timedeltas2.append((dfsingle.endtime[i] - dfsingle.starttime[j]).days * -1)
                else:
                    timedeltas2.append(999)
                k +=1 #move to the next itteration
        if len(arr) > 0:
            timedeltas =   np.array(timedeltas, dtype="object")
            timedeltas = np.resize(timedeltas, (len(timedeltas),))            
            timedeltas2 =   np.array(timedeltas2, dtype="object")
            timedeltas2 = np.resize(timedeltas2, (len(timedeltas2),))
            s1 = np.tile(i, (len(arr)))
            s2 = np.stack((s1, arr, timedeltas, timedeltas2), axis=1)
            finalpairs = np.vstack((finalpairs, s2))        
print("Creating the pairs took %s seconds" % (time.time() - start_time))
if len(np.unique(finalpairs)) == 1:
    finalpairs = np.array([])
else:
    finalpairs = finalpairs[1:]
    days = lag 
    k = (finalpairs[:,2] <= days) * 1
    k2 =( finalpairs[:,3] <= days  )* 1
    l = k + k2
    finalpairs = finalpairs[l >=1]  
    finalpairs = finalpairs[:, :2] 
 
#np.savetxt(path, finalpairs, delimiter=',', fmt="%s") #to save hazard pairs if needed

#%%Functions to make groups out of hazard pairs


def makegroups(sub, i):
    G = []
    unique4 = np.unique(sub)
    unique4 = unique4[unique4>i]
    for s in range(len(unique4)):
        j = unique4[s]
        f = np.where(sub == j)[0]
        n = np.unique(sub[f][sub[f]!=j])
        rowss = rows_in(n, sub)
        if len(f) == 1:
            coll = sub[f].tolist()[0]
            coll.append(i)
            if sorted(coll) not in G:
                G.append(sorted(coll))
        elif len(rowss) == (len(n)*(len(n)-1))/2:
            coll = n.tolist()
            coll.append(i)
            coll.append(j)
            if sorted(coll) not in G:
                G.append(sorted(coll))                        
        else:
            t4 = np.unique(sub[f])
            t4 = t4[t4!=j]
            rowss2 = rows_in(t4, sub)
            if len(np.where(sub[f][:,0] == j)[0]) == 0:
                continue                                
            else:
                for k in range(np.where(sub[f][:,0] == j)[0][0], len(f)):
                    ki = True
                    if k != 0:
                        t4 = np.append(t4, t4[0])
                        t4 = t4[1:]
                    if (np.unique(sub[f[k]][sub[f[k]]!=j]) < i) :
                        continue
                    n = np.unique(sub[f[k]][sub[f[k]]!=j])
                    f2 = np.where(rowss2 == n[0])[0]
                    if len(f2) == 0:
                        coll = list(n)
                        coll.append(i)
                        coll.append(j)
                        if sorted(coll) not in G:
                            G.append(sorted(coll))
                    else:
                        t5 = np.unique(rowss2[f2])
                        t5 = t5[t5!=n[0]]
                        t5 = t4[np.isin(t4,t5)]
                        for l in t5:                            
                            if k != l:
                                t = np.hstack((n,  l))
                                rowss = rows_in(t, sub)
                                if len(rowss) == (len(t)*(len(t)-1))/2:                                    
                                    if  (l > i) :
                                        n =  np.hstack((n, l))
                                    else:
                                        ki = False
                                        break
                        if ki == True:
                            coll = list(n)
                            coll.append(i)
                            coll.append(j)
                            if sorted(coll) not in G:
                                G.append(sorted(coll))
                    
    return G      

def rows_in(sub, finalpairs):
    testing = np.isin(finalpairs, sub)
    testing2 = np.sum(testing,axis=1)
    return finalpairs[testing2==2]

#%% making the Multi-hazard groups
G = [] #list for groups
unique, counts = np.unique(finalpairs[:,0], return_counts=True)  
start_time = time.time()
for j in range(0, len(unique)):
    i = unique[j]
    if j % 100000 == 0: #how often to print the time passed
        if j != 0: 
            print('j = {}'.format(j), 'seconds', (time.time() - start_time)) #tracking how fast the groups are being created
    f = np.where(finalpairs == i)[0]
    o = (np.where(finalpairs == i)[1] * -1) + 1
    if len(f) == 0:
        continue
    if len(f) ==1:
        T = finalpairs[f].tolist()
        G = G + T
        finalpairs = np.delete(finalpairs, f, 0)
    else:
        t = finalpairs[f,o]
        s1 = np.tile(i, (len(t)))
        s2 = np.stack((s1, t), axis=1)
        sub = rows_in(t, finalpairs)
        if len(sub) == 0:
            T = finalpairs[f].tolist()
            G = G + T
            finalpairs = np.delete(finalpairs,f, 0)
        else:
            unique2 = np.unique(sub)
            isinsub = np.isin(t, unique2)
            f1 = f[np.invert(isinsub)]
            if len(f1) > 0:
                T = finalpairs[f1].tolist()
                G = G + T
                finalpairs = np.delete(finalpairs, f1, 0)
            col = makegroups(sub, i)          
            G = G + col
   
print("Creating the groups took %s seconds" % (time.time() - start_time))  

G.sort()  
G = list(set(map(tuple, G)))
G = [list(row) for row in G]


#%%Making final dataset
path2= 'C:/Users/....csv' #path where final csv should be saved
dfmulti = pd.DataFrame(columns = dfsingle.columns)
for i in range(len(G)):
    j = G[i]
    dfmulti_i = dfsingle.loc[j].reset_index(drop=True)
    dfmulti_i['Event'] = 'event'+ str(i)
    dfmulti  = pd.concat([dfmulti, dfmulti_i]).reset_index(drop=True)  
print("--- %s seconds ---" % (time.time() - start_time))  


dfmulti['Hazard'] = dfmulti['code2'].replace(acronym_map)
dfmulti = dfmulti[['Event', 'Hazard' ,'index1', 'code', 'starttime', 'endtime', 'Intensity', 'Unit', 'Geometry']]
#uncommend following line to save dataset
#df5.to_csv(path2)

