#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 16:08:47 2024

@author: da
"""

import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
fname = '/media/da/DA5T/Cold_Pool/Datasets/SHIPS/lsdiaga_1982_2022_sat_ts_5day.txt'

cols = ['ID','Name','YYMMDD','UTC','Vmax','Lat','Lon','MSLP','SHRD','SHDC','SDDC']


num_lines = sum(1 for _ in open(fname))   
n_iter = int(num_lines/139)
df = pd.DataFrame(columns=cols,index = range(n_iter))
col0 = [str(i) for i in range(-12,126,6)] + ['Vars','temp']

# %%


tab = pd.read_fwf(fname,widths=[5]*25,header = None)
tab.columns = col0
tab.set_index('Vars',inplace=True) 

SHRD = tab.loc['SHRD','0'].astype(float)
SHDC = tab.loc['SHDC','0'].astype(float)
SDDC = tab.loc['SDDC','0'].astype(float)


data = pd.read_csv(fname,header = None)
ind = range(0,num_lines,139)

name0 = 'ALBE'
k = 0
for i,j in enumerate(ind):
    header = data.iloc[j,0].split()
    namei = header[0]
    if namei == name0:
        df.iloc[i,0] = k
    else:
        k = k + 1
        df.iloc[i,0] = k
        name0 = namei
    
    df.iloc[i,1:8] = header[0:7]
    df['SHRD'][i] = SHRD.iloc[i]
    df['SHDC'][i] = SHDC.iloc[i]
    df['SDDC'][i] = SDDC.iloc[i]
    
    
#%% save to csv
df.to_csv('/media/da/DA5T/Cold_Pool/Datasets/SHIPS/SHIPS_extracted.csv')
    
    

