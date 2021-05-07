#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 14:05:17 2020

@author: X Zhou @ UZH patriarchi lab
"""



import numpy as np
import pandas as pd
import glob
import scipy.io


data_dir = ''

my_files = glob.glob(data_dir + '/*.csv')

green_chan = np.zeros(shape=(512,len(my_files)))

red_chan = np.zeros(shape=(512,len(my_files)))




for i in range(len(my_files)):
    nn = my_files[i].split('_T');
    N = nn[1].split("_")
    n = int(N[0])
    data_F = pd.read_csv(my_files[i],usecols=(1,2),skiprows=8,header=None)
    data_f = data_F.to_numpy()
    green_chan[:,n] = data_f[:,0]
    red_chan[:,n] = data_f[:,1]

Tmax = green_chan.shape[1]
baslineFrame = 1000
GdFF = np.zeros(shape=(512,1))
for j in range(len(green_chan)):
    GdFF[j,0] = (np.mean(green_chan[j,Tmax-baslineFrame:Tmax])-np.mean(green_chan[j,0:baslineFrame]))/np.mean(green_chan[j,0:baslineFrame])
    
scipy.io.savemat('/*_RAW.mat',{'green_chan':green_chan,'red_chan':red_chan})


boolean = GdFF[:,0]>0.65*max(GdFF)
GFPPix = green_chan[boolean].astype(np.int32)
TexRedPix = red_chan[boolean].astype(np.int32)

scipy.io.savemat('/*_threshold65perMax.mat',{'GFP':GFPPix,'TexRedPix':TexRedPix})

