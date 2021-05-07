#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Nov 24 19:23:41 2020

@author: X Zhou @ UZH patriarchi lab
"""



import numpy as np
import pandas as pd
import glob
import scipy.io


data_dir = '*'
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
NumofFram = 1000
GdFF = np.zeros(shape=(512,1))
for j in range(len(green_chan)):
    GdFF[j,0] = (np.mean(green_chan[j,Tmax-NumofFram:Tmax])-np.mean(green_chan[j,0:NumofFram]))/np.mean(green_chan[j,0:NumofFram])

scipy.io.savemat('*_RAW.mat',{'green_chan':green_chan,'red_chan':red_chan})

red_grad = np.gradient(np.mean(red_chan,0))


ligandT = np.argmax(red_grad,0)

responseGdFF = np.zeros(shape=(512,1))
for k in range(len(green_chan)):
    responseGdFF[k,0] = (np.mean(green_chan[k,ligandT-100:ligandT])-np.mean(green_chan[k,0:NumofFram]))/np.mean(green_chan[k,0:NumofFram])
    
    
boolean = responseGdFF[:,0]>0.65*max(responseGdFF)
GFPPix = green_chan[boolean].astype(np.int32)
TexRedPix = red_chan[boolean].astype(np.int32)
scipy.io.savemat('*_threshold65perMax.mat',{'GFP':GFPPix,'TexRedPix':TexRedPix})

