# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 21:47:54 2022

@author: Fernando Torres Leal
"""

import numpy as np
import os, os.path
from BinaryStereographicHologramFunctions import HogelGeneration
#%% Initializing hologram plane parameters

directorySubimages= '/home/fer/Desktop/Programs/3DHolograms/V8/Sub_images/'
directoryTransfer = '/home/fer/Desktop/Programs/3DHolograms/V8/Transfer function'

N_in = 5; #Number of inches
res = 3600; # Printing resolution
Ni = N_in*res; Nj = N_in*res; # Number of points in the matrix [ij]
Ui = N_in*0.0254; Uj = N_in*0.0254; # Size in real units of the matrix

N_TransFunction = len(os.listdir(directorySubimages)) # Number of files in the directory
scale = np.array([0.45]) # Scale array for generating holograms with different scaling factors
N_cicles= 500 # Number of cicles of the Gerchberch-Saxton Algorithm
crop_percentage = np.array([0,0]) # Crop percentage for each hogel

for p in range(N_TransFunction):
    
    directory1 = directorySubimages + os.listdir(directorySubimages)[p] # Directory to access 
    
    # Calculating number of projections from reading the file. Just works if files are named as intended
    N_proji=max([max([int(string) for string in list(os.listdir(directory1)[k].split('Sub-image')[1].split('.bmp')[0])[0]]) for k in range(len(os.listdir(directory1)))])+1
    N_projj=max([max([int(string) for string in list(os.listdir(directory1)[k].split('Sub-image')[1].split('.bmp')[0])[1]]) for k in range(len(os.listdir(directory1)))])+1
    
    if N_proji >= N_projj:
        N_proj = N_proji
    else:
        N_proj = N_projj
    
    H_size_i = int(Ni / (N_proj)) # Indicates the size for each hogel and dpi
    H_size_j = int(Nj / (N_proj)) # Indicates the size for each hogel and dpi
    
    
    # Allocating hologram plane
    H_plane = np.zeros([int(H_size_i*N_proji),int(H_size_j*N_projj)],dtype = complex) # Allocating for the matrix
    
    #% Getting the hogels and filling the hologram plane
    HogelGeneration(H_plane,scale,H_size_i,H_size_j,N_proji,N_projj,N_cicles,crop_percentage,res,directory1,directoryTransfer)
    

    
            

        