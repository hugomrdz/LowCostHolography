# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 21:38:23 2022

@author: Fernando
"""

import numpy as np
import matplotlib.pyplot as plt
from Function_Propagator import propagator
from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg
import os

#%%  Establishing directory for accesing and saving files

directoryTransfer = r'/home/fer/Desktop/Programs/3DHolograms/V8/Transfer function/'
directorySubimages ='/home/fer/Desktop/Programs/3DHolograms/V8/Sub_images/' 
directoryTransfer_bmp = r'/home/fer/Desktop/Programs/3DHolograms/V8/Transfer function bmp/'

N_TransFunction = len(os.listdir(directoryTransfer)) # Number of files in the directory

#%% Initializing parameters for deviating from optical axis

lambd = 632.8e-9; # beam wavelength
k = 2.0*np.pi/lambd; # wave number

N_in =5; #Number of inches of the hologram
res = 3600; # Printing resolution

Ni = N_in*res; Nj = N_in*res; # Number of points in the matrix [ij]
Ui = N_in*0.0254; Uj = N_in*0.0254; # Size in real units of the matrix

# Deviating angles that we found out tend to work better for these resolutions
if res == 2400:
    invangz,invangxy = 140.0 , 4.5
elif res == 3600:
     invangz,invangxy = 200.0 , 4.0

thz =  np.pi/invangz # Angle betweeen the Z-axis and the optical axis (Z-axis)
thxy = np.pi/invangxy  # Transversal angle XY plane

# Generating the components of the propagation dephasing vector
kx =k*np.sin(thz)*np.cos(thxy)
ky =k*np.sin(thz)*np.sin(thxy)

#%% Accesing the files and applying the plane wave to deviate the optical axis

for p in range(N_TransFunction):
    directoryTf = directoryTransfer + os.listdir(directoryTransfer)[p]
    directory11= directorySubimages + os.listdir(directorySubimages)[0]
    #N_proj = int(directory1.split(', ')[-1].split('x')[0]) # Number of hogels row,column, squared

    scale = os.listdir(directoryTransfer)[p].split('scale =')[1].split(',')[0]
    
    N_proji=max([max([int(string) for string in list(os.listdir(directory11)[k].split('Sub-image')[1].split('.bmp')[0])[0]]) for k in range(len(os.listdir(directory11)))])+1
    N_projj=max([max([int(string) for string in list(os.listdir(directory11)[k].split('Sub-image')[1].split('.bmp')[0])[1]]) for k in range(len(os.listdir(directory11)))])+1
    
    if N_proji >= N_projj:
        N_proj = N_proji 
    else:
        N_proj = N_projj
     
    H_size_i = int(Ni / (N_proj)) # Indicates the size for each hogel and dpi
    H_size_j = int(Nj / (N_proj)) # Indicates the size for each hogel and dpi
     
    H_plane = np.load(directoryTf) # Load the .npy file generated with mainHogelGeneration.py
    
    x = np.arange(-(H_size_j * N_projj)/2.,(H_size_j * N_projj)/2.,1) * (0.5*Uj) * (2./(H_size_j * N_projj))
    y = np.arange(-(H_size_i * N_proji)/2.,(H_size_i * N_proji)/2.,1) * (0.5*Ui) * (2./(H_size_i * N_proji))
    X,Y = np.meshgrid(x,y)
    
    ph_h = np.exp(1j*(kx*X + ky*Y)) # Plane wave to deviate the optical axis
    
    U = np.flipud(H_plane * ph_h) 
    
    Phase = np.angle(U);
    TGS = np.exp(1j* Phase);
    
    #% Generating the binary transfer function
    
    A = Phase;
    B = np.arcsin(np.abs(TGS)) /(np.max(np.abs(TGS)));
    T = 0.5 + 0.5*(np.sign(np.cos(A) + np.cos(B)));
    
    #%% Uncomment to test the transfer function

    # data = np.log(abs(propagator(X,Y,T,0,0.0015,2,lambd,0,0))**2)

    # #GUI control object
    # app = QtWidgets.QApplication(sys.argv)
    # #Window creation
    # window = pg.GraphicsLayoutWidget()
    # #Image object creation&Set image
    # image = pg.ImageItem() 
    # image.setImage(data)
    # #Create a box to store images&Set image object
    # view_box = pg.ViewBox()
    # view_box.addItem(image)
    # #Plot object creation&View created above_set box
    # plot = pg.PlotItem(viewBox=view_box)
    # #Add plot to window
    # window.addItem(plot)
    # #Window display
    # window.show()
    # #The end of the program
    # sys.exit(app.exec_())
    
    # plt.pcolormesh(np.log(abs(propagator(X,Y,T,0,0.254,2,lambd)))

    plt.imsave(directoryTransfer_bmp + directoryTf.split(',')[-4].split('/')[-1] +  ', off-axis angz= pi_' +str(invangz)+', angxy=pi_ ' + str(invangxy)+ ', ' +str(res) + ' dpi, scale =' + scale+', '+str(N_proji)+'x'+str(N_projj) +  'binary.bmp',T,cmap='gray')


#%%



    
