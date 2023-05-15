#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Fernando Torres Leal
"""

#%% Needed libraries to run the program

import numpy as np
from numpy import exp, angle, abs,square, sqrt, exp, arctan2, pi,zeros, size, logical_and
from numpy.fft import fft2, fftshift,ifft2
from matplotlib.image import imread
from cv2 import resize

#%% Defining the functions

def rgb2gray(image): 
    """
    gray = rgb2gray(image) is a function that takes a 3-D array where each dimension corresponds to color in RGB and converts it to a grayscale.
    
    Parameters
    ----------
    image : 3-D array with shape (n,m,3)

    Returns
    -------
    gray : 2-D array with shape (n,m), corresponding to the scale conversion from rgb to gray

    """
    gray = 0.2990*image[:, :, 0] + 0.5870 * image[:, :, 1] + 0.1140*image[:, :, 2]
    return gray

def hogel_prep(hogel, scale, crop_percentage):
    """
    hogel = hogel_prep(hogel, scale, crop_prcentage)

    Parameters
    ----------
    hogel : 2-D array corresponding to a image hogel
    scale : Scaling factor for the size of every hogel in range [0,1]. This can be a 1-D array of different scaling factors to compute.
    crop_percentage : Parameter for croping the edges of an image hogel in range [0,1) where 0 represents 0% croping the edges and 1 croping 100% of the image. Value
    of 1 for this parameter is won't do anything to the hogel.

    Returns
    -------
    zero_hogel : 2-D array of the scaled and cropped image hogel

    """
    
    zero_hogel = zeros((size(hogel,0),size(hogel,1))) # Array where the cropped and scaled hogel will be placed
    mask = zeros((size(hogel,0),size(hogel,1))) # Mask array for cropping the hogel
    im = resize(hogel,(int(size(hogel,1)*scale),int(size(hogel,0)*scale))) # Scaled hogel
        
    # Indices for colocating the scaled hogel   
    i_min = int((size(hogel,0)-size(im,0))/2)
    i_max = int(i_min + size(im,0))
    j_min = int((size(hogel,1)-size(im,1))/2)
    j_max = int(j_min + size(im,1))
    
    zero_hogel[i_min : i_max , j_min : j_max ]= im
    
    condition = logical_and((crop_percentage != 0.0),(crop_percentage < 1.0))
    
    # Cropping the hogel according to the crop percentage
    if condition[0] and condition[1]:
    
        i_min_mask, i_max_mask = int((i_max+i_min)/2 - (1-crop_percentage[0])*size(im,0)/2),int((i_max+i_min)/2 + (1-crop_percentage[0])*size(im,0)/2)
        j_min_mask, j_max_mask = int((j_max+j_min)/2 - (1-crop_percentage[1])*size(im,1)/2), int((j_max+j_min)/2 + (1-crop_percentage[1])*size(im,1)/2)
    
        mask[i_min_mask : i_max_mask, j_min_mask :j_max_mask] = 1
    
        zero_hogel = zero_hogel * mask
    else:
        return zero_hogel
        
    
    return zero_hogel

    
def HogelGeneration(H_plane, scale, H_size_i, H_size_j, N_proji,N_projj, N_cicles, crop_percentage, res, directory_subimages,directory_transfer):
    """
    HogelGeneration(H_plane, scale, H_size_i, H_size_j, N_proji,N_projj, N_cicles, crop_percentage, position, res, directory1) is a function that generates and saves the arrays
    of the retrieved complex phase with the Gerchberch-Saxton Algorithm in a .npy file. It does not have any output since the arrays can be heavy.

    Parameters
    ----------
    H_plane : 2-D complex array where the hologram will fill with each hogel
    scale : Scaling factor for the size of every hogel in range [0,1]. This can be a 1-D array of different scaling factors to compute.
    H_size_i : Size of each hogel in i-index given a resolution
    H_size_j : Size of each hogel in j-index given a resolution
    N_proji :  Number of projections of the scene in i-index
    N_projj :  Number of projections of the scene in j-index
    N_cicles : Number of cicles to iterate in the Gerchberch-Saxton Algorithm
    crop_percentage : 
    res : Wanted resolution to print. For example 3600 dpi
    directory1 : Directory where 

    Returns
    -------
    Saves the files of the binary holograms as bmp in the given directory

    """
    #Sub_imag =  np.zeros([int(H_size_i*N_proji),int(H_size_j*N_projj)],dtype = float) 
    for h in range(len(scale)):

        ii, jj = 0, 0

        # Filling the holohram plane with hogels
        for i in iter(range(N_proji)):
            for j in iter(range(N_projj)):
                current_hogel = rgb2gray(imread(directory_subimages + '/Sub-image'+str(i)+str(j)+'.bmp'))
                current_hogel = hogel_prep(current_hogel, scale[h], crop_percentage)
                #Sub_imag[int(0 + ii*H_size_i): int(H_size_i  + ii*H_size_i) , int(0 + jj*H_size_j): int(H_size_j + jj*H_size_j)] = current_hogel
                H_plane[int(0 + ii*H_size_i): int(H_size_i + ii*H_size_i), int(0 + jj*H_size_j): int(H_size_j + jj*H_size_j)] = gerchberch_saxton(current_hogel, N_cicles, H_size_i, H_size_j)
                jj = jj + 1
            jj = 0
            ii = ii + 1

        np.save(directory_transfer + '/' + directory_subimages.split('_')[-1] + ',' + str(res) + ' dpi, scale =' + str(scale[h])+', '+str(N_proji)+'x'+str(N_projj) + 'binary.npy', H_plane)

def gerchberch_saxton(hogel, N_cicles,size_i,size_j):
    """
    A = gerchberch_saxton(hogel, N_cicles, size_i, size_j) is a function to the Gerchberch-Saxton Algorithm
    
    Parameters
    ----------
    hogel : 2-D array of shape (size_i,size_j) of one image hogel to apply the Gerchberch-Saxton Algorithm
    N_cicles : Number of cicles the algorithm is going to loop in order to recover the phase.
    size_i : Size in i-index of the hogel
    size_j : Size in j-index of the hogel

    Returns
    -------
    A : 2-D complex array of the recovered phase of the hogel 

    """
    
    A = ifft2(hogel,[size_i,size_j]) # Putting the image in the plane of the transfer function

    for k in range(N_cicles):
        B = exp(1j * angle(A)) # Obtain the field in the transfer function plane
        C = fft2(B,[size_i,size_j]) # Propagating into the image plane
        D = np.abs(hogel) * exp(1j * angle(C)) #Calculating the field in the image plane
        A = ifft2(D,[size_i,size_j])# Putting the field again in the transfer function
        print(k)

    return A

def propagator(X,Y,T,m,w0,zp,lambd,x0,y0):
    """
    Uz = propagator(X,Y,T,m,w0,zp,lambd,x0,y0) is a function that propagates a LG beam throught the transparency to
    retrieve the hologram. It is useful to verify the binary holograms before printing.
    
    Parameters
    ----------
    X : 2-D array spacial grid for de x discretization
    Y : 2-D array spacial grid for de y discretization
    T : 2-D array corresponding to the binary hologram or transfer function
    m : Topological charge of the LG, beam. Use value 0 for gaussian beam
    w0 : beam waist in meters (Small waist work better for high resolution binary holograms)
    zp : Propagation distance in meters
    lambd : Wavelength of the beam in meters
    x0 : Displacement in x of the center of the beam
    y0 : Displacement in x of the center of the beam

    Returns
    -------
    Uz : Propagation of the beam throught the transfer function

    """
    delta = X[0,1] -X[0,0] # Step size of the discretization
    
    k = 2.0*pi/lambd; # wave number
    
    E = lambda m,w0,x,y: (sqrt(square(x-x0)+ square(y-y0))/w0)**(abs(m)) * exp(-(square(x-x0)+ square(y-y0))/w0**2) * exp(1j * m * arctan2(y-y0,x-x0)) # LG beam function
    U0 = 20 * E(m,w0,X,Y) * T
    Uz = fftshift((1./(1j*lambd*zp)) * exp(1j*k*zp)* (delta * fftshift(fft2(U0))))
    
    return Uz



