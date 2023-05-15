"""
Blender Subimages generation

This code uses the Python API from Blender to generate different images from 
different angles of a scene and export its to a certain file.

@author: Fernando Torres Leal
"""

import bpy
import numpy as np
from math import radians


#%% Directories to access and save files

# Reccomended way to structure/organize the directories to run the program
current_object = "Lego"
file_format = ".stl"
directory_path = '/home/fer/Desktop/Programs/3DHolograms/DataBase/'
stl_path = directory_path + current_object + file_format
save_path_directory = "/home/fer/Desktop/Programs/3DHolograms/V8/Sub_images/"+"Sub_images"

#%% Deleting everything from scene and creating camera and illumination
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete()
bpy.ops.object.camera_add()
bpy.ops.object.light_add(type='POINT')

#%% Initializing sub images parameters

N_in = 5; #Number of inches
res = 3600;

Ni = N_in*res; Nj = N_in*res; # Number of points in the matrix [ij]
Ui = N_in*0.0254; Uj = N_in*0.0254; # Size in real units of the matrix

N_proj_i = 2 # N-1 of projections in i
N_proj_j = 7 # N-1 of projections in j

max_abs_azimut = radians(48);
max_abs_polar = radians(24)

H_size_i = Ni / (N_proj_i + 1) # Indicates the size for each hogel and dpi
H_size_j = Nj / (N_proj_j + 1) # Indicates the size for each hogel and dpi

if N_proj_i >= N_proj_j:
    N_proj = N_proj_i
else:
    N_proj = N_proj_j


azimut_range = np.arange(-N_proj_j/2,N_proj_j/2+1,1)* (2/N_proj_j) * (0.5*max_abs_azimut)
polar_range = np.linspace(np.pi/2 -max_abs_polar/2 ,np.pi/2 + max_abs_polar/2, N_proj_i+1)


#%% Establishing parameters of the scene

bpy.data.scenes["Scene"].render.resolution_x = H_size_j
bpy.data.scenes["Scene"].render.resolution_y = H_size_i

# Set world background to transparent
bpy.data.scenes["Scene"].render.film_transparent

#%% Importing scene, normalizing the size and centering the object

bpy.ops.import_mesh.stl(filepath= stl_path, directory= directory_path)
object = bpy.context.active_object

max_dim = np.max(np.abs(bpy.data.objects[current_object].dimensions))

bpy.data.objects[current_object].scale[0] = 1/max_dim
bpy.data.objects[current_object].scale[1] = 1/max_dim
bpy.data.objects[current_object].scale[2] = 1/max_dim

bpy.ops.object.origin_set(type="ORIGIN_CENTER_OF_VOLUME") 
#bpy.ops.object.origin_set(type="ORIGIN_GEOMETRY")

bpy.ops.object.align(align_mode='OPT_2', relative_to='OPT_1', align_axis={'Z'})
bpy.ops.object.align()

#%% Set camera parameters
camera_settings = bpy.data.objects["Camera"]
bpy.context.view_layer.objects.active = camera_settings
bpy.context.object.data.name = "Camera"
bpy.data.cameras["Camera"].type = "PERSP"
bpy.data.cameras["Camera"].sensor_width = 50
camera_settings.rotation_mode = "XYZ"
bpy.ops.object.select_all(action='DESELECT')

#%% Adding constraint to track the object with the camera
bpy.ops.object.constraint_add(type='TRACK_TO')
camera_settings.constraints["Track To"].target = bpy.data.objects[current_object]

camera_settings.rotation_euler[0] = radians(90)
camera_settings.rotation_euler[2] = radians(90)

#%% Establishing light parameters
bpy.context.view_layer.objects.active = bpy.data.objects["Point"]
bpy.context.object.data.name = "Point"
bpy.data.lights["Point"].type = "POINT"
bpy.data.lights["Point"].energy = 70 #(4e5)
bpy.data.lights["Point"].shadow_soft_size = 0.1 * np.max(bpy.data.objects[current_object].dimensions)
bpy.data.lights["Point"].use_shadow = False
bpy.data.objects["Point"].location[0] = bpy.data.objects[current_object].dimensions[0]*(1.8)
bpy.data.objects["Point"].location[1] = bpy.data.objects[current_object].dimensions[1]*(0.5)
bpy.data.objects["Point"].location[2] = bpy.data.objects[current_object].dimensions[2]*(1.3)

#%% Establishing render parameters
bpy.data.scenes["Scene"].render.image_settings.file_format="BMP"
bpy.data.scenes["Scene"].render.image_settings.color_mode="RGB"
bpy.data.scenes["Scene"].render.use_overwrite
bpy.data.scenes["Scene"].render.use_file_extension
bpy.data.scenes["Scene"].render.engine="BLENDER_EEVEE"
bpy.context.scene.camera = camera_settings


#%% Generating the different views

view_rad = 1.3* np.max(np.abs(bpy.data.objects[current_object].dimensions))

for i in range(N_proj_i +1):
    for j in range(N_proj_j +1):
        x_view = view_rad*np.cos(azimut_range[j])*np.sin(polar_range[i])
        y_view = view_rad*np.sin(azimut_range[j])*np.sin(polar_range[i])
        z_view = view_rad*np.cos(polar_range[i])
        
        camera_settings.location[0] = x_view
        camera_settings.location[1] = y_view
        camera_settings.location[2] = z_view      
             
        # Establishing dir for saving images
        save_path = save_path_directory + str(N_proj_i+1)+"x"+str(N_proj_j+1)+"_" + str(res)+"_"+ current_object + "/Sub-image"
        bpy.data.scenes["Scene"].render.filepath = save_path + str(i) + str(j)
        bpy.ops.render.render(write_still=True)