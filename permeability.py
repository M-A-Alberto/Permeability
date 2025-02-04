# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 15:12:58 2022

@author: A. Martín-Asensio
"""

import os

import numpy as np

import matplotlib.pyplot as plt

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

import pandas as pd

from skimage import  io, filters



from skimage.measure import regionprops, label



from scipy import stats



from scipy import ndimage as ndi
mpl.style.use("classic")
plt.rcParams['figure.facecolor']='white'


os.chdir(r"XXXXX") #Input directory containing the time-lapse tif file with only ONE channel


def linear_fit(x,a,b):
    y = a*x+b
    return  y


#Plotting functions

def plot_images(images):
    
    
    fig = plt.figure()
    plt.imshow(images[0,:,:],vmax = 1000,cmap="Reds_r")
    fig.patch.set_facecolor('white')
    plt.colorbar()
    plt.axis("off")
    plt.savefig("./Images/First_image_NPs.tif",bbox_inches = "tight",dpi = 600)
    plt.close()
    
    fig = plt.figure()
    plt.imshow(images[0,:,:],vmax = 1000)
    fig.patch.set_facecolor('white')
    plt.colorbar()
    plt.axis("off")
    plt.savefig("./Images/First_image_NPs_colors.tif",bbox_inches = "tight",dpi = 600)
    plt.close()
        
    fig = plt.figure()
    plt.imshow(images[-1,:,:],vmax = 1000,cmap="Reds_r")
    fig.patch.set_facecolor('white')
    plt.colorbar()
    plt.axis("off")
    plt.savefig("./Images/Last_image_NPs.tif",bbox_inches = "tight",dpi = 600)
    plt.close()
    
        
    fig = plt.figure()
    plt.imshow(images[-1,:,:],vmax = 1000)
    fig.patch.set_facecolor('white')
    plt.axis("off")
    plt.colorbar()
    plt.savefig("./Images/Last_image_NPs_colors.tif",bbox_inches = "tight",dpi = 600)
    plt.close()

    return

def plot_images_axis(images):
    
    
    fig = plt.figure()
    plt.imshow(images[0,:,:],vmax = 1000,cmap="Reds_r")
    fig.patch.set_facecolor('white')
    plt.colorbar()
    plt.savefig("./Images/First_image_NPs_axis.tif",bbox_inches = "tight",dpi = 600)
    plt.close()
    
    fig = plt.figure()
    plt.imshow(images[0,:,:],vmax = 1000)
    fig.patch.set_facecolor('white')
    plt.colorbar()
    plt.savefig("./Images/First_image_NPs_colors_axis.tif",bbox_inches = "tight",dpi = 600)
    plt.close()
        
    fig = plt.figure()
    plt.imshow(images[-1,:,:],vmax = 1000,cmap="Reds_r")
    fig.patch.set_facecolor('white')
    plt.colorbar()
    plt.savefig("./Images/Last_image_NPs_axis.tif",bbox_inches = "tight",dpi = 600)
    plt.close()
    
        
    fig = plt.figure()
    plt.imshow(images[-1,:,:],vmax = 1000)
    fig.patch.set_facecolor('white')
    plt.colorbar()
    plt.savefig("./Images/Last_image_NPs_colors_axis.tif",bbox_inches = "tight",dpi = 600)
    plt.close()

    return

def plot_crops_axis(nuclei_img,images,x_max_ecm,x_min_ecm,y_max,y_min,x_min_channel,x_max_channel,x_max_pillars,x_min_pillars,side):
        
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(nuclei_img,vmax = 1000,cmap = "Blues_r")
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Crops_"+side+"_cells_axis.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(images[-1,:,:],vmax = 1000,cmap = "Reds_r")
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Last_Crops_"+side+"_NPs_axis.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(images[-1,:,:],vmax = 1000)
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Last_Crops_"+side+"_NPs_colors_axis.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(images[0,:,:],vmax = 1000,cmap = "Reds_r")
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Last_Crops_"+side+"_NPs_axis.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(images[0,:,:],vmax = 1000)
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Initial_Crops_"+side+"_NPs_colors_axis.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    return

def plot_crops(nuclei_img,images,x_max_ecm,x_min_ecm,y_max,y_min,x_min_channel,x_max_channel,x_max_pillars,x_min_pillars,side):
        
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(nuclei_img,vmax = 1000,cmap = "Blues_r")
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    plt.axis("off")
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Crops_"+side+"_cells.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(images[-1,:,:],vmax = 1000,cmap = "Reds_r")
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    
    ax.set_ylabel("x axis (Pixels)")
    plt.axis("off")
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Last_Crops_"+side+"_NPs.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(images[-1,:,:],vmax = 1000)
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    plt.axis("off")
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Last_Crops_"+side+"_NPs_colors.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(images[0,:,:],vmax = 1000,cmap = "Reds_r")
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    plt.axis("off")
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Last_Crops_"+side+"_NPs.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    
    fig, ax = plt.subplots()
    
    mapable = ax.imshow(images[0,:,:],vmax = 1000)
    plt.colorbar(mapable)
    ax.set_xlabel("y axis (Pixels)")
    plt.axis("off")
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crops")
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1.5, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    fig.patch.set_facecolor('white')
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1.5, edgecolor='limegreen', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1.5, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    plt.savefig("./Images/Initial_Crops_"+side+"_NPs_colors.tif",bbox_inches="tight",dpi = 600)
    plt.show()
    plt.close()
    return

#def main():

carpetas = ["Images","plots","Times"]
working_directory = os.getcwd()

#Create folders
for carpeta in carpetas:
    
    if os.path.exists(carpeta)==False:
        os.makedirs(carpeta)

#Read image
for direct in os.listdir("."):
    if ".tif" in direct:
        image = io.imread(direct)


N = image.shape[0]
nuclei_img = image[int(N/2),:,:]

#Plot image

fig = plt.figure()
plt.imshow(nuclei_img,vmax = 1000)
plt.xlabel("y axis (Pixels)")

plt.ylabel("x axis (Pixels)")
fig.patch.set_facecolor('white')
plt.colorbar()
plt.title("Original image")
plt.savefig("./Images/Original_image_nuclei.tif",bbox_inches = "tight",dpi = 600)

plt.show()
plt.close()

#Ask for time interval between pictures in s and
delta_t = float(input("time interval between pictures (s): "))

#Ask for the magnification used and calculate the calibration factor (Values for OrcaFlash v3)
magnification = {"10x":0.65,"20x":6.5/20,"40x": 6.5/40} #microns/pixel
mag = input("Magnification (10x, 20x, 40x): ")
c = magnification[mag+"x"]

#plot_images(images)

#plot_images_axis(images)


#Ask for the side in which the gel is respective to the central channel
side = input("Gel side (l/r): ")

cont = "n"

#Delimit ROIs
while cont == "n":
        
    y_min = int(input("Choose crop x_min: "))
    
    y_max = int(input("Choose crop x_max: "))
    
    
    
    x_min_ecm = int(input("Choose crop ECM y_min: "))
    x_max_ecm = int(input("Choose crop ECM y_max: "))
    
    
    
    fig, ax = plt.subplots()
    
    ax.imshow(nuclei_img,vmax = 1000)
        
    ax.set_xlabel("y axis (Pixels)")
    
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crop")
    
    rect = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1, edgecolor='r', facecolor='none')
    ax.add_patch(rect)
    fig.patch.set_facecolor('white')
    
    plt.show()
    A_ecm = (y_max-y_min)*(x_max_ecm-x_min_ecm)*c**2
    print("A_ecm = ",A_ecm,"microns^2")
    cont = input("Would you like to continue? ")
    if cont!="n":
        ECM_region = image[:,y_min:y_max,x_min_ecm:x_max_ecm]
    plt.close()
    

    

cont = "n"


while cont == "n":
        
    x_min_channel = int(input("Choose crop channel y_min: "))
    x_max_channel = int(input("Choose crop channel y_max: "))
    
    fig, ax = plt.subplots()
    
    ax.imshow(nuclei_img,vmax = 1000)
        
    ax.set_xlabel("y axis (Pixels)")
    
    ax.set_ylabel("x axis (Pixels)")
    
    ax.set_title("Crop")
    rect = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1, edgecolor='g', facecolor='none')
    ax.add_patch(rect)
    
    
    plt.show()
    cont = input("Would you like to continue? ")
    if cont!="n":
        channel_region = image[:,y_min:y_max,x_min_channel:x_max_channel]
    plt.close()
if side == "l":
    
    x_min_pillars = x_max_ecm+1
    x_max_pillars = x_min_channel-1
if side == "r":
    x_min_pillars = x_max_channel+1
    x_max_pillars = x_min_ecm-1

#Calculate pillars length and ROI
w = (y_max-y_min)*c

pillars_region = image[:,y_min:y_max,x_min_pillars:x_max_pillars]


#Find position of the channel to calculate the channel displacement over time
threshold = filters.threshold_li(nuclei_img)

mask_0 = nuclei_img > 200

mask_0 = ndi.median_filter(mask_0,30)

channel, n_features = label(mask_0,return_num = True)
    
channel_detected = regionprops(channel)

miny, minx, maxy, maxx = channel_detected[0].bbox

pos0 = minx

time = []

I_ecm = []

displacement = np.zeros(N)
positions = []
#I_ecm_bad = []


for i in range(image.shape[0]):
    
    #Track central channel and select ROIs
    detection = image[i,:,:]
    threshold = filters.threshold_li(detection)
    mask = detection > 200
    
    
    mask = ndi.median_filter(mask,30)
    channel, n_features = label(mask,return_num = True)
    
    channel_detected = regionprops(channel)
    
    miny,pos,maxy,maxx = channel_detected[0].bbox
    
    displacement[i] = pos-pos0
    positions.append(pos)
    x_min_ecm_new = int(x_min_ecm+displacement[i])
    x_max_ecm_new = int(x_max_ecm+displacement[i])
    
    x_min_channel_new = int(x_min_channel+displacement[i])
    x_max_channel_new = int(x_max_channel+displacement[i])
    
    x_min_pillars_new = int(x_min_pillars+displacement[i])
    x_max_pillars_new = int(x_max_pillars+displacement[i])
    
    ECM_ROI = image[i,y_min:y_max,x_min_ecm_new:x_max_ecm_new]
    
    Channel_ROI = image[i,y_min:y_max,x_min_channel_new:x_max_channel_new]
    
    Pillars_ROI = image[i,y_min:y_max,x_min_pillars_new:x_max_pillars_new]
    
    
    A_ecm = (y_max-y_min)*(x_max_ecm_new-x_min_ecm_new)*c**2
    
    
    
    I_f = np.sum(ECM_ROI)
    
    if i==0:
        I_0 = I_f
        I_vascular = np.sum(Channel_ROI)
        I_pillars = np.sum(Pillars_ROI)
        
        
    #Calculate time and ECM fluorescence intensity
    time.append(i*delta_t)
    I_ecm.append(I_f-I_0)
    
    
    if i>0:
        #Calculate gradient
        dI = np.gradient(I_ecm,time)
    
    
    plt.figure(figsize=(10,8))
    
    ax = plt.subplot(321)
    mapable = ax.imshow(detection,vmax = np.mean(detection)+np.std(detection))
    rect1 = patches.Rectangle((x_min_ecm_new,y_min), x_max_ecm_new-x_min_ecm_new,y_max-y_min, linewidth=1, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    rect2 = patches.Rectangle((x_min_channel_new,y_min), x_max_channel_new-x_min_channel_new,y_max-y_min, linewidth=1, edgecolor='blue', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars_new,y_min), x_max_pillars_new-x_min_pillars_new,y_max-y_min, linewidth=1, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    ax.axis("Off")
    ax.set_title("Image")
    
    
    ax2 = plt.subplot(322)
    mapable2 = ax2.imshow(Channel_ROI,vmax = np.mean(Channel_ROI)+np.std(Channel_ROI))
    ax2.axis("Off")
    ax2.set_title("Channel ROI")
    
    
    ax3 = plt.subplot(323)
    mapable3 = ax3.imshow(Pillars_ROI,vmax = np.mean(Pillars_ROI)+np.std(Pillars_ROI))
    ax3.axis("Off")
    ax3.set_title("Pillars ROI")
    
    
    ax4 = plt.subplot(324)
    mapable4 = ax4.imshow(ECM_ROI,vmax = np.mean(ECM_ROI)+np.std(ECM_ROI))
    ax4.axis("Off")
    ax4.set_title("ECM ROI")
    
    ax6 = plt.subplot(326)
    ax6.plot(time,I_ecm)
    ax6.set_xlabel("Time (s)")
    ax6.set_ylabel("Fluorescence Intensity (a. u.)")
    ax6.set_title("ECM_ROI Intensity")
    ax6_2 = ax6.twinx()
    
    if i>0:
        ax6_2.plot(time,dI,"g")
    ax6_2.set_ylabel("Intensity gradient (s$^{-1}$)")
    #ax6_2.set_ylim([np.min(dI),np.max(dI)])
    ax5 = plt.subplot(325)
    
    ax5.plot([x_max_channel_new-x_min_channel_new,x_max_channel_new-x_min_channel_new],[0,detection[y_min:y_max,x_min_channel_new:x_max_ecm_new].shape[0]],color = "yellow",lw = 2)
    ax5.plot([x_min_ecm_new-x_min_channel_new,x_min_ecm_new-x_min_channel_new],[0,detection[y_min:y_max,x_min_channel_new:x_max_ecm_new].shape[0]],color = "red",lw = 2)
    mapable5 = ax5.imshow(detection[y_min:y_max,x_min_channel_new:x_max_ecm_new],vmax = np.mean(detection[y_min:y_max,x_min_channel_new:x_max_ecm_new])+np.std(detection[y_min:y_max,x_min_channel_new:x_max_ecm_new]))
    ax5.axis("Off")
    ax5.set_title("ROIs merged")
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.suptitle(f"Time = {i*delta_t} s",fontsize = 18)
    
    plt.colorbar(mapable,ax = ax,location='right')
    
    
    plt.savefig("./Times/"+str(i)+"_full.tif",bbox_inches="tight")
    plt.show()
    plt.close()
    
    if i == int(N/2):
        
        ax = plt.subplot(111)
        mapable = ax.imshow(detection,vmax = 1000)
        rect1 = patches.Rectangle((x_min_ecm_new,y_min), x_max_ecm_new-x_min_ecm_new,y_max-y_min, linewidth=3, edgecolor='r', facecolor='none')
        ax.add_patch(rect1)
        rect2 = patches.Rectangle((x_min_channel_new,y_min), x_max_channel_new-x_min_channel_new,y_max-y_min, linewidth=3, edgecolor='blue', facecolor='none')
        ax.add_patch(rect2)
        rect3 = patches.Rectangle((x_min_pillars_new,y_min), x_max_pillars_new-x_min_pillars_new,y_max-y_min, linewidth=3, edgecolor='yellow', facecolor='none')
        ax.add_patch(rect3)
        plt.colorbar(mapable)
        
        plt.savefig(r"./Images/NPs_crops.tif",bbox_inches = "tight",dpi = 300)
        
    
    
    
    print(f"Time: {i*delta_t} s")
    print(f"# Features detected: {n_features}")
    print(f"Chip displacement (pixels): {displacement[i]}")
    print(f"New ECM ROI limits: {x_min_ecm_new} - {x_max_ecm_new} (length = {x_max_ecm_new-x_min_ecm_new} pixels)")
    print(f"New ECM ROI Area (um^2): {A_ecm}")
    print(f"New Channel ROI limits: {x_min_channel_new} - {x_max_channel_new} (length = {x_max_channel_new-x_min_channel_new} pixels)")
    print(f"New Pillars ROI limits: {x_min_pillars_new} - {x_max_pillars_new} (length = {x_max_pillars_new-x_min_pillars_new} pixels)")
    
    
    
    """ 
    #Code section to compare between chip tracking and not tracking
    ECM_ROI_bad = images[i,y_min:y_max,x_min_ecm:x_max_ecm]
    
    Channel_ROI_bad = images[i,y_min:y_max,x_min_channel:x_max_channel]
    
    Pillars_ROI_bad = images[i,y_min:y_max,x_min_pillars:x_max_pillars]
    
    I_f_bad = np.sum(ECM_ROI_bad)
    
    if i==0:
        I_0_bad = np.sum(ECM_ROI_bad)
        
        
    I_ecm_bad.append(I_f_bad-I_0_bad)   
    ax = plt.subplot(321)
    ax.imshow(detection,vmax = 1000)
    rect1 = patches.Rectangle((x_min_ecm,y_min), x_max_ecm-x_min_ecm,y_max-y_min, linewidth=1, edgecolor='r', facecolor='none')
    ax.add_patch(rect1)
    rect2 = patches.Rectangle((x_min_channel,y_min), x_max_channel-x_min_channel,y_max-y_min, linewidth=1, edgecolor='blue', facecolor='none')
    ax.add_patch(rect2)
    rect3 = patches.Rectangle((x_min_pillars,y_min), x_max_pillars-x_min_pillars,y_max-y_min, linewidth=1, edgecolor='yellow', facecolor='none')
    ax.add_patch(rect3)
    ax.axis("Off")
    
    ax2 = plt.subplot(322)
    ax2.imshow(Channel_ROI_bad,vmax = 1000)
    ax2.axis("Off")
    
    ax3 = plt.subplot(323)
    ax3.imshow(Pillars_ROI_bad,vmax = 1000)
    ax3.axis("Off")
    
    ax4 = plt.subplot(324)
    ax4.imshow(ECM_ROI_bad,vmax = 1000)
    ax4.axis("Off")
    
    
    ax6 = plt.subplot(326)
    ax6.plot(time,I_ecm,label = "Good")
    ax6.plot(time,I_ecm_bad,label = "Bad")
    ax6.set_xlabel("Time (s)")
    ax6.set_ylabel("Fluorescence Intensity (a. u.)")
    ax6.set_title("ECM_ROI Intensity")
    #ax6.legend()
    
    ax5 = plt.subplot(325)
    
    ax5.plot([x_max_channel-x_min_channel,x_max_channel-x_min_channel],[0,detection[y_min:y_max,x_min_channel:x_max_ecm].shape[0]],color = "yellow",lw = 2)
    ax5.plot([x_min_ecm-x_min_channel,x_min_ecm-x_min_channel],[0,detection[y_min:y_max,x_min_channel:x_max_ecm].shape[0]],color = "red",lw = 2)
    ax5.imshow(detection[y_min:y_max,x_min_channel:x_max_ecm],vmax = 1000)
    ax5.axis("Off")
    plt.tight_layout()
    
    plt.suptitle(f"Time = {i*delta_t} s")
    
    plt.savefig("./Times/"+str(i)+"_full_bad.jpg",bbox_inches="tight")
    plt.show()
    plt.close()
    
    """
    
    
    

plt.plot(np.array(time)/delta_t,I_ecm)

for i,t in enumerate(time):
    if i%10 == 0:
        plt.vlines(i,0,np.max(I_ecm))

plt.xlabel("Image number")
plt.show()
plt.close()


#Select curve region to perform linear fit

print(r"Select curve region to perform linear fit")

initial_t = int(input("Choose first image: "))
last_t = int(input("Choose last image: "))

save_data = pd.DataFrame(columns = ["Variable","Value"])

time = np.array(time)

#Perform linear fit and calculate permeability
I_fit = np.array(I_ecm[initial_t:last_t+1])
T_fit = np.array(time[initial_t:last_t+1])


Pars = stats.linregress(T_fit,I_fit)
a = Pars[0] #Pendiente
b = Pars[1] #Ordenada
y_fit = linear_fit(T_fit,a,b)
P_fit = a*A_ecm/(w*(I_vascular-I_pillars))
P_err = Pars[-1]*A_ecm/(w*(I_vascular-I_pillars))
plt.plot(time/60,I_ecm,label = "Original")

plt.plot(T_fit/60,I_fit, label = "Fitted cuve")

plt.plot(T_fit/60,y_fit,label = "Fit")

plt.legend(loc = "upper left")

plt.xlim([-0.5,(np.max(time)/60+0.5)])

plt.ylabel("$\Delta$I (a.u.)")
plt.xlabel("Time (min)")

plt.savefig("./plots/Linear fit.tif",bbox_inches = "tight",dpi = 600)
plt.savefig("./plots/Linear fit_low_res.jpg",bbox_inches = "tight")

plt.show()
plt.close()

Delta_t = np.max(time)

I_final = I_ecm[last_t+1]

I_i = I_ecm[initial_t]

P = (A_ecm*(I_f-I_i)/Delta_t)/(w*(I_vascular-I_pillars))



variables = ["Magnification (microns/pixel)","x_min (pixels)","x_max (pixels)","y_min_ECM (pixels)","y_max_ECM (pixels)","y_min_channel (pixels)","y_max_channel (pixels)","y_min_pillars (pixels)","y_max_pillars (pixels)","A_ECM (microns^2)","w_pillars (microns)","I_f (a.u.)","I_i (a.u.)","I_vascular (a.u.)","I_pillars (a.u.)","Experiment time (Delta t, s) ","Permeability (microns/s)","Time interval between pictures (delta_t,s)","First image","Last image","Slope","err","Intercept","R^2","pvalue","Permeability_fit (microns_s)","Permeability_err"]
values = [c,y_min,y_max,x_min_ecm,x_max_ecm,x_min_channel,x_max_channel,x_min_pillars,x_max_pillars,A_ecm,w,I_final,I_0,I_vascular,I_pillars,Delta_t,P,delta_t,initial_t,last_t,a,Pars[-1],b,Pars[2],Pars[3],P_fit,P_err]
save_data.loc[:,"Variable"] = variables
save_data.loc[:,"Value"] = values
print(save_data)


#Save results
save_data.to_csv("Experiment log_"+side+".txt",sep = "\t",decimal = ".")

positions_data = pd.DataFrame()

positions_data.loc[:,"Positions"] = positions

positions_data.to_csv("Displacement.txt",sep = "\t",decimal = ".")

#Save Intensity data
curve_I = pd.DataFrame()

curve_I.loc[:,"t"] = time
curve_I.loc[:,"I"] = I_ecm

curve_I_fit = pd.DataFrame()


curve_I_fit.loc[:,"T_fit"] = T_fit
curve_I_fit.loc[:,"I_fit"] = I_fit
curve_I_fit.loc[:,"fit"] = y_fit


curve_I.to_csv("Curve.txt",sep = "\t",decimal = ".",index = None)
curve_I_fit.to_csv("Curve_fit.txt",sep = "\t",decimal = ".",index = None)


