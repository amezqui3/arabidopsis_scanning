#!/usr/bin/env python
# coding: utf-8

# # Get petioles

# The petiole is the part that connects the leaf to the main stem. This notebook attempts to segment out leaves while keeping the petiole intact.
# 
# Based on previous discussions, I assume that the pots will be transparent or below the level of the plant, and the leaves will not go outside of the pot. Therefore I do not attempt to address these problems in this notebook.

# To obtain the .tif file used in this notebook from a stack of scans, enter the command
# 
# python3 seqToStack.py Col-0_Y_Slices_day12/Col-0\ Y_*.tif col_0.tif

# In[1]:


import os
import math
import importlib
import pandas as pd

import numpy as np
from matplotlib import pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

import scipy.ndimage as ndimage

import tifffile as tf

import arabidopsis_utils as thales

from PIL import Image, ImageDraw

from itertools import groupby


# In[2]:


dst = 'petiole/'
tiff_file  = 'u112-3.tif'


# In[3]:


_ , fname = os.path.split(tiff_file)
bname = os.path.splitext(fname)[0]
img = tf.imread(tiff_file)


# ## Separate pots

# Since we are not thresholding out air or soil, we create a 2-d image by taking the maximum value of each 0th coordinate column. After the first cell, the strategy follows Erik's code

# In[5]:


# 180:500 for col-0
# 210:500 for u112-3
max_image = np.max(img[210:500], axis = 0)

fig, ax = plt.subplots(figsize = (24,36))
im = ax.imshow(max_image, cmap='inferno', origin='lower', vmax =75)

plt.show()


# In[6]:


test4 = np.max(img[180:500], axis = 0)


# In[7]:


pts = max_image.astype(np.uint8)


# In[8]:


# be aggressive here, choice is only to get pot outline for bounding box.
pts[pts < 30] = 0 # 30 for col-0 35 for u112-3
pts[pts > 0] = 1
border = np.array([[-1, -1, -1],
                   [-1, 8, -1],
                   [-1, -1, -1]])

surface = ndimage.convolve(pts, border, np.int8, 'constant', cval=0)
surface[ surface < 0 ] = 0
surface = surface.astype(np.uint8)
surface[ surface > 0 ] = 1


# In[9]:


# fig, ax = plt.subplots(figsize = (24,36))
# im = ax.imshow(pts, cmap='inferno', origin='lower')

# plt.show()


# In[10]:


# get bounding box of all 4 plants
corners = thales.four_corners(pts)
tpt, rpt, bpt, lpt = corners
center = np.mean(corners,axis=0)

rhull = np.vstack((lpt,tpt,rpt,bpt,lpt))


# In[11]:


# plt.figure(figsize=(10,10))
# plt.imshow(surface, cmap='inferno', origin='lower');
# plt.plot(rhull[:,1], rhull[:,0], 'y--', lw=3);


# In[12]:


# get coordinates midpoints, centers (shown below)
bcoefs = thales.four_borderlines(corners)
midpoints = thales.four_midpoints(corners)
mcoefs = thales.midline_splits(midpoints)
centers = thales.four_centers(corners, midpoints)


# In[13]:


plt.figure(figsize=(10,10))
plt.imshow(surface, cmap='inferno', origin='lower');
for i in range(len(bcoefs)):
    plt.axline(corners[i,:][::-1], slope=bcoefs[i,1], c='y', lw=3)
    plt.plot(corners[i,1], corners[i,0], 'yo', ms=20)
    plt.plot(centers[i,1], centers[i,0], 'cX', ms=15)
    plt.plot(midpoints[i,1], midpoints[i,0], 'gh', ms=15)
for i in range(len(mcoefs)):
    plt.axline(midpoints[i,:][::-1], slope=mcoefs[i,1], c='g', lw=3)

#plt.plot(center[1], center[0], 'mD', ms=15);
#plt.tight_layout()

#filename = dst + bname + '_split_pots.jpg'
#plt.savefig(filename, dpi=96, pil_kwargs={'optimize':True})


# For each pot, create a mask that will show only that pot. Each pot is drawn clockwise, starting from top.

# In[14]:


# need to flip these around to draw the polygon
midpoints_reverse = midpoints.astype(int)
midpoints_reverse = midpoints_reverse[:,[1,0]]


# In[15]:


center_reverse = center.astype(int)
center_reverse = center_reverse[[1,0]]


# In[16]:


corners_reverse = corners.astype(int)
corners_reverse = corners_reverse[:,[1,0]]


# In[17]:


width = img.shape[1]
height = img.shape[2]


# In[18]:


mask = {}


# In[19]:


# draw each pot, based on the four corners, going clockwise
polygon0 = [tuple(corners_reverse[0]),tuple(midpoints_reverse[0]),tuple(center_reverse),tuple(midpoints_reverse[3])]

img0 = Image.new('L', (width, height), 0)
ImageDraw.Draw(img0).polygon(polygon0, outline=1, fill=1)
mask[0] = np.array(img0)


# In[20]:


polygon1 = [tuple(midpoints_reverse[0]),tuple(corners_reverse[1]),tuple(midpoints_reverse[1]),tuple(center_reverse)]

img1 = Image.new('L', (width, height), 0)
ImageDraw.Draw(img1).polygon(polygon1, outline=1, fill=1)
mask[1] = np.array(img1)


# In[21]:


polygon2 = [tuple(center_reverse),tuple(midpoints_reverse[1]),tuple(corners_reverse[2]),tuple(midpoints_reverse[2])]

img2 = Image.new('L', (width, height), 0)
ImageDraw.Draw(img2).polygon(polygon2, outline=1, fill=1)
mask[2] = np.array(img2)


# In[22]:


polygon3 = [tuple(midpoints_reverse[3]),tuple(center_reverse),tuple(midpoints_reverse[2]),tuple(corners_reverse[3])]

img3 = Image.new('L', (width, height), 0)
ImageDraw.Draw(img3).polygon(polygon3, outline=1, fill=1)
mask[3] = np.array(img3)


# An area of improvement would be to figure out a less repetetive way to code this (using a dictionary and/or for loop)

# In[23]:


# isolate each pot, by multiplying all non-pot pixels by zero
pots = {}

for i in range(0, 4):
    pots[i] = np.multiply(mask[i],max_image)


# In[24]:


fig, ax = plt.subplots(2,2, figsize=(15,15), sharex=True, sharey= True)

ax[0,0].imshow(pots[0], cmap='inferno', origin='lower', vmin = 19, vmax = 75);
ax[0,1].imshow(pots[1], cmap='inferno', origin='lower', vmin = 19, vmax = 75);
ax[1,0].imshow(pots[2], cmap='inferno', origin='lower', vmin = 19, vmax = 75);
ax[1,1].imshow(pots[3], cmap='inferno', origin='lower', vmin = 19, vmax = 75);

fig.suptitle(bname, fontsize = 20)

fig.tight_layout()
filename = dst + bname +'petiole_separated_pots_4_binary.jpg'
plt.savefig(filename, dpi=96, format='jpg', pil_kwargs={'optimize':True}, bbox_inches='tight')


# # Get conical volume of "allowed" voxels

# Now we create a second mask for each pot. The idea here is that an arabidopsis plant is roughly conical (albeit a very flat cone), so if all voxels outside the cone are zero, that will take away most of the air and soil.

# ### Now get conical mask for each pot

# In[25]:


# crop each pot, and get new centers


# In[26]:


# get stack of slices, but only show one pot
# use 500 as top of all pots

pot_stacks = {}

for i in range(0, 4):
    pot_stacks[i] = img[0:500]*mask[i]
# pot0_stack = img[0:upcutoff]*mask0 #[190:450]
# pot1_stack = img[0:upcutoff]*mask1 #[195:450]
# pot2_stack = img[0:upcutoff]*mask2 #[230:450]
# pot3_stack = img[0:upcutoff]*mask3 #[190:450]


# In[27]:


# cut out most of the black area, make image quicker to work with

for i in range(0, 4):
    pot_stacks[i] = thales.clean_zeroes(pot_stacks[i])[0]


# In[28]:


# to get new center of each pot, pick a random nonzero slice and use function from before
# the center of each plant is close enough to the center of the pot, so we will just center the cone there
def get_new_center(stack):
    pot_stack = stack[0]
    corners2 = thales.four_corners(pot_stack) # the 2 is just an artifact from function writing
    tpt, rpt, bpt, lpt = corners2
    center2 = np.mean(corners2,axis=0)
    center2_reverse = center2.astype(int)
    center2_reverse = center2_reverse[[1,0]]
    return center2_reverse


# In[29]:


centers = {}

for i in range(0, 4):
    centers[i] = get_new_center(pot_stacks[i])


# In[30]:


# https://stackoverflow.com/questions/44865023/how-can-i-create-a-circular-mask-for-a-numpy-array
def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask


# In[31]:


# get a circular mask that just shows the small center circle. When there is a massive dropoff in intensity here,
# that is where the plant starts
circles = {}
for i in range(0, 4):
    circles[i] = create_circular_mask(pot_stacks[i][0].shape[0],pot_stacks[i][0].shape[1], centers[i], 50)


# In[32]:


cylinders = {}
for i in range(0, 4):
    cylinders[i] = pot_stacks[i]*circles[i]


# In[33]:


avgs = {}
for i in range(0,4):
    avgs[i] = np.mean(cylinders[i], axis=(1,2))


# In[34]:


def find_sub_list(sl,l):
    sll=len(sl)
    for ind in (i for i,e in enumerate(l) if e==sl[0]):
        if l[ind:ind+sll]==sl:
            return ind,ind+sll-1


# In[35]:


dwcutoffs = {}

for i in range(0, 4):
    moving_average = [i for i in list(np.convolve(avgs[i], np.ones(10)/10, mode='valid'))]
    diffs = [j-i for i, j in zip(moving_average[:-1], moving_average[1:])]
    descending_lists = [k for g, k in ((i,list(j)) for i,j in (groupby(diffs, key=lambda x: x < 0))) if g == True and len(k) > 1]
    
    k = 0.0
    for j in descending_lists:
        if -sum(j)>k:
            k = -sum(j)
            descending_sublist = j
        else:
            pass
        
    dwcutoffs[i] = find_sub_list(descending_sublist, diffs)[1]+6


# In[36]:


#dwcutoffs = {}

# dwcutoffs[0] = 190
# dwcutoffs[1] = 195
# dwcutoffs[2] = 230
# dwcutoffs[3] = 190

# pot0_stack = img[0:500]*mask0 #[190:450]
# pot1_stack = img[0:500]*mask1 #[195:450]
# pot2_stack = img[0:500]*mask2 #[230:450]
# pot3_stack = img[0:500]*mask3 #[190:450]


# In[37]:


for i in range(0,4):
    pot_stacks[i] = pot_stacks[i][dwcutoffs[i]:]


# In[38]:


# get cone of allowed pixels
def get_conical_stack(pot_stack, center, radius_coeff):
    # make a copy so the input does not get changed (in case we want to change coefficient)
    pot_stack_copy = pot_stack.copy()
    
    # Multiply each layer of the pot stack by the conical mask
    for i in range(0, len(pot_stack)):
        circular_mask = create_circular_mask(pot_stack[0].shape[0],pot_stack[0].shape[1], center, radius_coeff*i)
        pot_stack_copy[i] = pot_stack[i]*circular_mask
        
    return pot_stack_copy   


# In[39]:


# if radius is too small, plant leaves will get cut off. Think 3d

max_pots = {}

for i in range(0, 4):
    max_pots[i] = np.max(get_conical_stack(pot_stacks[i], centers[i], 6), axis = 0)


# In[40]:


# max_pots_sideways1 = {}

# for i in range(0, 4):
#     max_pots_sideways1[i] = np.max(pot_stacks[i][:,200:700,:], axis = 1)


# In[41]:


# max_pots_sideways2 = {}

# for i in range(0, 4):
#     max_pots_sideways2[i] = np.max(pot_stacks[i][:,:,200:700], axis = 2)


# In[42]:


# single_pot = np.max(pot_stacks[2], axis = 0)

# fig, ax = plt.subplots(figsize = (24,36))
# im = ax.imshow(max_pots_sideways2[2], cmap='inferno', origin='lower', vmin = 22)

# plt.show()


# In[43]:


fig, ax = plt.subplots(2,2, figsize=(15,15), sharex=True, sharey= True)

ax[0,0].imshow(max_pots[0], cmap='inferno', origin='lower', vmin = 22);
ax[0,1].imshow(max_pots[1], cmap='inferno', origin='lower', vmin = 22);
ax[1,0].imshow(max_pots[2], cmap='inferno', origin='lower', vmin = 22);
ax[1,1].imshow(max_pots[3], cmap='inferno', origin='lower', vmin = 22);

fig.suptitle(bname, fontsize = 20)

fig.tight_layout()
filename = dst + bname +'petiole_separated_pots_4.jpg'
plt.savefig(filename, dpi=96, format='jpg', pil_kwargs={'optimize':True}, bbox_inches='tight')


# # Isolate individual leaves
# Look at the numpy array to find connected components in 3d

# In[44]:


pots_3d = {}

for i in range(0,4):
    pots_3d[i] = get_conical_stack(pot_stacks[i], centers[i], 6)
    pots_3d[i][pots_3d[i]<24] = 0
    pots_3d[i][pots_3d[i]>0] = 1


# In[45]:


structure = np.array([[[1,1,1],[1,1,1],[1,1,1]],[[1,1,1],[1,1,1],[1,1,1]],[[1,1,1],[1,1,1],[1,1,1]]])


# In[ ]:


snaps = {}

for i in range(0,4):
    labeled, ncomponents = ndimage.measurements.label(pots_3d[i], structure)
    
    unique, counts = np.unique(labeled, return_counts=True)
    large_components = list(np.where(counts>500)[0])
            
    labels_dict = {}

    for j in large_components:
        labeled2 = labeled.copy()
        labeled2[labeled2 != j] = 0 
        labels_dict[j] = labeled2
    
    snaps[i] = [thales.clean_zeroes(np.sum(labels_dict[j], axis=0))[0] for j in labels_dict.keys()]


# In[ ]:


_, axs = plt.subplots(3, 4, figsize=(18, 12))
axs = axs.flatten()
for img, ax in zip(snaps[2][1:], axs):
    ax.imshow(img, cmap='inferno', origin='lower')    

plt.show() 


# In[ ]:


#leaf_group = thales.clean_zeroes(labels_dict[11])[0]


# In[ ]:


#leaf_group[:, 100:150, 25:125][0]


# In[ ]:


# single_pot = np.max(pot_stacks[2], axis = 0)

# fig, ax = plt.subplots(figsize = (24,36))
# #im2 = ax.imshow(np.sum(leaf_group[:, 50:150, 25:125], axis=0), alpha = 0.5, cmap='inferno', origin='lower')
# #im = ax.imshow(np.sum(watershed_test, axis=0), cmap='inferno', origin='lower')
# im = ax.imshow(leaf_group[18], cmap='inferno', origin='lower')
# im2 = ax.imshow(np.sum(leaf_group, axis=0), alpha = 0.5, cmap='inferno', origin='lower')

# plt.show()


# In[ ]:


# top leaf: 50,240, 160
# bottom leaf: 50, 40, 305
# middle leaf: 18, 110, 110


# # Watershed

# In[ ]:


#from skimage.segmentation import watershed


# In[ ]:





# In[ ]:


#np.where(watershed_test !=1)


# In[ ]:


#markers = np.zeros(shape=leaf_group.shape)


# In[ ]:


# markers[50,240, 160] = 1
# markers[50, 40, 305] = 2
# markers[18, 110, 110] = 3


# In[ ]:


#watershed_test = watershed(leaf_group, markers=markers, connectivity=structure, offset=None, mask=None, compactness=0, watershed_line=False)


# In[ ]:


# single_pot = np.max(pot_stacks[2], axis = 0)

# fig, ax = plt.subplots(figsize = (24,36))
# im = ax.imshow(np.sum(watershed_test, axis=0), cmap='inferno', origin='lower')
# im2 = ax.imshow(np.sum(leaf_group, axis=0), alpha = 0.2, cmap='inferno', origin='lower')

# plt.show()


# In[ ]:




