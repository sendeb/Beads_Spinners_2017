# Import packages
from __future__ import division, unicode_literals#, print_function
import numpy as np
from math import radians, sin, cos, floor
import matplotlib as mpl
import matplotlib.pyplot as plt
import trackpy as tp
import scipy
import pims
from scipy import interpolate, signal
import cv2
import time
import xml.etree.ElementTree
import tifffile
from datetime import datetime
import metamorph_timestamps
from scipy.ndimage.filters import median_filter
import sys
import os

mpl.rc('figure',  figsize=(16, 10))
mpl.rc('image', cmap='gray')

## Paths to our files. Update as data is recorded.

## TODO: notice how tif files now have prefix!!!
### We will have to change how code SAVES and OPENS tifs/params/kymos/etc so that pipeline runs smoothly.

concentrations = ['100nM', '1uM', '10uM', '100uM', '1mM', 'MotMed']

# Format strings as CONCENTRATION/CONCENTRATION_STREAMNAME.tif
paths = {
'1mM': [
    '1mM/1mM_Stream2_serine1milmol3msexposure_lampv7_55fps.tif',
    '1mM/1mM_Stream3_serine1milmol3msexposure_lampv7_55fps.tif', # move to correct folder
    '1mM/1mM_Stream133wssecondsafter_serine1milmol3msexposure_lampv7_55fps.tif'],
'1uM' : [
    '1uM/1uM_3ms_LV12_Stream0.tif',
    '1uM/1uM_3ms_LV12_Stream1.tif',
    '1uM/1uM_3ms_LV12_Stream2.tif'
    ],
'10uM' : [
    '10uM/10uM_3ms_LV12_Stream0.tif',
    '10uM/10uM_3ms_LV12_Stream1.tif',
    '10uM/10uM_3ms_LV12_Stream2.tif'
    ],
'100nM' : [
    '100nM/100nM_3ms_LV12_Stream0.tif',
    '100nM/100nM_3ms_LV12_Stream1.tif',
    '100nM/100nM_3ms_LV12_Stream2.tif'
    ],
'100um' : [
    '100uM/100uM_3ms_LV12_Stream0.tif',
    '100uM/100uM_3ms_LV12_Stream1.tif',
    '100uM/100uM_3ms_LV12_Stream2.tif'
    ],
'MotMed' : [
    'MotMed/MotMed_3ms_LV12_Stream0.tif',
    'MotMed/MotMed_Stream1_motmed3msexposure_lampv7_55fps.tif',
    'MotMed/MotMed_Stream2_motmed3msexposure_lampv7_55fps.tif', # actually in 1mM, move to correct folder
    'MotMed/MotMed_3ms_LV12_Stream1.tif',
    'MotMed/MotMed_3ms_LV12_Stream2.tif'
    ]
}

angs = []
# Set this!
for i in np.linspace(0, 360, 12): #select num. intervals per circle.
    angs.append(i)

def create_directories(list_of_directories):
    map(os.mkdir, filter(lambda dir : not os.path.isdir(dir), list_of_directories))

def get_video_path(args):
    # Arg1 = folder w/ concentration. Arg2 = stream number in that folder (see paths dict)
    path = paths[args[1]][int(args[2])]
    videos_dir = unicode.join('/', unicode.split(path, '/')[:-1])
    video_name = unicode.split(unicode.split(path, '/')[-1], '.')[0]
    return '/' + str(video_name), str(videos_dir)

# Path to tif stack 

def convert_to_8bit(image):
    # from int32
    im = image.astype(np.float64)
    im2 = (im - im.min())
    im2 = im2*255/im2.max()
    im2 = np.uint8(im2)
    return im2

def press(event):
  if event.key == 'n':
     plt.close()
  if event.key == 'escape':
    print('Exiting!')
    exit(0)
  else:
    plt.close()

def show(image):
    fig, ax = plt.subplots()
    fig.canvas.mpl_connect('key_press_event', press)
    ax.set_title('Is this a good choice of parameters? If yes, press \'y\', else press ESC.')
    plt.imshow(image)
    plt.show()

def create_rot_mask(orig_mask, deg):
    rotmymask = []
    ang = radians(deg)
    for i in orig_mask:
        rotMatrix = np.array([[cos(ang), -sin(ang)], [sin(ang),  cos(ang)]])
        rotmymask.append(list(rotMatrix.dot(i)))
    return np.array(rotmymask)

def adj_ctr_mask(mask, deg, cell_num, centers):
    adj_mask = []
    rotated_mask_about_00 = create_rot_mask(mask, deg)
    for mask_boundary in rotated_mask_about_00:
        # Move the origin of the mask to the center of the bacterium.
        # (Mask is originally at origin (0,0).)
        adj_mask.append(list(mask_boundary+centers[cell_num]))
    return np.array(adj_mask)


def angle_kym(ang, cell_num, frames, mymask, centers, Show=False):
    ang_ar=[]
    t0 = time.time()
    for i in range(frames.shape[0]):
        frame = frames[i].astype(np.uint8)
        box = np.int64(adj_ctr_mask(mymask, ang, cell_num, centers)) # this is the box rotated at ang deg.
        cv2.drawContours(frame,[box],0,(0,0,0),1)
	if Show:	
       	    if i == 0 and ang == 0: # shows the windows on top of 75th frame
       	        show(frame) # only showing filter do a 360 on first frame.
        mask = np.zeros(frame.shape,np.uint8)
        cv2.drawContours(mask,[box],0,1,-1) # cv2.drawContours(mask,[box],0,255,-1)
        ang_ar.append(cv2.mean(frame,mask=mask)[0])
    return ang_ar # for each frame, computes the pixel average for a window rot'd at a given theta 

def invert_colors(kymograph):
    kymograph -= kymograph.min()
    kymograph /= kymograph.max()
    return (1-kymograph) * 255


def build_kymograph(cell_num, frames, mask, centers, Show=False):
    kymograph = []
    for ang in angs:
        kymograph.append(angle_kym(ang, cell_num, frames, mask, centers, Show=Show))
    return np.array(kymograph)

def process_kymograph(kymograph):
    # Appends kymograph to kymograph_images and returns
    # computed trace from that saved kymograph image.
    trace = []
    # Remove background by subtracting median of each vertical column from itself.
    no_background=[]
    for i in range(kymograph.shape[1]):
        # black cells on white background (switch signs if reversed)
        no_background.append(kymograph[:,i]-np.median(kymograph,1)) 
    no_background=np.array(no_background).T
    
    # Change negative values to 0.
    clipped_background = no_background.clip(min=0)
    return clipped_background # the processed kymograph
    
def compute_trace(processed_kymograph):
    # Extract 1D signal using LA trick.
    eps = 1e-12
    def exp_func(x):
        return np.dot(np.arange(len(x)), np.power(x, 10))/(eps + np.sum(np.power(x, 10)))
    weighted_sum = np.apply_along_axis(exp_func, 0, processed_kymograph)
    
    # Derivative of 1D signal. Continuous parts show angular velocity of cell (not 100% sure on this.)
    # conv = np.convolve([-1.,1],weighted_sum, mode='full')[:-1]
    # median_filtered_conv = median_filter(conv, 7) #pick window size based on result. second arg and odd number.
    trace = weighted_sum
    return trace

def plot_kymograph(kymograph):
    plt.title('Kymograph', fontsize=20)
    plt.ylabel('Angles', fontsize=20)
    plt.xlabel('Frame', fontsize=20)
    plt.imshow(kymograph[:,:300])


#Returns the indices (frame locations) of when the sign switches.
def sign_switch(oneDarray):
    inds = []
    for ind in range(len(oneDarray)-1):
        if (oneDarray[ind]<0 and oneDarray[ind+1]>0) or (oneDarray[ind]>0 and oneDarray[ind+1]<0):
            inds.append(ind)
    return np.array(inds)

def get_date(path_to_tif):
    # Parse timestamp strings from the XML 
    # of the Metamorph metadata.
    # Returns an equivalent numpy array (in miliseconds).
    scope_times = []
    scope_tif_file = tifffile.TiffFile(path_to_tif)
    for t in range(len(scope_tif_file)):
        metadata = scope_tif_file[t].image_description
        root = xml.etree.ElementTree.fromstring(metadata).find('PlaneInfo')        
        for neighbor in root:
            if neighbor.attrib['type'] == 'time':
                if neighbor.attrib['id'] == 'acquisition-time-local':
                    first = neighbor.attrib['value']
                    return first


def print_to_csv(data, fname, meta, tifname):
    acquisition_time = get_date(tifname)
    with open(fname + ".csv", "wb") as f:
        f.write("Rotation Data (1D Kymograph Signal [radians vs time])," + '\n')
        f.write(tifname + ',\n')
        f.write(acquisition_time + ',\n,\n')
        header = "time (ms),"
        for i in range(len(data)):
            header += "cell" + str(i) + ","
        header += "\n"
        f.write(header)
        T = len(data[0])
        print "len of run " + str(T)
        Cells = len(data)
        print "len data " + str(Cells)
        for t in range(T):
            new_row_in_file = str(meta[t]) + ','
            for cell_data in range(Cells-1):
                new_row_in_file += (str(data[cell_data][t]) + ',')
            new_row_in_file += (str(data[Cells-1][t]))
            f.write(new_row_in_file + '\n')
    f.close()

# Goal:
# * Filter particles better
# * Organize code
# * Read data from .csv
