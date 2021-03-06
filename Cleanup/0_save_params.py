from utilities import *
# # 
#  Ex. python traces.py 1mM 2
# [Concentration] [File No. in paths dictionary]
# (See utilities.py)

## Setting up directories to save data in
## expecting the working directory to contain folders with data from each concentration. See utilities.py.
required_directories = ["./params", "./features", "./traces", "./kymographs"]
create_directories(required_directories)
create_directories(concentrations)
for D in required_directories:
	create_directories(map(lambda concentration : D + '/' + concentration, concentrations))

video_name, videos_dir = get_video_path(sys.argv)
fname = videos_dir + video_name
tifname = fname + '.tif'
raw_frames = pims.TiffStack(tifname, as_grey=False)
frames = [np.fromstring(f.data, dtype=np.int16) for f in raw_frames] # int16 may have to change depending on dtype
frames = [np.reshape(f, (-1, raw_frames[0].shape[0], raw_frames[0].shape[1]) )[0] for f in frames]

#*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^
#*^*^*^*^*^*^*^*^*^*^*^       Inverting Image          *^*^*^*^*^*^*^*^*^*^*^*^*^*^
#*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^

# If we have BLACK cells on WHITE background, UNCOMMENT me!
frames = [np.invert(frame, dtype=np.int16) for frame in frames]

#*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^
#*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^
#*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^

bit_frames = []
for i in range(len(frames)):
    bit_frames.append(convert_to_8bit(frames[i]))
frames = np.array(bit_frames)

#avg is the average along the z axis of the image stack aka average image
avg = np.mean(frames, axis = 0)

##################################################################################
#######################       SET Parameters      ################################
##################################################################################

diameter = 3 ## approximate size in pixels of object you're trying to locate.
ecc = 0.7 # 0 means circular
minmass = 100 ## min integral of brightness for a particle

##################################################################################
##################################################################################
##################################################################################

f = tp.locate(avg, diameter=diameter, invert=False, minmass=minmass)
f = f[(f['ecc'] < ecc)]

# Uncomment below to view distribution of a value.
if len(sys.argv) >= 4 and sys.argv[3] == '--s':
    fig, ax = plt.subplots()
    feature_to_view = 'mass'
    ax.hist(f[feature_to_view], bins=20)
    # Optionally, label the axes.
    ax.set(xlabel=feature_to_view, ylabel='count')
    fig.canvas.mpl_connect('key_press_event', press)
    ax.set_title('Is this a good choice of parameters? If yes, press \'y\', else press ESC.')
    plt.show()

### REMEMBER TO SAVE ADDITIONAL PARAMS IN DICTIONARY
#### ANDDDDD CHANGE CODE WHICH EXTRACTS PARAMS IN 1_CREATE_KYMOGRAPHS.py !!!!!!!!!!!!!

#check if mean image params are good: diameter, ecc, mimass
fig, ax = plt.subplots()
fig.canvas.mpl_connect('key_press_event', press)
ax.set_title('Is this a good choice of parameters? If yes, press \'y\', else press ESC.')
tp.annotate(f, avg)


# If got this far, save paramters to file: VIDEONAME + params .npy
# Parameters are the variables that result in the current
# annotated (circled) average image.
params = {'ecc' : ecc, 'diameter': diameter, 'minmass' : minmass}
np.save('params/' + videos_dir + video_name + '_params', params)
print('Sucessfully saved!')
exit(0)
