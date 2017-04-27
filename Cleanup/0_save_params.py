from utilities import *
# # 
#  Ex. python traces.py 1mM 2
# [Concentration] [File No. in paths dictionary]
# (See utilities.py)

video_name, videos_dir = get_video_path(sys.argv)
fname = videos_dir + video_name
tifname = fname + '.tif'
raw_frames = pims.TiffStack(tifname, as_grey=False)
frames = [np.fromstring(f.data, dtype=np.int16) for f in raw_frames] # int16 may have to change depending on dtype
frames = [np.reshape(f, (-1, raw_frames[0].shape[0], raw_frames[0].shape[1]) )[0] for f in frames]

bit_frames = []
for i in range(len(frames)):
    bit_frames.append(convert_to_8bit(frames[i]))
frames = np.array(bit_frames)

#avg is the average along the z axis of the image stack aka average image
avg = np.mean(frames, axis = 0)

# Set parameters.
diameter = 7 ## approximate size of object you're trying to locate.
ecc = 0.2
minmass = 150 ## if features same size, this is dimness
# topn = 20 # max number of cells
#possibly filter particles using ecc vals stationary cells will not look circular

f = tp.locate(avg, diameter=diameter, invert=False, minmass=minmass)
f = f[(f['ecc'] < ecc)]

# Uncomment below to view distribution of a value.
if len(sys.argv) >= 4 and sys.argv[3] == '--v':
    fig, ax = plt.subplots()
    ax.hist(f['mass'], bins=20)
    # Optionally, label the axes.
    ax.set(xlabel='mass', ylabel='count');
    fig.canvas.mpl_connect('key_press_event', press)
    ax.set_title('Is this a good choice of parameters? If yes, press \'y\', else press ESC.')
    plt.show()

### REMEMBER TO SAVE PARAMS IN DICTIONARY
#### ANDDDDD EXTRACT PARAMS IN CREATE_KYMOGRAPHS.py !!!!!!!!!!!!!

#check if mean image params are good: diameter, ecc, mimass
fig, ax = plt.subplots()
fig.canvas.mpl_connect('key_press_event', press)
ax.set_title('Is this a good choice of parameters? If yes, press \'y\', else press ESC.')
tp.annotate(f, avg)


# If got this far, save paramters to file: VIDEONAME + params .npy
# Parameters are the variables that result in the current
# annotated (circled) average image.
params = {'ecc' : ecc, 'diameter': diameter, 'minmass' : minmass}
np.save('params' + video_name + '_params', params)
print('Sucessfully saved!')
exit(0)
