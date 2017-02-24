from utilities import *
# # 
videos_dir = 'videos/'
video_name = '2016-11-04_spinners_switching_tz17_MM_met_OD_1_2'
fname = videos_dir + video_name
# fname = argv[1]
tifname = fname + '.tif'
meta = metamorph_timestamps.get(tifname)
raw_frames = pims.TiffStack(tifname, as_grey=False)
frames = [np.fromstring(f.data, dtype=np.int16) for f in raw_frames] # int16 may have to change depending on dtype
frames = [np.reshape(f, (-1, raw_frames[0].shape[0], raw_frames[0].shape[1]) )[0] for f in frames]

bit_frames = []
for i in range(len(frames)):
    bit_frames.append(convert_to_8bit(frames[i]))
frames = np.array(bit_frames)

avg = np.mean(frames, axis = 0)

# Set parameters.
diameter = 13
ecc = 0.3

#possibly filter particles using ecc vals stationary cells will not look circular
f = tp.locate(avg, diameter=diameter, invert=False)
f = f[(f['ecc'] < ecc)]
fig, ax = plt.subplots()
fig.canvas.mpl_connect('key_press_event', press)
ax.set_title('Is this a good choice of parameters? If yes, press \'y\', else press ESC.')
tp.annotate(f, avg)

# If got this far, save paramters to file: VIDEONAME + params .npy
# Parameters are the variables that result in the current
# annotated (circled) average image.
params = {'ecc' : ecc, 'diameter': diameter}
np.save('params/' + video_name + '_params', params)
print('Sucessfully saved!')
exit(0)