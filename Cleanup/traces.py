from utilities import *
# # 
fname = './videos/2016-11-04_spinners_switching_tz17_MM_met_OD_1_2'
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

#possibly filter particles using ecc vals stationary cells will not look circular
f = tp.locate(avg, diameter=13, invert=False) #change 15 later, need to tune
f = f[(f['ecc'] < 0.1)]
# f.head() # shows the first few rows of data
fig, ax = plt.subplots()
fig.canvas.mpl_connect('key_press_event', press)
ax.set_title('Is this a good choice of parameters? If yes, press \'y\', else press ESC.')
tp.annotate(f, avg)

centers = []
# num_elems = len(f.x) # number of particles detected
num_elems = min(10,len(f.x)) # number of particles detected
xs = list(f.x)
ys = list(f.y)
for i in range(num_elems):
    x = xs[i]
    y = ys[i]
    center = [x, y]
    centers.append(center)
radius = 20 # pixel radius

w, l = 2.45, radius # choose dimensions of mask
mymask = np.array([[w,0],[-w,0],[-w,l],[w,l]])

bacterial_traces = []
kymograph_images = []
for cell_num in range(num_elems): #### NOTE: For now only 10 cells until we get things working! ####
    t0 = time.time()
    print 'Percent complete:', cell_num*100./num_elems, '%'
    kymograph = build_kymograph(cell_num, frames, mymask, centers)
    print "step1", time.time() - t0
    # kymograph = invert_colors(kymograph) -- this line for black cells on white background
    trace = compute_trace(kymograph_images, kymograph)
    print "step2", time.time() - t0
    bacterial_traces.append(trace)

np.save('traces_test', bacterial_traces)
# traces = np.load('traces_test.npy')

print_to_csv(bacterial_traces, 'test_csv', meta, tifname)
