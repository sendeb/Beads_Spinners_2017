from utilities import *
#  Ex. python 1_create_kymographs.py 1mM [Concentration] 2 [File No. in paths dictionary]
# (See utilities.py)
from matplotlib.patches import Circle
from matplotlib import animation



video_name, videos_dir = get_video_path(sys.argv)
fname = videos_dir + video_name
# fname = argv[1]
tifname = fname + '.tif'
meta = metamorph_timestamps.get(tifname)
raw_frames = pims.TiffStack(tifname, as_grey=False)
frames = [np.fromstring(f.data, dtype=np.int16) for f in raw_frames] # int16 may have to change depending on dtype
frames = [np.reshape(f, (-1, raw_frames[0].shape[0], raw_frames[0].shape[1]) )[0] for f in frames]

#need these for cycling through cells
frame_0 = frames[0]
frame_1 = frames[1]

if len(sys.argv) >= 4 and sys.argv[3] == '--s':
    Show=True
    print 'Will be showing frames...'
else:
    Show = False

bit_frames = []
for i in range(len(frames)):
    bit_frames.append(convert_to_8bit(frames[i]))
frames = np.array(bit_frames)

params = np.load('params/' + videos_dir + video_name + '_params.npy')[()] # dict created in 0_save_params.py
diameter = params['diameter']
ecc = params['ecc']
minmass = params['minmass']

avg = np.mean(frames, axis = 0)

#possibly filter particles using ecc vals stationary cells will not look circular
f = tp.locate(avg, diameter=diameter, invert=False, minmass=minmass) #change 15 later, need to tune
f = f[(f['ecc'] < ecc)]

#cycle through cells and filter the ones we want to keep. We want to do this so that the centers array is accurate and can be processed directly.

centers = []
num_elems = len(f.x) # number of particles detected
xs = list(f.x)
ys = list(f.y)
for i in range(num_elems):
    x = xs[i]
    y = ys[i]
    center = [x, y]
    centers.append(center)

filtered_set = set()


#################################################################
#################################################################
#################################################################

fig, ax = plt.subplots()
im=ax.imshow(frame_0, aspect='equal')
F = ax.scatter(x=[c[0] for c in centers], y=[c[1] for c in centers], s=240, facecolors=len(centers)*["none"], color=len(centers)*["blue"], picker=5)  # 5 points tolerance

def on_pick(event):
    ind = event.ind[0]
    print centers[ind], "clicked"
    filtered_set.add(tuple(centers[ind])) # use set to avoid duplicates being stored. (use tuple because can hash.)
    F._facecolors[ind,:] = (1, 0, 0, 0)
    F._edgecolors[ind,:] = (1, 0, 0, 1)
    fig.canvas.draw()

fig.canvas.mpl_connect('pick_event', on_pick)

def init():
    im.set_data(frame_0)

on0 = False
def animate(i):
    print "eye", i
    global on0
    if on0:
        im.set_data(frame_0)
    else:
        im.set_data(frame_1)
    on0 = not on0
    return im

anim = animation.FuncAnimation(fig, animate, init_func=init, interval=100)
plt.show()

# after closing plot:
filtered_centers = list(filtered_set)
filtered_centers = [list(c) for c in filtered_centers]

num_filtered_centers = len(filtered_centers)

#################################################################
#################################################################
#################################################################

radius = 8 # pixel radius of cell == length of filter
w, l = 2.45, radius # choose dimensions of rotating window
mymask = np.array([[w,0],[-w,0],[-w,l],[w,l]])

kymograph_images = []
for cell_num in range(num_filtered_centers): #### NOTE: For now only 10 cells until we get things working! ####
    t0 = time.time()
    print 'Percent complete:', cell_num*100./num_filtered_centers, '%'
    unprocessed_kymograph = build_kymograph(cell_num, frames, mymask, filtered_centers, Show=Show)
    print "step1", time.time() - t0
    # kymograph = invert_colors(unprocessed_kymograph) -- this line for black cells on white background
    processed_kymograph = process_kymograph(unprocessed_kymograph)
    kymograph_images.append(processed_kymograph)
    print "step2", time.time() - t0

np.save('kymographs/' + videos_dir + video_name + '_kymographs', kymograph_images)
print('Sucessfully saved!')
