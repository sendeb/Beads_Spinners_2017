from utilities import *
#  Ex. python 1_create_kymographs.py 1mM [Concentration] 2 [File No. in paths dictionary]
# (See utilities.py)
from matplotlib.patches import Circle
from matplotlib import animation

mpl.rc('figure',  figsize=(8, 5))

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




###############################################################
save_cell = False

def press_for_cycle(event):
  global save_cell
  if event.key == 'n':
     save_cell = False
     plt.close()
  if event.key == 'c':
     cycle()
  if event.key == 'escape':
    print('Exiting!')
    exit(0)
  if event.key == 'y':
    save_cell = True
  # else:
  #   plt.close()

def cycle():
    freq = 0.5 # how fast to switch between frames
    T = 10./2 # duration of frame switches
    for i in range(int(T/freq)):
        ax1.imshow(frame_0, aspect='equal')
        fig1.canvas.draw()
        plt.pause(freq)
        print 'i', i

        ax1.imshow(frame_1, aspect='equal')
        fig1.canvas.draw()
        plt.pause(freq)

#################################################################

if len(sys.argv) >= 4 and sys.argv[3] == '--s':
    Show=True
    print 'Will be showing frames...'
else:
    Show = False

bit_frames = []
for i in range(len(frames)):
    bit_frames.append(convert_to_8bit(frames[i]))
frames = np.array(bit_frames)

params = np.load('params' + video_name + '_params.npy')[()] # dict created in 0_save_params.py
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
F = ax.scatter(x=[c[0] for c in centers], y=[c[1] for c in centers], facecolors=len(centers)*["none"], color=len(centers)*["blue"], picker = 5)  # 5 points tolerance

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

#################################################################
#################################################################
#################################################################

radius = 8 # pixel radius of cell == length of filter
w, l = 2.45, radius # choose dimensions of rotating window
mymask = np.array([[w,0],[-w,0],[-w,l],[w,l]])

kymograph_images = []
for cell_num in range(num_elems): #### NOTE: For now only 10 cells until we get things working! ####
    t0 = time.time()
    print 'Percent complete:', cell_num*100./num_elems, '%'
    unprocessed_kymograph = build_kymograph(cell_num, frames, mymask, filtered_centers, Show=Show)
    print "step1", time.time() - t0
    # kymograph = invert_colors(unprocessed_kymograph) -- this line for black cells on white background
    processed_kymograph = process_kymograph(unprocessed_kymograph)
    kymograph_images.append(processed_kymograph)
    print "step2", time.time() - t0

np.save('kymographs' + video_name + '_kymographs', kymograph_images)
kymographs = np.load('kymographs' + video_name + '_kymographs.npy')
print('Sucessfully saved!')
