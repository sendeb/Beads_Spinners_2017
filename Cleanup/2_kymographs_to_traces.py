from utilities import *
# # Make sure you ALREADY manually filtered the kymographs.
#  Ex. python traces.py 1mM 2
# [Concentration] [File No. in paths dictionary]
# (See utilities.py)
video_name, videos_dir = get_video_path(sys.argv)
fname = videos_dir + video_name
tifname = fname + '.tif'
meta = metamorph_timestamps.get(tifname)
ang_chunks = 72
kymograph_images = np.load('kymographs/' + videos_dir + video_name + '_kymographs.npy')
num_elems = len(kymograph_images)
bacterial_traces = []

for cell_num, processed_kymograph in enumerate(kymograph_images): #### NOTE: For now only 10 cells until we get things working! ####
	print 'Saving!'
	t0 = time.time()
	print 'Percent complete:', cell_num*100./num_elems, '%'
	print "step1", time.time() - t0
	trace = compute_trace2(processed_kymograph)
	print "step2", time.time() - t0
	bacterial_traces.append(trace * 360.0 / ang_chunks)

	if len(sys.argv) >= 4 and sys.argv[3] == '--s':
			# Shows first 300 frames with trace on top after pressing yes
			plt.figure()
			plt.xlabel('Frame', fontsize=20)
			plt.ylabel('Angle', fontsize=20)
			plt.title('Kymograph with Position Detected', fontsize=20)
			plt.imshow(processed_kymograph[:,:300], cmap='gray')
			plt.plot(trace[:300], 'r-', lw = 3)
			plt.show()

np.save('traces/' + videos_dir + video_name + '_traces', bacterial_traces)
print_to_csv(bacterial_traces, 'traces/' + videos_dir + video_name + '_traces_csv', meta, tifname)
print('Sucessfully saved CSV!')
