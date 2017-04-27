from utilities import *
# # Make sure you ALREADY manually filtered the kymographs.
#  Ex. python traces.py 1mM 2
# [Concentration] [File No. in paths dictionary]
# (See utilities.py)
video_name, videos_dir = get_video_path(sys.argv)
fname = videos_dir + video_name
tifname = fname + '.tif'
meta = metamorph_timestamps.get(tifname)
ang_chunks = 12

save_kymograph = False
def press2(event):
	global save_kymograph 
	if event.key == 'n':
		 save_kymograph = False
		 plt.close()
	if event.key == 'escape':
		print('Exiting!')
		exit(0)
	if event.key == 'y':
		save_kymograph = True 
		plt.close()

# fname = argv[1]
kymograph_images = np.load('kymographs' + video_name + '_kymographs.npy') # maybe have to do [()]
num_elems = len(kymograph_images)
bacterial_traces = []

for cell_num, processed_kymograph in enumerate(kymograph_images): #### NOTE: For now only 10 cells until we get things working! ####
	save_kymograph = False
	if len(sys.argv) >= 4 and sys.argv[3] == '--s':
		fig, ax = plt.subplots()
		fig.canvas.mpl_connect('key_press_event', press2)
		ax.set_title('Save kymograph? If yes, press \'y\', else press \'n\' to skip, or ESC to exit.')
	chunks = 4
	stacked = np.vstack(np.hsplit(processed_kymograph[:, : processed_kymograph.shape[1] - (processed_kymograph.shape[1] % chunks)], chunks))
	print processed_kymograph.shape
	print stacked.shape
	plt.imshow(stacked)
	plt.show()
		
	if save_kymograph: 
		print 'Saving!'
		t0 = time.time()
		print 'Percent complete:', cell_num*100./num_elems, '%'
		print "step1", time.time() - t0
		trace = compute_trace(processed_kymograph)
		print "step2", time.time() - t0
		bacterial_traces.append(trace * 360.0 / ang_chunks)

		# Shows first 500 frames with trace on top after pressing yes
		plt.figure()
		plt.xlabel('Frame', fontsize=20)
		plt.ylabel('Angle', fontsize=20)
		plt.title('Kymograph with Position Detected', fontsize=20)
		plt.imshow(processed_kymograph[:,:300], cmap='gray', extent=[0,300,2*np.pi,0], aspect="auto")
		# plt.plot(trace[:300], 'r-', lw = 3)
		plt.show()
	else:
		print 'Skipping!'

np.save('traces' + video_name + '_traces', bacterial_traces)
print_to_csv(bacterial_traces, 'traces' + video_name + '_traces_csv', meta, tifname) # do later
print('Sucessfully saved CSV!')
