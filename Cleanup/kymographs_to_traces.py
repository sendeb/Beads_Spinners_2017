from utilities import *
# # Make sure you ALREADY manually filtered the kymographs.
videos_dir = 'videos/'
video_name = '2016-11-04_spinners_switching_tz17_MM_met_OD_1_2'
fname = videos_dir + video_name
# fname = argv[1]:
kymograph_images = np.load('kymographs/' + video_name + '_kymographs.npy') # maybe have to do [()]
num_elems = len(kymograph_images)
bacterial_traces = []
for cell_num, processed_kymograph in enumerate(kymograph_images): #### NOTE: For now only 10 cells until we get things working! ####
    t0 = time.time()
    print 'Percent complete:', cell_num*100./num_elems, '%'
    print "step1", time.time() - t0
    trace = compute_trace(processed_kymograph)
    print "step2", time.time() - t0
    bacterial_traces.append(trace)

np.save('traces/' + video_name + '_traces', bacterial_traces)
# print_to_csv(bacterial_traces, 'test_csv', meta, tifname) # do later
