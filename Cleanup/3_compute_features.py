#from utilities import *
import numpy as np
import sys

# Input: some traces .npy (cols = time steps, rows = bact #)
#				2D arrays are a list of traces, each trace is one element in the list.
#
# Output: list of 2D arrays containing same number of rows
#			(one new time series for each bacterium (for each new engineered feature),
#			each array represents some feature, and
#			for each time step, the value of that feature for
#			the particular bacterium.

paths = {
'100nM' : ['100nM_3ms_LV12_Stream0_traces.npy',
'100nM_3ms_LV12_Stream1_traces.npy',
'100nM_3ms_LV12_Stream2_traces.npy'],
'100uM':['100uM_3ms_LV12_Stream0_traces.npy',
'100uM_3ms_LV12_Stream1_traces.npy',
'100uM_3ms_LV12_Stream2_traces.npy'],
'10uM' : ['10uM_3ms_LV12_Stream0_traces.npy',
'10uM_3ms_LV12_Stream1_traces.npy',
'10uM_3ms_LV12_Stream2_traces.npy'],
'1mM': ['1mM_Stream133wssecondsafter_serine1milmol3msexposure_lampv7_55fps_traces.npy',
'1mM_Stream2_serine1milmol3msexposure_lampv7_55fps_traces.npy',
'1mM_Stream3_serine1milmol3msexposure_lampv7_55fps_traces.npy'],
'1uM' : ['1uM_3ms_LV12_Stream0_traces.npy',
'1uM_3ms_LV12_Stream1_traces.npy',
'1uM_3ms_LV12_Stream2_traces.npy'],
'MotMed' : [ 'MotMed_3ms_LV12_Stream0_traces.npy',
'MotMed_3ms_LV12_Stream1_traces.npy',
'MotMed_3ms_LV12_Stream2_traces.npy',
'MotMed_Stream1_motmed3msexposure_lampv7_55fps_traces.npy',
'MotMed_Stream2_motmed3msexposure_lampv7_55fps_traces.npy'] }

def get_video_path(args):
    # Arg1 = folder w/ concentration. Arg2 = stream number in that folder (see paths dict)
    path = paths[args[1]][int(args[2])]
    videos_dir = unicode.join('/', unicode.split(path, '/')[:-1])
    video_name = unicode.split(unicode.split(path, '/')[-1], '.')[0]
    return '/' + str(video_name), str(videos_dir)

interval_bias = lambda s: np.sum((np.array(s)+1)/2)/len(s) # CCW / (CCW + CW); s is interval over which to compute bias, s is signs of rotation direction.

def compute_bias(trace, window=7):
  # Input: trace is single bacterium raw time series
  # 			window is size of widnow over which we (locally) compute bias.
  # Output: time series of bias at all overlapping intervals of window-length.  (len: len(trace) - window)

  bias = []
  
  # 1. Derivative of 1D signal. (Angular velocity) Use to get signs, which tell us CCW or CW. 
  conv = np.convolve([-1.,1], trace, mode='full')
  # Optionally: 
  #				median_filtered_conv = median_filter(conv, 7) # pick window size based on result. second arg is odd number.

  # 2. Get direction of rotation (CCW & CW)
  signs = np.sign(conv)
  # Optionally:
  #				filtered_signs = median_filter(signs, 5) # pick window size based on result. second arg is odd number.

  # 3. Compute bias over each window-length interval
  for i in range(len(trace) - window):
    interval = signs[i:i+window]
    bias.append(interval_bias(interval))
  return bias

def compute_mean_interswitch_interval(trace):
  # Output: real number representing average time between switches.
  conv = np.convolve([-1.,1], trace, mode='full')
  signs = np.sign(conv)

  # Compute the length of each run in one direction until switch.
  runs = []
  current_run = 0
  saved_most_recent_result = False
  for i in range(len(signs) - 1):
    s0 = signs[i]
    s1 = signs[i + 1]
    if s0 * s1 >= 0: # while signs same, increment the current run length
      # Doesn't this code just comptue the number of times the bacterium switches / number of time steps ?
      current_run += 1
      saved_most_recent_result = False
    else: # signs different ==> switched rotation direction! (s0 * s1 < 0)
      if current_run != 0: # might have started on a switch
        runs.append(current_run)
      saved_most_recent_result = True
      current_run = 0
  if not saved_most_recent_result:
    runs.append(current_run) # because the last interval may not be saved in runs (only saved if signs switch)
  return np.mean(runs)

def compute_cw_MISI(trace):
  # Output: real number representing average time between switches.
  conv = np.convolve([-1.,1], trace, mode='full')
  signs = np.sign(conv)

  # Compute the length of each run in one direction until switch.
  runs = []
  current_run = 0
  saved_most_recent_result = False
  for i in range(len(signs) - 1):
    s0 = signs[i]
    s1 = signs[i + 1]
    if s0 * s1 >= 0: # while signs same, increment the current run length
      # Doesn't this code just comptue the number of times the bacterium switches / number of time steps ?
      current_run += 1
      saved_most_recent_result = False
    else: # signs different ==> switched rotation direction! (s0 * s1 < 0)
      if current_run != 0: # might have started on a switch
        runs.append(current_run)
      saved_most_recent_result = True
      current_run = 0
  if not saved_most_recent_result:
    runs.append(current_run) # because the last interval may not be saved in runs (only saved if signs switch)
  return np.mean(runs)




def compute_features_for_each_trace():
  for trace in traces:
    biases = compute_bias(trace, window=7)
    mean_interswitch_intervals = compute_mean_interswitch_interval(trace)
    features['bias'].append(biases)
    features['mean_interswitch_interval'].append(mean_interswitch_intervals)

if __name__ == '__main__':
  traces = []
  concentrations = ['10uM']#,'100nM','1uM', '10uM', '100uM', '1mM', 'MotMed'] 
  for concentration in concentrations:
    video_name = '/' + concentration 
    for trace_path in paths[concentration]: # get all traces for this conc
      stream = np.load('traces/' + trace_path)
      for tr in stream:
        traces.append(tr)

    features = 	{
            'bias'								:	[],
            'mean_interswitch_interval'			:	[],
            'etc'								:	[]
          }

    # TODO:
    '''
    TODO: Remove paths dict and other code above and revert to:
    video_name, videos_dir = get_video_path(sys.argv)
    print video_name, videos_dir
    traces = np.load('traces' + video_name + '_traces.npy')
    '''

    compute_features_for_each_trace()
    # add prefix?
    np.save('features' + video_name + '_features', features)

    # To open, run:
    # features = np.load('features' + video_name + '_features.npy')[()] # Note: [()] allows us to load the dict() we saved.

    # plot distribution of bias and mean_is_interval for each concentration?
