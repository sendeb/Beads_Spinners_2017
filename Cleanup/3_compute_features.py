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
fps = 55.0 # NOTE: CHANGE ME!!!! :NOTE

def get_video_path(args):
    # Arg1 = folder w/ concentration. Arg2 = stream number in that folder (see paths dict)
    path = paths[args[1]][int(args[2])]
    videos_dir = unicode.join('/', unicode.split(path, '/')[:-1])
    video_name = unicode.split(unicode.split(path, '/')[-1], '.')[0]
    return '/' + str(video_name), str(videos_dir)

interval_bias = lambda s: np.sum((-np.array(s)+1)/2)/len(s) # CCW / (CCW + CW); s is interval over which to compute bias, s is signs of rotation direction. NOTE: correct if cw is positive, ccw is negative.

def compute_bias(trace, window=7, first=None):
  # Input: trace is single bacterium raw time series
  # 			window is size of widnow over which we (locally) compute bias.
  # Output: time series of bias at all overlapping intervals of window-length.  (len: len(trace) - window)
  # first : if not None, bias on first x frames.

  bias = []
  
  # 1. Derivative of 1D signal. (Angular velocity) Use to get signs, which tell us CCW or CW. 
  conv = np.convolve([-1.,1], trace, mode='full') # SHOULD WE UNWRAP?
  # Optionally: 
  #				median_filtered_conv = median_filter(conv, 7) # pick window size based on result. second arg is odd number.

  # 2. Get direction of rotation (CCW & CW)
  signs = np.sign(conv) #Positive values correspond to cw rotation. Negative = ccw rotation.

  # Optionally:
  #				filtered_signs = median_filter(signs, 5) # pick window size based on result. second arg is odd number.

  # 3. Compute bias over each window-length interval
  # no sliding window as here:
  if first is None: # use overlapping window over whole trace 
    for i in range(len(trace) - window):
      interval = signs[i:i+window]
      bias.append(interval_bias(interval))
    interval = signs
    bias = interval_bias(interval)
    return bias
  else: # use first 'first' frames to compute bias
    interval = signs[:first]
    bias = interval_bias(interval)
    return bias

def compute_MISI(trace, frames=1700):
  # Output: real number representing average time between switches.
  conv = np.convolve([-1.,1], trace, mode='full')
  signs = np.sign(conv)[:frames]

  ccw_lengths, cw_lengths = list(), list()
  ccw_cur_length, cw_cur_length = 0, 0

  on_cw = False
  on_ccw = False

  for i in range(len(signs) -1):
    # Started a run
    s = signs[i]
    if np.isnan(s):
      print' isn an'
    s_next = signs[i+1]

    if not on_cw and not on_ccw: #Positive values correspond to ccw rotation. Negative = cw rotation.
      if s == 1 or s == 0:
        on_cw = True
      else:
        on_ccw = True
    if on_cw:
      cw_cur_length += 1
    elif on_ccw:
      ccw_cur_length += 1

    # Switch
    if (s == 1 or s == 0) and s_next == -1:
      on_cw = False
      cw_lengths.append(cw_cur_length)
      cw_cur_length = 0
      on_ccw = True
    elif (s == -1) and (s_next == 1 or s_next == 0):
      on_ccw = False
      ccw_lengths.append(ccw_cur_length)
      ccw_cur_length = 0
      on_cw = True
  return np.mean(cw_lengths)/fps, np.mean(ccw_lengths)/fps # divide by fps because the units of np.mean(...) is frames.
    
def compute_features_for_each_trace():
  for trace in traces:
    biases = compute_bias(trace)
    ccw_MISI, cw_MISI = compute_MISI(trace)
    features['bias'].append(biases)
    features['ccw_MISI'].append(ccw_MISI)
    features['cw_MISI'].append(cw_MISI)

if __name__ == '__main__':
  traces = []
  concentrations = ['MotMed', '1uM', '10uM', '100uM', '100nM', '1mM']
  for concentration in concentrations:
    video_name = '/' + concentration 
    for trace_path in paths[concentration]: # get all traces for this conc (flatten into groups of conc.)
      stream = np.load('traces/' + trace_path)
      for tr in stream:
        traces.append(tr)

    features = 	{
            'bias'								:	[],
            'ccw_MISI'			:	[],
            'cw_MISI'			:	[],
          }

    # TODO:
    '''
    TODO: Remove paths dict and other code above and revert to:
    video_name, videos_dir = get_video_path(sys.argv)
    print video_name, videos_dir
    traces = np.load('traces' + video_name + '_traces.npy')
    '''

    # Modifies features dict
    compute_features_for_each_trace()
    # add prefix?
    np.save('features' + video_name + '_features', features)

    # To open, run:
    # features = np.load('features' + video_name + '_features.npy')[()] # Note: [()] allows us to load the dict() we saved.

    # plot distribution of bias and mean_is_interval for each concentration?
