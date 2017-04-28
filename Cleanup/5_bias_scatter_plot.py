from utilities import *

# Input: some traces .npy (cols = time steps, rows = bact #)
#       2D arrays are a list of traces, each trace is one element in the list.
#
# Output: list of 2D arrays containing same number of rows
#     (one new time series for each bacterium (for each new engineered feature),
#     each array represents some feature, and
#     for each time step, the value of that feature for
#     the particular bacterium.



paths2 = {
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


fps = 80.0 # NOTE: CHANGE ME!!!! :NOTE
features = {'bias' : [], 'ccw_MISI' : [], 'cw_MISI' : []}
conc_map = {'100nM' : 100*1e-9, '1uM' : 1*1e-6, '10uM' : 10*1e-6, '100uM' : 100*1e-6, '1mM' : 1*1e-3, 'MotMed' : 0}

interval_bias = lambda s: np.sum((-np.array(s)+1)/2)/len(s) # CCW / (CCW + CW); s is interval over which to compute bias, s is signs of rotation direction. NOTE: correct if cw is positive, ccw is negative.

def compute_bias(trace, window=7, first=None):
  bias = []
  conv = np.convolve([-1.,1], trace, mode='full') # SHOULD WE UNWRAP?
  signs = np.sign(conv) #Positive values correspond to cw rotation. Negative = ccw rotation.
  interval = signs[:first]
  bias = interval_bias(interval)
  return bias

def compute_features_for_each_trace(first=800):
  for trace in traces:
    bias = compute_bias(trace, first=first) # Set first and frames so that it's about 10 s of data.
    features['bias'].append(bias)

if __name__ == '__main__':

  # The first strain or chemoattractant.
  avg_biases = {}
  for concentration in concentrations: # aggregating ALL traces for one concentration
    if concentration == 'MotMed':
      continue
    print "On concentration:", concentration
    traces = [] # Reset traces because save features from all traces of a single contration into one file. 
    features = {'bias' : [], 'ccw_MISI' : [], 'cw_MISI' : []} # Reset features
    for trace_path in paths[concentration]: # get all traces for this conc (flatten into groups of conc.)
      all_traces_for_one_stream = np.load('traces/' + trace_path.split('.')[0] + "_traces.npy")
      for trace in all_traces_for_one_stream:
        traces.append(trace)
    # Modifies features dict
    compute_features_for_each_trace(first=800) # compute features using ALL traces from ONE concentration
    avg_bias = np.mean(features['bias'])
    avg_biases[conc_map[concentration]] = avg_bias
  
  # For a second strain or chemoattractant.
  avg_biases2 = {}
  for concentration in concentrations: # aggregating ALL traces for one concentration
    if concentration == 'MotMed':
      continue
    print "On concentration:", concentration
    traces = [] # Reset traces because save features from all traces of a single contration into one file. 
    features = {'bias' : [], 'ccw_MISI' : [], 'cw_MISI' : []} # Reset features
    for trace_path in paths2[concentration]: # get all traces for this conc (flatten into groups of conc.)
      all_traces_for_one_stream = np.load('traces2/' + trace_path.split('.')[0] + ".npy")
      for trace in all_traces_for_one_stream:
        traces.append(trace)
    # Modifies features dict
    compute_features_for_each_trace(first=1700) # compute features using ALL traces from ONE concentration
    avg_bias = np.mean(features['bias'])
    avg_biases2[conc_map[concentration]] = avg_bias

    # To load:
    # Note: Use [()], which allows us to load the dict() we saved as a .npy file.

fig = plt.figure()
ax = plt.gca()
plt.title('Bias vs. Concentration')
plt.xlabel('Concentration of chemoattractant (serine)')
plt.ylabel('Avergage bias')
plt.scatter(avg_biases.keys(), avg_biases.values(), c='b', label='Non-adapting')
plt.scatter(avg_biases2.keys(), avg_biases2.values(), c='g', label='Wild type')
plt.legend(loc='upper right')
plt.xscale('log')
plt.xlim(1*1e-9, 1)
plt.show()
