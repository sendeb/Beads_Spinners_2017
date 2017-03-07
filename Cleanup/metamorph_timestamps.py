# from __future__ import division, unicode_literals
# import numpy as np
# import scipy, xml.etree.ElementTree, tifffile
# from libtiff import TIFF
# from datetime import datetime
from utilities import *

def get(path_to_tif):
	# Parse timestamp strings from the XML 
	# of the Metamorph metadata.
	# Returns an equivalent numpy array (in miliseconds).
	scope_times = []
	scope_tif_file = tifffile.TiffFile(path_to_tif)
	for t in range(len(scope_tif_file)):
	    metadata = scope_tif_file[t].image_description
	    root = xml.etree.ElementTree.fromstring(metadata).find('PlaneInfo')
	    for neighbor in root:
	        if neighbor.attrib['type'] == 'time':
	            if neighbor.attrib['id'] == 'acquisition-time-local':
	                ts = neighbor.attrib['value']
	                scope_times += [ts]


	reference_time = datetime.utcfromtimestamp(0)

	def ms(ts):
	    # Convert a Metamorph timestamp string "ts" to miliseconds.
	    # Metamorph timestamp is zero padded to the RIGHT of the
	    # number of milliseconds.
	    # Example: 42.998, 43.66, 43.233, 43.555, 44.899, 45.12
	    #                   ^^^error (should be 43.066)
	    # This function gracefully handles that format.
	    if ts[-3:][0] == '.':
	        milisecs = int(ts[-2:])
	        new_ts = ts[:-3]
	    else:
	        milisecs = int(ts[-3:])
	        new_ts = ts[:-4]
	    dt = datetime.strptime(new_ts, '%Y%m%d %H:%M:%S')
	    return (dt - reference_time).total_seconds() * 1000.0 + milisecs

    # Convert timestamp strings to miliseconds.
	scope_ms = []
	for t in scope_times:
	    scope_ms += [ms(t)]
	scope_ms = np.array(scope_ms)

	# Make sure t=0 is first sample.
	scope_ms = scope_ms - scope_ms[0]

	return scope_ms

