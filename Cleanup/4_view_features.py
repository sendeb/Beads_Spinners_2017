from utilities import *

'''
Generate grid of plots:
  Columns are Concentrations
  Rows are:   1) bias
              2) ccw MISI
              3) cw MISI
'''

features = ['bias', 'ccw_MISI', 'cw_MISI']
feature_xlabels = {'bias' : 'bias', 'ccw_MISI' : 'ccw MISI [s]', 'cw_MISI' : 'cw MISI [s]'}
feature_xlims = {'bias': (0, 0.67), 'ccw_MISI' : (0, 0.2), 'cw_MISI' : (0.01, 0.06) }
feature_ylims = {'bias': (0, 0.5), 'ccw_MISI' : (0, 0.667), 'cw_MISI' : (0, 0.85) }
C = len(concentrations)
fig, axes = plt.subplots(nrows=3, ncols=C)#sharex=True, sharey=True)

def plot_data(ax, data, concentration, xlabel, fontsize=12):
  # normalize histogram (by number of samples)
  n = len(data)
  weights = np.ones_like(data)/n
  ax.hist(data, weights=weights)
  # format (sub)plot
  ax.locator_params(nbins=3)

  xleft, xright = feature_xlims[xlabel]
  ax.set_xlim(left=xleft,right=xright)

  ybottom, ytop = feature_ylims[xlabel]
  ax.set_ylim(bottom=ybottom,top=ytop)

  xlabel = feature_xlabels[xlabel]
  ax.set_xlabel(xlabel, fontsize=fontsize)
  ax.set_ylabel('proportion of cells (n=' + str(n) + ')')
  ax.set_title(concentration, fontsize=fontsize)

for c_idx, concentration in enumerate(concentrations):
  print "On concentration:", concentration
  features = np.load('features/' + concentration + "/" + concentration + '_features.npy')[()] # Note: [()] allows us to load the dict() we saved.
  for f_idx, feature in enumerate(features):
    f = features[feature] # one feature at a time
    f = np.array(filter(lambda x: np.isnan(x) == False, f))
    plot_data(axes[f_idx][c_idx], f, concentration, xlabel=feature)
plt.tight_layout()
plt.show()
