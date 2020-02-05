"""
========
Barchart (final results)
========

A bar plot with errorbars and height labels on individual bars
"""
import numpy as np
import matplotlib
matplotlib.use('Tkagg')
import matplotlib.pyplot as plt
import pickle
import sys
from scipy.stats import sem
import matplotlib.patches as mpatches




METRIC = 'F1'
def load_results(fname, metric='pr_macro', ont='MF', Ntrials=1000):
    perf = pickle.load(open(fname, 'rb'))
    scores = perf[metric]
    return np.mean(scores), sem(scores)

org = 'Human'
ont = 'BRENDA Top Level Labels'


fnames = sys.argv[1:]
micro_means = []
macro_means = []
F1_means = []
micro_stds = []
macro_stds = []
F1_stds = []
lambdas = []
possible_colors = ['red', 'orange', 'purple', 'cyan', 'green']
colors = []
prev_lamb = -1
color_ind = 0
for fname in fnames:
    mean, std = load_results(fname, metric='pr_micro', ont=ont)    
    micro_means.append(mean)
    micro_stds.append(std)
    mean, std = load_results(fname, metric='pr_macro', ont=ont)    
    macro_means.append(mean)
    macro_stds.append(std)
    mean, std = load_results(fname, metric='F1', ont=ont)    
    F1_means.append(mean)
    F1_stds.append(std)
    if 'lamb' in fname:
        lamb = fname.split('lamb')[1].split('_')[1]
    else:
        lamb = 0
    lambdas.append(float(lamb))
    if prev_lamb != lamb and prev_lamb != -1:
        color_ind += 1
    colors.append(possible_colors[color_ind])
    prev_lamb = lamb

sort_inds = np.argsort(lambdas)
lambdas = np.array(lambdas)[sort_inds]

micro_means = np.array(micro_means)[sort_inds]
macro_means = np.array(macro_means)[sort_inds]
F1_means = np.array(F1_means)[sort_inds]

micro_stds = np.array(micro_stds)[sort_inds]
macro_stds = np.array(macro_stds)[sort_inds]
F1_stds = np.array(F1_stds)[sort_inds]

colors = np.array(colors)[sort_inds]

print(macro_stds)
x_coords = np.arange(0, len(fnames))
fig, ax = plt.subplots()
width = 0.2
ax.bar(x_coords*width, macro_means, width=width, label='Macro AUPR', yerr=macro_stds, color=colors)
ax.bar(x_coords*width + 5*width, micro_means, width=width, label='Micro AUPR', yerr=micro_stds, color=colors)
ax.bar(x_coords*width + 10*width, F1_means, width=width, label='F1 Score', yerr=F1_stds, color=colors)

patch_list = []
for i in range(0, len(lambdas)):
    patch_list.append(mpatches.Patch(color=possible_colors[i], label='Similarity Reg. Strength: ' + str(lambdas[i])))
plt.legend(handles=patch_list)
#ax.set_ylabel(METRIC, fontsize=18)
ax.set_xlabel('', fontsize=18)
ax.set_xticks(x_coords*5*width + 1.5*width)
ax.set_xticklabels(['Macro AUPR', 'Micro AUPR', 'F1 Score'])
#ax.set_xticklabels(lambdas)
#ax.set_xticks(np.arange(1, len(fnames) + 1))
#ax.set_xticklabels(fnames)
plt.locator_params(nbins=20)
ax.set_title(org + ' ' + ont + ' Cross Validation Performance', fontsize=18)
ax.grid()
ax.set_axisbelow(True)
plt.show()
