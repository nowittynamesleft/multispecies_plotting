"""
========
Barchart
========

A bar plot with errorbars and height labels on individual bars
"""
import sys
import numpy as np
import matplotlib.pyplot as plt


def load_results(fname):
    # Load individual GO term performance
    goterms = []
    auprs = []
    errs = []
    f = open(fname, 'r')
    for line in f:
        splitted = line.strip().split()
        goterms.append(splitted[0])
        auprs.append(float(splitted[1]))
        errs.append(float(splitted[2]))
    f.close()

    return goterms, auprs, errs


def compute_diff(auprs_1, auprs_2):
    auprs_1 = np.array(auprs_1, dtype=np.float)
    auprs_2 = np.array(auprs_2, dtype=np.float)
    return auprs_2 - auprs_1


def goid2name(fname):
    goid2goname = {}
    goterms = []
    f = open(fname, 'r')
    for line in f:
        goid, goname = line.strip().split('\t')
        goid2goname[goid] = goname
        goterms.append(goid)
    f.close()

    return goid2goname, goterms


# Main()

# models
fn_1 = str(sys.argv[1])
fn_2 = str(sys.argv[2])


goterms, auprs_1, errs_1 = load_results(fn_1)
_, auprs_2, errs_2 = load_results(fn_2)
diff = compute_diff(auprs_1, auprs_2)


# children only
"""
diff = [diff[i] for i in range(0, 58)]
goterms = [goterms[i] for i in range(0, 58)]
auprs_1 = [auprs_1[i] for i in range(0, 58)]
auprs_2 = [auprs_2[i] for i in range(0, 58)]
errs_1 = [errs_1[i] for i in range(0, 58)]
errs_2 = [errs_2[i] for i in range(0, 58)]
"""

idx = sorted(range(len(diff)), key=lambda x: diff[x])
diff = [diff[i] for i in idx]
goterms = [goterms[i] for i in idx]
auprs_1 = [auprs_1[i] for i in idx]
auprs_2 = [auprs_2[i] for i in idx]
errs_1 = [errs_1[i] for i in idx]
errs_2 = [errs_2[i] for i in idx]


for goterm in goterms:
    print (goterm)

# # Plot # #
N = len(goterms)
ind = np.arange(N)  # the x locations for the groups
width = 0.5       # the width of the bars
names = ['CNN', 'CNN-multitask', 'Diff']

l = []
fig, axes = plt.subplots(nrows=1, ncols=3)
panel_1 = axes[0].barh(ind, auprs_1, width, align='center',
                       color='deepskyblue', edgecolor='white',
                       ecolor='black', yerr=errs_1)
l.append(panel_1)
axes[0].spines['right'].set_visible(False)
axes[0].spines['bottom'].set_visible(False)
axes[0].spines['top'].set_visible(False)
axes[0].set_xticks([])
axes[0].set_xlim([0, 1.0])
axes[0].set_ylim([-1, N + 0.5])
axes[0].set_yticks(ind)
axes[0].set_yticklabels(goterms, fontsize=12)
axes[0].yaxis.set_ticks_position('left')
axes[0].tick_params(axis='y')
for tick in axes[0].xaxis.get_majorticklabels():
    tick.set_horizontalalignment("right")

panel_2 = axes[1].barh(ind, auprs_2, width, align='center',
                       color='tomato', edgecolor='white',
                       ecolor='black', yerr=errs_2)
l.append(panel_2)
axes[1].spines['right'].set_visible(False)
axes[1].spines['bottom'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].set_xticks([])
axes[1].set_xlim([0, 1.0])
axes[1].set_ylim([-1, N + 0.5])
axes[1].set_yticks([])

panel_3 = axes[2].barh(ind, diff, width, align='center',
                       color='k', edgecolor='white', left=0)
l.append(panel_3)
axes[2].spines['right'].set_visible(False)
axes[2].spines['bottom'].set_visible(False)
axes[2].spines['top'].set_visible(False)
axes[2].set_xticks([])
axes[2].set_yticks([])
# axes[len(names)].set_xlim([min(diff) - 0.02, max(diff) + 0.02])
axes[2].set_xlim([-0.2, 0.2])
axes[2].set_ylim([-1, N+0.5])
axes[2].axvline(0, color='k')

plt.suptitle('Area under the Precision-Recall curve', fontsize=16)
plt.legend(l, names, bbox_to_anchor=(-0.5, -0.1), fontsize=16, ncol=5)
plt.show()
