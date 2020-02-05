"""
========
Barchart
========

A bar plot with errorbars and height labels on individual bars
"""
import sys
import pickle
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def load_results(fname):
    # Load individual brenda term performance
    perf = pickle.load(open(fname, "rb"))
    auprs = np.nanmean(perf, axis=1)
    errs = np.nanstd(perf, axis=1)

    auprs[np.isnan(auprs)] = 0
    errs[np.isnan(errs)] = 0

    return auprs, errs


def compute_diff(auprs_1, auprs_2):
    auprs_1 = np.array(auprs_1, dtype=np.float)
    auprs_2 = np.array(auprs_2, dtype=np.float)
    return auprs_2 - auprs_1


# Main()
def plot_diff(fn_1, fn_2, label_1, label_2, title_label):
    auprs_1, errs_1 = load_results(fn_1)
    auprs_2, errs_2 = load_results(fn_2)
    
    # auprs_3, errs_3 = load_results(fn_3)

    diff = compute_diff(auprs_1, auprs_2)
    brendas = [1, 2, 3, 4, 5, 6]
    brenda_names = ['Oxidoreductases', 'Transferases', 'Hydrolases', 'Lyases', 'Isomerases', 'Ligases']
    term_numbers = [344, 657, 844, 117, 93, 88]
    # children only
    idx = sorted(range(len(diff)), key=lambda x: diff[x])
    diff = [diff[i] for i in idx]
    better_brenda = 0
    for i in idx:
        if diff[i] >= 0:
            better_brenda += 1

    print ("### Number of better brenda terms: %d" % (better_brenda))
    print ("### Number of total brenda terms: %d" % (len(diff)))
    brendas = [brendas[i] for i in idx]
    print(brendas)
    brenda_names = [brenda_names[i] for i in idx]
    term_numbers = [term_numbers[i] for i in idx]
    print(term_numbers)
    auprs_1 = [auprs_1[i] for i in idx]
    auprs_2 = [auprs_2[i] for i in idx]
    # auprs_3 = [auprs_3[i] for i in idx]
    print('First three brendas names and go ids:')
    print(brendas[:3])

    errs_1 = [errs_1[i] for i in idx]
    errs_2 = [errs_2[i] for i in idx]
    # errs_3 = [errs_3[i] for i in idx]


    # # Plot # #
    N = len(brendas)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.5       # the width of the bars
    names = [label_1, label_2, 'AUPR diff']

    l = []
    fig, axes = plt.subplots(nrows=1, ncols=3)
    panel_1 = axes[0].barh(ind, auprs_1, width, align='center',
                           color='skyblue', edgecolor='white',
                           ecolor='black', yerr=errs_1)
    l.append(panel_1)
    axes[0].spines['right'].set_visible(False)
    axes[0].spines['bottom'].set_visible(False)
    axes[0].spines['top'].set_visible(False)
    # axes[0].set_xticks(fontsize=10)
    axes[0].set_xlim([0, 1.0])
    axes[0].set_ylim([-1, N + 0.5])
    axes[0].set_yticks(ind)
    axes[0].set_yticklabels(brenda_names, fontsize=11)
    axes[0].yaxis.set_ticks_position('left')
    axes[0].tick_params(axis='y')
    axes[0].xaxis.grid()
    for tick in axes[0].xaxis.get_majorticklabels():
        tick.set_horizontalalignment("right")

    panel_2 = axes[1].barh(ind, auprs_2, width, align='center',
                           color='purple', edgecolor='white',
                           ecolor='black', yerr=errs_2)
    l.append(panel_2)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['bottom'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    # axes[1].set_xticks([])
    axes[1].set_xlim([0, 1.0])
    axes[1].set_ylim([-1, N + 0.5])
    axes[1].set_yticks([])
    axes[1].xaxis.grid()

    panel_3 = axes[2].barh(ind, diff, width, align='center',
                           color='black', edgecolor='white')
    l.append(panel_3)
    axes[2].spines['right'].set_visible(False)
    axes[2].spines['bottom'].set_visible(False)
    axes[2].spines['top'].set_visible(False)
    # axes[2].set_xticks([])
    axes[2].set_yticks([])
    # axes[len(names)].set_xlim([min(diff) - 0.02, max(diff) + 0.02])
    axes[2].set_xlim([min(diff) - 0.02, max(diff) + 0.02])
    axes[2].set_ylim([-1, N + 0.5])
    axes[2].set_yticks([])
    axes[2].xaxis.grid()
    # axes[2].axvline(0, color='k')

    plt.suptitle(title_label, fontsize=16)
    plt.legend(l, names, loc='upper right', bbox_to_anchor=(0, 0, 1, 1.05), fontsize=12, ncol=5)
    plt.show()

# models
if __name__ == "__main__":
    fn_1 = str(sys.argv[1])
    fn_2 = str(sys.argv[2])
    label_1 = sys.argv[3]
    label_2 = sys.argv[4]
    title_label = sys.argv[5]

    plot_diff(fn_1, fn_2, label_1, label_2, title_label)
