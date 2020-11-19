import matplotlib
matplotlib.use('macosx')

import sys
import numpy as np
from scipy.stats import sem
import matplotlib.pyplot as plt
import pickle
plt.style.use('ggplot')


def load_macros(fnames):
    # fnames: filenames of performances for different alphas (0.0 to 1.0 in increments of 0.1)
    alphas = []
    macros = []
    std_errs = []
    for fname in fnames:
        print(fname)
        alpha = float(fname.split('_alpha_')[1][:3])
        print(alpha)
        alphas.append(alpha)
        perfs = pickle.load(open(fname, 'rb'))
        print(perfs.shape)
        trial_macros = np.nanmean(perfs, axis=0)
        macro = np.nanmean(trial_macros)
        macros.append(macro)
        std_errs.append(sem(trial_macros))
    return alphas, macros, std_errs


def plot_bars(x, y, stds, x_label, y_label, title):
    x_pos = np.arange(0, len(x))
    fig, ax = plt.subplots(1)
    #ax.bar(x_pos, y, yerr=stds, zorder=2)
    ax.bar(x_pos, y, yerr=stds)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x)
    ax.set_yticks(np.arange(0, 1.05*max(y), 0.02))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    #ax.grid(zorder=0)
    ax.set_title(title)
    plt.savefig(title + "_figure.eps")
    plt.show()


if __name__ == '__main__':
    fnames = sys.argv[2:]
    alphas, macros, stds = load_macros(fnames)
    title = sys.argv[1]
    plot_bars(alphas, macros, stds, 'Alpha', 'Macro AUPR', title)
