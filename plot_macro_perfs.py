import matplotlib
matplotlib.use('TkAgg')
import sys
import numpy as np
from scipy.stats import sem
import matplotlib.pyplot as plt
import pickle
import scipy.io as sio

plt.style.use('ggplot')


def load_macros(fnames):
    # fnames: filenames of performances for different alphas (0.0 to 1.0 in increments of 0.1)
    macros = []
    std_errs = []
    for fname in fnames:
        print(fname)
        if fname[-5:] == '.pckl':
            perfs = pickle.load(open(fname, "rb"))
            trial_macros = np.nanmean(perfs, axis=0)
        elif fname[-4:] == '.mat':
            perfs = sio.loadmat(fname)
            perfs = perfs['all_aups']
            trial_macros = np.nanmean(perfs, axis=0)
        elif fname[-4:] == '.txt':
            trial_macros = []
            for line in open(fname, 'r'):
                fields = line.split(' ') 
                if len(fields) > 3 and is_number(fields[0]):
                    trial_macros.append(float(fields[1]))
                print(trial_macros)
            trial_macros = np.array(trial_macros)

        macro = np.nanmean(trial_macros)
        macros.append(macro)
        std_errs.append(sem(trial_macros))
    return macros, std_errs

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def plot_bars(x, y, stds, x_label, y_label, title):
    x_pos = np.arange(0, len(x))
    fig, ax = plt.subplots(1)
    #ax.bar(x_pos, y, yerr=stds, zorder=2)
    ax.bar(x_pos, y, yerr=stds, capsize=5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x)
    ax.set_yticks(np.arange(0, 1.05*max(y), 0.02))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    #ax.grid(zorder=0)
    ax.set_title(title)
    plt.savefig(title + '_macro_perfs.eps')
    plt.show()


if __name__ == '__main__':
    labels = sys.argv[1].split(',')
    fnames = sys.argv[3:]
    macros, stds = load_macros(fnames)
    plot_bars(labels, macros, stds, 'Branch', 'Macro AUPR', sys.argv[2])
