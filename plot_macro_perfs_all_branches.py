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


def plot_bars(x, y, stds, x_label, y_label, ax):
    x_pos = np.arange(0, len(x))
    #ax.bar(x_pos, y, yerr=stds, zorder=2)
    ax.bar(x_pos, y, width=1/len(x), yerr=stds, capsize=5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x)
    ax.set_yticks(np.arange(0, 1.05*max(y), 0.05))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


if __name__ == '__main__':
    labels = sys.argv[1].split(',')
    title = sys.argv[2]
    fnames = sys.argv[3:]
    branch_fnames = {'MF': [], 'BP': [], 'CC': []}
    for fname in fnames:
        if 'molecular_function' in fname:
            branch_fnames['MF'].append(fname)
        elif 'cellular_component' in fname:
            branch_fnames['CC'].append(fname)
        elif 'biological_process' in fname:
            branch_fnames['BP'].append(fname)
        else:
            print('One of the filenames doesn\'t have the strings \'molecular_function\', \'cellular_component\', or \'biological_process\' in it. All performance filenames must have one of these to be included in the plots.')
    assert len(branch_fnames['MF']) == len(branch_fnames['CC']) == len(branch_fnames['BP'])
    # I want to input number of files and have it know to put all files in one plots, with the number of subplots being with the number of files
    fig, axes = plt.subplots(1, 3, constrained_layout=True)
    for i, branch in enumerate(['MF', 'BP', 'CC']):
        macros, stds = load_macros(branch_fnames[branch])
        plot_bars(labels, macros, stds, branch + ' Methods', 'Macro AUPR', axes[i])
    fig.suptitle(title, fontsize=16) 
    #plt.tight_layout()
    plt.show()
    plt.savefig(title + '_macro_perfs.eps')
