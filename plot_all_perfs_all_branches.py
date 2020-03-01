import matplotlib
matplotlib.use('TkAgg')
import sys
import numpy as np
from scipy.stats import sem
import matplotlib.pyplot as plt
import pickle
import scipy.io as sio
from sklearn.metrics import accuracy_score, f1_score

plt.style.use('ggplot')


def real_AUPR(label, score):
    """Computing real AUPR . By Vlad and Meet"""
    label = label.flatten()
    score = score.flatten()

    order = np.argsort(score)[::-1]
    label = label[order]

    P = np.count_nonzero(label)
    # N = len(label) - P

    TP = np.cumsum(label, dtype=float)
    PP = np.arange(1, len(label)+1, dtype=float)  # python

    x = np.divide(TP, P)  # recall
    y = np.divide(TP, PP)  # precision

    pr = np.trapz(y, x)
    f = np.divide(2*x*y, (x + y))
    idx = np.where((x + y) != 0)[0]
    if len(idx) != 0:
        f = np.max(f[idx])
    else:
        f = 0.0

    return pr, f


def evaluate_performance(y_test, y_score, y_pred):
    """Evaluate performance"""
    n_classes = y_test.shape[1]

    # Compute macro-averaged AUPR
    pr_macro = 0.0
    n = 0
    for i in range(n_classes):
        pr, _ = real_AUPR(y_test[:, i], y_score[:, i])
        if sum(y_test[:, i]) > 0:
            n += 1
            pr_macro += pr
    pr_macro /= n

    # Compute micro-averaged AUPR
    pr_micro, _ = real_AUPR(y_test, y_score)

    # Computes accuracy
    acc = accuracy_score(y_test, y_pred)

    # Computes F1-score
    alpha = 3
    y_new_pred = np.zeros_like(y_pred)
    for i in range(y_pred.shape[0]):
        top_alpha = np.argsort(y_score[i, :])[-alpha:]
        y_new_pred[i, top_alpha] = np.array(alpha*[1])
    F1 = f1_score(y_test, y_new_pred, average='micro')

    return pr_macro, pr_micro, acc, F1


def load_perfs(fnames):
    # fnames: filenames of performances for different alphas (0.0 to 1.0 in increments of 0.1)
    macros = []
    macro_std_errs = []
    micros = []
    micro_std_errs = []
    accs = []
    acc_std_errs = []
    f1s = []
    f1_std_errs = []
    for fname in fnames:
        print(fname)
        '''
        if fname[-5:] == '.pckl':
            perfs = pickle.load(open(fname, "rb"))
            trial_macros = np.nanmean(perfs, axis=0)
        '''
        if fname[-5:] == '.pckl': # now assuming pred file instead
            pred_file = pickle.load(open(fname, "rb"))
            num_trials = len(pred_file['trial_splits'])
            trial_macros = []
            trial_micros = []
            trial_accs = []
            trial_f1s = []
            for trial in range(0, num_trials):
                print('Number of go ids: ' + str(len(pred_file['GO_IDs'])))
                curr_trial_test_inds = pred_file['trial_splits'][trial][1]
                curr_trial_preds = pred_file['trial_preds'][trial][curr_trial_test_inds]
                curr_trial_labels = pred_file['true_labels'][curr_trial_test_inds]
                curr_macro, curr_micro, curr_acc, curr_f1 = evaluate_performance(curr_trial_labels, curr_trial_preds, curr_trial_preds > 0.5)
                trial_macros.append(curr_macro)
                trial_micros.append(curr_micro)
                trial_accs.append(curr_acc)
                trial_f1s.append(curr_f1)
                
        elif fname[-4:] == '.mat': # worry about this one later, this is for GeneMANIA
            perfs = sio.loadmat(fname)
            perfs = perfs['all_aups']
            trial_macros = np.nanmean(perfs, axis=0)
        elif fname[-4:] == '.txt':
            trial_macros = []
            trial_micros = []
            trial_accs = []
            trial_f1s = []
            for line in open(fname, 'r'):
                fields = line.split(' ') 
                if len(fields) > 3 and is_number(fields[0]):
                    trial_micros.append(float(fields[0]))
                    trial_macros.append(float(fields[1]))
                    trial_accs.append(float(fields[2]))
                    trial_f1s.append(float(fields[3]))
                print(trial_macros)
            trial_macros = np.array(trial_macros)
            trial_micros = np.array(trial_micros)
            trial_accs = np.array(trial_accs)
            trial_f1s = np.array(trial_f1s)

        macro = np.nanmean(trial_macros)
        macros.append(macro)
        macro_std_errs.append(sem(trial_macros))

        micro = np.nanmean(trial_micros)
        micros.append(micro)
        micro_std_errs.append(sem(trial_micros))

        acc = np.nanmean(trial_accs)
        accs.append(acc)
        acc_std_errs.append(sem(trial_accs))

        f1 = np.nanmean(trial_f1s)
        f1s.append(f1)
        f1_std_errs.append(sem(trial_f1s))
    assert len(macros) == len(macro_std_errs)
    assert len(micros) == len(micro_std_errs)
    assert len(accs) == len(acc_std_errs)
    assert len(f1s) == len(f1_std_errs)
    return macros, macro_std_errs, micros, micro_std_errs, accs, acc_std_errs, f1s, f1_std_errs

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def plot_bars(x, y, stds, x_label, y_label, ax, start_pos):
    x_pos = start_pos + np.arange(0, len(x))
    #ax.bar(x_pos, y, yerr=stds, zorder=2)
    ax.bar(x_pos, y, width=1/len(x), yerr=stds, capsize=5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x)
    ax.set_yticks(np.arange(0, 1, 0.05))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


def plot_bars_grouped_by_metric(method_names, metric_lists, metric_stds, x_label, ax):
    metric_names = ['Macro AUPR', 'Micro AUPR', 'Acc', 'F1 score']
    #color_hexes = ['#003f5c', '#58508d', '#bc5090', '#ff6361', '#ffa600']
    #color_hexes = ['#130fff', '#ff0073', '#ff8600', '#a6ff00']
    color_hexes = ['#5d63a6', '#c97199', '#f1a07d', '#dfe092']
    # switch the metrics and the method names
    for i in range(0, len(method_names)):
        x_pos = i + np.arange(0, len(metric_lists))*(len(method_names)+1)
        methods_curr_metric = np.array(metric_lists)[:,i].tolist()
        methods_curr_metric_stds = np.array(metric_stds)[:,i].tolist()
        ax.bar(x_pos, methods_curr_metric, width=1, yerr=methods_curr_metric_stds, capsize=5, label=method_names[i], color=color_hexes[i])
    ax.set_xticks(np.arange(0, len(metric_lists))*(len(method_names)+1) + len(method_names)/2)
    ax.set_xticklabels(metric_names)
    ax.set_yticks(np.arange(0, 1.05, 0.05))
    ax.set_ylim([0, 1.1])
    ax.set_xlabel(x_label)
    #ax.set_ylabel(y_label)
    ax.legend()


def plot_bars_all_metrics(method_names, metric_lists, metric_stds, x_label, ax):
    metric_names = ['Macro AUPR', 'Micro AUPR', 'Acc', 'F1 score']
    for i in range(0, len(metric_lists)):
        x_pos = i + np.arange(0, len(method_names))*(len(metric_lists)+1)
        ax.bar(x_pos, metric_lists[i], width=1, yerr=metric_stds[i], capsize=5, label=metric_names[i])
    ax.set_xticks(np.arange(0, len(method_names))*(len(metric_lists)+1) + len(metric_lists)/2)
    ax.set_xticklabels(method_names)
    ax.set_yticks(np.arange(0, 1.05, 0.05))
    ax.set_xlabel(x_label)
    #ax.set_ylabel(y_label)
    ax.legend()


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
    # I want to input number of files and have it know to put all files in one plot, with the number of subplots being with the number of files
    fig, axes = plt.subplots(1, 3, constrained_layout=True)
    #fig_2, axes_2 = plt.subplots(1, 3, constrained_layout=True)

    for i, branch in enumerate(['MF', 'BP', 'CC']):
        macros, macro_stds, micros, micro_stds, accs, acc_stds, f1s, f1_stds = load_perfs(branch_fnames[branch])
        metric_lists = [macros, micros, accs, f1s]
        metric_stds = [macro_stds, micro_stds, acc_stds, f1_stds]
        print(metric_lists)
        print(metric_stds)
        plot_bars_grouped_by_metric(labels, metric_lists, metric_stds, branch, axes[i])
        #plot_bars_all_metrics(labels, metric_lists, metric_stds, branch, axes_2[i])
        '''
        plot_bars(labels, macros, macro_stds, branch, 'Macro AUPR', axes[i])
        plot_bars(labels, micros, micro_stds, branch, 'Micro AUPR', axes[i])
        plot_bars(labels, accs, acc_stds, branch, 'ACC', axes[i])
        plot_bars(labels, f1s, f1_stds, branch, 'F1', axes[i])
        '''
    fig.suptitle(title, fontsize=16) 
    #plt.tight_layout()
    plt.show()
    plt.savefig(title + '_macro_perfs.eps')