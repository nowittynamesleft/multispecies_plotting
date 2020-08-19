import matplotlib
#matplotlib.use('TkAgg')
import sys
import numpy as np
from scipy.stats import sem
import matplotlib.pyplot as plt
import pickle
import scipy.io as sio
from sklearn.metrics import accuracy_score, f1_score

plt.style.use('ggplot')


'''
CURRENT USAGE:
python3 plot_all_perfs_all_branches.py [METHOD LIST (COMMA SEPARATED)] [TITLE] [setting: 'loso', 'alpha_testing', or 'noloso'] [LABEL FILENAME (EXPERIMENTAL EVIDENCE CODES)] [PRED FILES FOR MF,BP,CC FOR ALL METHODS TO PLOT]
'''


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


def load_perfs(fnames, branch_label_dict, loso=False):
    # fnames: filenames of performances for different alphas (0.0 to 1.0 in increments of 0.1)
    macros = []
    macro_std_errs = []
    micros = []
    micro_std_errs = []
    accs = []
    acc_std_errs = []
    f1s = []
    f1_std_errs = []
    test_prot_set = set()
    for fname in fnames:
        print(fname)
        '''
        if fname[-5:] == '.pckl':
            perfs = pickle.load(open(fname, "rb"))
            trial_macros = np.nanmean(perfs, axis=0)
        '''
        if loso: # leave one species out; only one trial
            if '_scores.pckl' in fname: # now assuming BLAST pred files
                pred_file = pickle.load(open(fname, 'rb'))
                num_trials = len(pred_file['Y_hat_test'])
                print('BLAST baseline')
                trial_macros = []
                trial_micros = []
                trial_accs = []
                trial_f1s = []

                new_test_prot_set = set(pred_file['prots'])

                print('Number of go ids: ' + str(len(pred_file['go_IDs'])))
                curr_trial_preds = pred_file['Y_hat_test']
                print('Num test: ' + str(curr_trial_preds.shape[0]))

                if len(test_prot_set) > 0:
                    try:
                        assert new_test_prot_set == test_prot_set
                    except AssertionError:
                        print('Intersect:')
                        print(len(new_test_prot_set & test_prot_set))
                        print('Difference (blast - other):')
                        print(new_test_prot_set - test_prot_set)
                        print('Difference (other - blast):')
                        print(test_prot_set - new_test_prot_set)
                        print('Removing the (blast - other) samples')
                        prots_to_remove = new_test_prot_set - test_prot_set
                        inds_to_remove = [list(pred_file['prots']).index(prot) for prot in prots_to_remove]
                        keep_inds = np.array([i for i in range(0, len(pred_file['Y_hat_test'])) if i not in inds_to_remove])
                        print('Before')
                        curr_trial_preds = curr_trial_preds[keep_inds,:]
                        pred_file['Y_hat_test'] = curr_trial_preds
                        pred_file['prots'] = pred_file['prots'][keep_inds]
                        print('After')
                        print(curr_trial_preds.shape)

                else:
                    test_prot_set = new_test_prot_set
                '''
                curr_macro, curr_micro, curr_acc, curr_f1 = evaluate_performance(curr_trial_labels, curr_trial_preds, curr_trial_preds > 0.5)
                
                trial_macros.append(curr_macro)
                trial_micros.append(curr_micro)
                trial_accs.append(curr_acc)
                trial_f1s.append(curr_f1)
                '''
                trial_macros, trial_micros, trial_accs, trial_f1s = align_pred_file_with_label_file_loso(pred_file, branch_label_dict, blast=True)

            elif fname[-5:] == '.pckl': # now assuming pred file instead
                pred_file = pickle.load(open(fname, "rb"))
                trial_macros = []
                trial_micros = []
                trial_accs = []
                trial_f1s = []
                print('Number of go ids: ' + str(len(pred_file['GO_IDs'])))
                curr_trial_preds = pred_file['preds']
                
                print('Num test: ' + str(curr_trial_preds.shape[0]))
                '''
                trial_macros.append(curr_macro)
                trial_micros.append(curr_micro)
                trial_accs.append(curr_acc)
                trial_f1s.append(curr_f1)
                '''
                new_test_prot_set = set(pred_file['prot_IDs'])
                if len(test_prot_set) > 0:
                    try:
                        assert new_test_prot_set == test_prot_set
                    except AssertionError:
                        print('Intersect:')
                        print(len(new_test_prot_set & test_prot_set))
                        print('Difference (maxout - blast):')
                        print(new_test_prot_set - test_prot_set)
                        print('Difference (blast - maxout):')
                        print(test_prot_set - new_test_prot_set)
                else:
                    test_prot_set = new_test_prot_set
                trial_macros, trial_micros, trial_accs, trial_f1s = align_pred_file_with_label_file_loso(pred_file, branch_label_dict, blast=False)
        else:
            '''
            if '_scores.pckl' in fname: # now assuming BLAST pred files
                pred_file = pickle.load(open(fname, 'rb'))
                num_trials = len(pred_file['Y_hat_test_list'])
                print('BLAST baseline')
                trial_macros = []
                trial_micros = []
                trial_accs = []
                trial_f1s = []
                for trial in range(0, num_trials):
                    print('Number of go ids: ' + str(len(pred_file['go_IDs'])))
                    curr_trial_preds = pred_file['Y_hat_test_list'][trial]
                    curr_trial_labels = pred_file['Y_test_list'][trial]


                    print('Num test: ' + str(curr_trial_preds.shape[0]))
                    curr_macro, curr_micro, curr_acc, curr_f1 = evaluate_performance(curr_trial_labels, curr_trial_preds, curr_trial_preds > 0.5)
                    trial_macros.append(curr_macro)
                    trial_micros.append(curr_micro)
                    trial_accs.append(curr_acc)
                    trial_f1s.append(curr_f1)
            '''
            if '_scores.pckl' in fname: # now assuming BLAST pred files
                pred_file = pickle.load(open(fname, "rb"))
                trial_macros, trial_micros, trial_accs, trial_f1s = align_pred_file_with_label_file(pred_file, branch_label_dict, blast=True)

            elif fname[-5:] == '.pckl': # now assuming pred file instead
                pred_file = pickle.load(open(fname, "rb"))
                trial_macros, trial_micros, trial_accs, trial_f1s = align_pred_file_with_label_file(pred_file, branch_label_dict, blast=False)
                    
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


def plot_bars_grouped_by_metric(method_names, metric_lists, metric_stds, x_label, ax, alpha_testing=False):
    metric_names = ['M-AUPR', 'm-AUPR', 'Acc', 'F1 score']
    #color_hexes = ['#003f5c', '#58508d', '#bc5090', '#ff6361', '#ffa600']
    #color_hexes = ['#130fff', '#ff0073', '#ff8600', '#a6ff00']
    #color_hexes = ['#5d63a6', '#c97199', '#f1a07d', '#dfe092']
    if not alpha_testing:
        color_hexes = ['#003f5c', '#58508d','#bc5090','#ff6361','#ffa600', '#dfe092', '#a6ff00', '#98bac2', '#b6d1d5', '#d4e8e9', '#f3ffff']
        capsize = 2
    else:
        #color_hexes = ['#00141d', '#072c39', '#193d46', '#004c6d', '#29617d', '#45778d', '#608d9e', '#7ca3b0', '#98bac2', '#b6d1d5', '#d4e8e9', '#f3ffff']
        color_hexes = ['#00141d', '#072c39', '#193d46', '#004c6d', '#45778d', '#7ca3b0', '#b6d1d5', '#d4e8e9', '#f3ffff']
        capsize = 2

    print('Number of colors:')
    print(len(color_hexes))
    print('Number of methods:')
    print(len(method_names))
    # switch the metrics and the method names
    print('Metric lists shape')
    print(np.array(metric_lists).shape)
    print('Method names shape')
    print(np.array(method_names).shape)
    for i in range(0, len(method_names)):
        x_pos = i + np.arange(0, len(metric_lists))*(len(method_names)+1)
        print(len(method_names))
        print(np.array(metric_lists).shape)
        methods_curr_metric = np.array(metric_lists)[:,i].tolist()
        methods_curr_metric_stds = np.array(metric_stds)[:,i].tolist()
        ax.bar(x_pos, methods_curr_metric, width=1, yerr=methods_curr_metric_stds, capsize=capsize, label=method_names[i], color=color_hexes[i])
    ax.set_xticks(np.arange(0, len(metric_lists))*(len(method_names)+1) + len(method_names)/2)
    ax.set_xticklabels(metric_names)
    maximum_perf = np.max(np.array(metric_lists))
    ax.set_yticks(np.arange(0, maximum_perf + 0.05, 0.05)) # take the max of all performances as the height
    ax.set_yticks(np.arange(0, 1.05, 0.05))
    ax.set_ylim([0, maximum_perf + 0.1])
    ax.set_xlabel(x_label)
    #ax.set_ylabel(y_label)
    handles, labels = ax.get_legend_handles_labels()
    return handles, labels


def plot_bars_all_metrics(method_names, metric_lists, metric_stds, x_label, ax):
    metric_names = ['M-AUPR', 'm-AUPR', 'Acc', 'F1 score']
    for i in range(0, len(metric_lists)):
        x_pos = i + np.arange(0, len(method_names))*(len(metric_lists)+1)
        ax.bar(x_pos, metric_lists[i], width=1, yerr=metric_stds[i], capsize=2, label=metric_names[i])
    ax.set_xticks(np.arange(0, len(method_names))*(len(metric_lists)+1) + len(metric_lists)/2)
    ax.set_xticklabels(method_names)
    ax.set_yticks(np.arange(0, 1.05, 0.05))
    ax.set_xlabel(x_label)
    #ax.set_ylabel(y_label)
    ax.legend()


def load_label_mats(label_fname):
    label_pickle = pickle.load(open(label_fname,'rb'))
    to_short_name = {'molecular_function': 'MF', 'cellular_component': 'CC', 'biological_process': 'BP'}
    label_dict = {}
    for branch in to_short_name.keys():
        annot = label_pickle['annot'][branch]
        label_dict[to_short_name[branch]] = {}
        label_dict[to_short_name[branch]]['annot'] = np.array(annot.todense())
        label_dict[to_short_name[branch]]['prot_IDs'] = label_pickle['prot_IDs']
        label_dict[to_short_name[branch]]['go_IDs'] = label_pickle['go_IDs'][branch]
    return label_dict


def get_common_indices(pred_ids, string_annot_ids):
    common_ids = sorted(list(set(string_annot_ids).intersection(pred_ids)))
    #print ("### Number of elements in intersection:", len(common_ids))
    pred_ids_idx = [pred_ids.index(prot) for prot in common_ids] # pred_ids_idx is the array of indices in the annotation protein list of each protein common to both prediction and string annot protein lists
    string_annot_ids_idx = [string_annot_ids.index(prot) for prot in common_ids] # same thing for string protein list
    return pred_ids_idx, string_annot_ids_idx 


def align_pred_file_with_label_file(pred_dict, label_dict, blast=False):
    trial_macros = []
    trial_micros = []
    trial_accs = []
    trial_f1s = []
    num_trials = len(pred_dict['trial_splits'])
    label_mat = label_dict['annot']
    label_prot_ids = label_dict['prot_IDs']
    for trial in range(0, num_trials):
        print('Trial ' + str(trial) + ':')
        print('Number of go ids before intersecting: ' + str(len(pred_dict['GO_IDs'])))
        curr_trial_test_inds = pred_dict['trial_splits'][trial][1]
        curr_trial_train_inds = pred_dict['trial_splits'][trial][0]
        if blast:
            curr_trial_preds = pred_dict['trial_preds'][trial]
        else:
            curr_trial_preds = pred_dict['trial_preds'][trial][curr_trial_test_inds]
        #curr_trial_labels = pred_dict['true_labels'][curr_trial_test_inds]
        curr_trial_test_prot_ids = list(np.array(pred_dict['prot_IDs'])[curr_trial_test_inds])
        print('Num test before removing IEA/other evidence codes: ' + str(len(curr_trial_test_inds)))
        pred_prots_idx, string_annot_prots_idx = get_common_indices(curr_trial_test_prot_ids, label_prot_ids)
        pred_go_terms = list(pred_dict['GO_IDs'])
        label_go_terms = list(label_dict['go_IDs'])
        pred_go_idx, string_annot_go_idx = get_common_indices(pred_go_terms, label_go_terms)
        print('Number test prot ids with experimental evidence codes: ' + str(len(pred_prots_idx)))
        print('Number intersection GO ids with experimental evidence codes: ' + str(len(pred_go_idx)))
        try:
            assert len(pred_go_idx) == len(pred_dict['GO_IDs'])
        except AssertionError:
            print('NOT ALL GO TERMS WERE IN EXPERIMENTAL ANNOTATION SET')
        #print('string annot prots idx' + str(len(string_annot_prots_idx)))
        curr_trial_labels = label_mat[string_annot_prots_idx, :]
        curr_trial_labels = curr_trial_labels[:, string_annot_go_idx]
        curr_trial_preds = curr_trial_preds[pred_prots_idx, :]
        curr_trial_preds = curr_trial_preds[:, pred_go_idx]
        #print('Label shape and pred shape:')
        print(curr_trial_labels.shape)
        print(curr_trial_preds.shape)

        curr_macro, curr_micro, curr_acc, curr_f1 = evaluate_performance(curr_trial_labels, curr_trial_preds, curr_trial_preds > 0.5)
        print('Num train: ' + str(len(curr_trial_train_inds)))
        print('Num test: ' + str(len(curr_trial_test_inds)))
        trial_macros.append(curr_macro)
        trial_micros.append(curr_micro)
        trial_accs.append(curr_acc)
        trial_f1s.append(curr_f1)
    return trial_macros, trial_micros, trial_accs, trial_f1s


def align_pred_file_with_label_file_loso(pred_dict, label_dict, blast=False):
    trial_macros = []
    trial_micros = []
    trial_accs = []
    trial_f1s = []
    label_mat = label_dict['annot']
    label_prot_ids = label_dict['prot_IDs']
    print(pred_dict.keys())
    if blast:
        pred_go_ids = pred_dict['go_IDs']
        preds = pred_dict['Y_hat_test']
        pred_prot_ids = list(pred_dict['prots'])
    else:
        pred_go_ids = pred_dict['GO_IDs']
        preds = pred_dict['preds']
        pred_prot_ids = list(pred_dict['prot_IDs'])
    print('Number of go ids before intersecting: ' + str(len(pred_go_ids)))
    #curr_trial_labels = pred_dict['true_labels'][curr_trial_test_inds]
    print('Num test before removing IEA/other evidence codes: ' + str(preds.shape[0]))
    assert len(pred_prot_ids) == preds.shape[0]
    pred_prots_idx, string_annot_prots_idx = get_common_indices(pred_prot_ids, label_prot_ids)
    pred_go_terms = list(pred_go_ids)
    label_go_terms = list(label_dict['go_IDs'])
    pred_go_idx, string_annot_go_idx = get_common_indices(pred_go_terms, label_go_terms)
    print('Number test prot ids with experimental evidence codes: ' + str(len(pred_prots_idx)))
    print('Number intersection GO ids with experimental evidence codes: ' + str(len(pred_go_idx)))
    try:
        assert len(pred_go_idx) == len(pred_go_ids)
    except AssertionError:
        print('NOT ALL GO TERMS WERE IN EXPERIMENTAL ANNOTATION SET')
    #print('string annot prots idx' + str(len(string_annot_prots_idx)))
    curr_trial_labels = label_mat[string_annot_prots_idx, :]
    curr_trial_labels = curr_trial_labels[:, string_annot_go_idx]
    preds = preds[pred_prots_idx, :]
    preds = preds[:, pred_go_idx]
    #print('Label shape and pred shape:')
    print(curr_trial_labels.shape)
    print(preds.shape)

    curr_macro, curr_micro, curr_acc, curr_f1 = evaluate_performance(curr_trial_labels, preds, preds > 0.5)
    print('Num test: ' + str(preds.shape[0]))
    trial_macros.append(curr_macro)
    trial_micros.append(curr_micro)
    trial_accs.append(curr_acc)
    trial_f1s.append(curr_f1)
    return trial_macros, trial_micros, trial_accs, trial_f1s


if __name__ == '__main__':
    labels = sys.argv[1].split(',')
    title = sys.argv[2]
    leave_one_species_out = False
    alpha_testing = False
    if sys.argv[3] == 'loso':
        leave_one_species_out = True
    elif sys.argv[3] == 'alpha_testing':
        alpha_testing = True
    print('leave one species out?' + str(leave_one_species_out))
    print('Alpha testing?' + str(alpha_testing))
    label_fname = sys.argv[4] # added to get only evidence coded proteins extracted
    print('label fname ' + label_fname)
    label_dict = load_label_mats(label_fname)
    all_branches = sys.argv[5] == 'all'
    fnames = sys.argv[6:]
    if all_branches:
        print(str(len(labels)*3) + ' ' + str(len(fnames)))
        assert len(labels)*3 == len(fnames)
    else:
        assert len(labels) == len(fnames)

    
    branch_fnames = {'MF': [], 'BP': [], 'CC': []}
    ex = False
    for fname in fnames:
        if 'molecular_function' in fname:
            branch_fnames['MF'].append(fname)
        elif 'cellular_component' in fname:
            branch_fnames['CC'].append(fname)
        elif 'biological_process' in fname:
            branch_fnames['BP'].append(fname)
        else:
            print('One of the filenames doesn\'t have the strings \'molecular_function\', \'cellular_component\', or \'biological_process\' in it. All performance filenames must have one of these to be included in the plots. Exiting after checking the rest.')
            print(fname)
            ex = True
    if ex:
        exit()
    if all_branches:
        try:
            assert len(branch_fnames['MF']) == len(branch_fnames['CC']) == len(branch_fnames['BP'])
            for branch in ['MF', 'BP','CC']:
                print(branch)
                print(len(branch_fnames[branch]))
        except AssertionError:
            print('Not all branches have same number of files associated')
            for branch in ['MF', 'BP','CC']:
                print(branch)
                print(branch_fnames[branch])
                print(len(branch_fnames[branch]))
            exit()
        # I want to input number of files and have it know to put all files in one plot, with the number of subplots being with the number of files
        if leave_one_species_out:
            fig, axes = plt.subplots(1, 3, constrained_layout=True, figsize=(8.67, 4.61))
        else:
            fig, axes = plt.subplots(1, 3, constrained_layout=True, figsize=(8.71,5.66))
        #fig_2, axes_2 = plt.subplots(1, 3, constrained_layout=True)
        for i, branch in enumerate(['MF', 'BP', 'CC']):
            macros, macro_stds, micros, micro_stds, accs, acc_stds, f1s, f1_stds = load_perfs(branch_fnames[branch], label_dict[branch], loso=leave_one_species_out)
            metric_lists = [macros, micros, accs, f1s]
            metric_stds = [macro_stds, micro_stds, acc_stds, f1_stds]
            print(np.array(metric_lists).shape)
            print(np.array(metric_stds).shape)
            handles, method_lab = plot_bars_grouped_by_metric(labels, metric_lists, metric_stds, branch, axes[i], alpha_testing=alpha_testing)
            #plot_bars_all_metrics(labels, metric_lists, metric_stds, branch, axes_2[i])
            '''
            plot_bars(labels, macros, macro_stds, branch, 'Macro AUPR', axes[i])
            plot_bars(labels, micros, micro_stds, branch, 'Micro AUPR', axes[i])
            plot_bars(labels, accs, acc_stds, branch, 'ACC', axes[i])
            plot_bars(labels, f1s, f1_stds, branch, 'F1', axes[i])
            '''
        if leave_one_species_out:
            plt.subplots_adjust(top=0.925, bottom=0.225, left=0.05, right=0.99, hspace=0.2, wspace=0.2)
        else:
            plt.subplots_adjust(top=0.945, bottom=0.34, left=0.05, right=0.99, hspace=0.2, wspace=0.2)
    else:
        #fig_2, axes_2 = plt.subplots(1, 3, constrained_layout=True)
        fig, ax = plt.subplots(1, 1, constrained_layout=True)
        for i, branch in enumerate(['MF', 'BP', 'CC']):
            if len(branch_fnames[branch]) > 0:
                print('Found file with branch ' + branch +', assuming no other branches are present.')
                macros, macro_stds, micros, micro_stds, accs, acc_stds, f1s, f1_stds = load_perfs(branch_fnames[branch], label_dict[branch], loso=leave_one_species_out)
                metric_lists = [macros, micros, accs, f1s]
                metric_stds = [macro_stds, micro_stds, acc_stds, f1_stds]
                print(np.array(metric_lists).shape)
                print(np.array(metric_stds).shape)
                handles, method_lab = plot_bars_grouped_by_metric(labels, metric_lists, metric_stds, branch, ax, alpha_testing=alpha_testing)
                break

    fig.legend(handles, labels, loc='lower center')
    fig.suptitle(title, fontsize=16) 
    #plt.tight_layout()
    plt.show()
    file_title = ''.join(title.split(' '))
    fig.savefig(file_title + '_macro_perfs.eps', format='eps')
