import numpy as np
import matplotlib.pyplot as plt
from plot_all_perfs_all_branches import real_AUPR, get_common_indices
from scipy.stats import sem
import sys
import pickle
#plt.style.use('ggplot')
import argparse


def get_go_term_perfs(labels, preds):
    go_term_perfs = []
    for go_ind in range(0, labels.shape[1]):
        pr, _ = real_AUPR(labels[:, go_ind], preds[:, go_ind])
        go_term_perfs.append(pr)
    return np.array(go_term_perfs)


def load_pred_file_msmnn(fname, loso=False):
    pred_file = pickle.load(open(fname, "rb"))
    go_ids = pred_file['GO_IDs']
    print('Number of go ids: ' + str(len(go_ids)))
    if not loso:
        num_trials = len(pred_file['trial_splits'])
        trial_perfs = np.zeros((num_trials, len(go_ids)))
        for trial in range(0, num_trials):
            curr_trial_test_inds = pred_file['trial_splits'][trial][1]
            curr_trial_train_inds = pred_file['trial_splits'][trial][0]
            curr_trial_preds = pred_file['trial_preds'][trial][curr_trial_test_inds]
            curr_trial_labels = pred_file['true_labels'][curr_trial_test_inds]
            trial_perfs[trial, :] = get_go_term_perfs(curr_trial_labels, curr_trial_preds)
    else: 
        curr_trial_preds = pred_file['preds']
        curr_trial_labels = pred_file['true_labels']
        trial_perfs = get_go_term_perfs(curr_trial_labels, curr_trial_preds)
    return trial_perfs, go_ids


def align_pred_file_with_label_file(pred_dict, label_dict, ont, blast=False):
    trial_macros = []
    trial_micros = []
    trial_accs = []
    trial_f1s = []
    num_trials = len(pred_dict['trial_splits'])
    label_mat = np.asarray(label_dict['annot'][ont].todense())
    label_prot_ids = label_dict['prot_IDs']
    label_go_terms = list(label_dict['go_IDs'][ont])
    trial_preds = []
    trial_labels = []
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
        pred_prots_idx = np.array(pred_prots_idx)
        string_annot_prots_idx = np.array(string_annot_prots_idx)
        pred_go_terms = list(pred_dict['GO_IDs'])
        pred_go_idx, string_annot_go_idx = get_common_indices(pred_go_terms, label_go_terms)
        pred_go_idx = np.array(pred_go_idx)
        print(pred_go_idx[:10])
        string_annot_go_idx = np.array(string_annot_go_idx)
        print(string_annot_go_idx[:10])
        print('Number test prot ids with experimental evidence codes: ' + str(len(pred_prots_idx)))
        print('Number intersection GO ids with experimental evidence codes: ' + str(len(pred_go_idx)))
        try:
            assert len(pred_go_idx) == len(pred_dict['GO_IDs'])
        except AssertionError:
            print('NOT ALL GO TERMS WERE IN EXPERIMENTAL ANNOTATION SET')
        #print('string annot prots idx' + str(len(string_annot_prots_idx)))
        curr_trial_labels = label_mat[np.array(string_annot_prots_idx), :]
        curr_trial_labels = curr_trial_labels[:, string_annot_go_idx]
        curr_trial_preds = curr_trial_preds[pred_prots_idx, :]
        curr_trial_preds = curr_trial_preds[:, pred_go_idx]
        trial_preds.append(curr_trial_preds)
        trial_labels.append(curr_trial_labels)
        go_terms = np.array(pred_go_terms)[pred_go_idx]
    return trial_preds, trial_labels, go_terms


def check_ont(fname):
    onts = ['molecular_function', 'cellular_component', 'biological_process']
    for ont in onts:
        if ont in fname:
            return ont # assume that the first ontology name that's found in fname is the ontology of the file


def load_pred_file_msmnn_w_label_fname(fname, label_fname, loso=False):
    ont = check_ont(fname)
    pred_file = pickle.load(open(fname, "rb"))
    label_file = pickle.load(open(label_fname, "rb"))
    go_ids = pred_file['GO_IDs']
    print('Number of go ids: ' + str(len(go_ids)))
    if not loso:
        num_trials = len(pred_file['trial_splits'])
        trial_preds, trial_labels, intersect_go_ids = align_pred_file_with_label_file(pred_file, label_file, ont)
        trial_perfs = np.zeros((num_trials, len(intersect_go_ids)))
        for trial in range(0, num_trials):
            curr_trial_labels = np.asarray(trial_labels[trial])
            curr_trial_preds = trial_preds[trial]
            trial_perfs[trial, :] = get_go_term_perfs(curr_trial_labels, curr_trial_preds)
    else: 
        '''
        curr_trial_preds = pred_file['preds']
        curr_trial_labels = pred_file['true_labels']
        trial_perfs = get_go_term_perfs(curr_trial_labels, curr_trial_preds)
        '''
        print('LOSO not implemented. Exiting.')
        exit()
    return trial_perfs, intersect_go_ids


def load_pred_file_blast(fname, loso=False):
    pred_file = pickle.load(open(fname, "rb"))
    go_ids = pred_file['go_IDs']
    print('Number of go ids: ' + str(len(go_ids)))
    if not loso:
        num_trials = len(pred_file['Y_hat_test_list'])
        trial_perfs = np.zeros((num_trials, len(go_ids)))
        for trial in range(0, num_trials):
            curr_trial_preds = pred_file['Y_hat_test_list'][trial]
            curr_trial_labels = pred_file['Y_test_list'][trial]
            trial_perfs[trial, :] = get_go_term_perfs(curr_trial_labels, curr_trial_preds)
    else: 
        curr_trial_preds = pred_file['Y_hat_test']
        curr_trial_labels = pred_file['Y_test']
        trial_perfs = get_go_term_perfs(curr_trial_labels, curr_trial_preds)
    return trial_perfs, go_ids


def load_format(fname, label_fname=None, loso=False):
    if '_scores.pckl' in fname:
        if label_fname == None:
            return load_pred_file_blast(fname, loso=loso)
        else:
            print('BLAST not supported with label file yet. Exiting.')
            exit()
    elif '.pckl' in fname:
        if label_fname == None:
            return load_pred_file_msmnn(fname, loso=loso)
        else:
            return load_pred_file_msmnn_w_label_fname(fname, label_fname=label_fname, loso=loso)
    else:
        print('unrecognized file ' + str(fname) + ' ; exiting')
        exit()


def inds_from_list(go_terms_1, go_terms_2):
    return [list(go_terms_2).index(go_terms_1[ind]) for ind in range(0, len(go_terms_1))]


def get_mean_and_sem(trial_perfs):
    assert len(trial_perfs.shape) > 1
    num_trials = trial_perfs.shape[0]
    mean = np.nanmean(trial_perfs, axis=0)
    err = sem(trial_perfs, axis=0, nan_policy='omit') 
    return mean, err

def first_n_last_n(l, n):
    new_l = l[:n]
    new_l.extend(l[-n:])
    return new_l


def remove_indices(*args, indices=None): # removes specified indices from all arrays given as positional arguments
    new_arrs = []
    for arr in args:
        mask = np.ones(np.array(arr).shape, dtype=bool)
        mask[indices] = False
        new_arrs.append(np.array(arr)[mask].copy())
    return new_arrs


def goid_to_name(fname):
    gonames = []
    goterms = []
    f = open(fname, 'r')
    for line in f:
        goid, goname = line.strip().split('\t')
        gonames.append(goname)
        goterms.append(goid)
    f.close()
    goid_name_dict = {goterms[i]: gonames[i] for i in range(0,len(goterms))}

    return goid_name_dict


def truncate_go_names(go_names):
    new_names = []
    for go_name in go_names:
        if len(go_name) < 50:
            new_names.append(go_name)
        else:
            new_names.append(go_name[:47] + '...')
    return new_names


def plot_diff(goterms, gonames, label_1, label_2, auprs_1, auprs_2, errs_1, errs_2, diff, ont, title_label):
    # # Plot # #
    N = len(goterms)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.5       # the width of the bars
    names = [label_1, label_2, 'AUPR diff (middle - first)']

    l = []
    fig, axes = plt.subplots(nrows=1, ncols=3)
    panel_1 = axes[0].barh(ind, auprs_1, width, align='center',
                           color='skyblue', edgecolor='white',
                           ecolor='black', xerr=errs_1)
    l.append(panel_1)
    axes[0].spines['right'].set_visible(False)
    axes[0].spines['bottom'].set_visible(False)
    axes[0].spines['top'].set_visible(False)
    # axes[0].set_xticks(fontsize=10)
    axes[0].set_xlim([0, 1.0])
    axes[0].set_ylim([-1, N + 0.5])
    axes[0].set_yticks(ind)
    axes[0].set_yticklabels(gonames, fontsize=11, rotation=20)
    #axes[0].set_yticklabels(goterms, fontsize=11)
    axes[0].yaxis.set_ticks_position('left')
    axes[0].tick_params(axis='y')
    axes[0].grid(b=True)
    for tick in axes[0].xaxis.get_majorticklabels():
        tick.set_horizontalalignment("right")

    panel_2 = axes[1].barh(ind, auprs_2, width, align='center',
                           color='purple', edgecolor='white',
                           ecolor='black', xerr=errs_2)
    l.append(panel_2)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['bottom'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    # axes[1].set_xticks([])
    axes[1].set_xlim([0, 1.0])
    axes[1].set_ylim([-1, N + 0.5])
    axes[1].set_yticks([])
    axes[1].grid(b=True)

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
    axes[2].grid(b=True)
    # axes[2].axvline(0, color='k')
    
    plt.suptitle(title_label, fontsize=16)
    plt.legend(l, names, loc='upper right', bbox_to_anchor=(0, 0, 1, 1.05), fontsize=12, ncol=5)

def main(fname_1, fname_2, combined_model, label_fname, go_name_fname, loso, label_1, label_2, title_label):
   
    auprs_1, go_ids_1 = load_format(fname_1, label_fname=label_fname, loso=loso) 
    auprs_2, go_ids_2 = load_format(fname_2, label_fname=label_fname, loso=loso) 
    auprs_combined, go_ids_combined = load_format(combined_model, label_fname=label_fname, loso=loso) 

    goid_name_dict = goid_to_name(go_name_fname)
    ont = check_ont(fname_1)
    if ont == 'molecular_function':
        short_ont = 'MF'
    elif ont == 'biological_process':
        short_ont = 'BP'
    elif ont == 'cellular_component':
        short_ont = 'CC'
    print('Assuming ontology is ' + ont)

    assert np.all(go_ids_1 == go_ids_2)
    '''
    try:
        assert np.all(sorted(go_ids_1) == sorted(go_ids_2))
        inds_from_1 = inds_from_list(go_ids_1, go_ids_2)
        auprs_2 = auprs_2[:, inds_from_1]
        go_ids_2 = go_ids_2[inds_from_1]
        try:
            assert np.all(go_ids_1 == go_ids_2)
        except AssertionError:
            print('Sorting from list failed.')
    except AssertionError:
        print(sorted(go_ids_1))
        print(sorted(go_ids_2))
    '''
    diff = auprs_2 - auprs_1
    diff_comb_1 = auprs_combined - auprs_1
    diff_comb_2 = auprs_combined - auprs_2
    if not loso:
        auprs_1, errs_1 = get_mean_and_sem(auprs_1)
        auprs_2, errs_2 = get_mean_and_sem(auprs_2)
        auprs_combined, errs_combined = get_mean_and_sem(auprs_combined)
        diff, diff_err = get_mean_and_sem(diff)
        diff_comb_1, diff_comb_err_1 = get_mean_and_sem(diff_comb_1)
        diff_comb_2, diff_comb_err_2 = get_mean_and_sem(diff_comb_2)
        '''
        remove_inds = list(np.where(auprs_1 == np.nan)[0])
        remove_inds.extend(list(np.where(auprs_2 == np.nan)[0]))
        print('where nan before?')
        print(np.where(auprs_1 == np.nan))
        print('remove inds:')
        print(remove_inds)
        auprs_1, errs_1, auprs_2, errs_2, diff, diff_err = remove_indices(auprs_1, errs_1, auprs_2, errs_2, diff, diff_err, indices=remove_inds)
        print('where nan?')
        print(np.where(auprs_1 == np.nan))
        '''


    idx = sorted(range(diff.shape[0]), key=lambda x: diff[x])
    #idx = np.where(errs_2 < 0.1)
    #print(errs_2)
    diff = [diff[i] for i in idx]
    better_go = 0
    for i in idx:
        if diff[i] >= 0:
            better_go += 1

    print ("### Number of better goterms: %d" % (better_go))
    print ("### Number of total goterms: %d" % (len(diff)))
    goterms = [go_ids_1[i] for i in idx]
    #gonames = [gonames[i] for i in idx]
    auprs_1 = [auprs_1[i] for i in idx]
    auprs_2 = [auprs_2[i] for i in idx]
    # auprs_3 = [auprs_3[i] for i in idx]
    n = 10
    print('Top ' + str(n) + ' better go ids for ' + label_1 + ' and performances:')
    print(goterms[:n])
    print('Perf for ' + label_1 + ': ')
    print(auprs_1[:n])
    print('Perf for ' + label_2 + ': ')
    print(auprs_2[:n])
    print('Top ' + str(n) + ' better go ids for ' + label_2 + ' and performances:')
    print(goterms[-n:])
    print('Perf for ' + label_1 + ': ')
    print(auprs_1[-n:])
    print('Perf for ' + label_2 + ': ')
    print(auprs_2[-n:])
    '''
    print('All better GO ids for ' + label_1)
    print(np.array(goterms)[np.array(auprs_1) > np.array(auprs_2)])
    print('All better GO ids for ' + label_2)
    print(np.array(goterms)[np.array(auprs_2) > np.array(auprs_1)])
    '''

    errs_1 = [errs_1[i] for i in idx]
    errs_2 = [errs_2[i] for i in idx]

    auprs_combined = [auprs_combined[i] for i in idx]
    errs_combined = [errs_combined[i] for i in idx]
    diff_comb_1 = [diff_comb_1[i] for i in idx]
    diff_comb_2 = [diff_comb_2[i] for i in idx]

    auprs_1 = first_n_last_n(auprs_1, n)
    auprs_2 = first_n_last_n(auprs_2, n)
    errs_1 = first_n_last_n(errs_1, n)
    errs_2 = first_n_last_n(errs_2, n)
    goterms = first_n_last_n(goterms, n)
    diff = first_n_last_n(diff, n)

    auprs_combined = first_n_last_n(auprs_combined, n)
    errs_combined = first_n_last_n(errs_combined, n)
    diff_comb_1 = first_n_last_n(diff_comb_1, n)
    diff_comb_2 = first_n_last_n(diff_comb_2, n)

    remove_inds = []
    for i, go_term in enumerate(goterms):
        if go_term not in goid_name_dict.keys(): # if it's not in the list, remove it
            remove_inds.append(i)
    (errs_1, errs_2, auprs_1, auprs_2, goterms, diff, 
            auprs_combined, errs_combined, diff_comb_1, diff_comb_2) = remove_indices(errs_1, 
                    errs_2, auprs_1, auprs_2, goterms, diff, auprs_combined, errs_combined, 
                    diff_comb_1, diff_comb_2, indices=remove_inds) 
    gonames = [goid_name_dict[go_term] for go_term in goterms]

    gonames = truncate_go_names(gonames)
    if title_label is None:
        title_label = 'Top ' + str(n) + ' GO Terms with Positive Performance Difference For Each Model (AUPR) -- ' + short_ont
    plot_diff(goterms, gonames, label_1, label_2, auprs_1, auprs_2, errs_1, errs_2, diff, ont, title_label)

    # okay, now I want the GO term performances when combining the models, compared to each individual one, for each of those GO terms that are in the top 10
    combined_title = 'AUPR Comparison of Combined Model with Single-source Model'
    plot_diff(goterms, gonames, label_1, 'Combined ' + label_1 + '/' + label_2, auprs_1, auprs_combined, errs_1, errs_combined, diff_comb_1, ont, combined_title)
    plot_diff(goterms, gonames, label_2, 'Combined ' + label_1 + '/' + label_2, auprs_2, auprs_combined, errs_2, errs_combined, diff_comb_2, ont, combined_title)
    
    plt.show()
     

if __name__ == "__main__":
    '''
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('fn_1', type=str)
    parser.add_argument('fn_2', type=str)
    parser.add_argument('combined_model', type=str)
    parser.add_argument('--loso', action='store_true', type=bool,
                    help='whether to do leave one species out (not implemented)
    '''

    fn_1 = str(sys.argv[1])
    fn_2 = str(sys.argv[2])
    combined_model = str(sys.argv[3])
    label_fname = str(sys.argv[4])
    goname_fname = str(sys.argv[5])
    loso = sys.argv[6] == 'loso'
    label_1 = sys.argv[7]
    label_2 = sys.argv[8]
    if len(sys.argv) > 9:
        title_label = sys.argv[9]
        main(fn_1, fn_2, combined_model, label_fname, goname_fname, loso, label_1, label_2, title_label)
    else:
        main(fn_1, fn_2, combined_model, label_fname, goname_fname, loso, label_1, label_2, None)
