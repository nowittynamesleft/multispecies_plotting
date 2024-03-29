"""
========
Barchart (final results)
========

A bar plot with errorbars and height labels on individual bars
"""
import numpy as np
import pickle
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


"""
# THIS PART MAKES A LIST OF GO IDS IN ALL THE TAXA IN THE ANNOTATION FILE (THAT PASS THE FIRST CRITERIA OF > 30 AND < 150 EXAMPLES)


single_pckl = str(sys.argv[1])
go_type = str(sys.argv[2])

single_annot = pickle.load(open(single_pckl, 'rb'))
single_Y = single_annot['annot'][go_type].toarray()
single_goids = single_annot['go_IDs'][go_type]
single_gonames = single_annot['go_names'][go_type]
single_prots = np.asarray(single_annot['prot_IDs'])
num_total_taxa = 6

test_funcs = np.where(np.logical_and(single_Y.sum(axis=0) > 30, single_Y.sum(axis=0) <= 150))[0] # test_func_inds

common_goids = []
common_gonames = []
for ii in test_funcs:
    go_id = single_goids[ii]
    go_name = single_gonames[ii] 
    indx = np.where(single_Y[:, ii] > 0)[0] # find all inds of prots with this go term
    tmp_prots = single_prots[indx] # names of those prots
    tax = set([prot.split('.')[0] for prot in tmp_prots]) # names of taxa for those prots
    if len(tax) == num_total_taxa: # if there are all the taxa in the set just made, this go_id is present in all the taxa
        common_goids.append(go_id)  # so add it to the common_goids list
        common_gonames.append(go_name)
        print (go_id, '\t', go_name)

pickle.dump(common_goids, open("bacteria_" + go_type + "_train_goids.pckl", "wb")) # these are the goids for the thing
"""

# This checks whether the single org annot list with the single org threshold has go terms that have at least 50 examples in the multi org list (i.e., more than 20 examples not in the single org annot list)

single_pckl = str(sys.argv[1])
multi_pckl = str(sys.argv[2])

single_annot = pickle.load(open(single_pckl, 'rb'))
multi_annot = pickle.load(open(multi_pckl, 'rb'))
go_type = 'molecular_function'

single_Y = single_annot['annot'][go_type].toarray()
single_goids = single_annot['go_IDs'][go_type]
single_gonames = single_annot['go_names'][go_type]

multi_Y = multi_annot['annot'][go_type].toarray()
multi_prots = multi_annot['prot_IDs']
idx = [ii for ii, prot in enumerate(multi_prots) if not prot.startswith("553174")] # checks whether it's ecoli strain K12 substrain whatever (highest annotated bacteria)
multi_prots = np.asarray(multi_prots)
multi_prots = multi_prots[idx]
multi_Y = multi_Y[idx]
multi_goids = multi_annot['go_IDs'][go_type]

test_funcs = np.where(np.logical_and(single_Y.sum(axis=0) > 30, single_Y.sum(axis=0) <= 100))[0]

common_goids = []
common_gonames = []
all_taxa = set()
for ii in test_funcs:
    go_id = single_goids[ii]
    go_name = single_gonames[ii]
    if go_id in multi_goids: # if the test func for the single org annot thing is in the multigoid list
        idx = multi_goids.index(go_id) 
        if sum(multi_Y[:, idx]) > 50: # if the number of examples in the multi org annot are more than 50
            indx = np.where(multi_Y[:, idx] > 0)[0]
            tmp_prots = multi_prots[indx] # 
            tax = set([prot.split('.')[0] for prot in tmp_prots]) # taxa that are included in the multiple taxa annotations
            all_taxa = all_taxa.union(tax)
            if len(tax) == 5:
                common_goids.append(go_id)
                common_gonames.append(go_name)
                print (go_id, '\t', go_name)


pickle.dump(common_goids, open("553174-model-org_" + go_type + "_train_goids.pckl", "wb"))

single_idx = [single_goids.index(go) for go in common_goids]
multi_idx = [multi_goids.index(go) for go in common_goids]

single_Y = single_Y[:, single_idx]
multi_Y = multi_Y[:, multi_idx]

single_count = single_Y.sum(axis=0)
multi_count = multi_Y.sum(axis=0)
multi_count += single_count

# # Plot # #
ind = np.arange(len(common_goids))  # the x locations for the groups
width = 0.20       # the width of the bars
N = max(ind)


cs = ['skyblue', 'purple']

fig, ax = plt.subplots(1, 1)

ax.bar(ind, single_count, width, color=cs[0], align='center')
ax.bar(ind + width, multi_count, width, color=cs[1], align='center')
ax.set_ylabel("# Proteins", fontsize=18)
# ax.set_title('Yeast: ' + yeast_annot2name[yeast_level], fontsize=18)
# axarr[0].set_xlim([-0.1, N])
# ax.set_ylim([0, max(yeast_deepnf_res[yeast_level])+0.05])
ax.legend(('H', 'H + Y + M + F + W + B'), fontsize=18, loc='upper left')
ax.set_xticks(ind+width)

ax.set_xticklabels(common_gonames, fontsize=8, rotation=30, rotation_mode="anchor")
for tick in ax.xaxis.get_majorticklabels():
    tick.set_horizontalalignment("right")

plt.savefig('-'.join(all_taxa) + go_type + '_annot_stats_mult_single.png', bbox_inches='tight')
# plt.show()
