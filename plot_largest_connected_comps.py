import networkx as nx
import pickle
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

annot_file = sys.argv[1] # annot file for all organisms in this set
title_keyword = sys.argv[2]
network_files = sys.argv[3:] # individual network pickle files for all organisms in the set

net_ratios = []
taxa = []
num_nodes = []
for network_path in network_files:
    network_file = network_path.split('/')[-1]
    taxon = network_file.split('_')[0]
    taxa.append(taxon)
    net_pickle = pickle.load(open(network_path, 'rb'))
    net = net_pickle['nets']['experimental']

    graph = nx.from_scipy_sparse_matrix(net)

    largest_con_comp = max(nx.connected_component_subgraphs(graph), key=len)
    nc = float(len(largest_con_comp))
    n = float(len(nx.nodes(graph)))
    num_nodes.append(n)
    nc_over_n = nc/n
    net_ratios.append(nc_over_n)
    print(nc_over_n)

annot = pickle.load(open(annot_file, 'rb'))
annot_mats = annot['annot']

branches = ['molecular_function', 'biological_process', 'cellular_component']
annot_mat_list = [annot_mats[branch].todense() for branch in branches]
print(annot_mat_list)
all_branch_annots = np.concatenate(annot_mat_list, axis=1)
prot_ids = annot['prot_IDs']
# Need the number of proteins in each organism, and separated annotation matrices for each
org_of_prots = np.array([prot.split('.')[0] for prot in prot_ids])
species_inds = {}

annot_ratios = []
for i, taxon in enumerate(taxa):
    species_inds[taxon] = np.where(org_of_prots == taxon)[0]
    num_annotated_prots = np.sum(all_branch_annots[species_inds[taxon]].any(1))
    annot_ratios.append(float(num_annotated_prots)/float(num_nodes[i]))

fig, ax = plt.subplots()
x_positions = np.arange(len(taxa))
width = 0.3
#rects_1 = ax.bar(x_positions-0.25, net_ratios, width, label='N_c/N')
#rects_2 = ax.bar(x_positions+0.25, annot_ratios, width, label='Proportion annotated')

sums = np.array(net_ratios) + np.array(annot_ratios)
sort_inds = np.argsort(sums)
sorted_taxa = np.array(taxa)[sort_inds]
sorted_net_ratios = np.array(net_ratios)[sort_inds]
sorted_annot_ratios = np.array(annot_ratios)[sort_inds]
rects_1 = ax.bar(x_positions-0.15, sorted_net_ratios, width, label='N_c/N')
rects_2 = ax.bar(x_positions+0.15, sorted_annot_ratios, width, label='Proportion annotated')

ax.set_ylabel('N_c/N')
ax.set_xlabel('Taxa - Total number of proteins: ' + str(int(sum(num_nodes))))
ax.set_title('Ratio of largest connected components to total size of graph - ' + title_keyword)
ax.set_xticks(x_positions)
ax.set_xticklabels(taxa)

ax.legend()

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for i, rect in enumerate(rects):
        height = rect.get_height()
        ax.annotate('{0:.2f}'.format(height) + '-' + str(int(num_nodes[i])),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(rects_1)
autolabel(rects_2)
fig.tight_layout()
plt.show()
