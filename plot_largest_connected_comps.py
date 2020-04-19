import networkx as nx
import pickle
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use('ggplot')

'''
big_font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **big_font)
'''
#tax_id_sci_name_df = pd.read_csv('taxonomy-filtered-reviewed_yes+AND+annotated_yes.tab', delimiter='\t')
#tax_id_to_sci_name = dict(zip(tax_id_sci_name_df['Taxon'].astype(str), tax_id_sci_name_df['Scientific name']))
#tax_id_to_sci_name['1148'] = 'Synechocystis sp. (strain PCC 6803)'


annot_file = sys.argv[1] # annot file for all organisms in this set
title_keyword = sys.argv[2]
network_files = sys.argv[3:] # individual network pickle files for all organisms in the set

taxon_to_name = {}
with open('taxon_to_names.txt') as f:
    for line in f:
        fields = line.split('\t')
        taxon_to_name[fields[0]] = fields[1].strip()

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

    #largest_con_comp = max(nx.connected_component_subgraphs(graph), key=len)
    largest_con_comp = max(nx.connected_components(graph), key=len)
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

sorted_taxa_names = np.array([taxon_to_name[taxon] for taxon in sorted_taxa]) # changing taxon list to names
#sorted_sci_names = [tax_id_to_sci_name[taxon] for taxon in sorted_taxa]
sorted_net_ratios = np.array(net_ratios)[sort_inds]
sorted_annot_ratios = np.array(annot_ratios)[sort_inds]
sorted_num_nodes = np.array(num_nodes)[sort_inds]

rects_1 = ax.bar(x_positions-0.15, sorted_net_ratios, width, label='N_c/N')
rects_2 = ax.bar(x_positions+0.15, sorted_annot_ratios, width, label='Proportion annotated')

ax.set_ylabel('N_c/N')
ax.set_xlabel('Taxa - Total number of proteins: ' + str(int(sum(num_nodes))))
ax.set_title('Ratio of largest connected components to total size of graph - ' + title_keyword)
ax.set_xticks(x_positions)
ax.set_yticks(np.arange(0,11)/10)
#ax.set_xticklabels(sorted_sci_names)
print(sorted_taxa_names)
ax.set_xticklabels(sorted_taxa_names, rotation=45, ha='right')

ax.legend(loc='upper left')

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for i, rect in enumerate(rects):
        height = rect.get_height()
        annot = ax.annotate(str(int(sorted_num_nodes[i])),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=12)
        annot.set_color('red')
        '''
        ax.annotate('{0:.2f}'.format(height) + '\n' + str(int(sorted_num_nodes[i])),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=14)
        '''

autolabel(rects_1)
autolabel(rects_2)
fig.tight_layout()
plt.show()
