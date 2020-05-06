import networkx as nx
import pickle
import sys
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

plt.style.use('ggplot')

big_font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 26}

matplotlib.rc('font', **big_font)
#tax_id_sci_name_df = pd.read_csv('taxonomy-filtered-reviewed_yes+AND+annotated_yes.tab', delimiter='\t')
#tax_id_to_sci_name = dict(zip(tax_id_sci_name_df['Taxon'].astype(str), tax_id_sci_name_df['Scientific name']))
#tax_id_to_sci_name['1148'] = 'Synechocystis sp. (strain PCC 6803)'


annot_file = sys.argv[1] # annot file for all organisms in this set
title_keyword = sys.argv[2]
network_files = sys.argv[3:] # individual network pickle files for all organisms in the set

taxon_to_name = {}
with open('taxon_to_names.dat') as f:
    for line in f:
        fields = line.split('\t')
        taxon_to_name[fields[0]] = fields[1].strip()

net_ratios = []
taxa = []
num_nodes = []
num_edges = []
for network_path in network_files:
    network_file = network_path.split('/')[-1]
    taxon = network_file.split('_')[0]
    taxa.append(taxon)
    net_pickle = pickle.load(open(network_path, 'rb'))
    net = net_pickle['nets']['experimental']
    del net_pickle

    graph = nx.from_scipy_sparse_matrix(net)
    del net

    #largest_con_comp = max(nx.connected_component_subgraphs(graph), key=len)
    largest_con_comp = max(nx.connected_components(graph), key=len)
    n = float(len(nx.nodes(graph)))
    ne = len(graph.edges)
    del graph
    nc = float(len(largest_con_comp))
    del largest_con_comp
    num_nodes.append(n)
    num_edges.append(ne)
    nc_over_n = nc/n
    net_ratios.append(nc_over_n)
    
taxa_names = [taxon_to_name[taxon] for taxon in taxa]
df = pd.DataFrame([taxa_names, taxa, num_nodes, num_edges, net_ratios]).transpose()
df.columns = ['Scientific name', 'Taxonomy ID', 'Number of Nodes', 'Number of Edges', 'Ratio of Largest Connected Component']
print(df.to_latex())
annot = pickle.load(open(annot_file, 'rb'))
annot_mats = annot['annot']

branches = ['molecular_function', 'biological_process', 'cellular_component']
annot_mat_list = [annot_mats[branch] for branch in branches]
print(annot_mat_list)
all_branch_annots = sparse.hstack(annot_mat_list).tocsr()
prot_ids = annot['prot_IDs']
# Need the number of proteins in each organism, and separated annotation matrices for each
org_of_prots = np.array([prot.split('.')[0] for prot in prot_ids])
species_inds = {}

annot_ratios = []
for i, taxon in enumerate(taxa):
    species_inds[taxon] = np.where(org_of_prots == taxon)[0]
    #num_annotated_prots = np.sum(all_branch_annots[species_inds[taxon]].any(1))
    nnz_rows, _  = all_branch_annots[species_inds[taxon], :].nonzero()
    num_annotated_prots = len(set(nnz_rows))
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

rects_1 = ax.bar(x_positions-0.15, sorted_net_ratios, width, label=r'Proportion of nodes in largest connected component, $\frac{N_c}{N}$')
rects_2 = ax.bar(x_positions+0.15, sorted_annot_ratios, width, label=r'Proportion annotated')

#ax.set_ylabel(r'$\frac{N_c}{N}$')
#ax.set_xlabel('Taxa - Total number of proteins: ' + str(int(sum(num_nodes))))
ax.set_xlabel('Taxa')
ax.set_title(r'\textbf{Network and Annotation Completeness - ' + title_keyword + '}', y=1.08)
ax.set_xticks(x_positions)
ax.set_yticks(np.arange(0,11)/10)
#ax.set_xticklabels(sorted_sci_names)
print(sorted_taxa_names)
sorted_taxa_names_formatted = [r'\textbf{' + name + r'}' for name in sorted_taxa_names]
ax.set_xticklabels(sorted_taxa_names_formatted, rotation=45, ha='right')

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.5))


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


'''
# TODO: calculate network/annotation measures for detailed table in supplement
def create_latex_table(taxon_to_name, taxons, num_nodes):
    with open(title_keyword + 'num_nodes_table.tex'):
'''


#autolabel(rects_1)
#autolabel(rects_2)
# no annotating with num prots for now, just make the plots simple
fig.tight_layout()
plt.show()
