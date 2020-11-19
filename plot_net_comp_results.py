import pickle
import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np

x = pickle.load(open(sys.argv[1], 'rb'))
fig, ax = plt.subplots(2, 3)

def shorten_names(names):
    new_names = []
    for name in names:
        name_split = name.split(' ')
        first = name_split[0][0] + '. '
        print(name_split)
        other = ' '.join(name_split[1:])
        if len(other) > 20:
            new_name = first + other[:17] + '...'
        else:
            new_name = first + other
        new_names.append(new_name)
    return new_names


taxon_to_name = {}
with open('taxon_to_names.dat') as f:
    for line in f:
        fields = line.split('\t')
        taxon_to_name[fields[0]] = fields[1].strip()
taxon_names = [taxon_to_name[tax] for tax in x['taxa']]
taxon_names = shorten_names(taxon_names)

im = ax[0,0].imshow(x['avg_blast_mat'])
ax[0,0].set_title('Average BLAST edge')
ax[0,0].set_xticks(np.arange(0, len(taxon_names)))
ax[0,0].set_yticks(np.arange(0, len(taxon_names)))
ax[0,0].set_xticklabels(taxon_names)
ax[0,0].set_yticklabels(taxon_names)
fig.colorbar(im, ax=ax[0,0])

im = ax[0,1].imshow(x['avg_isorank_mat'])
ax[0,1].set_title('Average IsoRank edge')
ax[0,1].set_xticks(np.arange(0, len(taxon_names)))
ax[0,1].set_yticks(np.arange(0, len(taxon_names)))
ax[0,1].set_xticklabels(taxon_names)
ax[0,1].set_yticklabels(taxon_names)
fig.colorbar(im, ax=ax[0,1])

im = ax[1,0].imshow(x['blast_density'])
ax[1,0].set_title('BLAST Matrix Density')
ax[1,0].set_xticks(np.arange(0, len(taxon_names)))
ax[1,0].set_yticks(np.arange(0, len(taxon_names)))
ax[1,0].set_xticklabels(taxon_names)
ax[1,0].set_yticklabels(taxon_names)
fig.colorbar(im, ax=ax[1,0])

im = ax[1,1].imshow(x['isorank_density'])
ax[1,1].set_title('IsoRank Matrix Density')
ax[1,1].set_xticks(np.arange(0, len(taxon_names)))
ax[1,1].set_yticks(np.arange(0, len(taxon_names)))
ax[1,1].set_xticklabels(taxon_names)
ax[1,1].set_yticklabels(taxon_names)
fig.colorbar(im, ax=ax[1,1])

im = ax[0,2].imshow(x['clust_dist_mat'])
ax[0,2].set_title('Average Clustering Coeff. Difference')
ax[0,2].set_xticks(np.arange(0, len(taxon_names)))
ax[0,2].set_yticks(np.arange(0, len(taxon_names)))
ax[0,2].set_xticklabels(taxon_names)
ax[0,2].set_yticklabels(taxon_names)
fig.colorbar(im, ax=ax[0,2])
ax[1,2].axis('off')

for subplot in ax.flatten():
    plt.setp(subplot.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

plt.suptitle(sys.argv[2] + ' Network Comparisons -- BLAST, IsoRank and Average Clustering Coefficient Difference')
plt.show()
