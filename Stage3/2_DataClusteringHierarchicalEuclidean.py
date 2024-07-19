import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, inconsistent
import matplotlib.pyplot as plt
import numpy as np

# Read data from TSV file
data = pd.read_csv('input_data50.tsv', sep='\t', index_col=0)

# Perform hierarchical clustering using correlation-based distance and average linkage
# Change the method and metric according to your preference
clusters = linkage(data.T, method='average', metric='euclidean')

# Plot dendrogram
plt.figure(figsize=(12, 6))
dendrogram(clusters, labels=data.T.index, leaf_rotation=90)
plt.title('Hierarchical Clustering Dendrogram for Phenotypes')
plt.xlabel('Phenotypes')
plt.ylabel('Distance')

# Automatically determine the number of clusters using the inconsistency criterion
depth = 50  # Depth to calculate the inconsistency statistic
incons = inconsistent(clusters, depth)
threshold = 0.8  # Adjust this threshold as needed
num_clusters = 0
for i in range(1, len(incons)):
    if incons[i][3] - incons[i - 1][3] > 0.05:
        num_clusters = i + 1

# Assign cluster labels to each phenotype based on the determined number of clusters
# cluster_labels = fcluster(clusters, num_clusters, criterion='maxclust')
cluster_labels = fcluster(clusters, t=0.864, criterion='distance')

# Save dendrogram as PNG file
plt.savefig('dendrogram.png', dpi=300)  # Change the file name and dpi as needed
with open('cluster_labels.npy', 'wb') as f:
    np.save(f, cluster_labels)
print(num_clusters)

with open('cluster_labels.npy', 'rb') as f:
    cluster_labels = np.load(f)

# Save cluster assignments and list of phenotypes for each cluster to a text file
with open('cluster_assignments.txt', 'w') as f:
    for cluster_num in range(1, num_clusters + 1):
        phenotypes_in_cluster = data.columns[cluster_labels == cluster_num]
        f.write(f'Cluster {cluster_num}:\n')
        for phenotype in phenotypes_in_cluster:
            f.write(str(phenotype) + '\n')
        f.write('\n')

plt.show()
