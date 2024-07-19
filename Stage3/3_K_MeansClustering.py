import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
data = pd.read_csv("input_data50.tsv", sep="\t")

# Transpose the data
data_t = data.T

# Determine the optimal number of clusters using the elbow method
inertia = []
for k in range(1, 15):
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(data_t)
    inertia.append(kmeans.inertia_)

# Plot the elbow method results and save the figure
plt.figure(figsize=(8, 5))
plt.plot(range(1, 15), inertia, marker='o')
plt.xlabel('Number of Clusters')
plt.ylabel('Inertia')
plt.title('Elbow Method for Optimal Number of Clusters')
plt.savefig('elbow_method.png')
plt.show()

# Based on the elbow plot, choose the optimal number of clusters, e.g., k=5
optimal_k = 5
kmeans = KMeans(n_clusters=optimal_k, random_state=42)
clusters = kmeans.fit_predict(data_t)

# Add cluster labels to the transposed data
data_t['Cluster'] = clusters

data_t.Cluster.sort_values().to_csv('K-means_clusters.csv')

# Save the clustered data
data_t.to_csv("clustered_data.csv")

# Visualize the clusters and save the figure
plt.figure(figsize=(10, 7))
sns.heatmap(data_t.sort_values('Cluster'), cmap='viridis', yticklabels=False)
plt.title('Phenotype Clusters')
plt.savefig('phenotype_clusters.png')
plt.show()

print("Clustering completed and results saved to clustered_data.csv")
