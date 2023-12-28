import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the phenotypic data from the CSV file
phenotypic_data = pd.read_csv('new_output.txt', header=None)

# Read the PCA eigenvectors data from the file
with open('pcaoneFreq50.eigvecs', 'r') as file:
    eigvecs = np.array([list(map(float, line.strip().split())) for line in file])

# Row numbers to consider
phenotype_rows = [1, 2, 3, 4, 5, 6]

# Loop through each phenotype row and create a separate scatter plot and save it as a PNG file
for row_number in phenotype_rows:
    # Extract specific row from phenotypic data
    selected_phenotype = phenotypic_data.iloc[row_number - 1]  # Row numbers are 1-based

    # Subset data to reduce the number of points or increase the step value
    step = 5
    subset_indices = np.arange(0, eigvecs.shape[0], step)
    subset_eigvecs = eigvecs[subset_indices, :]

    # Create a scatter plot with colored points based on the selected phenotype
    plt.figure(figsize=(8, 6))
    plt.scatter(subset_eigvecs[:, 0], subset_eigvecs[:, 1], c=selected_phenotype.values[subset_indices], cmap='viridis', alpha=0.6)  # Change cmap to the desired color map
    plt.xlabel('PCA Component 1')
    plt.ylabel('PCA Component 2')
    plt.colorbar(label='Phenotype Value')
    plt.title(f'PCA Scatter Plot Colored by Phenotype Row {row_number}')

    # Save the scatter plot as a PNG file
    plt.savefig(f'pca_scatter_plot_row_{row_number}.png')

    # Close the current plot to release memory for the next plot
    plt.close()
