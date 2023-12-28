import pandas as pd

# Read loadings data from the tab-separated file without header
loadings_df = pd.read_csv('pcaoneFreq50.loadings', delimiter='\t', header=None)

# Calculate the absolute values of the loadings for components 1-10
components = list(range(10))  # Components 1 to 10
for component in components:
    loadings_df[component] = abs(loadings_df[component])

# Sort the phenotypes based on their influence on each component
top_phenotypes = {}
for component in components:
    top_phenotypes[component + 1] = loadings_df[component].sort_values(ascending=False).index.tolist()

# Save the output to a file
with open('top_phenotypes_outputFreq50.txt', 'w') as f:
    for component, phenotypes in top_phenotypes.items():
        f.write(f"Top phenotypes for Component {component}:\n")
        f.write(str(phenotypes))
        f.write("\n\n")

print("Top phenotypes for components 1 to 10 saved to top_phenotypes_output.txt")
