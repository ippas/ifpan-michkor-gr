import numpy as np
import matplotlib.pyplot as plt

# Read eigenvalues from file
with open('pcaoneFreq50.eigvals', 'r') as file:
    eigenvalues = np.loadtxt(file)

# Plot the scree plot with dark red bars
plt.figure(figsize=(10, 6))
plt.bar(range(1, len(eigenvalues) + 1), eigenvalues, alpha=0.7, align='center', color='darkgreen', label='Eigenvalues')
plt.xlabel('Principal Components')
plt.ylabel('Eigenvalues')
plt.title('Scree Plot')
plt.legend()

# Save the plot as ScreePlot.png
plt.savefig('ScreePlot.png')

