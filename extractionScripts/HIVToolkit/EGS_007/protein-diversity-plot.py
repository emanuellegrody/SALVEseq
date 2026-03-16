import pandas as pd
import matplotlib.pyplot as plt

# Read the file
df = pd.read_csv('/Users/egy2296/PycharmProjects/SALVEseq/HIVToolkit_myEdits/GeneticDiversity_AminoAcids.txt', header=None, names=['Protein', 'Position', 'Diversity'])

# Filter for diversity >= 0.1
df_filtered = df[df['Diversity'] >= 0.1]

# Save filtered data to a new .txt file
output_file = '/Users/egy2296/PycharmProjects/SALVEseq/HIVToolkit_myEdits/EGS_007/high_diversity_positions.txt'
df_filtered.to_csv(output_file, header=False, index=False, sep=',')
print(f"Filtered data saved to {output_file}")


# Create a figure and axis
fig, ax = plt.subplots(figsize=(12, 6))

# Get unique proteins
proteins = df['Protein'].unique()

# Plot each protein end-to-end
start_position = 0
for i, protein in enumerate(proteins):
    protein_data = df_filtered[df_filtered['Protein'] == protein]

    if not protein_data.empty:
        # Adjust positions to be end-to-end
        adjusted_positions = protein_data['Position'] + start_position

        # Plot the data
        ax.scatter(adjusted_positions, protein_data['Diversity'], label=protein)

        # Add protein name annotation
        midpoint = start_position + (protein_data['Position'].max() - protein_data['Position'].min()) / 2
        ax.annotate(protein, (midpoint, -0.01), rotation=45, ha='right', va='top')

        # Update start position for the next protein
        start_position += protein_data['Position'].max()

# Set labels and title
ax.set_xlabel('Cumulative Position')
ax.set_ylabel('Diversity')
ax.set_title('High Diversity Positions (>=0.1) Across Proteins')

# Adjust y-axis to start from 0
ax.set_ylim(bottom=0)

# Add legend
ax.legend(title="Proteins", bbox_to_anchor=(1.05, 1), loc='upper left')

# Adjust layout to prevent cutting off annotations
plt.tight_layout()

# Show the plot
#plt.show()

# Optionally, save the plot
plt.savefig('/Users/egy2296/PycharmProjects/SALVEseq/HIVToolkit_myEdits/EGS_007/protein_diversity_plot.svg', format='svg', bbox_inches='tight', dpi=300)
