import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import os
import matplotlib.colors as mcolors


CNA_link = r"C:/Users/lyaziddy/Downloads/brca_tcga/brca_tcga/CNA_data_cleaned.tsv"

df1 = pd.read_csv(CNA_link, sep="\t", index_col=False)

#print(df1)

overlapping_genes = ['ABCG8', 'ADM', 'ADRB2', 'AGT', 'AMPD1', 'ANGPTL3', 'ANGPTL4', 'ANKRD1', 'APOA5', 'APOB', 'CD34', 'CLCNKB', 'CORIN', 'CYP17A1', 'CYP2A6', 'CYP7A1', 'EDN1', 'ENO2', 'EPO', 'FGF2', 'GJB6', 'HBA1', 'HOTAIR', 'HSD11B2', 'HSPB7', 'IL6', 'KLHL3', 'KNG1', 'LPL', 'MB', 'MME', 'MRAS', 'MYH6', 'MYH7', 'MYL3', 'NKX2-5', 'NPR1', 'NR3C2', 'PCSK9', 'PLN', 'PRKD1', 'SLC12A1', 'SLC12A3', 'SLC2A1', 'TAC1', 'TCF7L2', 'TEK', 'TNNI3', 'UMOD', 'UTS2', 'WNK4']
#print(overlapping_genes)
#print(f"Number of overlapping genes: {len(overlapping_genes)}")

all_genes_in_df1 = set(overlapping_genes).issubset(df1['Hugo_Symbol'].values)
missing_genes = set(overlapping_genes) - set(df1['Hugo_Symbol'].values)
if missing_genes:
    print(f"These genes are missing from df1: {missing_genes}")
else:
    print("All genes are found in df1.")


df1_relevant = df1[df1['Hugo_Symbol'].isin(overlapping_genes)]
#print(df1_relevant)


#Calculate the sum, the average, min and max for each gene accross all samples ...
#Make a plot where higher CNAs are on top ... hypothesis: on a wider cohort, higher sum CNA means higher individual CNA if variability is controlled and range is small so this approach is valid
#Calculate also the number of each unique value for each gene ... try also to plot in a heatmap the concatenated number ot the number of features with their CNAs ... Genes in y axis and x CNA



df1_relevant.set_index("Hugo_Symbol", inplace = True)
#print(df1_relevant)

df2 = pd.DataFrame(index=df1_relevant.index)

# Calculate statistics for each row (gene)
df2['Min'] = df1_relevant.min(axis=1)
df2['Max'] = df1_relevant.max(axis=1)
df2['Average'] = df1_relevant.mean(axis=1)
df2['Sum'] = df1_relevant.sum(axis=1)

# Calculate unique value counts for each row (gene)
df2['Unique_Counts'] = df1_relevant.apply(lambda row: row.value_counts().to_dict(), axis=1)
#print(df2)

numeric_df = df1_relevant.apply(pd.to_numeric, errors='coerce')
#print(numeric_df)


# Sort each row (gene) in ascending order
#df1_sorted = pd.DataFrame(numeric_df.apply(lambda x: sorted(x), axis=1).tolist(), index=numeric_df.index)


# Plot heatmap
#plt.figure(figsize=(20, 12))
##sns.heatmap(df1_sorted, cmap='viridis', annot=False, cbar=True, xticklabels=False)     #, linewidths=0.5
#plt.title('Heatmap of Sorted Gene Expression Values')
#plt.xlabel('All tumor samples')
#plt.ylabel('Genes')
#plt.show()

# Sort df2 by the 'Sum' column in descending order
df2_sorted = df2.sort_values(by='Sum', ascending=False)

# Create a custom colormap
cmap_neg = mcolors.LinearSegmentedColormap.from_list('red_palette', ['darkred', 'red'])
cmap_pos = mcolors.LinearSegmentedColormap.from_list('green_palette', ['lightgreen', 'darkgreen'])

# Combine the colormaps for a single colormap
def combined_cmap(val):
    if val < 0:
        return cmap_neg((val - min(df2_sorted['Sum'])) / (0 - min(df2_sorted['Sum'])))
    else:
        return cmap_pos((val - 0) / (max(df2_sorted['Sum']) - 0))

# Normalize data for colormap
norm = mcolors.Normalize(vmin=df2_sorted['Sum'].min(), vmax=df2_sorted['Sum'].max())

# Create a color list for the bars
colors = [combined_cmap(val) for val in df2_sorted['Sum']]

# Plot horizontal bar plot with custom colormap
plt.figure(figsize=(6, 10 # Adjust height for better space

# Plot horizontal bar plot
bars = plt.barh(df2_sorted.index, df2_sorted['Sum'], color=colors, height=0.6)  # Adjust height for more spacing

# Add plot title and axis labels
plt.title('Sum of CNAs in selected genes across all tumor samples', fontsize=12)
plt.xlabel('Sum of CNA values', fontsize=10)
plt.ylabel('Selected relevant Genes', fontsize=10)

# Further increase spacing between y-axis labels and reduce margins
plt.subplots_adjust(left=0.35, right=0.95, top=0.008, bottom=0.007)  # Reduce top and bottom margins

# Use tight layout to adjust plot elements
plt.tight_layout()

# Add grid for better readability
plt.grid(axis='x', linestyle='--', alpha=0.7)

# Display plot
plt.show()