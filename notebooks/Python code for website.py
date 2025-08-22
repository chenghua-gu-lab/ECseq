#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection


# In[2]:


import re
import seaborn as sns


# Importing the dataset from Albert (with expression data and some domain information) 

# In[3]:


from scipy import stats
from sklearn.preprocessing import PowerTransformer


# In[4]:


df = pd.read_csv('/Users/shamika.bhandarkar/Desktop/Expression website coding/240701_uniprot_membrane_proteins+capillary_expression.csv')
len(df)


# The betsholtz dataset columns were not normally distributed - I log transform the columns from the betsholtz dataset to make it normalized

# In[5]:


columns_to_log = ['bets_aEC', 'bets_capilEC', 'bets_EC3', 'bets_vEC', 'bets_EC2',
       'bets_EC1', 'bets_MG', 'bets_PC', 'bets_vSMC', 'bets_aSMC',
       'bets_aaSMC', 'bets_AC', 'bets_OL', 'bets_FB2', 'bets_FB1']

df[columns_to_log] = df[columns_to_log].apply(lambda x:np.log1p(x))


# In[6]:


file_path = r'/Users/shamika.bhandarkar/Desktop/Expression website coding/normalized RNA_seq.xlsx'  
df.to_excel(file_path, index=False)


# **Jeong**
# 
# CapA = capillary/pre-capillary arteriolar
# 
# CapV = capillary/post-capillary venular 
# 
# REV = reactive endothelial venules 

# In[7]:


df.columns


# Jeong = Adams, Kalucka = Carmeliet, Bjornholm = Vanlandewijck

# In[ ]:





# Now we want to make a graph of **all the other cell types** in the datasets

# In[8]:


df.columns


# We want to use the following columns - 
# kalucka interferon, kalucka choroid plexus, bets_MG, bets_PC, bets_vSMC, bets_aSMC, bets_aaSMC, bets_AC, bets_OL, bets_FB1, bets_FB2, bjornholm_Microglia,bjornholm_VSMC, bjornholm_Pericyte, bjornholm_Myeloid, bjornholm_Astrocyte, bjornholm_Fibroblast 

# In[ ]:





# **Single gene**

# 
# 

# In[11]:


import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from scipy.interpolate import make_interp_spline
import numpy as np

# Define colors and labels
author_colors = {
    'bets': '#E15759',
    'bjornholm': '#76B7B2'
}

legend_labels = {
    'bets': 'Betsholtz',
    'bjornholm': 'Vanlandewijck'
}

label_lookup = {
    'MG': 'Microglia', 'PC': 'Pericyte', 'vSMC': 'Venous SMC', 'aSMC': 'Arterial SMC',
    'aaSMC': 'Arterial/Arteriolar SMC', 'AC': 'Astrocyte', 'OL': 'Oligodendrocyte',
    'FB1': 'Fibroblast 1', 'FB2': 'Fibroblast 2'
}

# --- Custom unified order ---
combined_order = [
    'bets_MG', 'bjornholm_Microglia', 'bets_PC', 'bjornholm_Pericyte', 'bets_AC', 'bjornholm_Astrocyte',
    'bets_FB2', 'bets_FB1', 'bjornholm_Fibroblast', 'bets_vSMC', 'bjornholm_VSMC',
    'bets_aSMC', 'bets_aaSMC', 'bets_OL', 'bjornholm_Myeloid'
]

# --- Input ---
gene_name = 'Prom1'  # change this to the gene name you want
gene_rows = df[df['Gene Names'].str.contains(gene_name, case=False, na=False)]

if gene_rows.empty:
    raise ValueError(f"{gene_name} not found in DataFrame")

protein_index = gene_rows.index.tolist()

ref_proteins = {
    'ITM2B_MOUSE': {'color': 'dimgrey', 'label': 'Max Expression (ITM2B)', 'style': 'x'},
    'CC141_MOUSE': {'color': 'black', 'label': 'Median Expression (Ccdc141)', 'style': 'x'}
}
ref_indices = {key: df.index[df['Entry Name'] == key].tolist()[0] for key in ref_proteins}

# --- Determine max value for scaling ---
max_value = max(
    df.loc[protein_index[0], combined_order].max(),
    *(df.loc[idx, combined_order].max() for idx in ref_indices.values()),
    1
)

# --- Plotting setup ---
fig, ax = plt.subplots(figsize=(14, 10))
bar_width = 0.1
base_x = -0.4
patches_list = []
x_positions = []

# --- Reference markers and curves ---
for ref_name, ref_data in ref_proteins.items():
    ref_idx = ref_indices[ref_name]
    cluster_xs = []
    cluster_ys = []
    for i, subgroup in enumerate(combined_order):
        x = base_x + i * bar_width + bar_width / 2
        y = df.loc[ref_idx, subgroup]
        cluster_xs.append(x)
        cluster_ys.append(y)
        ax.scatter(x, y, marker=ref_data['style'], color=ref_data['color'],
                   s=80, label=ref_data['label'] if i == 0 else "", zorder=5)
    if len(cluster_xs) >= 2:
        xnew = np.linspace(min(cluster_xs), max(cluster_xs), 300)
        ynew = make_interp_spline(cluster_xs, cluster_ys, k=2)(xnew)
        ax.plot(xnew, ynew, color=ref_data['color'], linestyle='--', alpha=0.8, zorder=4)

# --- Main protein bars ---
for i, subgroup in enumerate(combined_order):
    prefix = subgroup.split('_')[0]
    color = author_colors.get(prefix, 'black')

    x = base_x + i * bar_width
    y = df.loc[protein_index[0], subgroup]

    rect = patches.Rectangle((x, 0), bar_width, y, facecolor=color, edgecolor='black')
    patches_list.append(rect)
    x_positions.append(x + bar_width / 2)

    # Labeling
    raw_label = subgroup.split('_')[-1]
    label = label_lookup.get(raw_label, raw_label)
    ax.text(x + bar_width / 1.5, -0.1, label, ha='right', va='top',
            fontname='Arial', fontsize=10, rotation=30)

# --- Add all bars ---
pc = PatchCollection(patches_list, match_original=True)
ax.add_collection(pc)

# --- Formatting ---
ax.set_xticks([])  # Remove Betsholtz / Vanlandewijck labels
ax.set_ylabel('Normalized Average Abundance', fontname='Arial', fontsize=15, weight='bold')
ax.set_title(f'{gene_name} Abundance', fontname='Arial', fontsize=15, weight='bold')

# --- Legend ---
bar_handles = [
    patches.Patch(color=author_colors[key], label=legend_labels[key]) for key in author_colors
]
ref_handles = [
    Line2D([0], [0], color='dimgrey', marker='x', linestyle='None', label='High Expression (Itm2b)'),
    Line2D([0], [0], color='black', marker='x', linestyle='None', label='Moderate Expression (Ccdc141)')
]
ax.legend(handles=bar_handles + ref_handles, loc='upper right')

# --- Axis limits ---
ax.set_xlim(-0.5, base_x + len(combined_order) * bar_width + 0.1)
ax.set_ylim(0, max_value + 1)

plt.tight_layout()

plt.savefig('/Users/shamika.bhandarkar/Desktop/Expression website coding/Prom1_nonEC.png', dpi=300, bbox_inches='tight')
plt.show()


# Finding the **median** expressing protein for the first graph (with EC cells) 

# In[130]:


columns_of_interest = [
    'jeong_CapA', 'jeong_CapV', 'kalucka_capillary venous', 'kalucka_capillary', 
    'kalucka_capillary arterial', 'bets_capilEC', 'bjornholm_CapEC', 'jeong_Arterial', 
    'kalucka_artery shear stress', 'kalucka_artery', 'kalucka_large artery', 'bets_aEC',
    'bjornholm_ArtlEC', 'bjornholm_ArtEC', 'jeong_Venous', 'jeong_REV', 
    'kalucka_large vein', 'bets_vEC', 'bjornholm_VenlEC', 'bjornholm_VenEC'
]

# Replace 0s with NaN to ignore them in calculations
df_no_zeros = df[columns_of_interest].replace(0, np.nan)

# Filter rows where at least one bets_ column is non-zero (i.e., not NaN after replacement)
bets_columns = [col for col in columns_of_interest if col.startswith('bets_')]
has_bets_expression = df[bets_columns].replace(0, np.nan).notna().any(axis=1)

# Compute mean expression across columns_of_interest (ignoring zeros/NaNs)
df['total_expression'] = df_no_zeros.mean(axis=1, skipna=True)
# Apply filter
filtered_df = df[has_bets_expression].copy()

# Sort by expression
sorted_df = filtered_df.sort_values('total_expression').reset_index()

# Find the median index
median_index = len(sorted_df) // 2
median_protein_name = sorted_df.loc[median_index, 'Entry Name']

print("Median expressing protein (with non-zero expression in bets columns):", median_protein_name)


# Finding the **median** expressing protein for the second graph (non-EC cell types)

# In[73]:


columns_of_interest = ['bets_MG', 'bets_PC', 'bets_vSMC', 'bets_aSMC',
        'bets_aaSMC', 'bets_AC', 'bets_OL', 'bets_FB2', 'bets_FB1', 'bjornholm_Microglia', 'bjornholm_VSMC', 'bjornholm_Pericyte',
        'bjornholm_Myeloid', 'bjornholm_Astrocyte', 'bjornholm_Fibroblast']

# Replace 0s with NaN to ignore them in calculations
df_no_zeros = df[columns_of_interest].replace(0, np.nan)

# Filter rows where at least one bets_ column is non-zero (i.e., not NaN after replacement)
bets_columns = [col for col in columns_of_interest if col.startswith('bets_')]
has_bets_expression = df[bets_columns].replace(0, np.nan).notna().any(axis=1)

# Compute mean expression across columns_of_interest (ignoring zeros/NaNs)
df['total_expression'] = df_no_zeros.mean(axis=1, skipna=True)
# Apply filter
filtered_df = df[has_bets_expression].copy()

# Sort by expression
sorted_df = filtered_df.sort_values('total_expression').reset_index()

# Find the median index
median_index = len(sorted_df) // 2
median_protein_name = sorted_df.loc[median_index, 'Entry Name']

print("Median expressing protein (with non-zero expression in bets columns):", median_protein_name)


# Finding the **highest expressing** protein for the EC cells graph

# In[123]:


columns_of_interest = [
    'jeong_CapA', 'jeong_CapV', 'kalucka_capillary venous', 'kalucka_capillary', 
    'kalucka_capillary arterial', 'bets_capilEC', 'bjornholm_CapEC', 'jeong_Arterial', 
    'kalucka_artery shear stress', 'kalucka_artery', 'kalucka_large artery', 'bets_aEC',
    'bjornholm_ArtlEC', 'bjornholm_ArtEC', 'jeong_Venous', 'jeong_REV', 
    'kalucka_large vein', 'bets_vEC', 'bjornholm_VenlEC', 'bjornholm_VenEC'
]

# Replace 0s with NaN to ignore them in calculations
df_no_zeros = df[columns_of_interest].replace(0, np.nan)

# Filter rows where at least one bets_ column is non-zero (i.e., not NaN after replacement)
bets_columns = [col for col in columns_of_interest if col.startswith('bets_')]
has_bets_expression = df[bets_columns].replace(0, np.nan).notna().any(axis=1)

# Compute mean expression across columns_of_interest (ignoring zeros/NaNs)
df['total_expression'] = df_no_zeros.mean(axis=1, skipna=True)

# Apply filter
filtered_df = df[has_bets_expression]

# Find the protein with the highest expression in the filtered DataFrame
highest_idx = filtered_df['total_expression'].idxmax()
highest_protein_name = df.loc[highest_idx, 'Entry Name']

print("Highest expressing protein (with non-zero expression in bets columns):", highest_protein_name)


# Finding the **highest expressing** protein for the non-ECs graph

# In[74]:


columns_of_interest = ['bets_MG', 'bets_PC', 'bets_vSMC', 'bets_aSMC',
        'bets_aaSMC', 'bets_AC', 'bets_OL', 'bets_FB2', 'bets_FB1', 'bjornholm_Microglia', 'bjornholm_VSMC', 'bjornholm_Pericyte',
        'bjornholm_Myeloid', 'bjornholm_Astrocyte', 'bjornholm_Fibroblast' ]

# Replace 0s with NaN to ignore them in calculations
df_no_zeros = df[columns_of_interest].replace(0, np.nan)

# Filter rows where at least one bets_ column is non-zero (i.e., not NaN after replacement)
bets_columns = [col for col in columns_of_interest if col.startswith('bets_')]
has_bets_expression = df[bets_columns].replace(0, np.nan).notna().any(axis=1)

# Compute mean expression across columns_of_interest (ignoring zeros/NaNs)
df['total_expression'] = df_no_zeros.mean(axis=1, skipna=True)

# Apply filter
filtered_df = df[has_bets_expression]

# Find the protein with the highest expression in the filtered DataFrame
highest_idx = filtered_df['total_expression'].idxmax()
highest_protein_name = df.loc[highest_idx, 'Entry Name']

print("Highest expressing protein (with non-zero expression in bets columns):", highest_protein_name)


# In[54]:


import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
from scipy.interpolate import make_interp_spline

# Define clusters and colors
cav_clusters = {
    'Capillary': [
        'jeong_CapA', 'jeong_CapV', 'kalucka_capillary venous', 'kalucka_capillary', 
        'kalucka_capillary arterial', 'bets_capilEC', 'bjornholm_CapEC'
    ],
    'Artery': [
        'jeong_Arterial', 'kalucka_artery shear stress', 'kalucka_artery', 
        'kalucka_large artery', 'bets_aEC','bjornholm_ArtlEC', 'bjornholm_ArtEC'
    ],
    'Vein': [
        'jeong_Venous', 'jeong_REV', 'kalucka_large vein', 'bets_vEC', 
        'bjornholm_VenlEC', 'bjornholm_VenEC'
    ]
}

cav_colors = {
    'jeong': '#4E79A7',
    'kalucka': '#F28E2B',
    'bets': '#E15759',
    'bjornholm': '#76B7B2'
}

# Protein of interest
protein_name = 'ITM2A_MOUSE'
reference_protein_1 = 'BASI_MOUSE'  # max
reference_protein_2 = 'CC141_MOUSE'  # median

# Get indices
protein_index = df.index[df['Entry Name'] == protein_name].tolist()
ref1_index = df.index[df['Entry Name'] == reference_protein_1].tolist()
ref2_index = df.index[df['Entry Name'] == reference_protein_2].tolist()

# Ensure all are present
if not protein_index or not ref1_index or not ref2_index:
    raise ValueError("One or more protein names not found in the DataFrame.")

# Get max y for scaling
all_columns = cav_clusters['Capillary'] + cav_clusters['Artery'] + cav_clusters['Vein']
max_value = df.loc[[protein_index[0], ref1_index[0], ref2_index[0]], all_columns].max().max()
max_value = max(max_value, 1)

# Plot setup
fig, ax = plt.subplots(figsize=(14, 10))
bar_width = 0.1
group_width = bar_width * max(len(v) for v in cav_clusters.values())
base_x = -0.4
x_positions = []
cluster_positions = []

# Loop through clusters
for cluster_name, subgroup_list in cav_clusters.items():
    cluster_x = base_x + (len(subgroup_list) * bar_width) / 2
    cluster_positions.append(cluster_x)
    for j, subgroup in enumerate(subgroup_list):
        x = base_x + j * bar_width
        y = df.loc[protein_index[0], subgroup]
        prefix = subgroup.split('_')[0]
        color = cav_colors.get(prefix, 'black')
        
        rect = patches.Rectangle((x, 0), bar_width, y, facecolor=color, edgecolor='black')
        ax.add_patch(rect)
        ax.text(x + bar_width / 1.5, -0.1, subgroup.split('_')[-1], ha='right', va='top', fontname='Arial', fontsize=10, rotation=30)
        x_positions.append(x + bar_width / 2)
    base_x += group_width + bar_width

# Plot reference protein curves
def plot_reference_curve(ref_index, color, label):
    xs = []
    ys = []
    xpos = -0.4
    for cluster_name, subgroup_list in cav_clusters.items():
        cluster_xs = []
        cluster_ys = []
        for j, subgroup in enumerate(subgroup_list):
            x = xpos + j * bar_width + bar_width / 2
            y = df.loc[ref_index[0], subgroup]
            if not np.isnan(y):
                cluster_xs.append(x)
                cluster_ys.append(y)
                ax.scatter(x, y, color=color, marker='x', s=80, label=label if (x == cluster_xs[0]) else "")
        if len(cluster_xs) >= 2:
            xnew = np.linspace(min(cluster_xs), max(cluster_xs), 300)
            ynew = make_interp_spline(cluster_xs, cluster_ys, k=2)(xnew)
            ax.plot(xnew, ynew, color=color, linestyle='--', alpha=0.8)
        xpos += group_width + bar_width

# Add the two reference lines
plot_reference_curve(ref1_index, 'dimgrey', 'Max Expression (BASI)')
plot_reference_curve(ref2_index, 'black', 'Median Expression (CC141)')

# Final formatting
ax.set_xticks(cluster_positions)
ax.set_xticklabels(cav_clusters.keys(), fontname='Arial', fontsize=14, weight='bold')
ax.xaxis.set_tick_params(pad=50)
ax.set_ylabel('Normalized Average Abundance', fontname='Arial', fontsize=15, weight='bold')
ax.set_title(f'{protein_name} Abundance with Reference Curves', fontname='Arial', fontsize=15, weight='bold')

# Create handles for main protein bar colors
# Mapping for display names in legend
legend_labels = {
    'jeong': 'Adams',
    'kalucka': 'Carmeliet',
    'bjornholm': 'Vanlandewijck',
    'bets': 'Betsholtz'
}

bar_handles = [
    patches.Patch(color=cav_colors[key], label=legend_labels[key]) for key in cav_colors
]

# Create handles for reference markers
ref_handles = [
    plt.Line2D([0], [0], color='dimgrey', marker='x', linestyle='None', label='High Expression (Basigin)'),
    plt.Line2D([0], [0], color='black', marker='x', linestyle='None', label='Moderate Expression (Ccdc141)')
]

# Combine and add to legend
all_handles = bar_handles + ref_handles
ax.legend(handles=all_handles, loc='upper right')

ax.set_xlim(-0.5, base_x)
ax.set_ylim(0, max_value + 1)

plt.tight_layout()
plt.show()


# **Loop** through all the proteins in the df

# In[ ]:


import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
from scipy.interpolate import make_interp_spline
import re  # to sanitize filenames

# Ensure save_path exists
save_path = "/Volumes/Neurobio/GU LAB/PERSONAL FILES/Shamika/RNAseq website/EC plots"  # Update this path
os.makedirs(save_path, exist_ok=True)

# Define clusters and colors
cav_clusters = {
    'Capillary': [
        'jeong_CapA', 'jeong_CapV', 'kalucka_capillary venous', 'kalucka_capillary', 
        'kalucka_capillary arterial', 'bets_capilEC', 'bjornholm_CapEC'
    ],
    'Artery': [
        'jeong_Arterial', 'kalucka_artery shear stress', 'kalucka_artery', 
        'kalucka_large artery', 'bets_aEC','bjornholm_ArtlEC', 'bjornholm_ArtEC'
    ],
    'Vein': [
        'jeong_Venous', 'jeong_REV', 'kalucka_large vein', 'bets_vEC', 
        'bjornholm_VenlEC', 'bjornholm_VenEC'
    ]
}

cav_colors = {
    'jeong': '#4E79A7',
    'kalucka': '#F28E2B',
    'bets': '#E15759',
    'bjornholm': '#76B7B2'
}

legend_labels = {
    'jeong': 'Adams',
    'kalucka': 'Carmeliet',
    'bjornholm': 'Vanlandewijck',
    'bets': 'Betsholtz'
}

reference_protein_1 = 'BASI_MOUSE'
reference_protein_2 = 'CC141_MOUSE'

ref1_index = df.index[df['Entry Name'] == reference_protein_1].tolist()
ref2_index = df.index[df['Entry Name'] == reference_protein_2].tolist()

if not ref1_index or not ref2_index:
    raise ValueError("Reference proteins not found in the DataFrame.")

# Loop through entire DataFrame
for idx, row in df.iterrows():
    if idx % 1000 == 0 and idx != 0:
        print(f"{idx} proteins processed and saved...")
    
    # Get gene name (fallback to Entry Name if missing)
    gene_field = row.get('Gene Names', '')
    gene_name = gene_field.split()[0] if isinstance(gene_field, str) and gene_field.strip() else row['Entry Name']
    
    # Sanitize gene_name for filename
    gene_name_clean = re.sub(r'[\\/*?:"<>|]', "_", gene_name)
    
    protein_index = [idx]
    all_columns = cav_clusters['Capillary'] + cav_clusters['Artery'] + cav_clusters['Vein']
    max_value = df.loc[[protein_index[0], ref1_index[0], ref2_index[0]], all_columns].max().max()
    max_value = max(max_value, 1)

    # Plot setup
    fig, ax = plt.subplots(figsize=(14, 10))
    bar_width = 0.1
    group_width = bar_width * max(len(v) for v in cav_clusters.values())
    base_x = -0.4
    x_positions = []
    cluster_positions = []

    for cluster_name, subgroup_list in cav_clusters.items():
        cluster_x = base_x + (len(subgroup_list) * bar_width) / 2
        cluster_positions.append(cluster_x)
        for j, subgroup in enumerate(subgroup_list):
            x = base_x + j * bar_width
            y = df.loc[protein_index[0], subgroup]
            prefix = subgroup.split('_')[0]
            color = cav_colors.get(prefix, 'black')

            rect = patches.Rectangle((x, 0), bar_width, y, facecolor=color, edgecolor='black')
            ax.add_patch(rect)
            ax.text(x + bar_width / 1.5, -0.1, subgroup.split('_')[-1], ha='right', va='top',
                    fontname='Arial', fontsize=10, rotation=30)
            x_positions.append(x + bar_width / 2)
        base_x += group_width + bar_width

    # Reference curves
    def plot_reference_curve(ref_index, color, label):
        xs, ys = [], []
        xpos = -0.4
        for cluster_name, subgroup_list in cav_clusters.items():
            cluster_xs, cluster_ys = [], []
            for j, subgroup in enumerate(subgroup_list):
                x = xpos + j * bar_width + bar_width / 2
                y = df.loc[ref_index[0], subgroup]
                if not np.isnan(y):
                    cluster_xs.append(x)
                    cluster_ys.append(y)
                    ax.scatter(x, y, color=color, marker='x', s=80, label=label if (x == cluster_xs[0]) else "")
            if len(cluster_xs) >= 2:
                xnew = np.linspace(min(cluster_xs), max(cluster_xs), 300)
                ynew = make_interp_spline(cluster_xs, cluster_ys, k=2)(xnew)
                ax.plot(xnew, ynew, color=color, linestyle='--', alpha=0.8)
            xpos += group_width + bar_width

    plot_reference_curve(ref1_index, 'dimgrey', 'High Expression (Basigin)')
    plot_reference_curve(ref2_index, 'black', 'Moderate Expression (Ccdc141)')

    ax.set_xticks(cluster_positions)
    ax.set_xticklabels(cav_clusters.keys(), fontname='Arial', fontsize=14, weight='bold')
    ax.xaxis.set_tick_params(pad=50)
    ax.set_ylabel('Normalized Average Abundance', fontname='Arial', fontsize=15, weight='bold')
    ax.set_title(f'{gene_name} Abundance', fontname='Arial', fontsize=15, weight='bold')

    bar_handles = [patches.Patch(color=cav_colors[key], label=legend_labels[key]) for key in cav_colors]
    ref_handles = [
        plt.Line2D([0], [0], color='dimgrey', marker='x', linestyle='None', label='High Expression (Basigin)'),
        plt.Line2D([0], [0], color='black', marker='x', linestyle='None', label='Moderate Expression (Ccdc141)')
    ]
    ax.legend(handles=bar_handles + ref_handles, loc='upper right')
    ax.set_xlim(-0.5, base_x)
    ax.set_ylim(0, max_value + 1)

    plt.tight_layout()
    file_name = os.path.join(save_path, f"{gene_name_clean}.png")
    plt.savefig(file_name, dpi=300, bbox_inches='tight')
    plt.close()


# Trying to **rank** the proteins in the df based on all the other proteins (i.e., remaining distribution) by summing all the clusters

# In[ ]:





# In[ ]:





# In[ ]:





# 

# In[ ]:





# In[ ]:





# 

# In[ ]:





# In[ ]:





# In[11]:





# In[ ]:





# In[ ]:





# 

# 

# In[ ]:





# 

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# 

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




