#!/usr/bin/env python
# coding: utf-8

# In[237]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ### Amount of bacterial genomes

# In[238]:


threshold_num_bac = 10.0 # a threshold on a number of bacterial genomes (try 1, 5 or 10)


# In[239]:


num_bac_genomes = pd.read_csv("average_day7.txt", sep=" ")
num_bac_genomes = num_bac_genomes[num_bac_genomes['count'] > threshold_num_bac]
num_bac_genomes.reset_index(drop=True, inplace=True)

# In[240]:


#num_bac_genomes['count'].hist(bins=50)
#plt.xlabel('Count of bacterial genomes')
#plt.title('Distribution of counts of bacterial genomes')
#plt.savefig('counts_distr.png', dpi=300)


# In[241]:


num_bac_genomes.shape


# In[242]:


#num_bac_genomes.head()


# In[243]:


multiplier_num_bac = 21 # multiplier for a total number of bacteria to be simulated 250 (300 K the best)
# 1 15 25 50 100 250 500
num_bac_genomes['count'] = [round(x) for x in num_bac_genomes['count'] * multiplier_num_bac]
total_num_bac = num_bac_genomes['count'].sum()


# In[244]:


total_num_bac


# In[245]:


#num_bac_genomes.head()


# In[152]:


num_bac_genomes.to_csv(f"NUM_BAC_MULTIPLIED_total_bac_{total_num_bac}_mult_{multiplier_num_bac}_thresh_{threshold_num_bac}.tsv", sep='\t', index=None)


# ### Functional distances

# In[246]:

"""
genome_distance_matrix = pd.read_csv("genome_functional.distances.txt", sep=" ")


# In[247]:


genome_distance_matrix


# In[248]:


genome_distance_matrix_selected = genome_distance_matrix.loc[list(num_bac_genomes['genome']), list(num_bac_genomes['genome'])]


# In[257]:


## Distribution of distances for initial genomes quantity (174)

upper_triangle_mask = np.triu(np.ones(genome_distance_matrix.shape), k=1).astype(bool)
result = genome_distance_matrix.where(upper_triangle_mask).stack()
result_df = result.reset_index()
result_df.columns = ['Row', 'Column', 'Value']
result_df['Value'].hist(bins=50)
plt.xlabel('Similarity')
plt.title('Distribution of similarities for 174 bacterial genomes')
plt.savefig('sim_all.png', dpi=300)


# In[258]:


result_df['Value'].median()


# In[254]:


## Distribution of distances for selected genomes quantity (48)

upper_triangle_mask = np.triu(np.ones(genome_distance_matrix_selected.shape), k=1).astype(bool)
result = genome_distance_matrix_selected.where(upper_triangle_mask).stack()
result_df = result.reset_index()
result_df.columns = ['Row', 'Column', 'Value']
result_df['Value'].hist(bins=50)
plt.xlabel('Similarity')
plt.title('Distribution of similarities for 34 bacterial genomes')
plt.savefig('sim_selected.png', dpi=300)
# I would expect 3 agregates of bacteria based on this plot


# In[256]:


result_df['Value'].median()


# In[158]:


genome_distance_matrix_selected


# In[159]:


particle_types = dict(zip(list(genome_distance_matrix_selected.index), list(range(3, genome_distance_matrix_selected.shape[0] + 3))))


# In[160]:


base_coef = 1 / result_df['Value'].min() / 1e6 # For the particles of identical types


# In[161]:


result_df['Value'].min()


# In[162]:


base_coef


# In[163]:


result_df['Row'] = list(map(lambda x: particle_types[x], result_df['Row']))
result_df['Column'] = list(map(lambda x: particle_types[x], result_df['Column']))
result_df['Value'] = list(map(lambda x: format(1 / x / 1e6, '.15f'), result_df['Value']))


# In[164]:


result_df


# In[165]:


import matplotlib.pyplot as plt


# In[166]:


plt.hist([float(x) for x in result_df['Value']], bins=100)


# ### Write coefficients

# In[167]:


# Bacteria-Villi coefficients (fixed)

energy_bac_vil = 10.2
size_bac_vil = 1.0
limit_bac_vil = 2.5

identical_pair_fmt = f"3*{genome_distance_matrix_selected.shape[0] + 2}"
print(f"pair_coeff 1*2 {identical_pair_fmt}\t{energy_bac_vil}\t{size_bac_vil}\t{limit_bac_vil}")


# In[168]:


# Identical coefficients

bacteria_size = 0.12 # At the moment all bacteria have the same size and mass
identical_limit = 5 # distance limit for identical particles

identical_pair_fmt = f"3*{genome_distance_matrix_selected.shape[0] + 2}"
print(f"pair_coeff {identical_pair_fmt} {identical_pair_fmt}\t{base_coef:.15f}\t{bacteria_size}\t{identical_limit}")


# In[169]:


# Pair coefficients

bacteria_size = 0.12 # At the moment all bacteria have the same size and mass
identical_limit = 2.5 # distance limit for identical particles

with open('PAIR_COEFS.txt', 'w') as f:
    for i, row in result_df.iterrows():
        f.write(f"pair_coeff {row['Row']} {row['Column']}\t{row['Value']}\t{bacteria_size}\t{identical_limit}\n")


# ### calculate bacteria volume

# In[190]:


import numpy as np


# In[214]:


r_bac = 0.06
height_bac_space = 20
num_bac = 298120 # 119247 298120

vol_space = 40 * 40 * height_bac_space
v_bac = (4 / 3) * np.pi * (r_bac ** 3)
v_bac_tot = v_bac * num_bac
v_bac_share = v_bac_tot / vol_space


# In[215]:


v_bac_share # 298120


# In[183]:


v_bac_share # 119247


# In[216]:


r_bac = 0.06
height_bac_space = 3
num_bac = 119247 # 119247 298120

vol_space = 40 * 40 * height_bac_space
v_bac = (4 / 3) * np.pi * (r_bac ** 3)
v_bac_tot = v_bac * num_bac
v_bac_share = v_bac_tot / vol_space


# In[206]:


v_bac_share


# In[217]:


v_bac_share


# In[218]:


3/0.06


# In[ ]:
"""



