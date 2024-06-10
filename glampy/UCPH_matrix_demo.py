#!/usr/bin/env python
# coding: utf-8

# In[1]:


from glampy.sampler import PoissonDiskSampler
from glampy.generator import ArchitectureGenerator

import subprocess
import matplotlib.pyplot as plt

fids_list = ["FINAL_CONFORM_Z_42.5_45.5_mult_1900_numbac_2.25M_bacsize_0.04_"]

for fid in fids_list:
    print(fid)

    # In[2]:


    #bac_mass = 0.00000644579089
    #bac_num = 1000

    WITH_VILLI = True
    if WITH_VILLI:
        num_shift = 3
    else:
        num_shift = 1

    mode = 'hybrid' # angle or full or sphere or hybrid
    units = 'lj'
    atom_types = 1 # == 1 -- all ONE size ; <> 1 -- different sizes
    localization = 'ileum'

    #seed = 34882 # seed
    seed = 42

    spacing = 1.5  # spacing number from 0 to 2 (see docs)

    x = (0, 40)  # (μm) -> 0.5 cm 0 40
    y = (0, 40)  # (μm) -> 0.5 cm 0 40
    z = (0, 40)  # (μm) -> 0.5 cm 0 40
    # z = None  # if None, then calculated based on the chain(s) length
    box_size = (x, y, z)

    if localization == 'ileum':
        villus_height = 40  # 1200 # μm
        villus_width = 4  # μm 1 or 4

        min_intervillus_width = 1.5  # μm 4 or 1
        # mean_intervillus_width = 42  # μm

        graft_size = 4  # 1 or 4

        crypt_bottom_height = 4  # 1 or 4
        crypt_bottom_width = 4 # 1 or 4

        crypt_top_height = 4  # 1 or 4
        crypt_top_width = 4  # 1 or 4

        # min_dist = crypt_bottom_width + min_intervillus_width  # 353 μm -- should be revised with UCPH
        min_dist = graft_size + min_intervillus_width
    elif localization == 'caecum':
        villus_height = 300  # μm
        villus_width = 100  # μm

        min_intervillus_width = 25  # μm
        # mean_intervillus_width = 42  # μm

        graft_size = 300  # 300

        crypt_bottom_height = 100  # 583 μm
        crypt_bottom_width = 100

        crypt_top_height = 200  # 200
        crypt_top_width = 200  # 200

        min_dist = crypt_bottom_width + min_intervillus_width  # 353 μm -- should be revised with UCPH

    graft_mass = 1000
    crypt_bottom_mass = 1000
    crypt_top_mass = 1000
    villus_mass = 1000

    n_anchors = 42

    bottom_padding = 1


    # ### Initialize architecture generator

    # In[3]:


    archgen = ArchitectureGenerator(box_size, seed, min_dist, n_anchors, units, bottom_padding=bottom_padding, mode=mode)


    # ### Sample villi positions

    # In[4]:


    if WITH_VILLI: 
        n_actual = archgen.generate_grafting_layer(n_anchors, 10**6)


        # In[5]:


        if n_actual != n_anchors:
            print(f'Warning: Grafting layer too dense. {n_actual} grafting points instead of {n_anchors}.')
        else:
            print(f'{n_actual} points')


        # In[6]:


        #plt.figure(figsize=(7,7))
        #plt.scatter(archgen.coordinates[:, 0], archgen.coordinates[:, 1])


        # ### Construct brushes

        # In[7]:


        # Construct grafting beads (fixed layer)

        archgen.build(spacing=spacing,
                        bead_type="graft",
                        bead_size=graft_size,
                        bead_mass=graft_mass,
                        chain_height=graft_size,
                        bond_style="fene",
                        pair_style="lj/cut")


        # In[8]:


        # Construct villi beads (moving layer)

        archgen.build(spacing=spacing,
                        bead_type="villi",
                        bead_size=villus_width,
                        bead_mass=villus_mass,
                        chain_height=villus_height,
                        bond_style="fene",
                        pair_style="lj/cut")


        # ### Sample bacteria positions

        # In[9]:


        box_size


    # #### Get bacteria data from the file:

    # In[10]:


    import pandas as pd



    # In[12]:


    # Change the filename according to different multipliers # 1 15 25 50 100 250 500
    bacteria_data_file = pd.read_csv(f"NUM_BAC_MULTIPLIED_total_bac_2265701_mult_1900_thresh_10.0.tsv", sep="\t")


    # In[13]:


    bacteria_data_file
    for i in range(bacteria_data_file.shape[0]):
        bacteria_data_file.loc[i, 'genome'] = i + num_shift


    # In[14]:


    bacteria_data_file.head()


    # In[15]:


    ####


    # In[ ]:


    bacteria_box_size = ((0, 40), (0, 40), (42.5, 45.5)) # 42.5 -- 45.5 -- 62.5
    # bacteria_box_size = ((0, 10), (0, 10), (0, 10)) # 42.5 -- 45.5 -- 62.5

    n_bacteria = bacteria_data_file['count'].sum()

    bsize = 0.04
    # bsize = 0.12
    # bmass = 0.001728
    bmass = 0.000064

    bacteria_data = [(idx, bsize, bmass, bcnt) for idx, bcnt in zip(bacteria_data_file['genome'], bacteria_data_file['count'])]
    bacteria_data_dict = dict(zip(list(bacteria_data_file['genome']), [row1 + [row2] for row1, row2 in zip([[bsize, bmass]] * bacteria_data_file.shape[0], list(bacteria_data_file['count']))]))
    min_bac_dist = 0.07
    # min_bac_dist = 0.20

    n_bacteria_actual = archgen.generate_bacteria(n_bacteria, min_bac_dist, bacteria_box_size)


    # In[17]:


    print(n_bacteria, n_bacteria_actual)


    # In[18]:


    #fig = plt.figure(figsize=(7,7))
    #ax = fig.add_subplot(projection='3d') 
    #ax.scatter(archgen.coordinates_bacteria[:, 0], archgen.coordinates_bacteria[:, 1], archgen.coordinates_bacteria[:, 2], s=min_bac_dist)
    #plt.show()


    # ### Construct bacteria

    # In[19]:


    archgen.build(spacing=spacing,
                    bead_type="bacteria",
                    bead_size=villus_width,
                    bead_mass=villus_mass,
                    chain_height=villus_height,
                    bond_style="fene",
                    pair_style="lj/cut", 
                    bacteria_data=bacteria_data, 
                    bacteria_data_dict=bacteria_data_dict)

    print(archgen.atoms.atom_type.value_counts())


    # ### Write architecture in the file

    # In[20]:


    from os.path import isdir
    from os import makedirs


    # In[21]:


    if not isdir(fid):
        makedirs(fid)


    # In[22]:


    archgen.write(f"{fid}/merged.pos")

