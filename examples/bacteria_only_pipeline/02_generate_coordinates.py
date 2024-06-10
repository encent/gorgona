#!/usr/bin/env python

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Import BacteriaVilliCoordinatesGenerator class from generator module
from glampy.generator import BacteriaVilliCoordinatesGenerator

# Define and initialize parameters
params = {
	"input_file_path": "../../data/01_bacteria_amount/bac_amount_518967_mult_435.2017151420611_thresh_10.tsv",
	"generation_type": "bacteria",
	"mode": "hybrid",
	"bottom_padding": 0,
	"xlim": (0, 10),
	"ylim": (0, 10),
	"zlim_villi": (0, 10),
	"zlim_bacteria": (0, 10),
	"bsize": 0.04, # 0.04 or 0.12 for test
	"bmass": 0.000064, # 0.000064 or 0.001728 for test
    # Add more parameters if needed
}

# Create an instance of BacteriaVilliCoordinatesGenerator
bacteria_villi_generator = BacteriaVilliCoordinatesGenerator(**params)

# Call the generate_coordinates method
bacteria_villi_generator.generate_coordinates(output_dir_path="../../data/02_coordinates/", plot_figs=False)
