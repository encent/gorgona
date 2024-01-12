import numpy as np

def compute_density(mass, diameter):
    return round(mass / ((4/3) * np.pi * (diameter / 2) ** 3), 10)  # ro = m / V
