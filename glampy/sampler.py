from itertools import product
from typing import Optional

import numpy as np
import pandas as pd


class PoissonDiskSampler:
	"""
	Generate a Poisson-disk 2D or 3D point set using cell list accelerated dart throwing.
	"""

	def __init__(self, 
					seed: Optional[int],
					n_dim: int = 2) -> None:
		"""
		:param seed:		Seed used in random sampling. If None, a random seed will be used.
		:param n_dim:		Number of dimensions (2 or 3) of the points to be sampled. Default is 2.
		"""
		self.rng = np.random.RandomState(seed=seed)
		self.n_dim = n_dim

		# Create a list of 20 (2D) or 116 (3D) neighbors (5x5 excluding corners and center)
		self.neighbor_matrix = []
		if self.n_dim == 2:
			for i in [-2, -1, 0, 1, 2]:
				for j in [-2, -1, 0, 1, 2]:
					if not ((abs(i) == 2 and abs(j) == 2) or (i == 0 and j == 0)):
						self.neighbor_matrix.append((i, j))
		elif self.n_dim == 3:
			for i in [-2, -1, 0, 1, 2]:
				for j in [-2, -1, 0, 1, 2]:
					for k in [-2, -1, 0, 1, 2]:
						if not ((abs(i) == 2 and abs(j) == 2 and abs(k) == 2) or (i == 0 and j == 0 and k == 0)):
							self.neighbor_matrix.append((i, j, k))

	def generate(self,
					n: int,
					min_dist: float, 
					box_size: tuple[tuple[float, float], tuple[float, float], Optional[tuple[float, float]]],
					max_iter: int = 1000000000) -> np.ndarray:
		"""
		Generate a random sample according to the PoissonDiskSampler class instance.
		:param n:				Number of points to be sampled.
		:param min_dist:		Minimum distance between points (same to be used in architecture builder).
		:param box_size:		2D tuple corresponding to the 2D or 3D domain size.
		:param max_iter:		Iteration limit.
		:return:				np.ndarray of shape (n, 2) or (n, 3) with all point coordinates.
		"""
		def check_overlap(cell: np.ndarray, point: np.ndarray) -> bool:
			"""
			Check for overlap in the neighborhood of the current cell.
			:param cell:		np.ndarray of shape (1, 2) or (1, 3) containing cell indices.
			:param point:		np.ndarray of shape (1, 2) or (1, 3) containing the point coordinates.
			:return:			True if point overlaps, False if not.
			"""
			if self.n_dim == 2:
				for i, j in self.neighbor_matrix:
					x = cell[0] + i
					y = cell[1] + j
					if 0 <= x <= (x_max - 1) and 0 <= y <= (y_max - 1):
						coord = grid_lookup[x, y]
						if coord is not None:
							dist = (point - coord)**2
							if dist[0] + dist[1] < min_dist ** 2:
								return True
			elif self.n_dim == 3:
				for i, j, k in self.neighbor_matrix:
					x = cell[0] + i
					y = cell[1] + j
					z = cell[2] + k
					if 0 <= x <= (x_max - 1) and 0 <= y <= (y_max - 1) and 0 <= z <= (z_max - 1):
						coord = grid_lookup[x, y, z]
						if coord is not None:
							dist = (point - coord)**2
							if dist[0] + dist[1] + dist[2] < min_dist ** 2:
								return True

			return False

		def check_oob(point: np.ndarray) -> bool:
			"""
			Check if the point is out-of-bounds, which can happen because the domain the algorithm runs on is
			slightly larger than the requested domain because of rasterising errors.
			:param point:		np.ndarray of shape (1, 2) or (1, 3) containing the point coordinates.
			:return:			True if point is out-of-bounds, False if not.
			"""
			if self.n_dim == 2:
				return (point > np.array((box_size[0][1], box_size[1][1]))).any() or (point < np.array((box_size[0][0], box_size[1][0]))).any()
			elif self.n_dim == 3:
				return (point > np.array((box_size[0][1], box_size[1][1], box_size[2][1])) - min_dist).any() or (point < np.array((box_size[0][0], box_size[1][0], box_size[2][0])) + min_dist).any()

		# Create a grid of square background cells. The cells should be as large as possible but while still being fully
		# covered by the size of a point within it.
		cell_size = min_dist / np.sqrt(self.n_dim)

		# Round desired domain size up to nearest multiple of cell size
		x_max = int(np.ceil((box_size[0][1] - box_size[0][0]) / cell_size))
		y_max = int(np.ceil((box_size[1][1] - box_size[1][0]) / cell_size))
		if self.n_dim == 3:
			z_max = int(np.ceil((box_size[2][1] - box_size[2][0]) / cell_size))

		# List of indices of cells. Coordinate of bottom left corner = (i*cell_size, j*cell_size)
		if self.n_dim == 2:
			active_cells = list(product(range(x_max), range(y_max)))
		elif self.n_dim == 3:
			active_cells = list(product(range(x_max), range(y_max), range(z_max)))

		# x*y(*z) empty list
		if self.n_dim == 2:
			grid_lookup = np.full((x_max, y_max), None, dtype=object)
		elif self.n_dim == 3:
			grid_lookup = np.full((x_max, y_max, z_max), None, dtype=object)

		coordinates = []
		for _ in range(0, max_iter):
			# Choose random active cell
			if len(active_cells) == 0:
				break
			cell_id = self.rng.randint(low=0, high=len(active_cells))

			# Throw a dart
			rnd = self.rng.random_sample(self.n_dim)
			x = active_cells[cell_id][0] + rnd[0]
			y = active_cells[cell_id][1] + rnd[1]
			if self.n_dim == 3:
				z = active_cells[cell_id][2] + rnd[2]

			if self.n_dim == 2:
				point = np.array([x, y])*cell_size + np.array([box_size[0][0], box_size[1][0]])
			elif self.n_dim == 3:
				point = np.array([x, y, z])*cell_size + np.array([box_size[0][0], box_size[1][0], box_size[2][0]])

			# Check if point overlaps in neighbouring cells
			if not check_overlap(active_cells[cell_id], point) and not check_oob(point):
				coordinates.append(point)
				if self.n_dim == 2:
					grid_lookup[active_cells[cell_id][0], active_cells[cell_id][1]] = point
				elif self.n_dim == 3:
					grid_lookup[active_cells[cell_id][0], active_cells[cell_id][1], active_cells[cell_id][2]] = point
				del active_cells[cell_id]

			if len(coordinates) >= n:
				break

		return np.array(coordinates)
