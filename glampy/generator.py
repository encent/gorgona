import gzip
from abc import ABC, abstractmethod
from enum import Enum
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import random
import itertools

from os import makedirs
from os.path import isdir

from glampy.sampler import PoissonDiskSampler
from glampy.utils import compute_density


class BrushGenerator(ABC):
	"""
	Generate a LAMMPS data file containing a coarse-grained polymer brush grafted to a planar wall in a rectangular box.
	See https://lammps.sandia.gov/doc/read_data.html
	"""

	AtomTypes: Enum = Enum('AtomTypes', [])
	BondTypes: Enum = Enum('BondTypes', [])
	AngleTypes: Enum = Enum('AngleTypes', [])
	DihedralTypes: Enum = Enum('DihedralTypes', [])

	# in case of atom style = angle or full:
	masses: dict = {}

	pair_ij_coeffs: dict = {}
	bond_coeffs: dict = {}
	angle_coeffs: dict = {}
	dihedral_coeffs: dict = {}

	styles: dict[str] = {
		'pair': '',
		'bond': '',
		'angle': '',
		'dihedral': '',
	}

	def __init__(self,
				 box_size: tuple[tuple[float, float], tuple[float, float], Optional[tuple[float, float]]],
				 rng_seed: Optional[int],
				 min_dist: float,
				 units: str,
				 bottom_padding: float = 0,
				 mode: str = "sphere"):
		"""
		:param box_size:		3-tuple of floats describing the dimensions of the rectangular box. If the third (z)
								value is None, it will be automatically sized to contain the longest chain.
		:param rng_seed:		Seed used to initialize the PRNG. May be None, in which case a random seed will be used.
		:param min_dist:		Minimum distance between two centers of atoms: used as minimum distance for the
								Poisson-disk point set generator.
		:param bottom_padding:	Distance between the bottom edge of the box and the grafting layer. Must be positive.
		"""
		self.box_size = [list(x) for x in box_size]
		self.rng_seed = rng_seed
		self.min_dist = min_dist
		self.units = units
		self.bottom_padding = bottom_padding
		self.mode = mode

		# Grafted beads coordinates storage
		self.coordinates: np.ndarray = np.array([])
		self.coordinates_bacteria: np.ndarray = np.array([])

		# User parameters defined in a build() method
		self.spacing = 0.0

		self.density = 0.0
		self.bead_mass = 0.0
		self.bead_size = 0.0
		self.chain_height = 0.0
		self.n_beads = 0
		self.bond_style = ""
		self.pair_style = ""
		self.bead_type = ""
		self.beads_counter = 0
		self.current_height = 0
		self.current_beads = 0
		self.ground_layer = 0
		self.prev_layer_n_beads = 0
		self.multiplier = 1
		self.z_max = 0.0
		self.bacteria_data = []

		# Non-final lists
		self._atoms_list = []
		self._bonds_list = []
		self._angles_list = []
		self._dihedrals_list = []

		# Final dataframes
		self.atoms: pd.DataFrame = pd.DataFrame()
		self.bonds: pd.DataFrame = pd.DataFrame()
		self.angles: pd.DataFrame = pd.DataFrame()
		self.dihedrals: pd.DataFrame = pd.DataFrame()

		self.bac_box_size = ()

		if self.box_size[2]:
			self.not_grow_z = True
		else:
			self.not_grow_z = False


	@abstractmethod
	def _build_bead(self, mol_id: int, graft_coord: np.ndarray, bead_id: int) -> float:
		"""
		Adds a bead to the instance's atom/bond/angle/dihedral lists.
		Override this and implement according to the polymer model used. Note that LAMMPS expects ids to be 1-indexed.
		:param mol_id:      0-indexed molecule (chain) id
		:param graft_coord: 2-element ndarray containing xy coordinate of the grafting point
		:param bead_id:     0-indexed bead id
		:return Maximum z height
		"""
		pass

	@abstractmethod
	def _build_bacterium(self, mol_id: int, graft_coord: np.ndarray) -> None:
		"""
		Adds a bead to the instance's atom/bond/angle/dihedral lists.
		Override this and implement according to the polymer model used. Note that LAMMPS expects ids to be 1-indexed.
		:param mol_id:      0-indexed molecule (chain) id
		:param graft_coord: 2-element ndarray containing xy coordinate of the grafting point
		:param bead_id:     0-indexed bead id
		:return Maximum z height
		"""
		pass

	def generate_grafting_layer(self, n_chains: int, max_overlap_iter: int = 10**3) -> int:
		"""
		Generate coordinates of the grafting layer using a Poisson-disk point set generator.
		:param n_chains:         Number of grafting points (chains).
		:param max_overlap_iter: Iteration limit for the Poisson-disk point set generator.
		:return Number of grafting points generated.
		"""
		# Generate grafting point coordinates
		pdg = PoissonDiskSampler(self.rng_seed, n_dim=2)
		self.coordinates = pdg.generate(n_chains, self.min_dist, self.box_size, max_overlap_iter)

		return len(self.coordinates)
	
	def generate_bacteria(self, n_bacteria: int, min_bac_dist: float, bac_box_size: tuple) -> int:
		"""
		Generate coordinates of the grafting layer using a Poisson-disk point set generator.
		:param n_chains:         Number of grafting points (chains).
		:param max_overlap_iter: Iteration limit for the Poisson-disk point set generator.
		:return Number of grafting points generated.
		"""
		# Generate grafting point coordinates
		pdg = PoissonDiskSampler(self.rng_seed, n_dim=3)
		self.coordinates_bacteria = pdg.generate(n_bacteria, min_bac_dist, bac_box_size)
		self.bac_box_size = bac_box_size

		return len(self.coordinates_bacteria)

	def _build_chain(self) -> float:
		"""
		Create chains by looping over chains and beads in each chain, calling _build_bead() for each bead.
		:return Maximum z height
		"""
		if self.units == 'micro' or self.units == 'lj':
			self.multiplier = 0
			if self.current_height != 0:
				self.current_height += self.bead_size / 2
		z_max = 0
		z = 0
		# Loop over chains
		for mol_id, i in enumerate(self.coordinates):
			# Loop over successive beads in chain
			for j in range(0, self.n_beads):
				z = self._build_bead(mol_id, i, j)
				if z > z_max:
					z_max = z

		if self.units == 'micro' or self.units == 'lj':
			self.current_beads = int(z_max // (self.bead_size - self.spacing * self.bead_size / 2))
			self.current_height = z_max - self.bead_size / 2
			self.ground_layer += 1
			self.prev_layer_n_beads = self.n_beads

		return z_max

	def _build_bacteria(self, bacteria_data, bacteria_data_dict) -> None:
		# rng = np.random.RandomState(seed=self.rng_seed)
		bacteria_types = list(zip(*bacteria_data))[0]
		bacteria_diameters = list(zip(*bacteria_data))[1]
		bacteria_masses = list(zip(*bacteria_data))[2]
		bacteria_counts = list(zip(*bacteria_data))[3]
		types = list(itertools.chain.from_iterable(itertools.repeat(elem, count) for elem, count in zip(list(bacteria_types), list(bacteria_counts))))
		random.shuffle(types)
		num_bacteria = len(self.coordinates_bacteria)
		# types = rng.choice(bacteria_types, num_bacteria)
		villi_number = self.atoms.shape[0]
		atom_ids = list(map(lambda x: [x], list(range(villi_number, villi_number + len(types)))))
		bac_info_1 = list(map(lambda x: [x], types))
		# coordinates x y z
		bac_info_2 = list(map(lambda x: [bacteria_data_dict[x][0], compute_density(bacteria_data_dict[x][1], bacteria_data_dict[x][0])], types))
		atoms_list = np.concatenate((np.asarray(bac_info_1), self.coordinates_bacteria, np.asarray(atom_ids), np.asarray(bac_info_2)), axis=1)
		atoms = pd.DataFrame(atoms_list, columns=['atom_type', 'x', 'y', 'z', 'mol_id', 'diameter', 'density'])
		atoms.index = atoms.mol_id
		atoms.index += 1
		atoms.mol_id += 1
		atoms.atom_type = atoms.atom_type.astype(int)
		atoms.diameter = atoms.diameter.astype(float)
		atoms.mol_id = atoms.mol_id.astype(int)

		self.atoms = pd.concat([self.atoms, atoms])
		self.atoms.index = self.atoms.index.astype(int)

	def build(self,
			  spacing: float = 0.0,
			  bead_type: str = "graft",
			  bead_size: float = 329,
			  bead_mass: float = 8487040,
			  chain_height: float = 1750,
			  bond_style: str = "fene",
			  pair_style: str = "lj/cut",
			  bacteria_data: list = [(3, 0.04, 0.000064, 1000), (4, 0.08, 0.000512, 100), (5, 0.12, 0.001728, 10)],
			  bacteria_data_dict: dict = {3: [0.04, 0.000064, 1000], 4: [0.08, 0.000512, 100], 5: [0.12, 0.001728, 10]}
			  ) -> None:
		"""
		Create atom positions and molecular topology for a randomly-grafted monodisperse AdG-brush.
		:param spacing		Number from 0 to 2 (not included!) indicating the degree of intersection between atoms.
							0 means no intersection (consecutive beads).
							1 means intersection of size of half of a bead size (1 = 1 radius)
							2 means full intersection (no chain!) of size of a whole bead (2 = 2 radius = diameter)
		:param bead_type	Type of the beads: 'graft', 'crypt', 'intermediate', 'villus'
		:param bead_size	Diameter of the atoms (micrometers).
		:param bead_mass	Mass of the single bead
		:param chain_height	Desired size of the entire chain of beads (micrometers).
							Number of bead will be calculated according to chain_height and spacing parameters.
		:param bond_style	Could be "fene" or "harmonic"
		:param pair_style	Style of pair interactions. By default lj/cut
		"""
		# Add some parameters that user set for this build
		self.spacing = spacing
		self.bead_type = bead_type
		self.bead_size = bead_size
		self.bead_mass = bead_mass
		self.density = round(bead_mass / ((4/3) * np.pi * (bead_size / 2) ** 3), 10)  # ro = m / V
		self.chain_height = chain_height
		if self.bead_type == 'graft':
			self.n_beads = 1
		else:
			self.n_beads = int(self.chain_height // (self.bead_size - self.spacing * self.bead_size / 2))
		self.bond_style = bond_style
		self.pair_style = pair_style
		self.bacteria_data = bacteria_data

		# Set z size to z_max from _build_chain() if not set

		if self.bead_type != 'bacteria':
			z_max = self._build_chain()
			self.z_max = z_max
			if not self.not_grow_z:
				if not self.box_size[2] or z_max > self.current_height:
					if self.spacing == 0:
						self.box_size[2][1] = z_max - self.bead_size / 2
					else:
						self.box_size[2][1] = z_max + self.bead_size - (self.bead_size / 2) * self.spacing
		else:
			self._build_bacteria(bacteria_data, bacteria_data_dict)

		if self.bead_type != "bacteria":
			# Make dataframes from the non-final lists created by _build_chain()
			if self.mode == "sphere":
				self.atoms = pd.DataFrame(self._atoms_list, columns=['mol_id', 'atom_type', 'diameter', 'density', 'x', 'y', 'z'])
			elif self.mode == "hybrid":
				self.atoms = pd.DataFrame(self._atoms_list, columns=['atom_type', 'x', 'y', 'z', 'mol_id', 'diameter', 'density'])
			elif self.mode == "full":
				self.atoms = pd.DataFrame(self._atoms_list, columns=['mol_id', 'atom_type', 'q', 'x', 'y', 'z'])
			elif self.mode == "angle":
				self.atoms = pd.DataFrame(self._atoms_list, columns=['mol_id', 'atom_type', 'x', 'y', 'z'])
			self.bonds = pd.DataFrame(self._bonds_list, columns=['bond_type', 'atom1', 'atom2'])

			# Angles list based on the bonds list
			if self.mode != "sphere":
				for i in range(len(self._bonds_list) - 1):
					if self._bonds_list[i]['atom2'] == self._bonds_list[i + 1]['atom1']:
						self._angles_list.append({'angle_type': 1,
												'atom1': self._bonds_list[i]['atom1'],
												'atom2': self._bonds_list[i]['atom2'],
												'atom3': self._bonds_list[i + 1]['atom2']
												})

			self.angles = pd.DataFrame(self._angles_list, columns=['angle_type', 'atom1', 'atom2', 'atom3'])
			self.dihedrals = pd.DataFrame(self._dihedrals_list, columns=['dihedral_type', 'atom1', 'atom2', 'atom3',
																		'atom4'])

			# LAMMPS ids start at 1
			self.atoms.index += 1
			self.bonds.index += 1
			self.angles.index += 1
			self.dihedrals.index += 1

	def write(self, filename: str, generation_type: str) -> None:
		"""
		Write the LAMMPS data file.
		:param filename: Filename for the output file.
		"""
		num_atoms = len(self.atoms)
		num_bonds = len(self.bonds)
		num_angles = len(self.angles)
		num_dihedrals = len(self.dihedrals)

		num_atom_types = len(set(self.atoms.atom_type))
		num_bond_types = len(self.BondTypes) if not self.bonds.empty else 0
		num_angle_types = 1 if not self.angles.empty else 0
		num_dihedral_types = len(self.DihedralTypes) if not self.dihedrals.empty else 0

		# Auto-detect gzipped files
		o = gzip.open if filename.endswith('.gz') else open


		with o(filename, 'xt', newline='\n') as f:
			# Header
			f.write("#Header\n")
			f.write(f"{num_atoms} atoms\n")
			if num_bonds > 0:     f.write(f"{num_bonds} bonds\n")
			if num_angles > 0:    f.write(f"{num_angles} angles\n")
			if num_dihedrals > 0: f.write(f"{num_dihedrals} dihedrals\n\n")

			f.write(f"{num_atom_types} atom types\n")
			if num_bond_types > 0:     f.write(f"{num_bond_types} bond types\n")
			if num_angle_types > 0:    f.write(f"{num_angle_types} angle types\n")
			if num_dihedral_types > 0: f.write(f"{num_dihedral_types} dihedral types\n\n")

			# Box geometry
			f.write(f"{min(self.box_size[0][0], self.bac_box_size[0][0])} {max(self.box_size[0][1], self.bac_box_size[0][1])} xlo xhi\n")
			f.write(f"{min(self.box_size[1][0], self.bac_box_size[1][0])} {max(self.box_size[1][1], self.bac_box_size[1][1])} ylo yhi\n")
			if generation_type == "both":
				f.write(f"{self.bottom_padding} {max(self.box_size[2][1], self.bac_box_size[2][1]) + 0.5} zlo zhi\n\n")
			else:
				f.write(f"{self.bottom_padding} {max(self.box_size[2][1], self.bac_box_size[2][1])} zlo zhi\n\n")

			if self.mode != "sphere":
				if self.units == "micro" or self.units == "lj":
					for k, v in self.masses.items():
						self.masses[k] = self.bead_mass
				# Force field coeffs
				f.write("Masses\n\n")
				if generation_type == "both" or generation_type == "villi":
					for k, v in self.masses.items():
						f.write(f"{k.value} {v}\n")
				for t, d, m, _ in self.bacteria_data:
					f.write(f"{t} {m}\n")
				f.write("\n")

			if len(self.pair_ij_coeffs) > 0:
				f.write("PairIJ Coeffs" + (f" # {self.styles['pair']}" if self.styles['pair'] else '') + "\n\n")
				for k, vs in self.pair_ij_coeffs.items():
					f.write(f"{k[0].value} {k[1].value} " + ' '.join([f'{v:.3g}' for v in vs]) + "\n")
				f.write("\n")

			if len(self.bond_coeffs) > 0:
				f.write("Bond Coeffs" + (f" # {self.styles['bond']}" if self.styles['bond'] else '') + "\n\n")
				for k, vs in self.bond_coeffs.items():
					f.write(f"{k.value} " + ' '.join([f'{v:.3g}' for v in vs]) + "\n")
				f.write("\n")

			if len(self.angle_coeffs) > 0:
				f.write("Angle Coeffs" + (f" # {self.styles['angle']}" if self.styles['angle'] else '') + "\n\n")
				for k, vs in self.angle_coeffs.items():
					f.write(f"{k.value} " + ' '.join([f'{v:.3g}' for v in vs]) + "\n")
				f.write("\n")

			if len(self.dihedral_coeffs) > 0:
				f.write("Dihedral Coeffs" + (f" # {self.styles['dihedral']}" if self.styles['dihedral'] else '') + "\n\n")
				for k, vs in self.dihedral_coeffs.items():
					f.write(f"{k.value} " + ' '.join([f'{v:.3g}' for v in vs]) + "\n")
				f.write("\n")

			# Atom properties
			f.write("Atoms # full\n\n")
			self.atoms.to_csv(f, sep=' ', header=False, index=True, line_terminator='\n', float_format='%.6g')
			f.write("\n")

			# Molecular topology
			if len(self.bonds) > 0:
				f.write(f"Bonds\n\n")
				self.bonds.to_csv(f, sep=' ', header=False, index=True, line_terminator='\n', float_format='%.3g')
				f.write("\n")
			if len(self.angles) > 0:
				f.write("Angles\n\n")
				self.angles.to_csv(f, sep=' ', header=False, index=True, line_terminator='\n', float_format='%.3g')
				f.write("\n")
			if len(self.dihedrals) > 0:
				f.write("Dihedrals\n\n")
				self.dihedrals.to_csv(f, sep=' ', header=False, index=True, line_terminator='\n', float_format='%.3g')
				f.write("\n")


class ArchitectureGenerator(BrushGenerator):
	"""
	Generate a LAMMPS data file containing a XXX polymer brush grafted to a planar wall in a rectangular box.
	Mainly based on the KremerGrestBrushGenerator class. Very demo! (not sure about the physics)

	Kremer, K.; Grest, G. S. Dynamics of Entangled Linear Polymer Melts: A Molecular‐dynamics Simulation. J. Chem. Phys.
	1990, 92 (8), 5057–5086. https://doi.org/10.1063/1.458541.
	"""

#    AtomTypes = Enum('AtomTypes', ['graft', 'crypt-bottom', 'crypt-top', 'villus', 'bacteria'])
	AtomTypes = Enum('AtomTypes', ['graft', 'villi', 'bacteria'])
	BondTypes = Enum('BondTypes', ['fene'])

	# in case of atom style = angle or full:
	masses = {
		AtomTypes['graft']  : 1,
		AtomTypes['villi']  : 1,
	}

	# in case of atom style = angle or full:
	styles = {
		'pair': 'lj/cut',
		'bond': 'fene',
	}

	def __init__(self,
				 box_size: tuple[tuple[float, float], tuple[float, float], Optional[tuple[float, float]]],
				 rng_seed: Optional[int],
				 min_dist: int,
				 n_anchors: int,
				 units: str,
				 bottom_padding: float = 0,
				 mode: str = 'sphere',
				 ):
		"""
		:param box_size:        3-tuple of floats describing the dimensions of the rectangular box. If the third (z)
								value is None, it will be automatically sized to contain the longest chain.
		:param rng_seed:        Seed used to initialize the PRNG. May be None, in which case a random seed will be used.
		:param min_dist:        Minimum distance between two centers of atoms: used as minimum distance for the
								Poisson-disk point set generator.
		:param n_anchors:       Number of grafted beads (anchors of the growed brush polymers).
		:param bottom_padding:	Distance between the bottom edge of the box and the grafting layer. Must be positive.
								Use grafted bead size value!
		:param graft:           Generates grafted brushes when True, and non-grafted films when False
		"""
		self.box_size = box_size
		self.rng_seed = rng_seed
		self.min_dist = min_dist
		self.n_anchors = n_anchors
		self.units = units
		self.bottom_padding = bottom_padding
		self.mode = mode

		super().__init__(self.box_size, self.rng_seed, self.min_dist, self.units, self.bottom_padding, self.mode)

	def _build_bead(self, mol_id: int, graft_coord: np.ndarray, bead_id: int) -> float:

		# four options: angle, full, sphere, hybrid

		atom_type = self.AtomTypes[self.bead_type].value

		if self.mode == "sphere":
			# self._atoms_list.append({'mol_id'   : self.current_beads + bead_id + 1,
			self._atoms_list.append({'mol_id'   : mol_id + 1,
									 'atom_type': atom_type,
									 'diameter' : self.bead_size,
									 'density'  : self.density,
									 'x'        : graft_coord[0],
									 'y'        : graft_coord[1],
									 'z'        : self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
									 })  # float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
			# print(atom_type, self.current_beads, bead_id, self.bead_size, self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2)))
		elif self.mode == "hybrid":
			# self._atoms_list.append({'mol_id'   : self.current_beads + bead_id + 1,
			self._atoms_list.append({'atom_type': atom_type,
									 'x'        : graft_coord[0],
									 'y'        : graft_coord[1],
									 'z'        : self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2)),
									 'mol_id'   : mol_id + 1,
									 'diameter' : self.bead_size,
									 'density'  : self.density,
									 })  # float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
			# print(atom_type, self.current_beads, bead_id, self.bead_size, self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2)))
		elif self.mode == "full":
			# self._atoms_list.append({'mol_id'   : self.current_beads + bead_id + 1,
			self._atoms_list.append({'mol_id'   : mol_id + 1,
									 'atom_type': atom_type,
									 'q'  : 0,
									 'x'        : graft_coord[0],
									 'y'        : graft_coord[1],
									 'z'        : self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
									 })  # float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
			# print(atom_type, self.current_beads, bead_id, self.bead_size, self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2)))
		elif self.mode == "angle":
			# self._atoms_list.append({'mol_id'   : self.current_beads + bead_id + 1,
			self._atoms_list.append({'mol_id'   : mol_id + 1,
									 'atom_type': atom_type,
									 'x'        : graft_coord[0],
									 'y'        : graft_coord[1],
									 'z'        : self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
									 })  # float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
			# print(atom_type, self.current_beads, bead_id, self.bead_size, self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2)))

		# Molecular topology
		# Bonds and Angles

		if self.current_beads == 0 and bead_id == 0:
			pass
		elif self.current_beads != 0 and bead_id == 0:
			atom_id = len(self._atoms_list)
			self._bonds_list.append({'bond_type': self.BondTypes[self.bond_style].value,
									 'atom1'    : atom_id - self.n_anchors - self.beads_counter + (self.prev_layer_n_beads - 1) * (self.multiplier + 1),
									 'atom2'    : atom_id
									 })
			self.multiplier += 1

		elif self.current_beads != 0 and bead_id != 0:
			self.beads_counter += 1
			atom_id = len(self._atoms_list)
			self._bonds_list.append({'bond_type': self.BondTypes[self.bond_style].value,
									 'atom1'    : atom_id - 1,
									 'atom2'    : atom_id
									 })

		return self.current_height + float((bead_id + 1) * (self.bead_size - self.spacing * self.bead_size / 2))

	def _build_bacterium(self, mol_id: int, graft_coord: np.ndarray) -> None:

		# four options: angle, full, sphere, hybrid

		atom_type = self.AtomTypes[self.bead_type].value

		if self.mode == "sphere":
			# self._atoms_list.append({'mol_id'   : self.current_beads + bead_id + 1,
			self._atoms_list.append({'mol_id'   : mol_id + 1,
									 'atom_type': atom_type,
									 'diameter' : self.bead_size,
									 'density'  : self.density,
									 'x'        : graft_coord[0],
									 'y'        : graft_coord[1],
									 'z'        : graft_coord[2]
									 })  # float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
			# print(atom_type, self.current_beads, bead_id, self.bead_size, self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2)))
		elif self.mode == "hybrid":
			# self._atoms_list.append({'mol_id'   : self.current_beads + bead_id + 1,
			self._atoms_list.append({'atom_type': atom_type,
									 'x'        : graft_coord[0],
									 'y'        : graft_coord[1],
									 'z'        : graft_coord[2],
									 'mol_id'   : mol_id + 1,
									 'diameter' : self.bead_size,
									 'density'  : self.density,
									 })  # float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
			# print(atom_type, self.current_beads, bead_id, self.bead_size, self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2)))
		elif self.mode == "full":
			# self._atoms_list.append({'mol_id'   : self.current_beads + bead_id + 1,
			self._atoms_list.append({'mol_id'   : mol_id + 1,
									 'atom_type': atom_type,
									 'q'  : 0,
									 'x'        : graft_coord[0],
									 'y'        : graft_coord[1],
									 'z'        : graft_coord[2]
									 })  # float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
			# print(atom_type, self.current_beads, bead_id, self.bead_size, self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2)))
		elif self.mode == "angle":
			# self._atoms_list.append({'mol_id'   : self.current_beads + bead_id + 1,
			self._atoms_list.append({'mol_id'   : mol_id + 1,
									 'atom_type': atom_type,
									 'x'        : graft_coord[0],
									 'y'        : graft_coord[1],
									 'z'        : graft_coord[2]
									 })  # float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2))
			# print(atom_type, self.current_beads, bead_id, self.bead_size, self.current_height + float(bead_id * (self.bead_size - self.spacing * self.bead_size / 2)))


# Class to generate bacteria amount per type
class BacteriaAmountGenerator():
	# create defaul __init__ method
	def __init__(self, 
			  input_file_path: str = "../data/input_data/average_day7.txt",
			  threshold_on_bacteria_amount: float = 10.0):
		self.threshold_on_bacteria_amount = threshold_on_bacteria_amount
		num_bac_genomes = pd.read_csv(input_file_path, sep=" ")
		num_bac_genomes = num_bac_genomes[num_bac_genomes['count'] > threshold_on_bacteria_amount]
		num_bac_genomes.reset_index(drop=True, inplace=True)
		self.num_bac_genomes = num_bac_genomes
		
	def generate_bacteria_amount(self, 
								 output_dir_path: str = "../data/01_bacteria_amount/", 
								 multiplier_on_bacteria_amount: int = 21):
		self.num_bac_genomes['count'] = [round(x) for x in self.num_bac_genomes['count'] * multiplier_on_bacteria_amount]
		total_num_bac = self.num_bac_genomes['count'].sum()
		if not isdir(output_dir_path):
			makedirs(output_dir_path)
		self.num_bac_genomes.to_csv(f"{output_dir_path}/bac_amount_{total_num_bac}_mult_{multiplier_on_bacteria_amount}_thresh_{self.threshold_on_bacteria_amount}.tsv", sep='\t', index=None)
		return self.num_bac_genomes.copy()

# Class to generate initial coordinates of bacteria and/or villi

class BacteriaVilliCoordinatesGenerator():
	# create defaul __init__ method
	def __init__(self, 
			  input_file_path: str = "../data/01_bacteria_amount/bac_amount_518967_mult_435.2017151420611_thresh_10.tsv",
			  generation_type: str = "bacteria", 
			  mode: str = "hybrid", 
			  units: str = "lj", 
			  seed: int = 42, 
			  spacing: float = 1.5, 
			  xlim: Tuple[float, float] = (0, 10), 
			  ylim: Tuple[float, float] = (0, 10), 
			  zlim_villi: Tuple[float, float] = (0, 10), 
			  zlim_bacteria: Tuple[float, float] = (0, 10), 
			  villus_height: int = 40, 
			  villus_width: int = 4, 
			  min_intervillus_width: float = 1.5, 
			  graft_size: int = 4, 
			  crypt_bottom_height: int = 4, 
			  crypt_bottom_width: int = 4, 
			  crypt_top_height: int = 4, 
			  crypt_top_width: int = 4, 
			  graft_mass: int = 1000, 
			  crypt_bottom_mass: int = 1000, 
			  crypt_top_mass: int = 1000, 
			  villus_mass: int = 1000, 
			  n_anchors: int = 42, 
			  bottom_padding: int = -1, 
			  bsize: float = 0.04, 
			  bmass: float = 0.000064):
		# generation_type: 'bacteria' or 'villi' (TBD) or 'both'
		# mode: 'sphere' or 'hybrid' or 'full' or 'angle'
		# units: 'lj' or 'micro'
		# bsize: 0.04 or 0.12
		# bmass: 0.000064 or 0.001728

		self.input_file_path = input_file_path
		self.generation_type = generation_type
		if self.generation_type == 'both':
			self.num_shift = 3
		elif self.generation_type == 'bacteria':
			self.num_shift = 1
		self.mode = mode
		self.units = units
		self.seed = seed
		self.spacing = spacing
		self.xlim = xlim
		self.ylim = ylim
		self.zlim_bacteria = zlim_bacteria
		self.zlim_villi = zlim_villi
		self.box_size = (xlim, ylim, zlim_villi)
		self.bacteria_box_size = (xlim, ylim, zlim_bacteria)

		##### illeum or caecum:
		self.villus_height = villus_height
		self.villus_width = villus_width
		self.min_intervillus_width = min_intervillus_width
		self.graft_size = graft_size
		self.crypt_bottom_height = crypt_bottom_height
		self.crypt_bottow_width = crypt_bottom_width
		self.crypt_top_height = crypt_top_height
		self.crypt_top_width = crypt_top_width
		self.min_dist = graft_size + min_intervillus_width
		##### end of illeum or caecum

		self.graft_mass = graft_mass
		self.crypt_bottom_mass = crypt_bottom_mass
		self.crypt_top_mass = crypt_top_mass
		self.villus_mass = villus_mass
		self.n_anchors = n_anchors
		self.bottom_padding = bottom_padding

		# bacteria data
		bacteria_data_file = pd.read_csv(input_file_path, sep="\t")
		for i in range(bacteria_data_file.shape[0]):
			bacteria_data_file.loc[i, 'genome'] = i + self.num_shift
		self.bacteria_data_file = bacteria_data_file.copy()
		self.n_bacteria = bacteria_data_file['count'].sum()
		self.bsize = bsize
		self.bmass = bmass
		self.bacteria_data = [(idx, bsize, bmass, bcnt) for idx, bcnt in zip(bacteria_data_file['genome'], bacteria_data_file['count'])]
		self.bacteria_data_dict = dict(zip(list(bacteria_data_file['genome']), [row1 + [row2] for row1, row2 in zip([[bsize, bmass]] * bacteria_data_file.shape[0], list(bacteria_data_file['count']))]))
		self.min_bac_dist = bsize * 1.5
		
	def generate_coordinates(self, 
						  output_dir_path: str = "../data/02_coordinates/", 
						  plot_figs: bool = False):
		archgen = ArchitectureGenerator(self.box_size, self.seed, self.min_dist, self.n_anchors, self.units, bottom_padding=self.bottom_padding, mode=self.mode)
		
		# generate villi coordinates
		if self.generation_type == 'both':
			n_actual = archgen.generate_grafting_layer(self.n_anchors, 10**6)
			if n_actual != self.n_anchors:
				print(f'Warning: Grafting layer too dense. {n_actual} grafting points instead of {self.n_anchors}.')
			else:
				print(f'{n_actual} points')
			if plot_figs:
				fig = plt.figure(figsize=(7,7))
				plt.scatter(archgen.coordinates[:, 0], archgen.coordinates[:, 1])
				plt.show()
			
			# Construct grafting beads (fixed layer)
			archgen.build(spacing=self.spacing,
							bead_type="graft",
							bead_size=self.graft_size,
							bead_mass=self.graft_mass,
							chain_height=self.graft_size,
							bond_style="fene",
							pair_style="lj/cut")

			# Construct villi beads (moving layer)
			archgen.build(spacing=self.spacing,
							bead_type="villi",
							bead_size=self.villus_width,
							bead_mass=self.villus_mass,
							chain_height=self.villus_height,
							bond_style="fene",
							pair_style="lj/cut")
			
		# generate bacteria coordinates
		n_bacteria_actual = archgen.generate_bacteria(self.n_bacteria, self.min_bac_dist, self.bacteria_box_size)
		print(self.n_bacteria, n_bacteria_actual)
		if plot_figs:
			fig = plt.figure(figsize=(7,7))
			ax = fig.add_subplot(projection='3d') 
			ax.scatter(archgen.coordinates_bacteria[:, 0], archgen.coordinates_bacteria[:, 1], archgen.coordinates_bacteria[:, 2], s=self.min_bac_dist)
			plt.show()
		archgen.build(spacing=self.spacing,
						bead_type="bacteria",
						bead_size=self.villus_width,
						bead_mass=self.villus_mass,
						chain_height=self.villus_height,
						bond_style="fene",
						pair_style="lj/cut", 
						bacteria_data=self.bacteria_data, 
						bacteria_data_dict=self.bacteria_data_dict)
		# print(archgen.atoms.atom_type.value_counts())

		# write to file
		if not isdir(output_dir_path):
			makedirs(output_dir_path)
		fid = f"{self.generation_type}_zlimbac_{self.zlim_bacteria[0]}_{self.zlim_bacteria[1]}_numbac_{n_bacteria_actual}_bacsize_{self.bsize}_bmass_{self.bmass}"
		archgen.write(f"{output_dir_path}/{fid}.pos", generation_type=self.generation_type)
