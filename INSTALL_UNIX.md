# Installation on macOS/Linux

1. Download and install **Ovito Basic**: https://www.ovito.org/
2. In Terminal app, check if you have Git: `git --version`. If you do have it, go to the next step. If not, let me know.
3. Install Miniforge:
	```bash
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
	# or
	wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
	bash Miniforge3-$(uname)-$(uname -m).sh
	```
4. After installation is done, close and open again the Terminal app. (base) should appear next to the command line prompt.
5. On your desktop create directory `Gorgona_workshop_13_11_2025`
6. In Terminal run these commands:
	```bash
	cd Desktop/Gorgona_workshop_13_11_2025
	git clone https://github.com/encent/gorgona # Clone gorgona repository
	cd gorgona
	mamba env create -f environment.yml # Create gorgona conda environment
	mamba activate gorgona # Activate the environment
	pip install -e . # Gorgona installation
	mamba install -c conda-forge lammps # LAMMPS installation
	git clone https://github.com/lammps/lammps # Clone lammps repository
	# Test (if lammps works properly)
	cd lammps/examples/crack
	mpirun -n 4 lmp_mpi -in in.crack
	```
7. If this worked correctly, congratulations! Everything is set for the workshop.

