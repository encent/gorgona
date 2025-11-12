# Installation on Mac (Intel chip)

1. Download and install Ovito Basic: https://www.ovito.org/
2. Download and install Git: 
3. You must have `Python 3` and `conda` installed on your system. Open Terminal app to check whether you have both:
	```bash
	python3 --version # to check that the installation went properly. Should print version of Python.
	conda --version # should print the version number. If not -- see the instructions below.
	```
4. If Python 3 is missing, follow the installation instructions from [the Python official website](https://www.python.org/)
5. If `conda` is missing, follow [these instructions](https://www.anaconda.com/docs/getting-started/miniconda/install) to install Miniconda.
6. After installation is done, close and open again the Terminal app. (base) should appear next to the command line prompt.
7. On your desktop create the directory `Gorgona_workshop_13_11_2025`
8. In the Anaconda Prompt window run this sequence of commands:
	```bash
	cd Desktop/Gorgona_workshop_13_11_2025
	git clone https://github.com/encent/gorgona # Clone gorgona repository
	cd gorgona
	conda env create -f environment.yml # Create gorgona conda environment
	conda activate gorgona # Activate the environment
	pip install -e . # Gorgona installation
	conda install lammps # LAMMPS installation
	conda install openkim-models
	```
9. Test. Go to the LAMMPS directory, then to the `Examples` folder. Once inside the `Examples`, go to the `crack` folder -- open it to the Terminal. In the command line prompt put this command: `mpirun -n 4 lmp_mpi -in in.crack.lmp`
11. If this worked correctly, congratulations! Everything is set for the workshop.
