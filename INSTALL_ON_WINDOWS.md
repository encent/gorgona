# Installation on Windows

1. Download and install Ovito Basic: https://www.ovito.org/
2. Download and install Git: https://github.com/git-for-windows/git/releases/download/v2.51.2.windows.1/Git-2.51.2-64-bit.exe
3. You must have `Python 3` and `conda` installed on your system. Open Windows Terminal or PowerShell application. Install Python 3 and conda:
	```bash
	python3 # should redirect you to the Microsoft Store Python 3 page, the click on Install button
	python3 --version # to check that the installation went properly. Should print version of Python.
	conda --version # should print the version number. If not -- see the instructions below.
	```
4. Install `conda` from here: https://www.anaconda.com/docs/getting-started/miniconda/install 
(Step 1: click on the PowerShell tab, copy the command, and paste into the Terminal; skip Step 2)
5. After installation is done, run Anaconda Prompt from the Windows Start Menu (right-click and select "Run as administrator").
6. On your desktop create the directory `Gorgona_workshop_13_11_2025`
7. In the Anaconda Prompt window run this sequence of commands:
	```bash
	cd Desktop\Gorgona_workshop_13_11_2025
	git clone https://github.com/encent/gorgona # Clone gorgona repository
	cd gorgona
	conda env create -f environment.yml # Create gorgona conda environment
	conda activate gorgona # Activate the environment
	pip install -e . # Gorgona installation
	```
8. Download and install Microsoft MPI (msmpisetup.exe only!): https://www.microsoft.com/en-us/download/details.aspx?id=100593
9. Download and install LAMMPS: https://rpm.lammps.org/windows/LAMMPS-64bit-latest-MSMPI.exe
10. Test. Go to the Windows Start Menu. From there, find the LAMMPS application (folder). Click on it and then click on the `Examples` folder. Once inside the `Examples`, go to the `crack` folder. Clear the address bar, put instead `cmd` and press Enter. In the command line prompt put this command: `mpiexec -n 4 lmp -in in.crack.lmp`
11. If this worked correctly, congratulations! Everything is set for the workshop.

