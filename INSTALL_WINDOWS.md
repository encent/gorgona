# Installation on Windows

1. Download and install **Ovito Basic**: https://www.ovito.org/
2. Download and install Git: https://github.com/git-for-windows/git/releases/download/v2.51.2.windows.1/Git-2.51.2-64-bit.exe
3. Download and install Miniforge: https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe
4. Once Miniforge installation is done, run Miniforge Prompt from the Windows Start Menu (right-click and select **Run as administrator**). You should see now the terminal window with the *(base)* prefix at the start of the prompt.
5. On your desktop create directory `Gorgona_workshop_13_11_2025`
6. In the Miniforge Prompt window run this sequence of commands:
	```bash
	cd Desktop\Gorgona_workshop_13_11_2025
	git clone https://github.com/encent/gorgona # Clone gorgona repository
	cd gorgona
	mamba env create -f environment.yml # Create gorgona mamba environment. It's important to use mamba, not conda, as conda is very slow.
	mamba activate gorgona # Activate the environment
	pip install -e . # Gorgona installation
	```
7. Download and install Microsoft MPI (msmpisetup.exe file only!): https://www.microsoft.com/en-us/download/details.aspx?id=100593
8. Download and install LAMMPS: https://rpm.lammps.org/windows/LAMMPS-64bit-latest-MSMPI.exe
9. Test. Go to the Windows Start Menu. From there, find the LAMMPS application folder. Click on it and then click on the `Examples` folder. Once inside the `Examples`, go to the `crack` folder. Clear the address bar, put instead `cmd` and press Enter. In the command line prompt put this command: `mpiexec -n 4 lmp -in in.crack.lmp`
10. If this worked correctly, congratulations! Everything is set for the workshop.
