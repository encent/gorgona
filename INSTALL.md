
# A guide on how to install and build LAMMPS on Mac with M1 chip

This guide is prepared in order to provide the instructions on how to install LAMMPS on Mac machine with M1 chip.

### Prerequisites

Be sure that you installed these software before:

1.	Install Xcode Command Line Tools. Type `xcode-select --install` in a command line
2.	[Install CMake](https://cmake.org/install/) (latest)
3.	[Install Homebrew](https://brew.sh/) (latest)
4.	[Install Anaconda](https://docs.anaconda.com/anaconda/install/mac-os/) (latest)
5.	[Install VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) (latest). Type `omz update` in a command line before to lunch the VMD for the first time.
6.	[Install Chimera](https://www.cgl.ucsf.edu/chimera/download.html) (latest) or [Ovito](https://www.ovito.org/macos-downloads/) (latest). These tools are usefool to visualize the initial conformation of atoms
7.	[Install ffmpeg](https://bbc.github.io/bbcat-orchestration-docs/installation-mac-manual/) library.
8.	[Install OpenMPI]: either `brew install openmpi` or `brew install mpich`


I would also suggest to create the separate Anaconda environment for the future LAMMPS experiments.

### LAMMPS installation

I would suggest first to install the LAMMPS as an executable (homebrew or anaconda). The reason is that it is simplier than to build LAMMPS from source, and you can anytime delete it and switch to the most 'powerful' way of working with LAMMPS, by building it from the source (keeping the recent updates). But if you plan to work with CLI and/or Python, than I would suggest to download LAMMPS as an executable for Anaconda:

Follow the instructions [here](https://docs.lammps.org/Install_conda.html). Make sure you install LAMMPS in a separate environment (it is highly recommended).

```bash
conda config --add channels conda-forge
conda activate my-lammps-env
conda install lammps
conda install openkim-models
```

### LAMMPS installation & build from source

Clone (or pull) the latest lammps release:

```bash
cd ~/Tools
git clone -b stable --depth 1 https://github.com/lammps/lammps.git lammps
```

Prepare build folder:

```bash
cd lammps
mkdir build
cd build
```

Build basic lammps pre-sets (with BPM package):

```bash
cmake -C ../cmake/presets/basic.cmake -C ../cmake/presets/nolib.cmake -D FFMPEG_EXECUTABLE=~/Tools/ffmpeg -D PKG_BPM=yes -D LAMMPS_MACHINE=mpi ../cmake
```

If you want to work with lammps by scripting on python, you may need to re-build with this: 

```bash
cmake -D PKG_PYTHON=yes -D BUILD_SHARED_LIBS=yes -D LAMMPS_MACHINE=name . # name = mpi, serial, mybox, titan, laptop, etc
```

If you want to have in your build some [auxiliary tool](https://docs.lammps.org/Tools.html), you may need to re-build with this:
```bash
cmake -D BUILD_TOOLS=yes -D BUILD_LAMMPS_SHELL=yes .
```

If you have some troubles with JPEG, PNG and/or MPEG, try to adjust the previous configuration with these options (where `path` should be the path to the corresponding `.h` and `.a` (`.so`) files):

```bash
cmake -D JPEG_INCLUDE_DIR=path -D JPEG_LIBRARY=path -D PNG_INCLUDE_DIR=path -D PNG_LIBRARY=path -D ZLIB_INCLUDE_DIR=path -D ZLIB_LIBRARY=path -D FFMPEG_EXECUTABLE=path .
```

After this, type:

```bash
make                        # perform make after CMake command
make install                # perform the installation into prefix
```

After that, symlinc lammps executable into the path:

```bash
sudo ln -s /Users/nbykov/Tools/lammps/build/lmp_mpi /usr/local/bin
```

### Test

Try this [example](https://github.com/lammps/lammps/tree/develop/examples/min) and compare it with the logs: `lmp_serial -in in.min`.
