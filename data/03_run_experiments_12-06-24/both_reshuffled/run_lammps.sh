#!/bin/bash

#SBATCH --job-name=lammps_both_reshuffled
#SBATCH --output=%x.%j.%N.out
#SBATCH --error=%x.%j.%N.err
#SBATCH --nodes=1
#SBATCH --ntasks=100
#SBATCH --time=07-00:00:00
#SBATCH --partition=genD
#SBATCH --qos=marathon
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1

module purge
module load Anaconda3
source ~/.bashrc
conda activate lammps
module load GCC/11.3.0
module load GCCcore/11.3.0
module load CMake/3.23.1-GCCcore-11.3.0
module load OpenMPI/4.1.4-GCC-11.3.0

n_steps=3000000
id_bac_coef=0.00013497191736
vil_bac_coef=0.00013497191736

bac_size=0.04
bac_nve_limit=2.5
id_bac_cf=2.5
diff_bac_cf=2.5
vil_bac_cf=2.5
global_cf=2.5
neigh_modify_size=55000
page_size=1000000
vil_bac_size=$(echo "2 + ${bac_size}/2" | bc -l)
bac_veloc=4e-09
temp_bac=4e-09

# export OMP_NUM_THREADS=100

mpirun lmp_mpi -in run.in -var n_steps ${n_steps} -var temp_bac ${temp_bac} -var bac_nve_limit ${bac_nve_limit} -var bac_veloc ${bac_veloc} -var neigh_modify_size ${neigh_modify_size} -var bac_size ${bac_size} -var id_bac_cf ${id_bac_cf} -var vil_bac_size ${vil_bac_size} -var vil_bac_cf ${vil_bac_cf} -var id_bac_coef ${id_bac_coef} -var vil_bac_coef ${vil_bac_coef} -var global_cf ${global_cf} -var diff_bac_cf ${diff_bac_cf} -var page_size ${page_size}
