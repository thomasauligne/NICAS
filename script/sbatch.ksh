#!/bin/ksh
#----------------------------------------------------------------------
# Korn shell script: sbatch
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2017 METEO-FRANCE
#----------------------------------------------------------------------

# Define root directory
rootdir=/home/gmap/mrpa/menetrie/codes/nicas

# Define model and xp
model=arp
xp=6B60

# Define resolution
resol=8
typeset -RZ3 resol

# Define mpicom
mpicom=2

# Define nproc
nproc="0008"

# Define data directory
datadir=${rootdir}/data/${model}/${xp}

# Define file name
filename=${model}_${xp}_resol-${resol}

# New working directory
workdir=${rootdir}/${filename}
rm -fr ${workdir}
mkdir ${workdir}

# Link to the distribution file
ln -sf ${datadir}/${model}_${xp}_distribution_${nproc}.nc ${datadir}/${model}_${xp}_resol-${resol}_distribution_${nproc}.nc

#----------------------------------------------------------------------
# Compute NICAS parameters
#----------------------------------------------------------------------

# Namelist
prefix=${filename}
sed -e "s|_DATADIR_|${datadir}|g" -e "s|_PREFIX_|${prefix}|g" -e "s|_RESOL_|${resol}|g" -e "s|_NPROC_|${nproc}|g" -e "s|_MPICOM_|${mpicom}|g" ${rootdir}/run/namelist_${model}_${xp}_bull > ${workdir}/namelist

# Job
#----------------------------------------------------------------------
cat<<EOFNAM >${workdir}/job_nicas.ksh
#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 40
#SBATCH -t 00:30:00
#SBATCH -p normal64,huge256
#SBATCH --exclusiv
#SBATCH -e ${workdir}/output
#SBATCH -o ${workdir}/output

if [ -n "\$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=\$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=\$omp_threads

ulimit -s unlimited
srun --mpi=pmi2 --ntasks 4 --ntasks-per-node 4 ${rootdir}/run/nicas < ${workdir}/namelist
EOFNAM

#----------------------------------------------------------------------

# Execute
sbatch ${workdir}/job_nicas.ksh
