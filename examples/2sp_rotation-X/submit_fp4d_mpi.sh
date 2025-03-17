#!/bin/bash

# queue request
#$ -l h="hnuc12|hnuc13"
# pe request
#$ -pe round_robin 96
# Job name 
#$ -N FP4D_NEO

#$ -S /bin/bash

#$ -V

#$ -cwd

# needs in 
#   $NSLOTS          
#       the number of tasks to be used
#   $TMPDIR/machines 
#       a valid machiche file to be passed to mpirun 
#   enables $TMPDIR/rsh to catch rsh calls if available

echo "Got $NSLOTS slots."
cat $TMPDIR/machines

#######################################################
### MPI JOB
#######################################################
#
#
# Remove all module environments.
 module purge

# Load MPI environments.
 module load sge
 module load netcdf/4.4.5-intel
 module load hdf5/1.12.2-intel
 cd $SGE_O_WORKDIR

export NMPI=2
export OMP_NUM_THREADS=48
export HDF5_MPI="OFF"
OUTFILE=out.fp4d.log
mpirun -machinefile $TMPDIR/machines -n $NMPI /home/hjun/Codes/FP4D_v1.1.a/FP4D >& ${OUTFILE}




