#!/bin/bash

# queue request
#$ -l h="hnuc34"

# pe request
#$ -pe mpi_64 64

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

export OMP_NUM_THREADS=1
OUTFILE=out.fp4d.log
#/home/youknow/Codes/FP4D/FP4D >& ${OUTFILE}
/home/youknow/Codes/FP4D_v1.1_mpi_oo/FP4D >& ${OUTFILE}





