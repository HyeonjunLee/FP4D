- Libraries for FP4D
    1) intel compiler
    2) NetCDF
    3) HDF5
    4) ....
  You can use the GNU compiler, but in that case,
  the rest of the library must also be compiled through the GNU compiler.

- How do you run it?
   qsub submit_fp4d_mpi.sh

- How do I check the result?
    USE THE PLOT TOOL THROUGH JUPYTER NOTEBOOK
      module load anaconda
      jupyter notebook

   /home/youknow/python/FP4D/v3_2/v3_2_example.ipynb

///FP4D///
  1) You can compile this code by running "make" in the top-level directory.
  2) The default compilation environments provided are 'HNUC' and 'KAIROS'.
     For other environments, you need to manually modify the 'Makefile'
     to compile successfully.

///eqmodule/// courtesy to J. Song
  BlaBlaBla

///Miscellaneous///
  1) "blacs_mod.F90" in src
    Not currently used.
  2) "interpol.f" in src
    Not currently used.
  3) "blacs_mod.F90"
    Not currently uesd.
