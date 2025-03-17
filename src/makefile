# Common commands
# -------------------------------
# To build:
#   make
# To delete .mod and .o files:
#   make clean

# First, make sure the environment variable FP4D_SYSTEM is set appropriately:

ifndef FP4D_SYSTEM
#$(error You must set the environment variable FP4D_SYSTEM to select which makefile to use. For example, to use makefiles/makefile.edison, set P1ICP_SYSTEM=edison. To set this environment variable, use the export or setenv command, depending on your shell)
FP4D_SYSTEM=hnuc
endif

SYSTEM_MAKEFILE = makefiles/makefile.$(FP4D_SYSTEM)
ifeq ($(wildcard $(SYSTEM_MAKEFILE)),)
$(error Based on the environment variable P1ICP_SYSTEM, make is trying to include the file $(SYSTEM_MAKEFILE), but this file does not exist)
endif

TARGET = FP4D

.PHONY: clean

all: $(TARGET)

include $(SYSTEM_MAKEFILE)

include makefile.depend

# Export variables that will be used in the quadpack subdirectory or for "make test":
export

#EQMODULE_DIRS = ../tools/eqmodule
vpath %.F90 $(EQMODULE_DIRS)
vpath %.c $(EQMODULE_DIRS)

%.o: %.F90
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c  $< -o $@

$(TARGET): $(OBJ_FILES)
	$(FLINKER) -o $(TARGET) $(OBJ_FILES) $(EXTRA_LINK_FLAGS) $(PETSC_LIB)
	@echo "**********************"
	@echo "FP4D Compiling is done"
	@echo "**********************"
	@echo "Compile completed at: $$(date '+%Y-%m-%d %H:%M:%S')" >> compile_log.txt
clean::
	rm -f *.o *.mod *.a *.nc *.dat *~ $(TARGET)
	@echo "**********************"
	@echo "make clean is done" 
	@echo "**********************"
