# Common commands
# -------------------------------
# To build:
#   make
# To delete .mod and .o files:
#   make clean

# First, make sure the environment variable FP4D_PLATFORM is set appropriately:
ifndef FP4D_PLATFORM
FP4D_PLATFORM=HNUC
endif

SYSTEM_MAKEFILE = ${FP4D_ROOT}/platforms/makefile.$(FP4D_PLATFORM)
ifeq ($(wildcard $(SYSTEM_MAKEFILE)),)
$(error Based on the environment variable FP4D_PLATFORM, make is trying to include the file $(SYSTEM_MAKEFILE), but this file does not exist)
endif

TARGET = FP4D

.PHONY: clean

all: $(TARGET)

include $(SYSTEM_MAKEFILE)

include makefile.depend

# Export variables that will be used in the quadpack subdirectory or for "make test":
#export

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
	rm -f *.o ${EQMODULE_DIRS}/*.o *.mod *.a *.nc *.dat compile_log.txt *~ $(TARGET)
	@echo "**********************"
	@echo "make clean is done" 
	@echo "**********************"
