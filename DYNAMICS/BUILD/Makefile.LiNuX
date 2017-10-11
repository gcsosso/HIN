FC_LINUX=gfortran
FFLAGS_LINUX=-ffree-line-length-0 -fbacktrace # -g -fbounds-check -Wall -fbacktrace -finit-real=nan # Debugging flags...
LDLIBS_LINUX=-L/home/gsosso/CODES/LIB/LIBXDR/xdrfile-1.1.1/BIN/lib -lxdrfile
EXE_LINUX=hin_dynamics_LiNuX.x

SRC_DIR=../SRC
XTC_INTERFACE_DIR=../../XTC_INTERFACE
BIN_DIR=../BIN

linux: xtc-interface.o module_read_input.o hin_dynamics.o
	$(FC_LINUX) $(FFLAGS_LINUX) xtc-interface.o module_read_input.o hin_dynamics.o $(LDLIBS_LINUX) -o $(BIN_DIR)/$(EXE_LINUX)
xtc-interface.mod: $(XTC_INTERFACE_DIR)/xtc-interface.o $(XTC_INTERFACE_DIR)/xtc-interface.f90
	$(FC_LINUX) $(FFLAGS_LINUX) -c $(XTC_INTERFACE_DIR)/xtc-interface.f90       
xtc-interface.o: $(XTC_INTERFACE_DIR)/xtc-interface.f90
	$(FC_LINUX) $(FFLAGS_LINUX) -c $(XTC_INTERFACE_DIR)/xtc-interface.f90
module_read_input.mod: $(SRC_DIR)/module_read_input.o $(SRC_DIR)/module_read_input.f90
	$(FC_LINUX) $(FFLAGS_LINUX) -c $(SRC_DIR)/module_read_input.f90
module_read_input.o: $(SRC_DIR)/module_read_input.f90
	$(FC_LINUX) $(FFLAGS_LINUX) -c $(SRC_DIR)/module_read_input.f90     
hin_dynamics.o: $(SRC_DIR)/hin_dynamics.f90
	$(FC_LINUX) $(FFLAGS_LINUX) -c $(SRC_DIR)/hin_dynamics.f90 

clean:
	rm -r -f *.o *.mod
