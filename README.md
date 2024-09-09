# HIN - (Heterogeneous) Ice Nucleation
This set of codes (written chiefly in Fortran 90) allows to calculate several properties of water and ice at the interface with molecules, biological material, or solid surfaces.

Please note that the latest version of the code as well as the relevant documentation refer to the `clathrate` branch of the code.

The code is divided into two main sections: `STRUCTURE` and `DYNAMICS`. 

The `DYNAMICS` section is functioning but outdated. It allows for the calculations of several dynamical properties of water, most notably the (self part of the) incoherent intermediate scattering function and the dynamical susceptibility, as discussed in , e.g., this paper: https://pubs.acs.org/doi/10.1021/jp507361f 

The `STRUCTURE` section is under constant development. It allows for the calculations of several structural properties of interfacial water and ice, as discussed in, e.g., this paper: https://doi.org/10.1039/D1CP05465A. Please refer to the Wiki for further details.

## Requirements
The code is intended to post-process molecular dynamics (MD) trajectories in `.xtc` format, coupled wit configuration files in `.gro` format - i.e., typical output formats of a GROMACS (https://www.gromacs.org) MD simulation. 
Both sections of the code rely on the `libxdrfile` library to efficiently process these MD trajectories. `libxdrfile` is provided within the `HIN` code (see the relevant LIBXDR directory) and should be compiled before building either the `STRUCTURE` or `DYNAMICS` sections of the code.
Makefiles for both Linux and MacOS are provided in the `BUILD` directory within either the `STRUCTURE` or `DYNAMICS` sections of the code. Any standard Fortran compiler should work. 

The R.I.N.G.S. code (https://sourceforge.net/projects/rings-code/) is optional - it is used to compute rings statistics as well as other topological features of the water/ice network, such as the "cages" discussed in, e.g., this paper: https://doi.org/10.1073/pnas.1509267112. Note that the R.I.N.G.S. code has now evolved into the I.S.A.A.C. code (http://people.se.cmich.edu/petko1vg/isaacs/phys/rings.html), which however we have not tested just yet.

The PLUMED code (https://www.plumed.org) is also optional - it is used to identify ice clusters within the liquid water phase.

## Documentation
Please refer to the Wiki: https://github.com/gcsosso/HIN/wiki, bearing in mind this documentation refers to the `clathrate` branch of the code.
