set(CMAKE_Fortran_COMPILER "/glade/u/apps/ch/opt/ncarcompilers/0.4.1/mpi/mpifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_PLATFORM_ID "")

set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "netcdff;netcdf;hdf5_hl;hdf5;sz;z;gfortran;m;dl;rt;dl;mpi_usempif08;mpi_usempi_ignore_tkr;mpi_mpifh;mpi;open-rte;open-pal;gfortran;m;pciaccess;osmcomp;ibverbs;rdmacm;pbs;gpfs;psm_infinipath;util;gfortran;m;quadmath;m;pthread;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/glade/u/apps/ch/os/usr/lib64;/glade/u/apps/ch/os/lib64;/glade/u/apps/ch/os/usr/lib;/glade/u/apps/ch/opt/netcdf/4.4.1.1/gnu/6.3.0/lib;/glade/u/apps/ch/os/lib;/glade/u/apps/ch/opt/openmpi/2.1.0/gnu/6.3.0/lib;/glade/u/apps/ch/opt/pbs_copy/lib;/glade/u/apps/ch/opt/gnu/6.3.0/lib/gcc/x86_64-pc-linux-gnu/6.3.0;/glade/u/apps/ch/opt/gnu/6.3.0/lib64;/lib64;/usr/lib64;/glade/u/apps/ch/opt/gnu/6.3.0/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")