# ------------------------------------------------------------------------------
# GNU gmake file.

# Change these lines to fit your system. fortranint should be 4 or 8 bytes. This
# was established when you compiled LAPACK and BLAS.
fortranint = 4
BLASLIB = -lblas
LAPACKLIB = -llapack
FORTRANLIB = -lgfortran

CPP = g++
MPICPP = mpic++
FORTRAN = gfortran

# Set the optimization level.
#opt = 
#opt = -g
opt = -O3

# Choose serial or OpenMP-parallelized versions. dc3dm does not currently
# support MPI, but hmmvp does. So if you compile with mode = mpi, libhmmvp_mpi.a
# is built, but bin/dc3dm is not.
#mode = s
mode = omp
#mode = mpi

# Probably does not need to be changed:
DC3DM = dc3dm
ifeq ($(mode),s)
	MODE_FLAGS =
	ext = omp
endif
ifeq ($(mode),omp)
	MODE_FLAGS = -fopenmp -DUTIL_OMP
	ext = omp
endif
ifeq ($(mode),mpi)
	CPP = $(MPICPP)
	MODE_FLAGS = -DUTIL_MPI
	ext = mpi
	DC3DM =
endif

# ------------------------------------------------------------------------------
# The rest should not have to be changed.

INCLUDE = -I .
LIBS = $(LAPACKLIB) $(BLASLIB) $(FORTRANLIB)
LIBDIRS =
OPTFLAGS = $(opt)
CPPFLAGS = $(OPTFLAGS) $(MODE_FLAGS) -DFORTRAN_INT_$(fortranint)
FFLAGS = $(OPTFLAGS) $(MODE_FLAGS) -Wall
LDFLAGS = $(MODE_FLAGS)

.SUFFIXES:
.SUFFIXES: .cpp .f .f90 .o

# dc3dm-specific files.
DCPPSRCS = src/RectMeshUnstruct.cpp src/BrickMeshBemBuilder.cpp \
src/Triangulation.cpp src/MeshAnalyzer.cpp src/Elastostatics.cpp src/dc3dm.cpp \
src/ValueSetter.cpp
DFSRCS = external/dc3omp.f
DOBJECTS = $(patsubst %.cpp,%.o,$(DCPPSRCS)) $(patsubst %.f,%.o,$(DFSRCS))

# hmmvp files.
HCPPSRCS = src/Hd.cpp src/Compress.cpp src/Hmat.cpp src/HmatIo.cpp \
src/KeyValueFile.cpp src/CodeAnalysis.cpp src/Mpi.cpp src/CHmat.cpp \
src/SFHmat.cpp
HOBJECTS = $(patsubst %.cpp,%.o,$(HCPPSRCS))

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.f
	$(FORTRAN) $(FFLAGS) -c $< -o $@

%.o : %.f90
	$(FORTRAN) $(FFLAGS) -c $< -o $@

all: libhmmvp $(DC3DM) mvp cmvp fmvp

# Library to compress and apply an H-matrix.
libhmmvp: $(HOBJECTS)
	ar rucs lib/libhmmvp_$(mode).a $(HOBJECTS)
	rm -f src/*.o

# dc3dm is a program to form and compress the DDM operator.
dc3dm: $(DOBJECTS) libhmmvp
ifneq ($(mode),mpi)
	$(CPP) src/dc3dmMain.cpp $(DOBJECTS) $(INCLUDE) $(CPPFLAGS) $(LDFLAGS) lib/libhmmvp_$(mode).a $(LIBS) -o bin/dc3dm
	rm -f src/*.o
endif

# C++ examples of matrix-vector product.
mvp:
	$(CPP) examples/mvp_$(ext).cpp $(INCLUDE) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) lib/libhmmvp_$(mode).a $(LIBS) -o examples/mvp_$(mode)

# C example of MVP.
cmvp:
ifeq ($(ext),omp)
	$(CPP) examples/cmvp_$(ext).c $(INCLUDE) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) lib/libhmmvp_$(mode).a $(LIBS) -o examples/cmvp_$(mode)
endif

# Fortran 90 example.
fmvp:
ifeq ($(ext),omp)
	$(FORTRAN) examples/fmvp_$(ext).f90 $(INCLUDE) $(LDFLAGS) lib/libhmmvp_$(mode).a $(LIBS) -lstdc++ -o examples/fmvp_$(mode)
endif

clean:
	rm -f src/*.o external/*.o lib/*.a bin/*
