#eqrupt03
OPTFLAGS =  -O3 -qopenmp -xCORE-AVX2 -ip -g -traceback
#ofp
#OPTFLAGS =  -O3 -qopenmp -axMIC-AVX512 -ip -g -traceback
#gpusolid
#OPTFLAGS =  -O3 -xCORE-AVX2 -ip -g -traceback

F90=mpiifort
F90FLAGS = $(OPTFLAGS) -fpp
LDFLAGS = -mkl=parallel

LINK=$(F90)

#Lattice H rectangle
#OBJS= m_const.o TDstressFS.o HACApK_lib.o m_HACApK_calc_entry_ij.o m_HACApK_base_LH.o m_HACApK_solve_LH.o m_HACApK_use_LH.o main_LHrec.o \

#Normal H
OBJS= m_const.o TDstressFS.o HACApK_lib.o m_HACApK_calc_entry_ij.o m_HACApK_base.o m_HACApK_solve.o m_HACApK_use.o main_new.o \

#Lattice H square
#OBJS= m_const.o TDstressFS.o HACApK_lib.o m_HACApK_calc_entry_ij.o m_HACApK_base_LH.o m_HACApK_solve_LH.o m_HACApK_use_LH.o main_LH2.o \

TARGET=hbiem
#TARGET=lhbiem

.SUFFIXES: .o .f90

$(TARGET): $(OBJS)
			$(LINK) -o $@ $(OBJS) $(LDFLAGS)

.f90.o: *.f90
			$(F90) -c $< $(F90FLAGS)
clean:
	rm -f *.o *.mod $(TARGET)

rmod:
	rm -f m_*.o *.mod
