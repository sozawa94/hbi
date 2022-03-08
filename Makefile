#eqrupt03
#OPTFLAGS =  -O3 -qopenmp -xCORE-AVX2 -ip -g -traceback
#F90= mpiifort

#ofp
#OPTFLAGS=-O3 -axMIC-AVX512 -g -traceback
#F90= mpiifort

#wisteria
OPTFLAGS = -Kfast -Kopenmp
F90=mpifrtpx -Kfast

F90FLAGS = $(OPTFLAG)
LDFLAGS = -mkl=parallel

LINK=$(F90)

#Lattice H rectangle
#OBJS= m_const.o TDstressFS.o HACApK_lib.o m_HACApK_calc_entry_ij.o m_HACApK_base_LH.o m_HACApK_solve_LH.o m_HACApK_use_LH.o main_new.o \

#Lattice H square
OBJS= m_const.o okada.o TDstressFS.o HACApK_lib.o m_HACApK_calc_entry_ij.o m_HACApK_base.o m_HACApK_solve.o m_HACApK_use.o main_LH.o \

TARGET=lhbiem

.SUFFIXES: .o .f90

$(TARGET): $(OBJS)
			$(LINK) -o $@ $(OBJS) $(LDFLAGS)

.f90.o: *.f90
			$(F90) -c $< $(F90FLAGS)
clean:
	rm -f *.o *.mod $(TARGET)

rmod:
	rm -f m_*.o *.mod
