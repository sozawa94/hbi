#gnu compiler
#F90= mpifort -fbacktrace -g -fallow-argument-mismatch -ffree-form -ffree-line-length-none -O3 -march=native -fopenmp

#Sherlock (intel MPI)
#F90 = mpiifort -O3 -qopenmp -heap-arrays -ip -g -traceback

#wisteria (Fujitsu compiler)
F90=mpifrtpx -Kfast -Kopenmp

F90FLAGS = $(OPTFLAG)
#LDFLAGS = -mkl=parallel

LINK=$(F90)

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
