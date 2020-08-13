#
# Build PIDOSITS
#
OPT = yes

FORTRAN = mpif90
FC = gfortran
#FC = ifort

FBLAS = -lopenblas #-llapack

ifeq ($(FC),gfortran)
  TARGET = pbj
  FFLAGS = -O3
  ifeq ($(OPT),no)
    TARGET = pbj_db
    FFLAGS = -O0 $(OMPFLAG) -g -fbounds-check -fbacktrace
  endif
endif

ifeq ($(FC),ifort)
  TARGET = ipbj
  FFLAGS = -O3
  ifeq ($(OPT),no)
    TARGET = ipbj_db
    FFLAGS = -O0 $(OMPFLAG) -g -check bounds -traceback -warn unused
  endif
endif

.SUFFIXES:	.f90 .o

OBJS = afcm.o conk.o invar.o readbc.o angle.o echo.o jima.o \
       readmt.o gammas.o main.o readsrc.o totvar.o \
       matsweep.o readxs.o version.o \
       idomats.o output.o idot.o solvar.o \
       procgrid.o solve.o check.o input.o psds.o \
       flxgthr.o idocomb.o pgsblk.o pgsred.o

SRCS = afcm.f90 conk.f90 invar.f90 readbc.f90 angle.f90 echo.f90 \
       jima.f90 readmt.f90 gammas.f90 main.f90 \
       readsrc.f90 totvar.f90 matsweep.f90 \
       readxs.f90 version.f90 idomats.f90 output.f90 \
       idot.f90 solvar.f90 \
       procgrid.f90 solve.f90 check.f90 input.f90 psds.f90 \
       flxgthr.f90 idocomb.f90 pgsblk.f90 pgsred.f90

$(TARGET) :	$(OBJS)
		$(FORTRAN)  $(FFLAGS) $(FBLAS) -o $@ $(OBJS) -llapack

afcm.o:		invar.o solvar.o
angle.o:	invar.o
check.o:	invar.o totvar.o
conk.o:		invar.o solvar.o
echo.o:		invar.o totvar.o
flxgthr.o:	invar.o solvar.o
gammas.o:	solvar.o
idocomb.o:	invar.o solvar.o
idomats.o:	invar.o solvar.o
idot.o:		invar.o solvar.o
input.o:	invar.o totvar.o
jima.o:		invar.o solvar.o
main.o:		invar.o solvar.o totvar.o
matsweep.o:	invar.o solvar.o
output.o:	invar.o solvar.o
pgsblk.o:	invar.o solvar.o
pgsred.o:	invar.o solvar.o
procgrid.o:	invar.o
psds.o:		invar.o totvar.o
readbc.o:	invar.o totvar.o
readmt.o:	invar.o totvar.o
readsrc.o:	invar.o totvar.o
readxs.o:	invar.o
solve.o:	invar.o solvar.o

#
# Fortran rules
#
%.o:	%.f90
	$(FORTRAN) $(FFLAGS) -c $<

#
# Cleanup
#
clean:
	rm -f *.o *.mod

#
# Link compiled files only. No recompile.
#
link:
	$(FORTRAN) $(FFLAGS) $(FBLAS) -o $(TARGET) $(OBJS)
