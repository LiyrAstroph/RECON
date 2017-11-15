SHELL=/bin/bash

CC       = mpicc
OPTIMIZE = -O2 -Wall


#------------target system---------
#SYSTEM="Darwin"
SYSTEM="Linux"
#SYSTEM="Cluster"
#SYSTEM="TianheII"

ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
LAPACK_INCL = -I/usr/include/lapacke
LAPACK_LIBS = -L/usr/lib64 -llapacke -llapack -lblas
DNEST_INCL  = -I /home/liyropt/Projects/GIT/DNest/
DNEST_LIBS  = -L /home/liyropt/Projects/GIT/DNest -ldnest
FFTW_INCL   = $(shell pkg-config --cflags fftw3) 
FFTW_LIBS   = $(shell pkg-config --libs fftw3) 

MPICHINCL     = $(shell pkg-config --cflags mpich) 
MPICHLIB    = $(shell pkg-config --libs mpich) 
endif

ifeq ($(SYSTEM), "Cluster")
GSL_INCL = -I/mbh/mbhd01/user/liyanrong/soft/gsl/include
GSL_LIBS = -L/mbh/mbhd01/user/liyanrong/soft/gsl/lib  -lgsl -lgslcblas -lm
MPICHLIB = -L/mbh/mbhd01/user/liyanrong/soft/mpich3/lib -lmpich
MPIINCL  = -I/mbh/mbhd01/user/liyanrong/soft/mpich3/include
LAPACK_INCL = -I/mbh/mbhd01/user/liyanrong/soft/lapack/include
LAPACK_LIBS = -L/mbh/mbhd01/user/liyanrong/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
FFTW_INCL = -I/mbh/mbhd01/user/liyanrong/soft/fftw/include
FFTW_LIBS = -L/mbh/mbhd01/user/liyanrong/soft/fftw/lib -lfftw3

DNEST_INCL  = -I /mbh/mbhd01/user/liyanrong/GIT/DNest/
DNEST_LIBS  = -L /mbh/mbhd01/user/liyanrong/GIT/DNest -ldnest
endif

ifeq ($(SYSTEM), "TianheII")
GSL_INCL =
GSL_LIBS = -lgsl -lgslcblas -lm
MPICHLIB = -lmpich
MPIINCL  =
LAPACK_INCL = -I/HOME/ihep_yrli_1/BIGDATA/soft/lapack/include
LAPACK_LIBS = -L/HOME/ihep_yrli_1/BIGDATA/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
FFTW_INCL =
FFTW_LIBS = -lfftw3

DNEST_INCL  = -I /HOME/ihep_yrli_1/BIGDATA/soft/DNest/
DNEST_LIBS  = -L /HOME/ihep_yrli_1/BIGDATA/soft/DNest -ldnest
endif


EXEC     = recon
SRC      = ./
OBJS     = $(SRC)/main.o $(SRC)/allvars.o $(SRC)/recon.o $(SRC)/system.o $(SRC)/help.o \
           $(SRC)/read.o $(SRC)/psd.o
 
INCL     = Makefile proto.h allvars.h      

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(MPICHINCL) $(DNEST_INCL) $(FFTW_INCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(MPICHLIB) $(DNEST_LIBS) $(FFTW_LIBS)

$(EXEC):$(OBJS)
	cd $(SRC)
	$(CC) $(OBJS) $(LIBS) -o $@

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)