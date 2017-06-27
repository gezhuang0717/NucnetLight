#///////////////////////////////////////////////////////////////////////////////
# <file type="public">
#
#   <description>
#     <abstract>
#       Makefile for compiling the example codes.
#     </abstract>
#   </description>
#   <license>
#     See the src/README.txt file in this distribution for copyright and license
#     information.
#   </license>
#
# </file>
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
# Edit the following lines to give the path to the wn_matrix codes, to choose
# the compiler (gcc or g++), the flag for compiling fortran with g++,
# and to choose whether you will check executable with Valgrind:
#///////////////////////////////////////////////////////////////////////////////

GC=gcc
FC=gfortran

AR= ar -rc

SPARSKITDIR = ../SPARSKIT2/
SPARSESOLVEDIR = ../src/
MATRIXDIR = ../../../wn_matrix/0.14/src/

VALGRIND=yes
PROFILE=no

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

#===============================================================================
# Include.
#===============================================================================

include Makefile.inc

#===============================================================================
# Some flags.
#===============================================================================

VPATH = $(MATRIXDIR):$(SPARSESOLVEDIR)

CINCLUDE= `xml2-config --cflags` `gsl-config --cflags` -I$(MATRIXDIR) -I$(SPARSESOLVEDIR)
CLIBS= `xml2-config --libs` `gsl-config --libs`
FLIBS= -L$(SPARSKITDIR) -lskit

#===============================================================================
# Compiler flags.
#===============================================================================

CFLAGS= -ansi -Werror -Wall -pedantic\
         -Wconversion -Wshadow \
         -Wpointer-arith -Wcast-qual \
         -Wwrite-strings \
         -fshort-enums -fno-common -Dinline= -g \

FFLAGS= -g -pedantic -Wall

ifeq ($(GC), gcc) 
	CFLAGS+= -W -Wcast-align -Wmissing-prototypes -Wstrict-prototypes \
                 -Wnested-externs
endif

ifeq ($(GC), g++) 
	CFLAGS+= -W -Wcast-align
	FFLAGS+= -lstdc++
endif

ifeq ($(GC), icc)
	CFLAGS+= -wd9 -wd981 -wd1292 -wd1419 -wd10148 -wd10156
endif

ifeq ($(GC), icpc)
	CFLAGS+= -wd9 -wd981 -wd1292 -wd1419 -wd10148 -wd10156
	FFLAGS+= -lstdc++
endif

ifeq ($(FC), ifort)
	FFLAGS+= -nofor-main
endif

ifeq ($(VALGRIND), yes)
	CFLAGS+= -O0
else
	CFLAGS+= -O2
endif

ifeq ($(PROFILE), yes)
	CFLAGS+= -pg
endif

CC=$(GC) $(CFLAGS) $(CINCLUDE)
FF=$(FC) $(FFLAGS)

#===============================================================================
# TMPDIR is the temporary directory for codes compilation, this is where
# object files are created. 
#===============================================================================

TMPDIR = ./tmp/

TMP_DIR := $(shell mkdir -p $(TMPDIR))

#===============================================================================
# Compile auxiliary C codes.
#===============================================================================

COBJ =  $(TMPDIR)WnMatrix.o                 \
	$(TMPDIR)WnSparseSolve.o            \
	$(TMPDIR)ilu_solvers.o              \
	$(TMPDIR)my_convergence_tester.o    \

$(COBJ): $(TMPDIR)%.o: %.c
	$(CC) -c -o $@ $<

#===============================================================================
# Compile unsupported sparse kit codes as modified.
#===============================================================================

FOBJ =  $(TMPDIR)blas1.o      \
	$(TMPDIR)itaux.o      \
	$(TMPDIR)exppro.o     \
	$(TMPDIR)phipro.o     \

$(FOBJ): $(TMPDIR)%.o: %.f
	$(FF) -c -o $@ $<

#===============================================================================
# Compile example codes.
#===============================================================================

EXECS = basic_solve_sparse                                 \
	solve_sparse_with_user_supplied_convergence_tester \
	solve_sparse_with_user_supplied_preconditioner     \
	solve_sparse_with_ilu_preconditioner               \
	solve_sparse_with_exppro                           \
	solve_sparse_with_phipro                           \

$(EXECS) : $(COBJ) $(FOBJ) lskit
	$(CC) -c -o $(TMPDIR)$@.o $@.c
	$(FF) $(COBJ) $(FOBJ) $(TMPDIR)$@.o $(CLIBS) $(FLIBS) -o $@

.PHONY all: $(EXECS)

#===============================================================================
# Clean up.
#===============================================================================

.PHONY: clean cleanall clean_sparskit

clean:
	rm -fr $(TMPDIR)

cleanall: clean clean_lskit
	rm -f $(EXECS) *.exe

clean_sparskit:
	rm -fr $(SPARSKITDIR)

