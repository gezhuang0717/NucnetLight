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

#===============================================================================
# Edit the following lines to give the path to the wn_matrix codes, to choose
# the compiler (gcc or g++), and to choose whether you will check
# executable with Valgrind:
#===============================================================================

STATMECHSRCDIR = ../src/
GC=gcc
VALGRIND=yes

#===============================================================================
# End of editing
#===============================================================================

VPATH = $(STATMECHSRCDIR)

CINCLUDE= -I$(STATMECHSRCDIR) `xml2-config --cflags` `gsl-config --cflags`
CLIBS= `xml2-config --libs` `gsl-config --libs`

#===============================================================================
# Compiler flags.
#===============================================================================

CFLAGS1= -ansi -pedantic -Werror -Wall -W \
         -Wconversion -Wshadow \
         -Wpointer-arith -Wcast-qual -Wcast-align \
         -Wwrite-strings \
         -fshort-enums -fno-common -Dinline= -g

ifeq ($(GC), gcc)
  CFLAGS2= $(CFLAGS1) -Wmissing-prototypes -Wstrict-prototypes -Wnested-externs
else
  CFLAGS2= $(CFLAGS1)
endif

ifeq ($(VALGRIND), yes)
  CFLAGS= $(CFLAGS2) -O0
else
  CFLAGS= $(CFLAGS2) -O2
endif

CC=$(GC) $(CFLAGS) $(CINCLUDE)

#===============================================================================
# TMPDIR is the temporary directory for codes compilation, this is where
# object files are created.
#===============================================================================

TMPDIR = ./tmp/

TMP_DIR := $(shell mkdir -p $(TMPDIR))

#===============================================================================
# Objects.
#===============================================================================

OBJS =  $(TMPDIR)Libstatmech.o \
	$(TMPDIR)boson_ground_state_functions.o

$(OBJS): $(TMPDIR)%.o: %.c
	$(CC) -c -o $@ $<

#===============================================================================
# Compile all codes.
#===============================================================================

EXEC = \
  create_fermion \
  fermion_accuracy \
  pressure_integrand \
  thermodynamic \
  fermion_derivatives \
  heat_capacity_vs_rho \
  thomas_fermi \
  create_boson \
  boson_function \
  boson_derivatives \
  photons \
  bec_n \
  bec_T \

.PHONY all : $(EXEC)

#===============================================================================
# Compile codes.
#===============================================================================

$(EXEC) : $(OBJS)
	$(CC) -c -o $(TMPDIR)$@.o $@.c
	$(CC) $(OBJS) $(TMPDIR)$@.o $(CLIBS) -o $@

#===============================================================================
# Clean up.
#===============================================================================

.PHONY: clean cleanall

clean: 
	rm -fr $(TMPDIR)

cleanall: clean
	rm -f $(EXEC) *.exe
