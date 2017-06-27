#--------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
# <file type="public">
#
#   <license>
#     This is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
# 
#     This software is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     Please see the src/README.txt file in this distribution for more
#     information.
#   </license>
#
#   <description>
#     <abstract>
#       Makefile for binding code.
#     </abstract>
#   </description>
#
#   <authors>
#     <current>
#       <author userid="mbradle" start_date="2010/12/07" />
#     </current>
#     <previous>
#     </previous>
#   </authors>
#
#   <compatibility>
#     Tested with GNU Make version 3.79.
#   </compatibility>
#
# </file>
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
# The following lines must be edited according to where you have
# placed your libnucnet and wn_matrix modules.  You may also
# edit to choose a different compiler (e.g., g++), to use Valgrind or
# not, or to profile:
#///////////////////////////////////////////////////////////////////////////////

GC=gcc

MODULES = ../../../../../../modules

LIBNUCNETDIR = $(MODULES)/libnucnet/0.8/src/
MATRIXSRCDIR = $(MODULES)/wn_matrix/0.13/src/
VALGRIND= yes
PROFILE= no

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

CINCLUDE= `xml2-config --cflags` `gsl-config --cflags` -I$(LIBNUCNETDIR) -I$(MATRIXSRCDIR)
CLIBS= `xml2-config --libs` `gsl-config --libs`

#===============================================================================
# Compiler flags.
#===============================================================================

CFLAGS1= -ansi -Werror -Wall -W -pedantic\
         -Wconversion -Wshadow \
         -Wpointer-arith -Wcast-qual -Wcast-align \
         -Wwrite-strings \
         -fshort-enums -fno-common -Dinline= -g \

ifeq ($(GC), gcc) 
	CFLAGS2= $(CFLAGS1) -Wmissing-prototypes -Wstrict-prototypes -Wnested-externs
else
	CFLAGS2= $(CFLAGS1)
endif

ifeq ($(VALGRIND), yes)
	CFLAGS3= $(CFLAGS2) -O0
else
	CFLAGS3= $(CFLAGS2) -O2
endif

ifeq ($(PROFILE), yes)
	CFLAGS= $(CFLAGS3) -pg
else
	CFLAGS= $(CFLAGS3)
endif

CC=$(GC) $(CFLAGS) $(CINCLUDE)

#===============================================================================
# TMPDIR is the temporary directory for codes compilation, this is where
# object files are created. 
#===============================================================================

TMPDIR = ./tmp/

MKTMP:
	mkdir -p tmp

#===============================================================================
# Compile matrix codes.
#===============================================================================

$(TMPDIR)WnMatrix.o: $(MATRIXSRCDIR)WnMatrix.c MKTMP
	$(CC) -c -o $(TMPDIR)WnMatrix.o $(MATRIXSRCDIR)WnMatrix.c

#===============================================================================
# Compile Libnucnet code.
#===============================================================================
$(TMPDIR)Libnucnet__Nuc.o: $(LIBNUCNETDIR)Libnucnet__Nuc.c MKTMP
	$(CC) -c -o $(TMPDIR)Libnucnet__Nuc.o $(LIBNUCNETDIR)Libnucnet__Nuc.c

#===============================================================================
# Compile main codes.
#===============================================================================
$(TMPDIR)binding.o: binding.c MKTMP
	$(CC) -c -o $(TMPDIR)binding.o binding.c

#--------------------------------------------------------------------------

OBJ =  $(TMPDIR)Libnucnet__Nuc.o   \
        $(TMPDIR)WnMatrix.o

#===============================================================================
# Compile code.
#===============================================================================

binding : $(TMPDIR)binding.o $(OBJ)
	$(CC) $(TMPDIR)binding.o $(OBJ) $(CLIBS) -o binding

#===============================================================================
# Clean up.
#===============================================================================

.PHONY: clean cleanall

clean: 
	rm -fr $(TMPDIR)

cleanall: clean
	rm -f binding

