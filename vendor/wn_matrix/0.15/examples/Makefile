#///////////////////////////////////////////////////////////////////////////////
# <file type="public">
#
#   <description>
#     <abstract>
#       Makefile for compiling the example codes.
#     </abstract>
#   </description>
#   <license>
#     This file was originally written by Bradley S. Meyer.
#
#     This is free software; you can redistribute it and#or modify it
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
#     copyright and license information.
#   </license>
#
# </file>
#///////////////////////////////////////////////////////////////////////////////

#===============================================================================
# Edit the following lines to give the path to the wn_matrix codes, to choose
# the compiler (gcc or g++), and to choose whether you will check
# executable with Valgrind:
#===============================================================================

MATRIXSRCDIR = ../src/
GC=gcc
VALGRIND=yes

#===============================================================================
# End of editing
#===============================================================================

CINCLUDE= -I$(MATRIXSRCDIR) `xml2-config --cflags` `gsl-config --cflags`
CLIBS= `xml2-config --libs` `gsl-config --libs`

#===============================================================================
# Compiler flags.
#===============================================================================

CFLAGS= -ansi -Werror -Wall -pedantic\
         -Wconversion -Wshadow \
         -Wpointer-arith -Wcast-qual \
         -Wwrite-strings \
         -fshort-enums -fno-common -Dinline= -g \
         -DWN_USE_MATRIX_LOOKUP

ifeq ($(GC), gcc) 
	CFLAGS+= -W -Wcast-align -Wmissing-prototypes -Wstrict-prototypes \
                 -Wnested-externs
endif

ifeq ($(GC), g++) 
	CFLAGS+= -W -Wcast-align
endif

ifeq ($(GC), icc)
	CFLAGS+= -wd9 -wd981 -wd1292 -wd1419 -wd10148 -wd10156
endif

ifeq ($(GC), icpc)
	CFLAGS+= -wd9 -wd981 -wd1292 -wd1419 -wd10148 -wd10156
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

#===============================================================================
# TMPDIR is the temporary directory for codes compilation, this is where
# object files are created.
#===============================================================================

TMPDIR = ./tmp/

MKTMP:
	mkdir -p $(TMPDIR)

#===============================================================================
# Compile matrix codes.
#===============================================================================

$(TMPDIR)WnMatrix.o: $(MATRIXSRCDIR)WnMatrix.c MKTMP
	$(CC) -c -o $(TMPDIR)WnMatrix.o $(MATRIXSRCDIR)WnMatrix.c

#===============================================================================
# Compile executables.
#===============================================================================

EXEC = create_matrix \
       print_matrix_to_file \
       convert_to_csr \
       convert_to_yale \
       insert_and_extract \
       ascii_matrix_to_xml \
       print_xml_matrix \
       print_row_or_column \
       insert_and_delete \
       get_transfer \
       get_copy_and_transpose \
       vector_to_xml \
       print_xml_vector \
       convert_coo_to_csr \
       convert_coo_to_yale \
       matrix_times_vector \
       solve_matrix_equation \
       get_arrow_matrix \
       solve_by_arrow

#===============================================================================
# Compile codes.
#===============================================================================

$(EXEC) : $(TMPDIR)WnMatrix.o
	$(CC) -c -o $(TMPDIR)$@.o $@.c
	$(CC) $(TMPDIR)WnMatrix.o $(TMPDIR)$@.o $(CLIBS) -o $@

.PHONY all : $(EXEC)
     
#===============================================================================
# Clean up.
#===============================================================================

.PHONY: clean cleanall

clean: 
	rm -fr $(TMPDIR)

cleanall: clean
	rm -fr $(EXEC) *.exe
