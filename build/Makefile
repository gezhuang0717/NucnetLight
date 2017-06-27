#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2011-2015 Clemson University.
# 
#  This file was originally written by Bradley S. Meyer.
# 
#  This is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA
# 
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to generate dependencies.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef BUILD_DEF

ifndef GC
  GC=g++
endif

ifndef GF
  GF=gfortran
endif

LIBNUCNET_VERSION = 0.23
LIBNUCEQ_VERSION = 0.5
LIBSTATMECH_VERSION = 0.8
WNMATRIX_VERSION = 0.15
SPARSESOLVE_VERSION = 0.5

SOURCEFORGE_URL = http://sourceforge.net/projects/

LIBNUCNET_URL = $(SOURCEFORGE_URL)/libnucnet/files/libnucnet
LIBNUCEQ_URL = $(SOURCEFORGE_URL)/libnuceq/files/libnuceq
LIBSTATMECH_URL = $(SOURCEFORGE_URL)/libstatmech/files/libstatmech
WNMATRIX_URL = $(SOURCEFORGE_URL)/wnmatrix/files/wn_matrix
SPARSESOLVE_URL = $(SOURCEFORGE_URL)/wnsparsesolve/files/wn_sparse_solve

LIBNUCNETDIR = $(VENDORDIR)/libnucnet/$(LIBNUCNET_VERSION)/src
LIBNUCEQDIR = $(VENDORDIR)/libnuceq/$(LIBNUCEQ_VERSION)/src
LIBSTATMECHDIR = $(VENDORDIR)/libstatmech/$(LIBSTATMECH_VERSION)/src
WNMATRIXDIR = $(VENDORDIR)/wn_matrix/$(WNMATRIX_VERSION)/src
SPARSESOLVEDIR = $(VENDORDIR)/wn_sparse_solve/$(SPARSESOLVE_VERSION)/src

VALGRIND= yes
PROFILE= no
OPENMP= yes
NODEPRECATED= yes

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

include $(BUILD_DIR)/Makefile_data.inc

VPATH= .:$(WNMATRIXDIR):$(LIBNUCNETDIR):$(LIBNUCEQDIR):$(LIBSTATMECHDIR)

#===============================================================================
# C includes.  Add the Boost directory if desired.  This Makefile assumes
# BOOST_DIR is defined as an environment variable.  Otherwise you can set
# it above as BOOST_DIR = path_to_boost, where path_to_boost is the
# directory in which your local copy of Boost resides.
#===============================================================================

CINCLUDE= `xml2-config --cflags` \
          `gsl-config --cflags`  \
          -I$(LIBNUCNETDIR)      \
          -I$(WNMATRIXDIR)       \
          -I$(LIBSTATMECHDIR)    \
          -I$(LIBNUCEQDIR)       \
          -I$(NUCNET_TARGET)     \

ifdef BOOST_DIR
  CINCLUDE += -I$(BOOST_DIR)
endif

#===============================================================================
# Set the environment variable WN_DEBUG to turn on webnucleo code debugging.
#===============================================================================

ifdef WN_DEBUG
  CFLAGS += -DWN_DEBUG
endif

#===============================================================================
# C libraries.
#===============================================================================

CLIBS= `xml2-config --libs` `gsl-config --cflags --libs`

#===============================================================================
# Check for clang.
#===============================================================================

CLANG_TEST := $(shell $(GC) -v 2>&1 > /dev/null | grep clang)
ifneq "$(CLANG_TEST)" ""
    GC=clang++
endif

#===============================================================================
# Compiler flags.
#===============================================================================

CFLAGS= -Werror -Wall \
         -Wpointer-arith -Wcast-qual \
         -Wwrite-strings \
         -fshort-enums -fno-common -g \

ifneq ($(GC), clang++)
  CFLAGS+= -ansi -Wconversion
endif

FFLAGS= -g

FLIBS=

ifeq ($(GC), icpc)
	CFLAGS+= -wd9 -wd279 -wd981 -wd1292 -wd1418 -wd1419 -wd2259 -wd10148 -wd10156
endif

ifeq ($(VALGRIND), yes)
	CFLAGS+= -O0
else
	CFLAGS+= -O2
endif

ifeq ($(PROFILE), yes)
	CFLAGS+= -pg
	FFLAGS+= -pg
endif

ifeq ($(OPENMP), yes)
  ifeq ($(GC), h5c++)
        CFLAGS+= -fopenmp
        FFLAGS+= -fopenmp
  endif
  ifeq ($(GC), icpc)
	CFLAGS+= -openmp
	FFLAGS+= -openmp
  endif
  ifeq ($(GC), clang++)
	CFLAGS+= -openmp
	FFLAGS+= -openmp
  else
	CFLAGS+= -fopenmp
	FFLAGS+= -fopenmp
  endif
endif

ifeq ($(NODEPRECATED), yes)
	CFLAGS+= -Wno-deprecated
endif

ifeq ($(GF), ifort)
        GF+= -nofor-main
endif

CC=$(GC) $(CFLAGS) $(CINCLUDE)

FF=$(GF) $(CINCLUDE) $(FFLAGS)

#===============================================================================
# OBJDIR is the temporary directory for codes compilation, this is where
# object files are created. 
#===============================================================================

ifndef OBJDIR
OBJDIR = ../../obj
endif

OBJ_DIR := $(shell mkdir -p $(OBJDIR))

#===============================================================================
# BINDIR is the temporary directory for code executables.
#===============================================================================

ifndef BINDIR
BINDIR = .
endif

BIN_DIR := $(shell mkdir -p $(BINDIR))

#===============================================================================
# VENDORDIR is the directory for storing webnucleo codes.
#===============================================================================

ifndef VENDORDIR
VENDORDIR = ../../vendor
endif

VENDOR_DIR := $(shell mkdir -p $(VENDORDIR))

#===============================================================================
# SPARSKITDIR is the directory for storing Sparskit.
#===============================================================================

ifndef SPARSKITDIR
SPARSKITDIR = $(VENDORDIR)/SPARSKIT2
endif

SPARSKIT_DIR := $(shell mkdir -p $(SPARSKITDIR))

#===============================================================================
# Build dependencies.
#===============================================================================

$(VENDORDIR)/wn_matrix.tar.gz:
	wget $(WNMATRIX_URL)/wn_matrix_$(WNMATRIX_VERSION).tar.gz -O $@

$(OBJDIR)/WnMatrix.o: $(VENDORDIR)/wn_matrix.tar.gz 
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(WNMATRIXDIR)/WnMatrix.c -o $@

$(VENDORDIR)/libstatmech.tar.gz:
	wget $(LIBSTATMECH_URL)/libstatmech_$(LIBSTATMECH_VERSION).tar.gz -O $@

$(OBJDIR)/Libstatmech.o: $(VENDORDIR)/libstatmech.tar.gz 
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBSTATMECHDIR)/Libstatmech.c -o $@

$(VENDORDIR)/libnuceq.tar.gz:
	wget $(LIBNUCEQ_URL)/libnuceq_$(LIBNUCEQ_VERSION).tar.gz -O $@

$(OBJDIR)/Libnuceq.o: $(VENDORDIR)/libnuceq.tar.gz $(OBJDIR)/Libnucnet.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCEQDIR)/Libnuceq.c -o $@

$(VENDORDIR)/libnucnet.tar.gz:
	wget $(LIBNUCNET_URL)/libnucnet_$(LIBNUCNET_VERSION).tar.gz -O $@

$(OBJDIR)/Libnucnet__Nuc.o: $(VENDORDIR)/libnucnet.tar.gz $(OBJDIR)/WnMatrix.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCNETDIR)/Libnucnet__Nuc.c -o $@

$(OBJDIR)/Libnucnet__Reac.o: $(VENDORDIR)/libnucnet.tar.gz $(OBJDIR)/WnMatrix.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCNETDIR)/Libnucnet__Reac.c -o $@

$(OBJDIR)/Libnucnet.o: $(VENDORDIR)/libnucnet.tar.gz $(OBJDIR)/WnMatrix.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCNETDIR)/Libnucnet.c -o $@

#--------------------------------------------------------------------------

WN_OBJ =$(OBJDIR)/WnMatrix.o		\
	$(OBJDIR)/Libnucnet__Nuc.o	\
	$(OBJDIR)/Libnucnet__Reac.o	\
	$(OBJDIR)/Libnucnet.o		\
	$(OBJDIR)/Libstatmech.o		\
	$(OBJDIR)/Libnuceq.o

#--------------------------------------------------------------------------
# nnt and other codes
#--------------------------------------------------------------------------

NNT_OBJ = $(OBJDIR)/auxiliary.o			\
          $(OBJDIR)/iter.o			\
          $(OBJDIR)/math.o			\
          $(OBJDIR)/two_d_weak_rates.o		\
          $(OBJDIR)/weak_detailed_balance.o	\
          $(OBJDIR)/write_output_xml.o          \
          $(OBJDIR)/wrappers.o

$(NNT_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean clean_all clean_dist

clean:
	rm -fr $(OBJDIR)

clean_dist:
	rm -fr $(VENDORDIR)

clean_all: clean clean_dist

clean_data:
	rm -fr $(DATA_DIR)

#===============================================================================
# End include.
#===============================================================================

BUILD_DEF = yes

endif
