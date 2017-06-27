#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2013-2014 Clemson University.
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
#//! \file
#//! \brief A makefile to generate hdf5 examples.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef NUCNET_TARGET
NUCNET_TARGET = ../..
endif

NNT_DIR = $(NUCNET_TARGET)/nnt
USER_DIR = $(NUCNET_TARGET)/user
BUILD_DIR = $(NUCNET_TARGET)/build
MY_USER_DIR = $(NUCNET_TARGET)/my_user

GC=h5c++

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

include $(BUILD_DIR)/Makefile

include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR)

#===============================================================================
# HDF5
#===============================================================================

CC += -DHDF5

#===============================================================================
# Set the hydro.  NNT_HYDRO_CODE is an environment variable.  In a bash shell,
# set this by typing at the command line, for example,
# 'export NNT_HYDRO_CODE=wind'.  This approach replaces the deprecated
# approach in which you set, for example, HYDRO_CODE=wind within this Makefile.
#===============================================================================

ifdef NNT_HYDRO_CODE
  ifeq ($(NNT_HYDRO_CODE), exponential_expansion)
        CFLAGS += -DHYDRO_EXP_EXPANSION
  else ifeq ($(NNT_HYDRO_CODE), traj)
        CFLAGS += -DHYDRO_TRAJ
  else ifeq ($(NNT_HYDRO_CODE), wind)
        CFLAGS += -DHYDRO_WIND
  else ifeq ($(NNT_HYDRO_CODE), shock)
        CFLAGS += -DHYDRO_SHOCK
  endif
else
        CFLAGS += -DHYDRO_EXP_EXPANSION
endif

#===============================================================================
# Objects.
#===============================================================================

HDF5_OBJ =  $(WN_OBJ)		\
            $(NNT_OBJ)		\
            $(HD5_OBJ)		\
            $(HYDRO_OBJ)	\
            $(SOLVE_OBJ)	\
            $(USER_OBJ)

#===============================================================================
# my_user routines, if desired.
#===============================================================================

ifdef NNT_MY_USER
  include $(MY_USER_DIR)/Makefile.inc
  VPATH += :$(MY_USER_DIR)
  ANALYSIS_OBJ += $(MY_USER_OBJ)
  CC += -DMY_USER
endif

#===============================================================================
# Executables.
#===============================================================================

HDF5_EXEC = write_zones_to_hdf5 			\
           combine_xml_to_hdf5			\
           print_zone_mass_fractions_in_steps	\
           print_zone_mass_fractions_for_step	\
           print_step_mass_fractions_in_zones	\
           print_step_properties_in_zones	\
           print_zone_properties_in_steps	\
           single_zone_hdf5			\

$(HDF5_EXEC): $(HDF5_OBJ)
	$(CC) -c -o $(OBJDIR)/$@.o $@.cpp
	$(CC) $(HDF5_OBJ) -o $(BINDIR)/$@ $(OBJDIR)/$@.o $(CLIBS)

.PHONY all_hdf5: $(HDF5_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_hdf5 cleanall_hdf5

clean_hdf5:
	rm -f $(HDF5_OBJ)

cleanall_hdf5: clean_hdf5
	rm -f $(BINDIR)/$(HDF5_EXEC) $(BINDIR)/$(HDF5_EXEC).exe
