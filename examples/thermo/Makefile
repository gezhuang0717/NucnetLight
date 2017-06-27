#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2011-2014 Clemson University.
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
#/////////////////////////////////////////////////////////////////////////////*/

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to generate thermo examples.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef NUCNET_TARGET
NUCNET_TARGET = ../..
endif

NNT_DIR = $(NUCNET_TARGET)/nnt
USER_DIR = $(NUCNET_TARGET)/user
BUILD_DIR = $(NUCNET_TARGET)/build

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

include $(BUILD_DIR)/Makefile

include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR)

CC_THERMO = $(CC) -I $(NUCNET_TARGET)

#===============================================================================
# Objects.
#===============================================================================

THERMO_OBJ = $(WN_OBJ)        \
             $(USER_OBJ)      \
             $(NNT_OBJ) 

#===============================================================================
# Executables.
#===============================================================================

THERMO_EXEC = compute_temperature_from_entropy 	         \
              compute_density_from_entropy	         \
              compute_energy_generation_rate_per_nucleon \
              compute_sound_speed			 \
              compute_thermo_quantity			 \
              compute_sdot

$(THERMO_EXEC): $(THERMO_OBJ)
	$(CC_THERMO) $(THERMO_OBJ) -o $(BINDIR)/$@ $@.cpp $(CLIBS)

.PHONY all_thermo: $(THERMO_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_thermo cleanall_thermo

clean_thermo:
	rm -f $(THERMO_OBJ)

cleanall_thermo: clean_thermo
	rm -f $(BINDIR)/$(THERMO_EXEC) $(BINDIR)/$(THERMO_EXEC).exe
