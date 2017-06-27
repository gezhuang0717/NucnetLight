#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2013 Clemson University.
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
#//! \brief A makefile to generate codes for the multi-zone project.
#//!
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
# Here are lines to be edited, if desired.
#///////////////////////////////////////////////////////////////////////////////

NUCNET_TARGET = ../..
VENDORDIR = $(NUCNET_TARGET)/vendor
OBJDIR = $(NUCNET_TARGET)/obj

NNT_DIR = $(NUCNET_TARGET)/nnt
USER_DIR = $(NUCNET_TARGET)/user
BUILD_DIR = $(NUCNET_TARGET)/build

GC = h5c++

FC = gfortran

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

#===============================================================================
# Includes.
#===============================================================================

include $(BUILD_DIR)/Makefile
include $(BUILD_DIR)/Makefile.sparse
include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR)

#===============================================================================
# Debugging, if desired.
#===============================================================================

ifdef DEBUG
  CFLAGS += -DDEBUG
  FFLAGS += -DDEBUG
endif

#===============================================================================
# Objects.
#===============================================================================

MULTI_ZONE_OBJ = $(WN_OBJ)          \
                 $(MY_USER_OBJ)     \
                 $(NNT_OBJ)         \
                 $(HD5_OBJ)         \
                 $(SOLVE_OBJ)       \
                 $(USER_OBJ)        \
                 $(ILU_OBJ)         \
                 $(MULTI_OBJ)       \
                 $(SP_OBJ)          \

CFLAGS += -DSPARSKIT2 -lgfortran
MULTI_ZONE_DEP= sparse
FLIBS= -L$(SPARSKITDIR) -lskit

MULTI_ZONE_DEP += $(MULTI_ZONE_OBJ)

#===============================================================================
# Executables.
#===============================================================================

MULTI_ZONE_EXEC = multi_zone_network multi_mass exp_multi_mass

$(MULTI_ZONE_EXEC): $(MULTI_ZONE_DEP)
	$(CC) -c -o $(OBJDIR)/$@.o $@.cpp
	$(CC) $(MULTI_ZONE_OBJ) -o $(BINDIR)/$@ $(OBJDIR)/$@.o $(CLIBS) $(FLIBS)

.PHONY all_multi_zone: $(MULTI_ZONE_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_multi_zone cleanall_multi_zone

clean_multi_zone:
	rm -f $(MULTI_ZONE_OBJ)

cleanall_multi_zone: clean_multi_zone
	rm -f $(BINDIR)/$(MULTI_ZONE_EXEC) $(BINDIR)/$(MULTI_ZONE_EXEC).exe
