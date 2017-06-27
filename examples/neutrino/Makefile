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
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to generate neutrino rate examples.
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

CC_NEUTRINO = $(CC) -I $(NUCNET_TARGET)

#===============================================================================
# Objects.
#===============================================================================

NEUTRINO_OBJ = $(WN_OBJ)        \
               $(NNT_OBJ)       \
               $(USER_OBJ)

#===============================================================================
# Executables.
#===============================================================================

NEUTRINO_EXEC = print_neutrino_cross_sections \
                print_neutrino_rates

$(NEUTRINO_EXEC): $(NEUTRINO_OBJ)
	$(CC_NEUTRINO) $(NEUTRINO_OBJ) -o $(BINDIR)/$@ $@.cpp $(CLIBS)

.PHONY all_neutrino : $(NEUTRINO_EXEC)

#===============================================================================
# Clean up.
#===============================================================================

.PHONY: clean_neutrino cleanall_neutrino

clean_neutrino: 
	rm -f $(NEUTRINO_OBJ)

cleanall_neutrino: clean_neutrino
	rm -f $(BINDIR)/$(NEUTRINO_EXEC) $(BINDIR)/$(NEUTRINO_EXEC).exe

#===============================================================================
# End include.
#===============================================================================

NEUTRINO_DEF = yes
