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
#//! \file
#//! \brief A makefile to generate weak rate examples.
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

#===============================================================================
# Includes.
#===============================================================================

include $(BUILD_DIR)/Makefile

CC_WEAK = $(CC) -I $(NUCNET_TARGET)

#===============================================================================
# Objects.
#===============================================================================

WEAK_OBJ = $(WN_OBJ)        \
           $(USER_OBJ)       \
           $(NNT_OBJ)

#===============================================================================
# Executables.
#===============================================================================

WEAK_EXEC = print_log_ft                             		\
            print_weak_forward_and_reverse           		\
            print_neutrino_energy_flows         		\
            print_neutrino_entropy_loss_rates         		\
            print_average_neutrino_energy         		\
            forward_and_reverse_vs_rhoe              		\
            forward_and_reverse_vs_t9                		\
            compute_yedot                            		\
            print_weak_approximate_rates             		\
            print_weak_approximate_average_neutrino_energies	\
            print_approximate_neutrino_entropy_loss_rates       \
            approximate_forward_and_reverse_vs_rhoe  		\
            approximate_forward_and_reverse_vs_t9    		\
            weak_xml_from_text					\
            weak_xml_from_text_with_nu

$(WEAK_EXEC): $(WEAK_OBJ)
	$(CC_WEAK) $(WEAK_OBJ) -o $(BINDIR)/$@ $@.cpp $(CLIBS)

.PHONY all_weak : $(WEAK_EXEC)

#===============================================================================
# Clean up.
#===============================================================================

.PHONY: clean_weak cleanall_weak

clean_weak: 
	rm -f $(WEAK_OBJ)

cleanall_weak: clean_weak
	rm -f $(BINDIR)/$(WEAK_EXEC) $(BINDIR)/$(WEAK_EXEC).exe

