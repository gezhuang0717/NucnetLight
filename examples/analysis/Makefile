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
#//! \brief A makefile to generate analysis examples.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef NUCNET_TARGET
NUCNET_TARGET = ../..
endif

NNT_DIR = $(NUCNET_TARGET)/nnt
USER_DIR = $(NUCNET_TARGET)/user
BUILD_DIR = $(NUCNET_TARGET)/build
MY_USER_DIR = $(NUCNET_TARGET)/my_user

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

include $(BUILD_DIR)/Makefile

include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR)

#===============================================================================
# Objects.
#===============================================================================

ANALYSIS_OBJ = $(WN_OBJ)        \
               $(NNT_OBJ)	\
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

ANALYSIS_EXEC = compute_flows 			\
        compute_abundance_moment_in_zones 	\
        compute_Ycdot				\
        compute_zone_mu 			\
        compare_equil				\
        compute_species_flows 			\
        compute_nucleon_number_flows 	        \
        print_abundances_for_element		\
        print_mass_fractions_in_zones		\
        print_abundances_vs_nucleon_number	\
        print_largest_mass_fractions_in_zones	\
        print_wse_abundances_in_zones           \
        print_properties			\
        print_species_currents			\
        print_zone_abundances			\
        print_zone_properties			\
        print_muh				\
        compare_nucleon_equil			\

$(ANALYSIS_EXEC): $(ANALYSIS_OBJ)
	$(CC) $(ANALYSIS_OBJ) -o $(BINDIR)/$@ $@.cpp $(CLIBS)

.PHONY all_analysis: $(ANALYSIS_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_analysis cleanall_analysis

clean_analysis:
	rm -f $(ANALYSIS_OBJ)

cleanall_analysis: clean_analysis
	rm -f $(BINDIR)/$(ANALYSIS_EXEC) $(BINDIR)/$(ANALYSIS_EXEC).exe
