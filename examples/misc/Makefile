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
#//! \brief A makefile to generate miscellaneous examples.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef NUCNET_TARGET
NUCNET_TARGET = ../..
endif

NNT_DIR = $(NUCNET_TARGET)/nnt
BUILD_DIR = $(NUCNET_TARGET)/build
USER_DIR = $(NUCNET_TARGET)/user
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

MISC_OBJ = $(WN_OBJ)        \
           $(USER_OBJ)      \
           $(NNT_OBJ)

#===============================================================================
# my_user routines, if desired.
#===============================================================================

ifdef NNT_MY_USER
  include $(MY_USER_DIR)/Makefile.inc
  VPATH += :$(MY_USER_DIR)
  MISC_OBJ += $(MY_USER_OBJ)
  CC += -DMY_USER
endif

#===============================================================================
# Executables.
#===============================================================================

MISC_EXEC = create_zone_xml_from_text 		\
            create_equil_zones_from_network     \
            compare_reaction_files_at_t9 	\
            extract_zones_to_xml		\
            select_network 			\
            merge_nuc 				\
            diff_current 			\
            create_multi_zone_xml_from_text	\
            test_versions 			\
            print_rates_by_xpath		\
            compute_separation_energies		\
            reactions_latex			\
            add_properties_to_zone_xml 		\
            extend_partf    			\
            update_mass_excesses_from_text 	\

$(MISC_EXEC): $(MISC_OBJ)
	$(CC) $(MISC_OBJ) -o $(BINDIR)/$@ $@.cpp $(CLIBS)

.PHONY all_misc: $(MISC_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_misc cleanall_misc

clean_misc:
	rm -f $(MISC_OBJ)

cleanall_misc: clean_misc
	rm -f $(BINDIR)/$(MISC_EXEC) $(BINDIR)/$(MISC_EXEC).exe
