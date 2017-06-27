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
#//! \brief A makefile to generate graph examples.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef NUCNET_TARGET
NUCNET_TARGET = ../..
endif

NNT_DIR = $(NUCNET_TARGET)/nnt
MY_USER_DIR = $(NUCNET_TARGET)/my_user
USER_DIR = $(NUCNET_TARGET)/user
BUILD_DIR = $(NUCNET_TARGET)/build

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

include $(BUILD_DIR)/Makefile

include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR)

#===============================================================================
# Objects.
#===============================================================================

HELPER_OBJ = $(OBJDIR)/graph_helper.o \
             $(OBJDIR)/color.o        \
             $(OBJDIR)/linestyle.o    \
             $(OBJDIR)/scaling.o      \

$(HELPER_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

GRAPH_OBJ = $(WN_OBJ)        \
            $(NNT_OBJ)	     \
            $(HELPER_OBJ)    \
            $(USER_OBJ)

ifdef NNT_MY_USER
  include $(MY_USER_DIR)/Makefile.inc
  VPATH += :$(MY_USER_DIR)
  GRAPH_OBJ += $(MY_USER_OBJ)
  CC += -DMY_USER
endif

#===============================================================================
# Add C++ 11 extensions for clang++.
#===============================================================================

ifeq ($(GC), clang++)
  CC += -std=c++11
endif

#===============================================================================
# Executables.
#===============================================================================

GRAPH_EXEC = net_view_graph			\
             zone_abundance_graph		\
             zone_current_graph			\
             zone_flow_graph			\
             zone_mu_graph  			\

$(GRAPH_EXEC): $(GRAPH_OBJ)
	$(CC) -c -o $(OBJDIR)/$@.o $@.cpp
	$(CC) $(GRAPH_OBJ) $(OBJDIR)/$@.o -o $(BINDIR)/$@ $(CLIBS)

.PHONY all_graph: $(GRAPH_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_graph cleanall_graph

clean_graph:
	rm -f $(GRAPH_OBJ)

cleanall_graph: clean_graph
	rm -f $(BINDIR)/$(GRAPH_EXEC) $(BINDIR)/$(GRAPH_EXEC).exe
