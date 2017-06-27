#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2011-2013 Clemson University.
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
#//! \brief A makefile to generate user routine objects.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef USER_DEF

HD5_OBJ = $(OBJDIR)/hdf5_routines.o                \

HYDRO_OBJ = $(OBJDIR)/hydro.o		           \

ILU_OBJ = $(OBJDIR)/ilu_solvers.o                  \

MULTI_OBJ = $(OBJDIR)/multi_zone_utilities.o       \

SOLVE_OBJ = $(OBJDIR)/matrix_solver.o              \
            $(OBJDIR)/evolve.o			   \
            $(OBJDIR)/network_limiter.o  	   \

USER_OBJ = $(OBJDIR)/user_rate_functions.o         \
           $(OBJDIR)/flow_utilities.o	           \
           $(OBJDIR)/aa522a25.o			   \
           $(OBJDIR)/neutrino_rate_functions.o     \
           $(OBJDIR)/network_utilities.o	   \
           $(OBJDIR)/rate_modifiers.o              \
           $(OBJDIR)/screen.o                      \
           $(OBJDIR)/thermo.o                      \
           $(OBJDIR)/coul_corr.o                   \
           $(OBJDIR)/weak_utilities.o              \
           $(OBJDIR)/remove_duplicate.o

$(HYDRO_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

$(SOLVE_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

$(USER_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

$(HD5_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

$(ILU_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

$(MULTI_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

#===============================================================================
# Clean up.
#===============================================================================

.phony: clean_user clean_solver

clean_user:
	rm -f $(USER_OBJ)

clean_solver:
	rm -f $(OBJDIR)/matrix_solver.o

clean_hydro:
	rm -f $(OBJDIR)/hydro.o

#===============================================================================
# End include.
#===============================================================================

USER_DEF = yes

endif

