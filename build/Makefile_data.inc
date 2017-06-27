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
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to download data.  It is included from the main
#//!        makefile.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef DATA_DEF

SF_DATA_URL=http://nucnet-tools.sourceforge.net/data_pub

#===============================================================================
# Get data.
#===============================================================================

ifndef DATA_DIR
DATA_DIR=$(NUCNET_TARGET)/data_pub
endif

DATA_URL=$(SF_DATA_URL)/2013-04-19
DATA=nnt_data

data:
	wget $(DATA_URL)/$(DATA).tar.gz
	mkdir -p $(DATA_DIR)
	tar xz -C ${DATA_DIR} -f ${DATA}.tar.gz
	rm ${DATA}.tar.gz

#===============================================================================
# Get multi-zone data.
#===============================================================================

MULTI_ZONE_URL = http://www.webnucleo.org/experimental/mbradle/download/special
MULTI_ZONE_DATA = multi_zone_data
MULTI_ZONE_DIR = $(NUCNET_TARGET)/data_pub

multi_zone_data:
	wget ${MULTI_ZONE_URL}/${MULTI_ZONE_DATA}.tar.gz
	mkdir -p ${MULTI_ZONE_DIR}
	tar xz -C ${MULTI_ZONE_DIR}/ -f ${MULTI_ZONE_DATA}.tar.gz
	rm ${MULTI_ZONE_DATA}.tar.gz

#===============================================================================
# Get neutrino data.
#===============================================================================

NEUTRINO_URL = $(SF_DATA_URL)/2013-06-20
NEUTRINO_DATA = neutrino
NEUTRINO_DIR = $(DATA_DIR)

neutrino_data:
	wget ${NEUTRINO_URL}/${NEUTRINO_DATA}.tar.gz
	mkdir -p ${NEUTRINO_DIR}
	tar xz -C ${NEUTRINO_DIR} -f ${NEUTRINO_DATA}.tar.gz
	rm ${NEUTRINO_DATA}.tar.gz
 
#===============================================================================
# Get weak data.
#===============================================================================

WEAK_URL = http://www.webnucleo.org/experimental/mbradle/download/special
WEAK_DATA = ffn_data
WEAK_DIR = $(NUCNET_TARGET)/data_pub

weak_data:
	wget ${WEAK_URL}/${WEAK_DATA}.tar.gz
	mkdir -p ${WEAK_DIR}
	tar xz -C ${WEAK_DIR}/../ -f ${WEAK_DATA}.tar.gz
	rm ${WEAK_DATA}.tar.gz

#===============================================================================
# End include.
#===============================================================================

DATA_DEF = yes

endif
