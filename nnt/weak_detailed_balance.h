////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2012 Clemson University.
//
// This file was originally written by Bradley S. Meyer.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file for routines to compute weak detailed balance.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_WEAK_DETAILED_BALANCE_HPP
#define NNT_WEAK_DETAILED_BALANCE_HPP

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"

namespace nnt
{

double
compute_reverse_weak_rate_for_reaction(
  Libnucnet__Nuc *,
  Libnucnet__Reaction *,
  double,
  double,
  double,
  double,
  double
);

} // namespace nnt
    
#endif // NNT_WEAK_DETAILED_BALANCE_HPP
