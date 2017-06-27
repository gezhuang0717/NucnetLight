//////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file for the network limiter code.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_NETWORK_LIMITER_H
#define USER_NETWORK_LIMITER_H

#include <iostream>
#include <string>
#include <set>

#include <boost/unordered_set.hpp>

#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

//##############################################################################
// Prototypes.
//##############################################################################

void
limit_evolution_network( nnt::Zone&, double );

void
limit_evolution_network( nnt::Zone& );

void
zero_out_small_abundances( nnt::Zone&, double );

void
zero_out_small_rates( nnt::Zone&, double );

std::set<std::string>
get_isolated_species(
  Libnucnet__Net *,
  std::string,
  std::string
);

} // namespace user

#endif // USER_NETWORK_LIMITER_H
