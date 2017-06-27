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
//! \brief A header file to define routines and structures for user-defined
//!        rate functions.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_RATE_MODIFIERS_H
#define USER_RATE_MODIFIERS_H

#include <iostream>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>

#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

//##############################################################################
// Structure.
//##############################################################################

struct view
{ 
  std::string id; 
  std::string nuc_xpath;
  std::string reac_xpath;
  double factor;

  view( 
    std::string &s_id,
    std::string &s_nuc_xpath,
    std::string &s_reac_xpath,
    double d_factor
  ) :
      id( s_id ),
      nuc_xpath( s_nuc_xpath ),
      reac_xpath( s_reac_xpath ),
      factor( d_factor ) {} 
}; 

struct change_nuc_xpath
{
  change_nuc_xpath(
    const std::string& new_xpath
  ):new_xpath(new_xpath){}

  void operator()(view& v)
  {
    v.nuc_xpath = new_xpath;
  }

  private:
    std::string new_xpath;

};

struct change_reac_xpath
{
  change_reac_xpath(
    const std::string& new_xpath
  ):new_xpath(new_xpath){}

  void operator()(view& v)
  {
    v.reac_xpath = new_xpath;
  }

  private:
    std::string new_xpath;

};

struct change_factor
{
  change_factor(
    double new_factor
  ):new_factor(new_factor){}

  void operator()(view& v)
  {
    v.factor = new_factor;
  }

  private:
    double new_factor;

};

typedef boost::multi_index::multi_index_container< 
  view, 
  boost::multi_index::indexed_by< 
    boost::multi_index::hashed_unique< 
      boost::multi_index::member< 
        view, std::string, &view::id
      > 
    >
  > 
> view_multi; 

//##############################################################################
// Prototypes.
//##############################################################################

void
modify_rates_in_view(
  nnt::Zone&,
  const char *,
  const char *,
  double
);

void
modify_rates( nnt::Zone& );

void
print_modified_reactions( nnt::Zone& );

int
modify_rates_views(
  const char *,
  const char *,
  const char *,
  const char *,
  view_multi *
);

void
modify_rates_for_reaction(
  nnt::Zone&,
  Libnucnet__Reaction *,
  double&,
  double&
);

} // namespace user

#endif // USER_RATE_MODIFIERS_H
