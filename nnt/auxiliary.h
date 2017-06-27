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
//! \brief A header file for a variety of useful auxiliary utilities.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_AUX_H
#define NNT_AUX_H

#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_set.hpp>

#include <Libnucnet.h>
#include <Libstatmech.h>
#include <Libnuceq.h>

#include "nnt/string_defs.h"
#include "nnt/param_defs.h"
#include "nnt/math.h"
#include "nnt/wrappers.h"
#include "nnt/iter.h"

namespace nnt
{

//############################################################################
// Prototypes.
//############################################################################

int
print_abundance(
  Libnucnet__Species *,
  Libnucnet__Zone *
);

int
compare_reactions_by_string(
  const Libnucnet__Reaction *,
  const Libnucnet__Reaction *
);

double
compute_fermi_dirac_factor(
  double,
  double
);

double
compute_one_minus_fermi_dirac_factor(
  double,
  double
);

double
compute_kT_in_MeV(
  double
);

double
compute_electron_chemical_potential_kT(
  double,
  double,
  double
);

int
zone_compare_by_first_label(
  const Libnucnet__Zone *,
  const Libnucnet__Zone *
);

int
species_sort_by_z_then_a(
  const Libnucnet__Species *,
  const Libnucnet__Species *
);

int
species_sort_function(
  const Libnucnet__Species *,
  const Libnucnet__Species *
);

double
compute_reaction_nuclear_Qvalue(
  Libnucnet__Net *,
  Libnucnet__Reaction *,
  double
);

int
is_electron_capture_reaction(
  Libnucnet__Reaction *
);

int
is_positron_capture_reaction(
  Libnucnet__Reaction *
);

int
is_beta_plus_reaction(
  Libnucnet__Reaction *
);

int
is_beta_minus_reaction(
  Libnucnet__Reaction *
);

std::string *
create_nuc_xpath_from_set(
  Libnucnet__Net *,
  boost::unordered_set<std::string>
);

std::string *
create_reac_xpath_from_list(
  Libnucnet__Net *,
  boost::ptr_list<std::string>
);

size_t
reaction_element_count(
  Reaction&,
  std::string,
  std::string
);

gsl_vector *
get_new_gsl_vector_from_std_vector( const std::vector<double>& );

void
get_property_array_size(
  const char *, const char *, const char *, const char *, size_t *
);

void
get_property_array(
  const char *, const char *, const char *, const char *, gsl_vector *
);

void
get_property_matrix(
  const char *, const char *, const char *, const char *, gsl_matrix *
);

void
set_zone_abundances_to_equilibrium( nnt::Zone& );

Libnuceq *
get_zone_equilibrium( nnt::Zone& );

void
update_equilibrium_with_zone_data( Libnuceq *, nnt::Zone& );

int
get_clusters(
  Libnuceq__Cluster *,
  std::vector<Libnuceq__Cluster *> *
);

void
set_zone_cluster_data( nnt::Zone&, Libnuceq * );

gsl_vector *
get_rate_function_property_gsl_vector(
  Libnucnet__Reaction *,
  const char *
);

std::vector<Libnucnet__Species *>
get_nuc_anchors( Libnucnet__Nuc * );

std::vector<std::string> get_stable_species();

/**
 * \brief Return a vector giving the property in each zone.
 *
 * \param s_name A string giving the name of the property.
 * \param s_tag1 A string giving an 
 *                 extra tag associated with the property (optional).
 * \param s_tag2 A string giving a second extra tag associated with the
 *		      property (optional).
 * \return p_vector A vector containing the properties.
 */

template <typename T>
std::vector<T>
get_property_vector_from_zone_list(
  zone_list_t& zone_list,
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2
)
{

  std::vector<T> new_vector;

  BOOST_FOREACH( Zone zone, zone_list )
  {

    if(
        !Libnucnet__Zone__getProperty(
           zone.getNucnetZone(),
           s_name,
           s_tag1,
           s_tag2
        )
    )
    {

      if( s_tag1 && s_tag2 )
      {
        std::cerr <<
          "Not all zones have the property " <<
          s_name << " " << s_tag1 << " " << s_tag2 <<
          ", so can't create vector." <<
          std::endl;
      }
      else if( s_tag1 && !s_tag2 )
      {
        std::cerr <<
          "Not all zones have the property " <<
          s_name << " " << s_tag1 <<
          ", so can't create vector." <<
          std::endl;
      }
      else
      {
        std::cerr <<
          "Not all zones have the property " <<
          s_name <<
          ", so can't create vector." <<
          std::endl;
      }

      exit( EXIT_FAILURE );

    }

    try{
      new_vector.push_back(
        boost::lexical_cast<T>(
          Libnucnet__Zone__getProperty(
            zone.getNucnetZone(),
            s_name,
            s_tag1,
            s_tag2
          )
        )
      );
    }
    catch( const boost::bad_lexical_cast& e )
    {
      std::cerr << "Can't cast this property properly." << std::endl;
      throw e;
    }

  }

  return new_vector;

}
      
template <typename T>
std::vector<T>
get_property_vector_from_zone_list(
  zone_list_t& zone_list,
  const char * s_name,
  const char * s_tag1
)
{

  return
    get_property_vector_from_zone_list<T>(
      zone_list, s_name, s_tag1, NULL
    );

}

template <typename T>
std::vector<T>
get_property_vector_from_zone_list(
  zone_list_t& zone_list,
  const char * s_name
)
{

  return
    get_property_vector_from_zone_list<T>(
      zone_list, s_name, NULL, NULL
    );

}

} // namespace nnt

#endif // NNT_AUX_H
