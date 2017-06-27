////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2013-2014 Clemson University.
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
//
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file multi_zone_utilities.h
//! \brief A header file to define useful multi-zone utilities.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_ZONE_UTILITIES_H
#define MULTI_ZONE_UTILITIES_H

//##############################################################################
// Includes.
//##############################################################################

#include <omp.h>

#include <gsl/gsl_blas.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/evolve.h"
#include "user/network_limiter.h"

namespace user
{

//##############################################################################
// Defines.
//##############################################################################

#define D_DT_E_CMP   1.e-8
#define D_DT_E_STEP  2.
#define D_DT_E_REG   1.15

#define IT_MAX  10

//##############################################################################
// Types and classes.
//##############################################################################

struct nucnet_t { typedef boost::graph_property_tag kind; };
typedef boost::property< nucnet_t, Libnucnet * > graph_p;

class Zone_Vertex {
  public:
    Zone_Vertex(){};
    ~Zone_Vertex(){};
    void setNucnetZone( Libnucnet__Zone * p_zone ) { pZone = p_zone; }
    Libnucnet__Zone * getNucnetZone() { return pZone; }

  protected:
    Libnucnet__Zone * pZone;

};

class Zone_Link {
  public:
    Zone_Link(){};
    ~Zone_Link(){};
    void setWeight( double d_weight ){ dWeight = d_weight; }
    void setExtraString(
      std::string s_extra_string
    ){ sExtraString = s_extra_string; }
    double getWeight() const { return dWeight; }
    std::string getExtraString() const { return sExtraString; }
    boost::unordered_map<std::string, double> speciesYieldMap;

  protected:
    double dWeight;
    std::string sExtraString;

};

/**
 * A Boost graph type representing a zone link graph.
 */

typedef
boost::adjacency_list<
  boost::listS,
  boost::vecS,
  boost::bidirectionalS,
  Zone_Vertex,
  Zone_Link,
  graph_p
> zone_link_graph_t;

/**
 * A vertex hash entry.
 */

class Vertex_Entry {
  public:
    Vertex_Entry(
      const std::string &_s1,
      const std::string &_s2,
      const std::string &_s3,
      zone_link_graph_t::vertex_descriptor &_v
    ) : s1( _s1 ), s2( _s2 ), s3( _s3), v( _v ) {}
    zone_link_graph_t::vertex_descriptor getVertexDescriptor()
      const { return v; }
    std::string get_label1() const { return s1; };
    std::string get_label2() const { return s2; };
    std::string get_label3() const { return s3; };

  private:
    std::string s1;
    std::string s2;
    std::string s3;
    zone_link_graph_t::vertex_descriptor v;

};

typedef boost::multi_index_container<
  Vertex_Entry,
  boost::multi_index::indexed_by<
    boost::multi_index::hashed_unique< 
      boost::multi_index::composite_key<
        Vertex_Entry,
        boost::multi_index::const_mem_fun<
          Vertex_Entry,
          std::string,
          &Vertex_Entry::get_label1
        >,
        boost::multi_index::const_mem_fun<
          Vertex_Entry,
          std::string,
          &Vertex_Entry::get_label2
        >,
        boost::multi_index::const_mem_fun<
          Vertex_Entry,
          std::string,
          &Vertex_Entry::get_label3
        >
      >
    >
  >
> vertex_multi_index;

//##############################################################################
// Prototypes.
//##############################################################################

zone_link_graph_t
make_zone_link_graph( std::vector<nnt::Zone>& );

void
fill_multi_zone_vertex_hash( zone_link_graph_t&, vertex_multi_index& );

gsl_vector *
create_full_vector( std::vector<nnt::Zone>& );

void
update_from_full_vector(
  std::vector<nnt::Zone>&,
  gsl_vector *,
  std::string
);

std::pair< std::vector<WnMatrix *>, std::vector<gsl_vector *> >
get_zone_jacobians_and_rhs_vectors(
  std::vector<nnt::Zone>&
);

std::vector<WnMatrix *>
get_zone_jacobian_matrices(
  std::vector<nnt::Zone>&
);

double
get_new_timestep_from_zones(
  std::vector<nnt::Zone>&, double, double, double, double
);

void
limit_zone_networks( std::vector<nnt::Zone>&, double );

void
limit_zone_networks( std::vector<nnt::Zone>& );

int
check_multi_zone_mass_fractions(
  std::vector<nnt::Zone>&,
  double
);

void
normalize_multi_zone_abundances(
  std::vector<nnt::Zone>&
);

void
add_views_to_zones(
  std::vector<nnt::Zone>&,
  std::vector<std::pair<std::string,std::string> >&
);

gsl_vector *
get_mass_numbers_array( Libnucnet__Nuc * );

void
print_links( zone_link_graph_t& );

int
zone_compare_by_labels(
  const Libnucnet__Zone *,
  const Libnucnet__Zone *
);

std::vector<nnt::Zone> get_vector_of_zones( Libnucnet * );

void
printout_abundances_in_zones(
  std::vector<nnt::Zone>&,
  double,
  double
);

int
multi_mass_calculation_check( std::vector<nnt::Zone>& );

zone_link_graph_t::vertex_descriptor
get_vertex_from_multi_zone_hash(
  vertex_multi_index& vm,
  std::string,
  std::string,
  std::string
);

void
multi_zone_zero_small_abundances(
  std::vector<nnt::Zone>&,
  double
);

} // namespace user

#include "user/multi_zone_utilities_impl.hpp"

#endif // MULTI_ZONE_UTILITIES_H
