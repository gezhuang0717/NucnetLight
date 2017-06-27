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
//! \brief Useful multi-zone utilities.
//!
////////////////////////////////////////////////////////////////////////////////

#include "user/multi_zone_utilities.h"

namespace user
{

//##############################################################################
// make_zone_link_graph().
//##############################################################################

/**
 * \brief Routine to create a zone_link graph from a set of zones.  The
 *        user must separately set the links.
 *
 * \param zones The vector of zones.
 * \return A zone_link graph (just the vertices).
 * 
*/

zone_link_graph_t
make_zone_link_graph( std::vector<nnt::Zone>& zones )
{

  zone_link_graph_t g;
  zone_link_graph_t::vertex_descriptor v;

  for( size_t i = 0; i < zones.size(); i++ )
  {

    v = boost::add_vertex( g );

    g[v].setNucnetZone( zones[i].getNucnetZone() );

  }

  return g;

}

//##############################################################################
// fill_multi_zone_vertex_hash().
//##############################################################################

void
fill_multi_zone_vertex_hash(
  zone_link_graph_t& g,
  vertex_multi_index& my_vertices
)
{

  zone_link_graph_t::vertex_iterator vi, vi_end;

  boost::tie( vi, vi_end ) = boost::vertices( g );

  for( ; vi != vi_end; vi++ )
  {
    zone_link_graph_t::vertex_descriptor v = *vi;
    my_vertices.insert(
      Vertex_Entry(
        std::string( Libnucnet__Zone__getLabel( g[v].getNucnetZone(), 1 ) ),
        std::string( Libnucnet__Zone__getLabel( g[v].getNucnetZone(), 2 ) ),
        std::string( Libnucnet__Zone__getLabel( g[v].getNucnetZone(), 3 ) ),
        v
      )
    );
  }

}

//##############################################################################
// get_new_timestep_from_zones().
//##############################################################################

/**
 * \brief Routine to determine the next timestep for a set of zones
 *        based on changes to abundances.
 *
 * \param zones The vector of zones.
 * \param d_dt_old The previous timestep.
 * \param d_regt The time change regulator.
 * \param d_regy The abundance change regulator.
 * \param d_ymin The minimum abundance to check for change.
 * \return The new timestep.
 * 
*/

double
get_new_timestep_from_zones(
  std::vector<nnt::Zone>& zones,
  double d_dt_old,
  double d_regt,
  double d_regy,
  double d_ymin
)
{

  double d_dt, d_dt_check;

  d_dt = ( 1. + d_regt ) * d_dt_old;

  BOOST_FOREACH( nnt::Zone zone, zones )
  {

    d_dt_check = d_dt_old;

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt_check,
      d_regt,
      d_regy,
      d_ymin
    );

    if( d_dt_check < d_dt ) d_dt = d_dt_check;

  }

  return d_dt;

}

//##############################################################################
// limit_zone_networks().
//##############################################################################

void
limit_zone_networks(
  std::vector<nnt::Zone>& zones,
  double d_cutoff
)
{

  #pragma omp parallel for schedule( dynamic, 1 )
    for( size_t i = 0; i < zones.size(); i++ )
    {
      limit_evolution_network( zones[i], d_cutoff );
    }

}

void
limit_zone_networks(
  std::vector<nnt::Zone>& zones
)
{

  #pragma omp parallel for schedule( dynamic, 1 )
    for( size_t i = 0; i < zones.size(); i++ )
    {
      limit_evolution_network( zones[i] );
    }

}

//##############################################################################
// check_multi_zone_mass_fractions().
//##############################################################################

/**
 * \brief Routine to compare mass fractions in zones to threshold value.
 *
 * \param zones The vector of zones.
 * \param d_threshold The threshold value.
 * \return Routine returns 1 (true) if all zones have a sum of mass fractions
 *         within the threshold and 0 (false) if not.
 * 
*/


int
check_multi_zone_mass_fractions(
  std::vector<nnt::Zone>& zones,
  double d_threshold
)
{

  std::vector<double> dxsum( zones.size() );

  #pragma omp parallel for schedule( dynamic, 1 )
    for( size_t i = 0; i < zones.size(); i++ )
    {
      dxsum[i] =
        1. - Libnucnet__Zone__computeAMoment( zones[i].getNucnetZone(), 1 );
    }

  if( *(std::max_element( dxsum.begin(), dxsum.end() ) ) > d_threshold )
    return 0;
  else
    return 1;

}

//##############################################################################
// normalize_multi_zone_abundances().
//##############################################################################

/**
 * \brief Routine to normalize abundances in zones.
 *
 * \param zones The vector of zones.
 * 
*/

void
normalize_multi_zone_abundances(
  std::vector<nnt::Zone>& zones
)
{

  std::vector<double> dxsum( zones.size() );

  #pragma omp parallel for schedule( dynamic, 1 )
    for( size_t i = 0; i < zones.size(); i++ )
    {
      zones[i].normalizeAbundances();
    }

}

//##############################################################################
// create_full_vector().
//##############################################################################

gsl_vector *
create_full_vector( std::vector<nnt::Zone>& zones )
{

  gsl_vector_view view;
  gsl_vector * p_mass_numbers = NULL;
  size_t i_species, i_offset, i_delta;
  int is_multi_mass;

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  is_multi_mass = multi_mass_calculation_check( zones );

  if( is_multi_mass )
  {
    i_offset = i_species + 1;
    i_delta = 1;
    p_mass_numbers =
      get_mass_numbers_array(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
        )
      );
  }
  else
  {
    i_offset = i_species;
    i_delta = 0;
  }

  gsl_vector * p_vector = gsl_vector_alloc( zones.size() * i_offset );

  for( size_t i = 0; i < zones.size(); i++ )
  {

    if( is_multi_mass )
    {
      gsl_vector_set(
        p_vector,
        i * i_offset,
        boost::lexical_cast<double>( zones[i].getProperty( nnt::s_ZONE_MASS ) )
      );
    }

    gsl_vector * p_abunds =
      Libnucnet__Zone__getAbundances( zones[i].getNucnetZone() );

    view =
      gsl_vector_subvector(
        p_vector, 
        i * i_offset + i_delta,
        i_species
      );

    if( is_multi_mass )
    {
      if( !p_mass_numbers )
      {
        std::cerr << "Missing mass numbers vector because a zone is" <<
           " missing its zone mass." << std::endl;
        exit( EXIT_FAILURE );
      }
      gsl_vector_scale(
        p_abunds,
        boost::lexical_cast<double>( zones[i].getProperty( nnt::s_ZONE_MASS ) )
      );
      gsl_vector_mul( p_abunds, p_mass_numbers );
    }

    gsl_vector_memcpy( &view.vector, p_abunds );

    gsl_vector_free( p_abunds );

  }

  if( p_mass_numbers) gsl_vector_free( p_mass_numbers );

  return p_vector;

}
    
//##############################################################################
// update_from_full_vector().
//##############################################################################

void
update_from_full_vector(
  std::vector<nnt::Zone>& zones,
  gsl_vector * p_vector,
  std::string s_type
)
{

  gsl_vector_view view;
  gsl_vector * p_mass_numbers = NULL, * p_abunds;
  size_t i_species, i_offset, i_delta;
  int is_multi_mass;

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  p_abunds = gsl_vector_alloc( i_species );

  is_multi_mass = multi_mass_calculation_check( zones );

  if( is_multi_mass )
  {
    i_offset = i_species + 1;
    i_delta = 1;
    p_mass_numbers =
      get_mass_numbers_array(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
        )
      );
  }
  else
  {
    i_offset = i_species;
    i_delta = 0;
  }

  for( size_t i = 0; i < zones.size(); i++ )
  {

    view =
      gsl_vector_subvector(
        p_vector, 
        i * i_offset + i_delta,
        i_species
      );

    gsl_vector_memcpy( p_abunds, &view.vector ); 

    if( is_multi_mass )
    {
      if( s_type == "abundances" )
      {
	zones[i].updateProperty(
	  nnt::s_ZONE_MASS,
	  boost::lexical_cast<std::string>(
	    gsl_vector_get(
	      p_vector,
              i * i_offset
	    )
	  )
	);
      }
      else if( s_type == "abundance changes" )
      {
	zones[i].updateProperty(
	  nnt::s_ZONE_MASS_CHANGE,
	  boost::lexical_cast<std::string>(
	    gsl_vector_get(
	      p_vector,
              i * i_offset
	    )
	  )
	);
      }

      double d_mass =
        boost::lexical_cast<double>( zones[i].getProperty( nnt::s_ZONE_MASS ) );

      if( WnMatrix__value_is_zero( d_mass ) )
      {
        gsl_vector_scale( p_abunds, 0. );
      }
      else
      {
        gsl_vector_scale( p_abunds, 1. / d_mass );
        gsl_vector_div( p_abunds, p_mass_numbers );
      }

    }
      
    if( s_type == "abundances" )
    {
      Libnucnet__Zone__updateAbundances(
        zones[i].getNucnetZone(),
        p_abunds
      );
    }
    else if( s_type == "abundance changes" )
    {
      Libnucnet__Zone__updateAbundanceChanges(
        zones[i].getNucnetZone(),
        p_abunds
      );
    }
    else
    {
      std::cerr << "No such vector type." << std::endl;
      exit( EXIT_FAILURE );
    }

  }

  if( p_mass_numbers) gsl_vector_free( p_mass_numbers );
  gsl_vector_free( p_abunds );

}

//##############################################################################
// get_zone_jacobians_and_rhs_vectors().
//##############################################################################

std::pair< std::vector<WnMatrix *>, std::vector<gsl_vector *> >
get_zone_jacobians_and_rhs_vectors( std::vector<nnt::Zone>& zones )
{

  std::pair< std::vector<WnMatrix *>, std::vector<gsl_vector *> > my_pair;
  gsl_vector * p_mass_numbers = NULL;
  int is_multi_mass;

  my_pair.first.assign( zones.size(), NULL );
  my_pair.second.assign( zones.size(), NULL );

  is_multi_mass = multi_mass_calculation_check( zones );

  if( is_multi_mass )
  {
    p_mass_numbers =
      get_mass_numbers_array(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
        )
      );
  }

  #pragma omp parallel for schedule( dynamic, 1 )
    for( size_t i = 0; i < zones.size(); i++ )
    {

      boost::tie( my_pair.first[i], my_pair.second[i] ) =
	get_evolution_matrix_and_vector( zones[i] );

      Libnucnet__Zone__clearRates( zones[i].getNucnetZone() );

      if( is_multi_mass )
      {
        gsl_vector_scale(
          my_pair.second[i],
          boost::lexical_cast<double>(
            zones[i].getProperty( nnt::s_ZONE_MASS )
          )
        );
        gsl_vector_mul( my_pair.second[i], p_mass_numbers );
      }

    }

  if( p_mass_numbers ) gsl_vector_free( p_mass_numbers );

  return my_pair;

}

//##############################################################################
// get_zone_jacobian_matrices().
//##############################################################################

std::vector<WnMatrix *>
get_zone_jacobian_matrices( std::vector<nnt::Zone>& zones )
{

  std::vector<WnMatrix *> matrices;

  matrices.assign( zones.size(), NULL );

  #pragma omp parallel for schedule( dynamic, 1 )
    for( size_t i = 0; i < zones.size(); i++ )
    {

      matrices[i] = get_evolution_matrix( zones[i] );

      Libnucnet__Zone__clearRates( zones[i].getNucnetZone() );

      WnMatrix__scaleMatrix( matrices[i], -1. );

    }

  return matrices;

}

//##############################################################################
// get_mass_numbers_array().
//##############################################################################

gsl_vector *
get_mass_numbers_array( Libnucnet__Nuc * p_nuc )
{

  gsl_vector * p_mass_numbers_array =
    gsl_vector_alloc( Libnucnet__Nuc__getNumberOfSpecies( p_nuc ) );

  nnt::species_list_t species_list = nnt::make_species_list( p_nuc );

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    gsl_vector_set(
      p_mass_numbers_array,
      Libnucnet__Species__getIndex( species.getNucnetSpecies() ),
      (double) Libnucnet__Species__getA( species.getNucnetSpecies() )
    );
  }

  return p_mass_numbers_array;

}

//##############################################################################
// add_views_to_zones().
//##############################################################################

void
add_views_to_zones(
  std::vector<nnt::Zone>& zones,
  std::vector<std::pair<std::string,std::string> >& views
)
{

  typedef std::pair<std::string,std::string> view_t;

  BOOST_FOREACH( view_t v, views )
  {
    zones[0].getNetView( v.first.c_str(), v.second.c_str() );
  }

  #pragma omp parallel for schedule( dynamic, 1 )
    for( size_t i = 1; i < zones.size(); i++ )
    {

      Libnucnet__Zone__copy_net_views(
	zones[i].getNucnetZone(),
	zones[0].getNucnetZone()
      );

    }

}

//##############################################################################
// print_links().
//##############################################################################

void
print_links( zone_link_graph_t& g )
{

  boost::property_map<zone_link_graph_t, boost::vertex_index_t>::type
    id = boost::get( boost::vertex_index, g );

  zone_link_graph_t::edge_iterator ei, ei_end;

  boost::tie( ei, ei_end ) = boost::edges( g );

  for( ; ei != ei_end; ei++ )
  {
    
    zone_link_graph_t::vertex_descriptor u = boost::source( *ei, g );
    zone_link_graph_t::vertex_descriptor v = boost::target( *ei, g );

    std::cout <<
      "Zone " <<
      id[u] <<
      " (" << 
      Libnucnet__Zone__getLabel( g[u].getNucnetZone(), 1 ) <<
      ", " <<
      Libnucnet__Zone__getLabel( g[u].getNucnetZone(), 2 ) <<
      ", " <<
      Libnucnet__Zone__getLabel( g[u].getNucnetZone(), 3 ) <<
      ") --> " <<
      "Zone " <<
      id[v] <<
      " (" << 
      Libnucnet__Zone__getLabel( g[v].getNucnetZone(), 1 ) <<
      ", " <<
      Libnucnet__Zone__getLabel( g[v].getNucnetZone(), 2 ) <<
      ", " <<
      Libnucnet__Zone__getLabel( g[v].getNucnetZone(), 3 ) <<
      "): Weight = " <<
      g[*ei].getWeight() << std::endl;

  }

}

//############################################################################
// zone_compare_by_labels()
//############################################################################

/**
  \brief Zone comparison function which sorts the zones by
	  the alphabetic value associated with their first label,
          then by the numerical value of the second.

  \param p_zone1 The first zone.
  \param p_zone2 The second zone.
  \return 1 if zone 1 is after zone 2, -1 if zone 1 is before
	  zone 2, and 0 if the two zones are the same.
*/

int
zone_compare_by_labels(
  const Libnucnet__Zone *p_zone1,
  const Libnucnet__Zone *p_zone2
)
{

  int i_comp =
    strcmp(
      Libnucnet__Zone__getLabel( p_zone1, 1 ),
      Libnucnet__Zone__getLabel( p_zone2, 1 )
    );

  if( i_comp < 0 )
     return -1;
  else if( i_comp > 0 )
     return 1;
  else
    if( 
       atof( Libnucnet__Zone__getLabel( p_zone1, 2 ) ) <
       atof( Libnucnet__Zone__getLabel( p_zone2, 2 ) )
    )
      return -1;
    else
      return 1;
	  
}

//##############################################################################
// get_vertex_from_multi_zone_hash().
//##############################################################################

/**
  \brief Routine to retrieve a vertex from a vertex hash.

  \param vm The hash.
  \param s1 The first zone label.
  \param s2 The second zone label.
  \param s3 The third zone label.
  \return Routine returns a vertex descriptor corresponding to the zone.
*/

zone_link_graph_t::vertex_descriptor
get_vertex_from_multi_zone_hash(
  vertex_multi_index& vm,
  std::string s1,
  std::string s2,
  std::string s3
)
{

  vertex_multi_index::iterator it = vm.find( boost::make_tuple( s1, s2, s3 ) );

  if( it == vm.end() )
  {
    return zone_link_graph_t::null_vertex();
  }
  else
  {
    return it->getVertexDescriptor();
  }

}

//############################################################################
// get_vector_of_zones().
//############################################################################

std::vector<nnt::Zone>
get_vector_of_zones( Libnucnet * p_nucnet )
{

  std::vector<nnt::Zone> zones;

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_nucnet );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
    zones.push_back( zone );
  }

  return zones;

}

//##############################################################################
// printout_abundances_in_zones().
//##############################################################################

void
printout_abundances_in_zones(
  std::vector<nnt::Zone>& zones,
  double d_t,
  double d_dt
)
{

/**
 * \brief Routine to print abundances in zones to the screen.
 *
 * \param zones The vector of zones.
 * \param d_t The current time (s).
 * \param d_dt The current timestep (s).
 * 
*/

  double d_abund, d_change;

  for( size_t i = 0; i < zones.size(); i++ )
  {

    std::cout << "Zone: (" <<
      Libnucnet__Zone__getLabel( zones[i].getNucnetZone(), 1 ) << 
      ", " <<
      Libnucnet__Zone__getLabel( zones[i].getNucnetZone(), 2 ) << 
      ", " <<
      Libnucnet__Zone__getLabel( zones[i].getNucnetZone(), 3 ) << 
      ")" << std::endl;

    nnt::species_list_t species_list =
      nnt::make_species_list(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zones[i].getNucnetZone() )
        )
      );

    BOOST_FOREACH( nnt::Species species, species_list )
    {

      d_abund =
        Libnucnet__Zone__getSpeciesAbundance(
          zones[i].getNucnetZone(),
          species.getNucnetSpecies()
        );

      if( d_abund > 1.e-20 )
      {

	d_change =
	  Libnucnet__Zone__getSpeciesAbundanceChange(
	    zones[i].getNucnetZone(),
	    species.getNucnetSpecies()
	  );

        std::cout <<
          Libnucnet__Species__getName( species.getNucnetSpecies() ) <<
          "   " <<
          d_abund * Libnucnet__Species__getA( species.getNucnetSpecies() ) <<
          "   " <<
          d_change * Libnucnet__Species__getA( species.getNucnetSpecies() ) <<
          "   " <<
          fabs( d_change / d_abund ) <<
          std::endl;

       }

     }

     std::cout << " 1 - Xsum = " <<
       1 -
       Libnucnet__Zone__computeAMoment(
         zones[i].getNucnetZone(), 1
       ) << std::endl;

     std::cout << std::endl;

  }

}

//##############################################################################
// multi_mass_calculation_check().
//##############################################################################

int
multi_mass_calculation_check(
  std::vector<nnt::Zone>& zones
)
{

  int is_multi_mass = 0;

  BOOST_FOREACH( nnt::Zone zone, zones )
  {
    if( zone.hasProperty( nnt::s_ZONE_MASS ) )
    {
      is_multi_mass = 1;
    }
    else
    {
      if( is_multi_mass == 1 )
      {
        std::cerr <<
          "Not all zones for multi-mass calculation have zone mass." <<
          std::endl;
        exit( EXIT_FAILURE );
      }
    }
  }

  return is_multi_mass;

}

//##############################################################################
// multi_zone_zero_small_abundances().
//##############################################################################

void
multi_zone_zero_small_abundances(
  std::vector<nnt::Zone>& zones,
  double d_cutoff
)
{

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  #pragma omp parallel for schedule( dynamic, 1 )
    for( size_t i = 0; i < zones.size(); i++ )
    {

      zero_out_small_abundances( zones[i], d_cutoff );

    }

}

} // namespace user
