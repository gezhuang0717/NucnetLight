////////////////////////////////////////////////////////////////////////////////
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
//! \file
//! \brief Example code to compute (n,gamma)-(gamma,n) equilibrium.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>
#include <Libnuceq.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

//##############################################################################
// print_species_abundances().
//##############################################################################

int
print_species_abundances(
  Libnuceq__Species * p_eq_species,
  nnt::Zone * p_zone
)
{

  Libnucnet__Species * p_species =
    Libnuceq__Species__getNucSpecies( p_eq_species );

  fprintf(
    stdout,
    "%d  %.8e  %.8e\n",
    Libnucnet__Species__getA( p_species ) -
      Libnucnet__Species__getZ( p_species ),
    Libnucnet__Zone__getSpeciesAbundance( p_zone->getNucnetZone(), p_species ),
    Libnuceq__Species__getAbundance(
      p_eq_species
    )
  );

  return 1;

}

//##############################################################################
// print_zone_abundances().
//##############################################################################

void
print_zone_abundances(
  nnt::Zone zone,
  Libnuceq * p_equil,
  unsigned int i_Z
)
{

  Libnuceq__Cluster *p_cluster;
  size_t i;
  gsl_vector * p_z_vector;
  char s_xpath[nnt::i_BUF_SIZE];

  //============================================================================
  // Set Ye.
  //============================================================================
  
  Libnuceq__setYe(
    p_equil,
    Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
  );

  //============================================================================
  // Create (n,g)-(g,n) equilibrium clusters.
  //============================================================================
  
  p_z_vector =
    Libnucnet__Zone__getSummedAbundances( zone.getNucnetZone(), "z" ); 

  for( i = 2; i < WnMatrix__get_gsl_vector_size( p_z_vector ); i++ )
  {
    sprintf( s_xpath, "[z = %lu]", (unsigned long) i );
    p_cluster = Libnuceq__getCluster( p_equil, s_xpath );
    if( !p_cluster )
      p_cluster = Libnuceq__newCluster( p_equil, s_xpath );
    Libnuceq__Cluster__updateConstraint(
      p_cluster,
      gsl_vector_get( p_z_vector, i )
    );
  }

  gsl_vector_free( p_z_vector );

  //============================================================================
  // Compute equilibrium.
  //============================================================================

  Libnuceq__computeEquilibrium(
    p_equil,
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
  );

  //============================================================================
  // Print diagnostics.
  //============================================================================
  
  sprintf( s_xpath, "[z = %d]", i_Z );
  
  fprintf(
    stdout,
    "For step number %s, T9 = %s, rho (g/cc) = %s\n",
    Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ),
    zone.getProperty( nnt::s_T9 ).c_str(),
    zone.getProperty( nnt::s_RHO ).c_str()
  );

  fprintf(
    stdout,
    "Z = %d, mu/kT = %g\n",
    i_Z,
    Libnuceq__Cluster__getMukT(
      Libnuceq__getCluster( p_equil, s_xpath )
    )
  );

  //============================================================================
  // Print abundances and return.
  //============================================================================

  Libnuceq__Cluster__iterateSpecies(
    Libnuceq__getCluster( p_equil, s_xpath ),
    (Libnuceq__Species__iterateFunction) print_species_abundances,
    &zone
  );

  fprintf( stdout, "\n" );

}

//##############################################################################
// main()
//##############################################################################

int main( int argc, char **argv ) {

  Libnucnet * p_my_nucnet;
  Libnuceq * p_equil;
  unsigned int i_Z;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 4 )
  {
      fprintf(
        stderr, "\nUsage: %s input_file zone_xpath Z\n", argv[0]
      );
      fprintf(
        stderr, "\n  input_file: name of input nuc xml file\n"
      );
      fprintf(
        stderr,
        "\n  zone_xpath: XPATH expression to get zones\n"
      );
      fprintf(
        stderr, "\n  Z: Z for abundances\n\n"
      );
      return EXIT_FAILURE;
   }

  //============================================================================
  // Get the input.
  //============================================================================

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, argv[2] );

  //============================================================================
  // Create the equilibrium.
  //============================================================================

  p_equil =
    Libnuceq__new( Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) );

  i_Z = (unsigned int) atoi( argv[3] );

  //============================================================================
  // Iterate the zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  boost::ptr_list<nnt::Zone> zone_list = nnt::make_zone_list( p_my_nucnet );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    print_zone_abundances(
      zone,
      p_equil,
      i_Z
    );

  }

  //============================================================================
  // Clean up.
  //============================================================================

  Libnuceq__free( p_equil );
  Libnucnet__free( p_my_nucnet );

  //============================================================================
  // Done!
  //============================================================================

  return EXIT_SUCCESS;

}

