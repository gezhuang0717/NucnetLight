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
//! \brief Example code to print out the heavy nuclide chemical potential
//!     in zones from a network xml file.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>
#include <Libnuceq.h>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/iter.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " my_output.xml \"[position() >= last() - 10]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 3 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " computes quantities related to heavy nuclei for selected zones" <<
      std::endl;
    fprintf(
      stderr,
      "\nUsage: %s xml_file zone_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zones of interest\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );

  }

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__NucView * p_my_nuc_view;
  Libnuceq * p_my_equil;
  Libnuceq__Cluster * p_my_cluster;
  double d_yh;
  const char S_XPATH[] = "[z >= 6]";

  //============================================================================
  // Check input.
  //============================================================================

  check_input( argc, argv );

  //============================================================================
  // Read input data.
  //============================================================================

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      argv[2]
    );

  //============================================================================
  // Create equilibrium and cluster.
  //============================================================================

  p_my_equil =
    Libnuceq__new(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
    );

  p_my_cluster =
    Libnuceq__newCluster(
      p_my_equil,
      S_XPATH
    );

  //============================================================================
  // Get the species list.  Use a view to select only the relevant species.
  //============================================================================

  p_my_nuc_view =
    Libnucnet__NucView__new(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      S_XPATH
    );

  boost::ptr_list<nnt::Species> my_species_list =
    nnt::make_species_list( Libnucnet__NucView__getNuc( p_my_nuc_view ) );

  Libnucnet__NucView__free( p_my_nuc_view );

  //============================================================================
  // Iterate zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  boost::ptr_list<nnt::Zone> my_zone_list = nnt::make_zone_list( p_my_nucnet );

  BOOST_FOREACH( nnt::Zone zone, my_zone_list )
  {

    //--------------------------------------------------------------------------
    // Set the Ye for the equilibrium.
    //--------------------------------------------------------------------------

    Libnuceq__setYe(
      p_my_equil,
      Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
    );

    //--------------------------------------------------------------------------
    // Compute yh and update the cluster constraint.
    //--------------------------------------------------------------------------

    d_yh = 0;

    BOOST_FOREACH( nnt::Species species, my_species_list )
    {
      d_yh +=
	Libnucnet__Zone__getSpeciesAbundance(
	  zone.getNucnetZone(),
	  species.getNucnetSpecies()
	);
    }

    Libnuceq__Cluster__updateConstraint(
      p_my_cluster,
      d_yh
    );

    //--------------------------------------------------------------------------
    // Compute equilibrium and print out.
    //--------------------------------------------------------------------------

    Libnuceq__computeEquilibrium(
      p_my_equil,
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
    );

    std::cout <<
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) << "  " <<
      zone.getProperty( nnt::s_TIME ) << "  " <<
      zone.getProperty( nnt::s_T9 ) << "  " <<
      zone.getProperty( nnt::s_RHO ) << "  " <<
      d_yh << "  " <<
      Libnuceq__Cluster__getMukT( p_my_cluster ) <<
      std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnuceq__free( p_my_equil );

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
