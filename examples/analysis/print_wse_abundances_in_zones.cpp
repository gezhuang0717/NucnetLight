////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu. 
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
//! \brief Example code to print out wse abundances and mass fractions in
//!    zones selected by XPath.  Output is Z, A, abundance, mass fraction.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>
#include <Libnucnet.h>
#include <Libnuceq.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"
#include "user/aa522a25.h"
#include "user/weak_utilities.h"

#define D_SMALL   1.e-20        // The smallest abundance to print out.

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet *p_my_nucnet;
  Libnuceq *p_equil;
  Libnucnet__Zone * p_work_zone;
  gsl_vector * p_abundances;
  double d_abund;  

  if( argc != 3 && argc != 4 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file zone_xpath nuc_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zones\n\n"
    );
    fprintf(
      stderr,
      "  nuc_xpath = XPath expression to select species (optional)\n\n"
    );
    return EXIT_FAILURE;
  }

  p_my_nucnet = Libnucnet__new();

  if( argc == 3 )
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[1],
      NULL
    );
  else
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[1],
      argv[3]
    );

  Libnucnet__assignZoneDataFromXml(
    p_my_nucnet,
    argv[1],
    argv[2]
  );

  //============================================================================
  // Check zones.
  //============================================================================

  if( Libnucnet__getNumberOfZones( p_my_nucnet ) == 0 )
  {
    std::cerr << "No zones." << std::endl;
    return EXIT_FAILURE;
  }
 
  //============================================================================
  // Update with approximate rates.
  //============================================================================
  
  user::aa522a25__update_net(
    Libnucnet__getNet( p_my_nucnet )
  );
 
  //============================================================================
  // Set weak views in zones.
  //============================================================================
  
  user::set_weak_views_in_zones( p_my_nucnet );

  //============================================================================
  // Create a zone to store abundances temporarily.
  //============================================================================

  p_work_zone =
    Libnucnet__Zone__new( 
      Libnucnet__getNet( p_my_nucnet ),
      "work",
      NULL,
      NULL
    );

  //============================================================================
  // Create the species list.
  //============================================================================

  boost::ptr_list<nnt::Species> species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
    );

  //============================================================================
  // Create the zone list.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  boost::ptr_list<nnt::Zone> zone_list = nnt::make_zone_list( p_my_nucnet );

  //============================================================================
  // Create an equilibrium.
  //============================================================================

  p_equil =
    Libnuceq__new(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
    );

  //============================================================================
  // Iterate zones.
  //============================================================================

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    fprintf(
      stdout,
      "time = %.4e  t9 = %.4e  rho = %.4e\n",
      boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ),
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
    );

    //==========================================================================
    // Compute WSE. 
    //==========================================================================

    if( zone.getProperty( nnt::s_MU_NUE_KT ) != "-inf" )
    {

      Libnuceq__computeEquilibrium(
        p_equil, 
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
      );

    }
    else
    {

      user::register_my_rate_functions(
        Libnucnet__Net__getReac(
          Libnucnet__Zone__getNet(
            zone.getNucnetZone()
          )
        )
      );

      Libnuceq__computeEquilibrium(
        p_equil, 
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
      );

    }

    p_abundances = Libnuceq__getAbundances( p_equil );

    Libnucnet__Zone__updateAbundances( p_work_zone, p_abundances );

    gsl_vector_free( p_abundances );

    BOOST_FOREACH( nnt::Species species, species_list )
    {

      d_abund = 
         Libnucnet__Zone__getSpeciesAbundance(
           p_work_zone, 
           species.getNucnetSpecies()
         );

      if( d_abund > D_SMALL )
      {
        fprintf(
          stdout,
          "%5d  %5d  %.4e  %.4e\n",
          Libnucnet__Species__getZ( species.getNucnetSpecies() ),
          Libnucnet__Species__getA( species.getNucnetSpecies() ),
          d_abund,
          Libnucnet__Species__getA( species.getNucnetSpecies() ) * 
            d_abund
        );
      }

    }

    std::cout << std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnuceq__free( p_equil );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
