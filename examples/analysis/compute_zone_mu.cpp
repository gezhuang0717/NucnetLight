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
//! \brief Example code to compute mu/kT for species in zones.
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Libnucnet.h>

#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"

#include "user/coul_corr.h"
#include "user/weak_utilities.h"

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

  if( argc < 3 || argc > 4 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " compares network abundances to equilibria for selected zones" <<
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
    fprintf(
      stderr,
      "  cutoff = abundance cutoff for output (optional--default = 1.e-25)\n\n"
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
  typedef
    boost::tuple<std::string, unsigned int, unsigned int, double, double> tup;
  //user::coul_corr_data my_coul_corr_data;
  double d_cutoff = 1.e-25;

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
  // Check zones.
  //============================================================================

  if( Libnucnet__getNumberOfZones( p_my_nucnet ) == 0 )
  {
    std::cerr << "No zones." << std::endl;
    return EXIT_FAILURE;
  }
 
  //============================================================================
  // Abundance cutoff.
  //============================================================================

  if( argc == 4 ) d_cutoff = atof( argv[3] );

  //============================================================================
  // Iterate zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

    std::vector<tup> my_tuples;

    double d_t9 = boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) );

    double d_rho =
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );

    if( zone.hasProperty( nnt::s_TIME ) )
      std::cout <<
        boost::format( 
          "time(s) = %.5e  t9 = %.5e  rho(g/cc) = %.5e  Ye = %.4f\n"
        ) %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ) %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) %
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );
    else
      std::cout <<
        boost::format( 
          "t9 = %.5e  rho(g/cc) = %.5e  Ye = %.4f\n"
        ) %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) %
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );

    nnt::species_list_t species_list =
      nnt::make_species_list(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
      );
    
    double d_mass_p =
      Libnucnet__Species__getMassExcess(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          "h1"
        )
      );

    double d_mass_n =
      Libnucnet__Species__getMassExcess(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          "n"
        )
      );

    double d_mup =
      log(
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc(
              Libnucnet__Zone__getNet( zone.getNucnetZone() )
            ),
            "h1"
          )
        ) /
        Libnucnet__Species__computeQuantumAbundance(
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc(
              Libnucnet__Zone__getNet( zone.getNucnetZone() )
            ),
            "h1"
          ),
          d_t9,
          d_rho
        )
      );

    double d_mun =
      log(
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc(
              Libnucnet__Zone__getNet( zone.getNucnetZone() )
            ),
            "n"
          )
        ) /
        Libnucnet__Species__computeQuantumAbundance(
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc(
              Libnucnet__Zone__getNet( zone.getNucnetZone() )
            ),
            "n"
          ),
          d_t9,
          d_rho
        )
      );

    BOOST_FOREACH( nnt::Species species, species_list )
    {

      double d_abund =
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          species.getNucnetSpecies()
        );

      if( d_abund > d_cutoff )
      {

        my_tuples.push_back(
          boost::make_tuple(
	    Libnucnet__Species__getName( species.getNucnetSpecies() ),
	    Libnucnet__Species__getZ( species.getNucnetSpecies() ),
	    Libnucnet__Species__getA( species.getNucnetSpecies() ),
	    d_abund,
            log(
              d_abund /
              Libnucnet__Species__computeQuantumAbundance(
                species.getNucnetSpecies(),
                d_t9,
                d_rho
              )
            )
            -
            Libnucnet__Species__getZ( species.getNucnetSpecies() ) *
            (
              d_mup + ( d_mass_p / nnt::compute_kT_in_MeV( d_t9 ) )
            )
            -
            (
              Libnucnet__Species__getA( species.getNucnetSpecies() )
              -
              Libnucnet__Species__getZ( species.getNucnetSpecies() )
            ) *
            (
              d_mun + ( d_mass_n / nnt::compute_kT_in_MeV( d_t9 ) )
            )
            + Libnucnet__Species__getMassExcess( species.getNucnetSpecies() ) /
              nnt::compute_kT_in_MeV( d_t9 )
          )
        );

      }

    }

    std::cout << my_tuples.size() << std::endl;

    BOOST_FOREACH( tup my_tuple, my_tuples )
    {
      fprintf(
        stdout,
        "%5s  %3d  %3d  %.5e  %.5e\n",
        my_tuple.get<0>().c_str(),
        my_tuple.get<1>(),
        my_tuple.get<2>(),
        my_tuple.get<3>(),
        my_tuple.get<4>()
      );

    }

    std::cout << std::endl;

  }
            
  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
