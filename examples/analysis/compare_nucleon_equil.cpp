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
//! \brief Example code to compare (n,g) or (p,g) equilibrium abundances to
//!        those in a network file.
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <map>
#include <omp.h>

#include <boost/format.hpp>

#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"

#include "user/network_utilities.h"

#define S_CLUSTER_XPATH  "[z > 2]"

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__NucView * p_nuc_view;
  Libnucnet__Species * p_nucleon;
  double d_min_abund;
  std::string s_nuc, s_cluster, s_format;

  //============================================================================
  // Check input.
  //============================================================================

   if ( argc < 4 || argc > 6 )
   {
      fprintf(
        stderr,
        "\nUsage: %s xml_filename nucleon zone_xpath min_abund nuc_xpath\n\n",
        argv[0]
      );
      fprintf(
        stderr, "  xml_filename = input xml filename\n\n"
      );
      fprintf(
        stderr, "  nucleon = n for (n,g) or p for (p,g) equilibrium\n\n"
      );
      fprintf(
        stderr, "  zone_xpath = XPath to select zones for flows\n\n"
      );
      fprintf(
        stderr, "  min_abund = minimum abundance for print out (optional but required if nuc_xpath present--default = 1.e-25)\n\n"
      );
      fprintf(
        stderr, "  nuc_xpath = XPath to select species for print out (optional)\n\n"
      );
      return EXIT_FAILURE;
   }

  //============================================================================
  // Check nucleon input.
  //============================================================================

  if( std::string( argv[2] ) != "n" && std::string( argv[2] ) != "p" )
  {
    std::cerr << "Wrong nucleon type (should be n or p)." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Read file and exit if not present.
  //============================================================================

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, argv[3] );

  if( !p_my_nucnet ) {
    fprintf( stderr, "Input data not read!\n" );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Get species for printout.
  //============================================================================

  if( argc != 6 )
    p_nuc_view =
      Libnucnet__NucView__new(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
        ""
      );
  else
    p_nuc_view =
      Libnucnet__NucView__new(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
        argv[5]
      );

  nnt::species_list_t species_list =
    nnt::make_species_list( Libnucnet__NucView__getNuc( p_nuc_view ) );

  //============================================================================
  // Minimum abundance for print out.
  //============================================================================
  
  if( argc < 5 )
    d_min_abund = 1.e-25;
  else
    d_min_abund = atof( argv[4] );
  
  //============================================================================
  // Printout formats
  //============================================================================
  
  boost::format fmt1( "time(s) = %12.5e t9 = %12.5e  rho(g/cc) = %12.5e" );

  boost::format fmt2( "Ye = %10.4e" );

  boost::format fmt3( "Y(nucleon) = %10.4e" );

  boost::format fmt4( "Y_eq(nucleon) = %10.4e" );

  boost::format fmt5( "Yh = %10.4e" );

  boost::format fmt6( "mu / kT = %10.4e" );

  boost::format fmt7( "mu_eq / kT = %10.4e" );

  boost::format fmt8( "%5d %5d %12.4e %12.4e" );
  
  if( std::string( argv[2] ) == "n" )
  {
    s_nuc = "n";
    s_cluster = "z";
    s_format = "[z = %d]";
  }
  else
  {
    s_nuc = "h1";
    s_cluster = "n";
    s_format = "[a - z = %d]";
  }

  boost::format fmt_nuc( s_format );

  //============================================================================
  // Get the nucleon
  //============================================================================

  p_nucleon =
    Libnucnet__Nuc__getSpeciesByName(
      Libnucnet__Net__getNuc(
        Libnucnet__getNet( p_my_nucnet )
      ),
      s_nuc.c_str()
    );

  //===========================================================================
  // Set zones and equilibrium vector.
  //==========================================================================*/

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  std::vector<nnt::Zone> zones;
  std::vector<Libnuceq *> equils;

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

    zones.push_back( zone );

    equils.push_back(
      Libnuceq__new(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
      )
    );

  }

  //============================================================================
  // Set the equilibria and zone view.
  //============================================================================

  for(
     unsigned int i_nuc = 2;
     i_nuc <=
       Libnucnet__Nuc__getLargestNucleonNumber(
         Libnucnet__Net__getNuc(
           Libnucnet__getNet( p_my_nucnet )
         ),
         s_cluster.c_str()
       );
     i_nuc++
  )
  {
    Libnuceq__newCluster( equils[0], (fmt_nuc % i_nuc).str().c_str() );
    zones[0].getNetView( S_CLUSTER_XPATH, "" );
  }

  for( size_t i = 1; i < equils.size(); i++ )
  {
    Libnuceq__copy_clusters( equils[i], equils[0] );
    Libnucnet__Zone__copy_net_views(
      zones[i].getNucnetZone(),
      zones[0].getNucnetZone()
    );
  }

  for( size_t i = 0; i < equils.size(); i++ )
  {

    gsl_vector * p_abundances =
      Libnucnet__Zone__getSummedAbundances(
        zones[i].getNucnetZone(),
        s_cluster.c_str()
      );

    //==========================================================================
    // Set clusters.
    //==========================================================================

    for(
      unsigned int i_nuc = 2;
      i_nuc <=
        Libnucnet__Nuc__getLargestNucleonNumber(
          Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
          s_cluster.c_str()
        );
      i_nuc++
    )
    {
      Libnuceq__Cluster__updateConstraint(
        Libnuceq__getCluster( equils[i], (fmt_nuc % i_nuc).str().c_str() ),
        gsl_vector_get( p_abundances, i_nuc )
      );
    }

    gsl_vector_free( p_abundances );

    //==========================================================================
    // Set Ye.
    //==========================================================================

    Libnuceq__setYe(
      equils[i],
      Libnucnet__Zone__computeZMoment(
        zones[i].getNucnetZone(),
        1
      )
    );

  }

  //============================================================================
  // Loop on zones to compute the equilibrium.
  //============================================================================

  #pragma omp parallel for schedule( dynamic, 1 )
  for( size_t i = 0; i < zones.size(); i++ )
  {

    Libnuceq__computeEquilibrium(
      equils[i],
      boost::lexical_cast<double>( zones[i].getProperty( nnt::s_T9 ) ),
      boost::lexical_cast<double>( zones[i].getProperty( nnt::s_RHO ) )
    );

  }

  //============================================================================
  // Print conditions.
  //============================================================================

  for( size_t i = 0; i < zones.size(); i++ )
  {
    if( zones[i].hasProperty( nnt::s_TIME ) )
    {

      fmt1 %
        boost::lexical_cast<double>( zones[i].getProperty( nnt::s_TIME ) ) %
        boost::lexical_cast<double>( zones[i].getProperty( nnt::s_T9 ) ) %
        boost::lexical_cast<double>( zones[i].getProperty( nnt::s_RHO ) );

      std::cout << fmt1.str() << std::endl;

      fmt2 %
        Libnucnet__Zone__computeZMoment(
          zones[i].getNucnetZone(),
          1
        );

      fmt3 %
        Libnucnet__Zone__getSpeciesAbundance(
          zones[i].getNucnetZone(),
          p_nucleon
        );

      fmt4 %
        Libnuceq__Species__getAbundance(
          Libnuceq__getSpeciesByName(
            equils[i],
            Libnucnet__Species__getName( p_nucleon )
          )
        );

      fmt5 %
        user::compute_cluster_abundance_moment(
          zones[i],
          S_CLUSTER_XPATH,
          "a",
          0
        );

      fmt6 %
        log(
          Libnucnet__Zone__getSpeciesAbundance(
            zones[i].getNucnetZone(),
            p_nucleon
          ) /
          Libnucnet__Species__computeQuantumAbundance(
            p_nucleon,
            boost::lexical_cast<double>( zones[i].getProperty( nnt::s_T9 ) ),
            boost::lexical_cast<double>( zones[i].getProperty( nnt::s_RHO ) )
          )
        );

      if( s_nuc == "n" )
        fmt7 % Libnuceq__getMunkT( equils[i] );
      else
        fmt7 % Libnuceq__getMupkT( equils[i] );

      std::cout << fmt2.str() << std::endl;

      std::cout << fmt3.str() << std::endl;

      std::cout << fmt4.str() << std::endl;

      std::cout << fmt5.str() << std::endl;

      std::cout << fmt6.str() << std::endl;

      std::cout << fmt7.str() << std::endl;

    }
    else
    {
      std::cerr << "Invalid network file." << std::endl;
      return EXIT_FAILURE;
    }

    //==========================================================================
    // Set vector.
    //==========================================================================

    std::vector<
      boost::tuple< Libnucnet__Species *, double, double >
    > abunds;

    BOOST_FOREACH( nnt::Species species, species_list )
    {

      Libnucnet__Species * p_species = species.getNucnetSpecies();

      double d_abund =
        Libnucnet__Zone__getSpeciesAbundance(
          zones[i].getNucnetZone(),
          p_species
        );

      double d_eq_abund =
        Libnuceq__Species__getAbundance(
          Libnuceq__getSpeciesByName(
            equils[i],
            Libnucnet__Species__getName( species.getNucnetSpecies() )
          )
        );

      if( d_abund >= d_min_abund || d_eq_abund >= d_min_abund )
      {

        abunds.push_back(
          boost::make_tuple( p_species, d_abund, d_eq_abund )
        );

      }

    }

    //==========================================================================
    // Print species abundances.
    //==========================================================================

    std::cout <<
      "N = " << abunds.size() << std::endl;

    fprintf( stdout, "  Z     A     Network      equil\n" );
    printf( "===== =====   ==========   =================\n");

    for( size_t j = 0; j < abunds.size(); j++ )
    {

      fmt8 % 
        Libnucnet__Species__getZ( abunds[j].get<0>() ) %
        Libnucnet__Species__getA( abunds[j].get<0>() ) %
        abunds[j].get<1>() %
        abunds[j].get<2>();

        std::cout << fmt8.str() << std::endl;

    }
          
    std::cout << std::endl;

  }

  //===========================================================================
  // Clean up and exit.
  //==========================================================================*/

  BOOST_FOREACH( Libnuceq * p_equil, equils ) {Libnuceq__free( p_equil );}
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
