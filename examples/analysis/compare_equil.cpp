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
//! \brief Example code to compare various equilibria to a network calculation.
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Libnucnet.h>
#include <Libnuceq.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"
#include "user/aa522a25.h"
#include "user/coul_corr.h"

#include "user/weak_utilities.h"
#include "user/network_utilities.h"

#define S_XPATH   "[z > 2]"

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
      "\nUsage: %s xml_file zone_xpath wse\n\n",
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
      "  wse = also compute wse (enter any value)\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );

  }

}

//##############################################################################
// coul_wse_function().
//##############################################################################

double
coul_wse_function(
  Libnuceq * p_equil,
  user::coul_corr_data * p_data
)
{

  double d_result =
    user::my_coulomb_correction(
      Libnucnet__Nuc__getSpeciesByName(
        Libnuceq__getNuc( p_equil ),
        "h1"
      ),
      Libnuceq__getT9( p_equil ),
      Libnuceq__getRho( p_equil ),
      Libnuceq__computeZMoment( p_equil, 1 ),
      p_data
    );

  return d_result;

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnuceq *p_equil;
  Libnuceq__Cluster *p_cluster;
  nnt::Zone work_zone;
  gsl_vector * p_abundances, * p_z_wse = NULL, * p_z_nse, *p_z_qse, * p_z_net;
  size_t i;
  double d_ye_root, d_ye, d_yh, d_yh_nse, d_yh_wse;
  user::coul_corr_data my_coul_corr_data;

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
  // Set weak views in zones.
  //============================================================================
  
  user::set_weak_views_in_zones( p_my_nucnet );

  //============================================================================
  // Create a zone to store abundances temporarily.  Get a net view for
  // copying to other zones.
  //============================================================================

  work_zone.setNucnetZone(
    Libnucnet__Zone__new( 
      Libnucnet__getNet( p_my_nucnet ),
      "work",
      NULL,
      NULL
    )
  );

  work_zone.getNetView( S_XPATH, "" );

  //============================================================================
  // Create an equilibrium.
  //============================================================================

  p_equil =
    Libnuceq__new(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
    );

  //============================================================================
  // Output largest Z.
  //============================================================================

  fprintf(
    stdout,
    "Largest Z = %d\n",
    Libnucnet__Nuc__getLargestNucleonNumber(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      "z"
    )
  );

  //============================================================================
  // Iterate zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

    Libnucnet__Zone__copy_net_views(
      zone.getNucnetZone(),
      work_zone.getNucnetZone()
    );

    p_z_net = Libnucnet__Zone__getSummedAbundances( zone.getNucnetZone(), "z" );

    d_ye = Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );

    d_yh =
      user::compute_cluster_abundance_moment(
        zone,
        S_XPATH,
        "a",
        0
      );

    //==========================================================================
    // Add Coulomb correction if present.  Otherwise clear it.
    //==========================================================================
    
    if( zone.hasProperty( nnt::s_USE_SCREENING) )
    {

      my_coul_corr_data = user::get_coulomb_corr_data();

      Libnuceq__setNseCorrectionFactorFunction(
        p_equil,
        (Libnucnet__Species__nseCorrectionFactorFunction)
          user::my_coulomb_correction,
        &my_coul_corr_data
      );

    }
    else
    {
      Libnuceq__clearNseCorrectionFactorFunction( p_equil );
    }

    //==========================================================================
    // Compute NSE.
    //==========================================================================

    Libnuceq__setYe(
      p_equil,
      Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
    );

    Libnuceq__computeEquilibrium(
      p_equil, 
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
    );

    p_abundances = Libnuceq__getAbundances( p_equil );

    Libnucnet__Zone__updateAbundances(
      work_zone.getNucnetZone(),
      p_abundances
    );

    p_z_nse =
      Libnucnet__Zone__getSummedAbundances( work_zone.getNucnetZone(), "z" );

    gsl_vector_free( p_abundances );

    d_yh_nse =
      user::compute_cluster_abundance_moment(
        work_zone,
        S_XPATH,
        "a",
        0
      );

    //==========================================================================
    // Compute QSE.
    //==========================================================================

    p_cluster = Libnuceq__newCluster( p_equil, S_XPATH );

    Libnuceq__Cluster__updateConstraint(
      p_cluster,
      user::compute_cluster_abundance_moment(
        zone,
        S_XPATH,
        "a",
        0
      )
    );

    Libnuceq__computeEquilibrium(
      p_equil, 
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
    );

    p_abundances = Libnuceq__getAbundances( p_equil );

    Libnucnet__Zone__updateAbundances(
      work_zone.getNucnetZone(),
      p_abundances
    );

    p_z_qse =
      Libnucnet__Zone__getSummedAbundances(
        work_zone.getNucnetZone(),
        "z"
      );

    gsl_vector_free( p_abundances );

    Libnuceq__removeCluster( p_equil, p_cluster );

    //==========================================================================
    // Compute WSE.  Do it last for a zone because it modifies Ye.
    //==========================================================================

    if( argc == 4 )
    {

      if( zone.hasProperty( nnt::s_USE_SCREENING) )
      {

	zone.updateNseCorrectionFactorFunction(
	  (Libnucnet__Species__nseCorrectionFactorFunction)
	    user::my_coulomb_correction,
	  &my_coul_corr_data
	);

      }

      if(
	gsl_isinf(
	  zone.getElectronNeutrinoChemicalPotential()
	) != -1
      )
      {

	Libnuceq__clearYe( p_equil );

	Libnuceq__computeEquilibrium(
	  p_equil, 
	  boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
	  boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
	);

	p_abundances = Libnuceq__getAbundances( p_equil );

	d_ye_root = Libnuceq__computeZMoment( p_equil, 1 );

	Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_abundances );

	d_yh_wse = 
	  user::compute_cluster_abundance_moment(
	    zone,
	    S_XPATH,
	    "a",
	    0
	  );

	p_z_wse =
	  Libnucnet__Zone__getSummedAbundances( zone.getNucnetZone(), "z" );

	gsl_vector_free( p_abundances );

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

	user::update_my_rate_functions_data( zone );

	d_ye_root =
	  zone.computeRootFromQuantity(
	    (nnt::quantityFunction) user::my_yedot_root,
	    &zone
	  );

	d_yh_wse =
	  user::compute_cluster_abundance_moment(
	    zone,
	    S_XPATH,
	    "a",
	    0
	  );

	p_z_wse =
	  Libnucnet__Zone__getSummedAbundances( zone.getNucnetZone(), "z" );

      }

    }

    //==========================================================================
    // Print out.
    //==========================================================================

    if( zone.hasProperty( nnt::s_TIME ) )
      std::cout << "time = " << zone.getProperty( nnt::s_TIME ) << std::endl;

    std::cout << "t9 = " << zone.getProperty( nnt::s_T9 ) << std::endl;

    std::cout << "rho = " << zone.getProperty( nnt::s_RHO ) << std::endl;

    std::cout << "Ye = " << d_ye << std::endl;

    if( argc == 4 ) std::cout << "Ye_wse = " << d_ye_root << std::endl;

    std::cout << "Yh = " << d_yh << std::endl;

    std::cout << "Yh_nse = " << d_yh_nse << std::endl;

    if( argc == 4 ) std::cout << "Yh_wse = " << d_yh_wse << std::endl;

    for(
      i = 0;
      i <= Libnucnet__Nuc__getLargestNucleonNumber(
	     Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
	     "z"
	   );
      i++
    )
    {

      if( argc == 4 )
      {
	fprintf(
	  stdout,
	  "%lu %e %e %e %e\n",
	  (unsigned long) i,
	  gsl_vector_get( p_z_net, i ),
	  gsl_vector_get( p_z_nse, i ),
	  gsl_vector_get( p_z_qse, i ),
	  gsl_vector_get( p_z_wse, i )
	);
      }
      else
      {
	fprintf(
	  stdout,
	  "%lu %e %e %e\n",
	  (unsigned long) i,
	  gsl_vector_get( p_z_net, i ),
	  gsl_vector_get( p_z_nse, i ),
	  gsl_vector_get( p_z_qse, i )
	);
      }

    }

    gsl_vector_free( p_z_net );
    gsl_vector_free( p_z_nse );
    gsl_vector_free( p_z_qse );
    gsl_vector_free( p_z_wse );

    std::cout << std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Zone__free( work_zone.getNucnetZone() );
  Libnuceq__free( p_equil );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
