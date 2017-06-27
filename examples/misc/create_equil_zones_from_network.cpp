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
//! \brief Example code to create an output xml file with equilibrium
//!        abundances in place of the network abundances.
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"
#include "user/coul_corr.h"
#include "user/network_utilities.h"

#define B_USE_SCREENING  false

#define D_D_EPS       0.01

#define S_MY_CLUSTER  "my cluster"
#define S_FLAG        "equil flag"
#define S_DLT9         "equil t9 log derivative"
#define S_DLRHO        "equil rho log derivative"
#define S_DLYE         "equil ye log derivative"
#define S_NO_DERIVS    "no derivatives"

typedef struct{
  Libnuceq * pEquil;
  nnt::Zone * pZone;
} eq_data;

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
      " my_output.xml \"[z > 2]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 4 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " create new file with equilibrium abundances in place of network" <<
      " abundances for selected zones" <<
      std::endl;
    fprintf(
      stderr,
      "\nUsage: %s xml_file zone_xpath file_name deriv_xpath ...\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml file\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zones of interest\n\n"
    );
    fprintf(
      stderr,
      "  file_name = output xml file name\n\n"
    );
    fprintf(
      stderr,
      "  deriv_xpath = xpath for derivatives (or \"no derivatives\" for no calculation of derivatives)\n\n"
    );
    fprintf(
      stderr,
      "  constraints = XPath expression giving constraints (enter as many as you want)\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );

  }

}

//##############################################################################
// deriv().
//##############################################################################

double
deriv( double lx, void * p_data )
{

  eq_data * p_eq_data = (eq_data *) p_data;

  double x = exp( lx );

  if( p_eq_data->pZone->getProperty( S_FLAG ) == nnt::s_T9 )
  {
    Libnuceq__setYe(
      p_eq_data->pEquil,
      boost::lexical_cast<double>( p_eq_data->pZone->getProperty( nnt::s_YE ) )
    );
    Libnuceq__computeEquilibrium(
      p_eq_data->pEquil,
      x,
      boost::lexical_cast<double>( p_eq_data->pZone->getProperty( nnt::s_RHO ) )
    );
  }
  else if( p_eq_data->pZone->getProperty( S_FLAG ) == nnt::s_RHO )
  {
    Libnuceq__setYe(
      p_eq_data->pEquil,
      boost::lexical_cast<double>( p_eq_data->pZone->getProperty( nnt::s_YE ) )
    );
    Libnuceq__computeEquilibrium(
      p_eq_data->pEquil,
      boost::lexical_cast<double>( p_eq_data->pZone->getProperty( nnt::s_T9 ) ),
      x
    );
  }
  else if( p_eq_data->pZone->getProperty( S_FLAG ) == nnt::s_YE )
  {
    Libnuceq__setYe( p_eq_data->pEquil, x ); 
    Libnuceq__computeEquilibrium(
      p_eq_data->pEquil,
      boost::lexical_cast<double>( p_eq_data->pZone->getProperty( nnt::s_T9 ) ),
      boost::lexical_cast<double>( p_eq_data->pZone->getProperty( nnt::s_RHO ) )
    );
  }

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc(
        Libnucnet__NetView__getNet(
          p_eq_data->pZone->getNetView(
            p_eq_data->pZone->getProperty( S_MY_CLUSTER ).c_str(),
            ""
          )
        )
      )
    );

  double d_result = 0;

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    d_result +=
      Libnuceq__Species__getAbundance(
        Libnuceq__getSpeciesByName(
          p_eq_data->pEquil,
          Libnucnet__Species__getName( species.getNucnetSpecies() )
        )
      );

  }

  return log( d_result );
      
}

//##############################################################################
// compute_derivatives().
//##############################################################################

void
compute_derivatives(
  eq_data& my_data
)
{

  gsl_function F;
  double x;

  F.function = &deriv;
  F.params = &my_data;

  double d_dlt9, d_dlrho, d_dlye, d_abserr;

  my_data.pZone->updateProperty( S_FLAG, nnt::s_T9 );
  x = 
    log(
      boost::lexical_cast<double>( my_data.pZone->getProperty( nnt::s_T9 ) )
    );
  gsl_deriv_central( &F, x, x * D_D_EPS, &d_dlt9, &d_abserr );
  my_data.pZone->updateProperty(
    S_DLT9,
    boost::lexical_cast<std::string>( d_dlt9 )
  );
  
  my_data.pZone->updateProperty( S_FLAG, nnt::s_RHO );
  x =
    log(
      boost::lexical_cast<double>( my_data.pZone->getProperty( nnt::s_RHO ) )
    );
  gsl_deriv_central( &F, x, x * D_D_EPS, &d_dlrho, &d_abserr );
  my_data.pZone->updateProperty(
    S_DLRHO,
    boost::lexical_cast<std::string>( d_dlrho )
  );
  
  my_data.pZone->updateProperty( S_FLAG, nnt::s_YE );
  x =
    log(
      boost::lexical_cast<double>( my_data.pZone->getProperty( nnt::s_YE ) )
    );
  gsl_deriv_central( &F, x, x * D_D_EPS, &d_dlye, &d_abserr );
  my_data.pZone->updateProperty(
    S_DLYE,
    boost::lexical_cast<std::string>( d_dlye )
  );
  
}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  user::coul_corr_data my_coul_corr_data;
  eq_data my_data;

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
  // Create the zone list.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  //============================================================================
  // Get the primary views.
  //============================================================================

  zone_list.begin()->getNetView(
    "",
    ""
  );

  if( argc > 4 && std::string( argv[4] ) != S_NO_DERIVS )
  {
    zone_list.begin()->getNetView(
      argv[4],
      ""
    );
  }

  //============================================================================
  // Get the cluster views.
  //============================================================================

  size_t i_views = 0;

  for( int i = 5; i < argc; i++ )
  {
     zone_list.begin()->getNetView(
       argv[i],
       ""
     );
     i_views++;
  }

  //============================================================================
  // Equilibrium.
  //============================================================================

  my_data.pEquil =
    Libnuceq__new(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
    );

  //============================================================================
  // Iterate zones.
  //============================================================================

  for(
    nnt::zone_list_t::iterator it = zone_list.begin();
    it != zone_list.end();
    it++
  )
  {

      std::cout <<
        Libnucnet__Zone__getLabel( it->getNucnetZone(), 1 ) << std::endl;

    //==========================================================================
    // Screening.
    //==========================================================================

    if( B_USE_SCREENING )
    {
      my_coul_corr_data = user::get_coulomb_corr_data();

      it->updateNseCorrectionFactorFunction(
        (Libnucnet__Species__nseCorrectionFactorFunction)
          user::my_coulomb_correction,
        &my_coul_corr_data
      );
    }

    //==========================================================================
    // Clusters.
    //==========================================================================

    if( argc > 4 && std::string( argv[4] ) != S_NO_DERIVS )
      it->updateProperty( S_MY_CLUSTER, argv[4] );

    if( it != zone_list.begin() )
    {
      Libnucnet__Zone__copy_net_views(
        it->getNucnetZone(),
        zone_list.begin()->getNucnetZone()
      );
    }

    it->updateProperty(
      nnt::s_NUMBER_CLUSTERS,
      boost::lexical_cast<std::string>( i_views )
    );

    for( size_t i = 0; i < i_views; i++ )
    {

      it->updateProperty(
        nnt::s_CLUSTER_XPATH,
        boost::lexical_cast<std::string>( i ),
        argv[i+5]
      );

      it->updateProperty(
        nnt::s_CLUSTER_CONSTRAINT,
        boost::lexical_cast<std::string>( i ),
        boost::lexical_cast<std::string>(
          user::compute_cluster_abundance_moment(
            *it,
            argv[i+5],
            "a",
            0
          )
        )
      );
  
    }

    //==========================================================================
    // Set Ye.
    //==========================================================================
    
    it->updateProperty(
      nnt::s_YE,
      boost::lexical_cast<std::string>(
        user::compute_cluster_abundance_moment(
          *it,
          "",
          "z",
          1
        )
      )
    );

    //==========================================================================
    // Compute equilibrium.
    //==========================================================================

    nnt::update_equilibrium_with_zone_data( my_data.pEquil, *it );

    Libnuceq__computeEquilibrium(
      my_data.pEquil,
      boost::lexical_cast<double>( it->getProperty( nnt::s_T9 ) ),
      boost::lexical_cast<double>( it->getProperty( nnt::s_RHO ) )
    );

    gsl_vector * p_abunds = Libnuceq__getAbundances( my_data.pEquil );

    Libnucnet__Zone__updateAbundances(
      it->getNucnetZone(),
      p_abunds
    );

    gsl_vector_free( p_abunds );

    if( argc > 4 && std::string( argv[4] ) != S_NO_DERIVS )
    {
      it->updateProperty(
        nnt::s_YC,
        boost::lexical_cast<std::string>(
          user::compute_cluster_abundance_moment(
            *it,
            it->getProperty( S_MY_CLUSTER ).c_str(),
            "z",
            0
          )
        )
      );
      my_data.pZone = &(*it);
      compute_derivatives( my_data );
    }

  }

  //============================================================================
  // Output.
  //============================================================================

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_nucnet,
    "%.15e"
  );

  Libnucnet__writeToXmlFile( p_my_nucnet, argv[3] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnuceq__free( my_data.pEquil );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
