//////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code to convert text ffn rates file to xml file, storing
//!     both two-d weak rates and log10 ft values. 
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>
 
#include <Libnucnet.h>
#include <Libstatmech.h>

#include "nnt/two_d_weak_rates.h"
#include "nnt/string_defs.h"
#include "nnt/param_defs.h"

#define D_ALPHA    1.
#define D_MU_NUE_KT  GSL_NEGINF

//##############################################################################
// Prototypes.
//##############################################################################

double
convert_to_log10_ft(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  Libstatmech__Fermion *,
  double,
  double,
  double,
  const char *
);

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet__Net *p_my_net;
  Libnucnet__Reaction *p_reaction_bp, *p_reaction_ec, *p_reaction_bm;
  Libnucnet__Reaction *p_reaction_pc;
  Libnucnet__Species *p_reactant, *p_target;
  Libstatmech__Fermion *p_electron;

  std::ifstream input_file;
  std::string s_source;
  int i_row, i_col, i_number_rows, i_number_columns;
  unsigned int i_z, i_a;
  std::string s_t9, s_log10_rhoe, s_bp, s_ec, s_bm, s_pc;
  std::string s_row, s_col, s_new_value;
  std::string s_nu, s_nubar, s_av_e_nu_e, s_av_e_nubar_e;
  std::string s_line;

  //============================================================================
  // Set some pairs for extracting substrings.
  //============================================================================

  std::pair<int,int> pair_z = std::make_pair( 23, 2 );
  std::pair<int,int> pair_a = std::make_pair( 36, 3 );

  std::pair<int,int> pair_t9 = std::make_pair( 0, 6 );
  std::pair<int,int> pair_lrhoe = std::make_pair( 6, 4 );
  std::pair<int,int> pair_bp = std::make_pair( 17, 7 );
  std::pair<int,int> pair_ec = std::make_pair( 24, 8 );
  std::pair<int,int> pair_bm = std::make_pair( 48, 8 );
  std::pair<int,int> pair_pc = std::make_pair( 56, 8 );
  
  std::pair<int,int> pair_nu = std::make_pair( 40, 8 );
  std::pair<int,int> pair_nubar = std::make_pair( 72, 8 );

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 4 ) {
      fprintf(
        stderr, "\nUsage: %s nuc_file weak_file output_file \n", argv[0]
      );
      fprintf(
        stderr, "\n  nuc_file: name of input nuclear datal xml file\n"
      );
      fprintf(
        stderr, "\n  weak_file: name of input text file\n"
      );
      fprintf(
        stderr, "\n  output_file: name of output xml file\n\n"
      );
      return EXIT_FAILURE;
  }

  //============================================================================
  // Create net.
  //============================================================================

  p_my_net = Libnucnet__Net__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( p_my_net ),
    argv[1],
    NULL
  );

  //============================================================================
  // Open input file.
  //============================================================================

  input_file.open( argv[2] );
  if( !input_file )
  {
    fprintf( stderr, "Could not open file.\n" );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Read source and number of rows and columns.
  //============================================================================

  std::getline( input_file, s_source );
  boost::trim_if( s_source, boost::is_any_of( "\r" ) );
  boost::trim( s_source );

  std::getline( input_file, s_line );

  i_number_rows = atoi( s_line.substr( 0, 2 ).c_str() );

  i_number_columns = atoi( s_line.substr( 2, 3 ).c_str() );

  //============================================================================
  // Get electron.
  //============================================================================

  p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  //============================================================================
  // Read in input.
  //============================================================================

  while( !input_file.eof() )
  {

  //----------------------------------------------------------------------------
  // Read in and get target.
  //----------------------------------------------------------------------------

    std::getline( input_file, s_line );

    if( !input_file.good() ) break;

    i_z = atoi( s_line.substr( pair_z.first, pair_z.second ).c_str() );
    i_a = atoi( s_line.substr( pair_a.first, pair_a.second ).c_str() );

    p_reactant =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ),
        i_z,
        i_a,
        NULL
      );

    std::getline( input_file, s_line );
    std::getline( input_file, s_line );
    std::getline( input_file, s_line );

  //----------------------------------------------------------------------------
  // Create reactions and add reactants and products.
  //----------------------------------------------------------------------------

    //..........................................................................
    // beta plus.
    //..........................................................................

    p_reaction_bp = Libnucnet__Reaction__new();

    Libnucnet__Reaction__updateSource( p_reaction_bp, s_source.c_str() );

    Libnucnet__Reaction__addReactant(
      p_reaction_bp,
      Libnucnet__Species__getName( p_reactant )
    );

    p_target =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), i_z - 1, i_a, NULL
      );

    Libnucnet__Reaction__addProduct(
      p_reaction_bp,
      Libnucnet__Species__getName( p_target )
    );

    Libnucnet__Reaction__addProduct( p_reaction_bp, "positron" );
    Libnucnet__Reaction__addProduct( p_reaction_bp, "neutrino_e" );

    Libnucnet__Reaction__setUserRateFunctionKey(
      p_reaction_bp,
      nnt::s_TWO_D_WEAK_RATES
    );

    //..........................................................................
    // electron capture.
    //..........................................................................

    p_reaction_ec = Libnucnet__Reaction__new();

    Libnucnet__Reaction__updateSource( p_reaction_ec, s_source.c_str() );

    Libnucnet__Reaction__addReactant(
      p_reaction_ec,
      Libnucnet__Species__getName( p_reactant )
    );

    Libnucnet__Reaction__addReactant( p_reaction_ec, "electron" );

    p_target =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), i_z - 1, i_a, NULL
      );

    Libnucnet__Reaction__addProduct(
      p_reaction_ec,
      Libnucnet__Species__getName( p_target )
    );

    Libnucnet__Reaction__addProduct( p_reaction_ec, "neutrino_e" );

    Libnucnet__Reaction__setUserRateFunctionKey(
      p_reaction_ec,
      nnt::s_TWO_D_WEAK_RATES_LOG10_FT
    );

  //----------------------------------------------------------------------------
  // Reverse reactions. Decrease Z.
  //----------------------------------------------------------------------------

    p_reactant =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), --i_z, i_a, NULL
      );

    //..........................................................................
    // beta minus.
    //..........................................................................

    p_reaction_bm = Libnucnet__Reaction__new();

    Libnucnet__Reaction__updateSource( p_reaction_bm, s_source.c_str() );

    Libnucnet__Reaction__addReactant(
      p_reaction_bm,
      Libnucnet__Species__getName( p_reactant )
    );

    p_target =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), i_z + 1, i_a, NULL
      );

    Libnucnet__Reaction__addProduct(
      p_reaction_bm,
      Libnucnet__Species__getName( p_target )
    );

    Libnucnet__Reaction__addProduct( p_reaction_bm, "electron" );
    Libnucnet__Reaction__addProduct( p_reaction_bm, "anti-neutrino_e" );

    Libnucnet__Reaction__setUserRateFunctionKey(
      p_reaction_bm,
      nnt::s_TWO_D_WEAK_RATES
    );

    //..........................................................................
    // positron capture.
    //..........................................................................

    p_reaction_pc = Libnucnet__Reaction__new();

    Libnucnet__Reaction__updateSource( p_reaction_pc, s_source.c_str() );

    Libnucnet__Reaction__addReactant(
      p_reaction_pc,
      Libnucnet__Species__getName( p_reactant )
    );

    Libnucnet__Reaction__addReactant( p_reaction_pc, "positron" );

    p_target =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), i_z + 1, i_a, NULL
      );

    Libnucnet__Reaction__addProduct(
      p_reaction_pc,
      Libnucnet__Species__getName( p_target )
    );

    Libnucnet__Reaction__addProduct( p_reaction_pc, "anti-neutrino_e" );

    Libnucnet__Reaction__setUserRateFunctionKey(
      p_reaction_pc,
      nnt::s_TWO_D_WEAK_RATES_LOG10_FT
    );

  //----------------------------------------------------------------------------
  // Read in data.
  //----------------------------------------------------------------------------

    for( i_col = 0; i_col < i_number_columns; i_col++ )
    {

      s_col = boost::lexical_cast<std::string>( i_col );

      for( i_row = 0; i_row < i_number_rows; i_row++ )
      {

        s_row = boost::lexical_cast<std::string>( i_row );

        std::getline( input_file, s_line );

        s_t9 = s_line.substr( pair_t9.first, pair_t9.second );
        s_log10_rhoe = s_line.substr( pair_lrhoe.first, pair_lrhoe.second );

        s_bp = s_line.substr( pair_bp.first, pair_bp.second );
        s_ec = s_line.substr( pair_ec.first, pair_ec.second );
        s_bm = s_line.substr( pair_bm.first, pair_bm.second );
        s_pc = s_line.substr( pair_pc.first, pair_pc.second );

        s_nu = s_line.substr( pair_nu.first, pair_nu.second );
        s_nubar = s_line.substr( pair_nubar.first, pair_nubar.second );

        boost::trim( s_t9 );
        boost::trim( s_log10_rhoe );
        boost::trim( s_bp );
        boost::trim( s_ec );
        boost::trim( s_bm );
        boost::trim( s_pc );

        boost::trim( s_nu );
        boost::trim( s_nubar );

    //..........................................................................
    // Update t9 tables for each reaction.
    //..........................................................................

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bp,
          nnt::s_T9,
          s_row.c_str(),
          NULL,
          s_t9.c_str()
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_ec,
          nnt::s_T9,
          s_row.c_str(),
          NULL,
          s_t9.c_str()
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bm,
          nnt::s_T9,
          s_row.c_str(),
          NULL,
          s_t9.c_str()
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_pc,
          nnt::s_T9,
          s_row.c_str(),
          NULL,
          s_t9.c_str()
        ); 

    //..........................................................................
    // Update rhoe tables for each reaction.
    //..........................................................................

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bp,
          nnt::s_LOG10_RHOE,
          s_col.c_str(),
          NULL,
          s_log10_rhoe.c_str()
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_ec,
          nnt::s_LOG10_RHOE,
          s_col.c_str(),
          NULL,
          s_log10_rhoe.c_str()
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bm,
          nnt::s_LOG10_RHOE,
          s_col.c_str(),
          NULL,
          s_log10_rhoe.c_str()
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_pc,
          nnt::s_LOG10_RHOE,
          s_col.c_str(),
          NULL,
          s_log10_rhoe.c_str()
        ); 

    //..........................................................................
    // Update rate tables for each reaction.
    //..........................................................................

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bp,
          nnt::s_LOG10_RATE,
          s_row.c_str(),
          s_col.c_str(),
          s_bp.c_str()
        );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_ec,
          nnt::s_LOG10_RATE,
          s_row.c_str(),
          s_col.c_str(),
          s_ec.c_str()
        );

        s_new_value =
          boost::lexical_cast<std::string>(
            convert_to_log10_ft(
              p_reaction_ec,
              p_my_net,
              p_electron,
              atof( s_t9.c_str() ),
              pow( 10., atof( s_log10_rhoe.c_str() ) ), 
              D_MU_NUE_KT,
              s_ec.c_str()
            )
          );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_ec,
          nnt::s_LOG10_FT,
          s_row.c_str(),
          s_col.c_str(),
          s_new_value.c_str()
        );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bm,
          nnt::s_LOG10_RATE,
          s_row.c_str(),
          s_col.c_str(),
          s_bm.c_str()
        );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_pc,
          nnt::s_LOG10_RATE,
          s_row.c_str(),
          s_col.c_str(),
          s_pc.c_str()
        );

        s_new_value =
          boost::lexical_cast<std::string>(
            convert_to_log10_ft(
              p_reaction_pc,
              p_my_net,
              p_electron,
              atof( s_t9.c_str() ),
              pow( 10., atof( s_log10_rhoe.c_str() ) ), 
              D_MU_NUE_KT,
              s_pc.c_str()
            )
          );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_pc,
          nnt::s_LOG10_FT,
          s_row.c_str(),
          s_col.c_str(),
          s_new_value.c_str()
        );

    //..........................................................................
    // Update average nu and nubar energies.
    //..........................................................................

        s_av_e_nu_e =
          boost::lexical_cast<std::string>(
            pow(
              10.,
              atof( s_nu.c_str() )
            ) /
            (
              pow(
                10.,
                atof( s_bp.c_str() )
              )
              +
              pow(
                10.,
                atof( s_ec.c_str() )
              )
            )
          );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bp,
          nnt::s_AVERAGE_ENERGY_NU_E,
          s_row.c_str(),
          s_col.c_str(),
          s_av_e_nu_e.c_str()
        );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_ec,
          nnt::s_AVERAGE_ENERGY_NU_E,
          s_row.c_str(),
          s_col.c_str(),
          s_av_e_nu_e.c_str()
        );

        s_av_e_nubar_e =
          boost::lexical_cast<std::string>(
            pow(
              10.,
              atof( s_nubar.c_str() )
            ) /
            (
              pow(
                10.,
                atof( s_bm.c_str() )
              )
              +
              pow(
                10.,
                atof( s_pc.c_str() )
              )
            )
          );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bm,
          nnt::s_AVERAGE_ENERGY_NUBAR_E,
          s_row.c_str(),
          s_col.c_str(),
          s_av_e_nubar_e.c_str()
        );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_pc,
          nnt::s_AVERAGE_ENERGY_NUBAR_E,
          s_row.c_str(),
          s_col.c_str(),
          s_av_e_nubar_e.c_str()
        );

      }

    }

  //----------------------------------------------------------------------------
  // Add reactions.
  //----------------------------------------------------------------------------

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_my_net ),
      p_reaction_bp
    );

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_my_net ),
      p_reaction_ec
    );

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_my_net ),
      p_reaction_bm
    );

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_my_net ),
      p_reaction_pc
    );

  //----------------------------------------------------------------------------
  // Read end.
  //----------------------------------------------------------------------------

    std::getline( input_file, s_line );

  }     
  
  //============================================================================
  // Close input file.
  //============================================================================

  input_file.close();

  //============================================================================
  // Write to XML.
  //============================================================================

  Libnucnet__Reac__writeToXmlFile(
    Libnucnet__Net__getReac( p_my_net ),
    argv[3]
  );

  //============================================================================
  // Clean up.
  //============================================================================

  Libstatmech__Fermion__free( p_electron );
  Libnucnet__Net__free( p_my_net );

  //============================================================================
  // Done!
  //============================================================================

  return EXIT_SUCCESS;

}

//##############################################################################
// convert_to_log10_ft().
//##############################################################################

double
convert_to_log10_ft(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Net *p_net,
  Libstatmech__Fermion *p_electron,
  double d_t9,
  double d_rhoe,
  double d_mu_nue_kT,
  const char *s_value
)
{

  double d_I, d_eta_F;

  d_eta_F = 
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      d_t9 * GSL_CONST_NUM_GIGA,
      d_rhoe * GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    ) +
    Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV( d_t9 );

  d_I = 
    nnt::ffnIV__compute_Ie(
      p_reaction,
      p_net, 
      Libstatmech__Fermion__getRestMass( p_electron ),
      d_t9,
      d_eta_F,
      d_mu_nue_kT
    );

  if( d_I > 0. )
    return
      log10(
        M_LN2
      ) +
      log10( d_I ) -
      atof( s_value ) +
      log10( D_ALPHA );
  else
    return 0;

}
