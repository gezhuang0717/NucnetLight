////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu and Bradley S. Meyer.
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
//! \brief Example code to extend partition function data with the
//!        high-temperature data from Rauscher.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <fstream>
#include <iostream>

#include <nnt/auxiliary.h>

#include <Libnucnet.h>

#define I_RAUSCHER  48

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] <<
      " ../../data_pub/my_net.xml ../../data_pub/partf.txt updated_net.xml" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 4 )
  {
    fprintf(
      stderr,
      "\nUsage: %s net_xml partf_txt output_xml\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_xml = network xml filename\n\n"
    );
    fprintf(
      stderr, "  partf_txt = Rauscher high-temperature partf file\n\n"
    );
    fprintf(
      stderr, "  output_xml = output xml file\n\n"
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

  Libnucnet__Net * p_net;
  Libnucnet__Species * p_species;
  size_t i, i_size_old, i_size_new;
  std::string line, name, source;
  std::ifstream infile;
  unsigned int i_z, i_a;
  double d_j0, d_tmp;
  gsl_vector * p_t9_old, * p_log10_partf_old;
  gsl_vector * p_t9_new, * p_log10_partf_new;
  gsl_vector * p_t9_update, * p_log10_partf_update;
  gsl_vector_view view;

  double d_t9[] =
    { 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 35., 40.,
      45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95.,100.,
     105.,110.,115.,120.,125.,130.,135.,140.,145.,150.,155.,160.,
     165.,170.,175.,180.,190.,200.,210.,220.,230.,240.,250.,275.
    };

  //============================================================================
  // Check input.
  //============================================================================

  check_input( argc, argv );

  //============================================================================
  // Read and store input.
  //============================================================================

  p_net = Libnucnet__Net__new_from_xml( argv[1], NULL, NULL );

  p_t9_new = gsl_vector_alloc( I_RAUSCHER );
  p_t9_new->data = d_t9;
  p_log10_partf_new = gsl_vector_alloc( I_RAUSCHER );

  //============================================================================
  // Open input partf file.
  //============================================================================

  infile.open( argv[2] );

  //============================================================================
  // Skip information lines.
  //============================================================================

  for( i = 0; i < 63; i++ ) std::getline( infile, line, '\n' );

  //============================================================================
  // Loop over lines to end.
  //============================================================================

  while( infile >> name >> i_z >> i_a >> d_j0 )
  {

    for( i = 0; i < 48; i++ )
    {
      infile >> d_tmp;
      gsl_vector_set( p_log10_partf_new, i, log10( d_tmp ) );
    }

    p_species =
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc( p_net ),
        name.c_str()
      );

    if( p_species )
    {

      p_t9_old = Libnucnet__Species__getPartitionFunctionT9( p_species );

      p_log10_partf_old =
        Libnucnet__Species__getPartitionFunctionLog10( p_species );

      if( p_t9_old && p_log10_partf_old )
      {

	i_size_old = WnMatrix__get_gsl_vector_size( p_t9_old );

	i_size_new = i_size_old + I_RAUSCHER;
	
	p_t9_update = gsl_vector_alloc( i_size_new );

	p_log10_partf_update = gsl_vector_alloc( i_size_new );

	view =
	  gsl_vector_subvector( p_t9_update, 0, i_size_old );

	gsl_vector_memcpy( &view.vector, p_t9_old );

	view =
	  gsl_vector_subvector( p_t9_update, i_size_old, I_RAUSCHER );

	gsl_vector_memcpy( &view.vector, p_t9_new );

	view =
	  gsl_vector_subvector( p_log10_partf_update, 0, i_size_old );

	gsl_vector_memcpy( &view.vector, p_log10_partf_old );

	view =
	  gsl_vector_subvector( p_log10_partf_update, i_size_old, I_RAUSCHER );

	gsl_vector_memcpy( &view.vector, p_log10_partf_new );

	Libnucnet__Species__updatePartitionFunctionData(
	  p_species,
	  p_t9_update,
	  p_log10_partf_update
	);

        gsl_vector_free( p_t9_update );
        gsl_vector_free( p_log10_partf_update );

      }
      else
      {

	Libnucnet__Species__updatePartitionFunctionData(
	  p_species,
	  p_t9_new,
	  p_log10_partf_new
	);

      }

      source = std::string( Libnucnet__Species__getSource( p_species ) );

      source += " + Rauscher high-temperature partition function data";

      Libnucnet__Species__updateSource(
        p_species,
        source.c_str()
      );

    }

  }

  //============================================================================
  // Write to output.
  //============================================================================

  Libnucnet__Net__writeToXmlFile(
    p_net,
    argv[3]
  );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  infile.close();
  gsl_vector_free( p_t9_new );
  gsl_vector_free( p_log10_partf_new );
  Libnucnet__Net__free( p_net );
  return EXIT_SUCCESS;

}

