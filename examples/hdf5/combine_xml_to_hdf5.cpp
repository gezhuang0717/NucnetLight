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
//! \brief Example code to convert xml to hdf5.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>
#include <fstream>
#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/hdf5_routines.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " my_output_files.txt my_output.h5" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 3 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " combines xml to hdf5." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s file_list hdf5_file\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  file_list = file containing list of network xml output files\n\n"
    );
    fprintf(
      stderr,
      "  hdf5_file = output hdf5 file\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  std::ifstream file_list;
  std::string file_name;
  Libnucnet *p_my_nucnet;
  
  check_input( argc, argv );

  file_list.open( argv[1] );

  if( !( std::getline( file_list, file_name ) ) )
  {
    std::cerr << "No files in file list." << std::endl;
    exit( EXIT_FAILURE );
  }

  p_my_nucnet =
    Libnucnet__new_from_xml(
      file_name.c_str(),
      NULL,
      NULL,
      NULL
  );

  user::hdf5::create_output( argv[2], p_my_nucnet );

  user::hdf5::append_zones( argv[2], p_my_nucnet );

  Libnucnet__free( p_my_nucnet );

  while( std::getline( file_list, file_name ) )
  { 
    if( file_name == "" ) break;
    std::cout << file_name << std::endl;
    p_my_nucnet =
      Libnucnet__new_from_xml(
        file_name.c_str(),
        NULL,
        NULL,
        NULL
    );
    user::hdf5::append_zones( argv[2], p_my_nucnet );
    Libnucnet__free( p_my_nucnet );
  }

  file_list.close();

  return EXIT_SUCCESS;

}
