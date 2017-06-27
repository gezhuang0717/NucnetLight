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
//! \brief Example code to test versions of various libraries.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <stdlib.h>
#include <iostream>
#include <libxml/xmlversion.h>
#include <gsl/gsl_version.h>
#include <boost/version.hpp>

int main( )
{

  std::cout <<
    std::endl <<
    "Codes are compiled against:" << std::endl << std::endl;

  std::cout << 
     "libxml version " <<
     LIBXML_DOTTED_VERSION  << std::endl;

  std::cout <<
     "gsl version " <<
     GSL_VERSION << std::endl;

  std::cout <<
     "Boost version " <<
     BOOST_VERSION / 100000 << "." <<
     BOOST_VERSION / 100 % 1000 << "." <<
     BOOST_VERSION % 100  << std::endl;

  std::cout << std::endl;

  LIBXML_TEST_VERSION

  return EXIT_SUCCESS;

}
