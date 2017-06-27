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
//! \brief Example code for simple hydrodynamics.
////////////////////////////////////////////////////////////////////////////////

#include "hydro.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

//##############################################################################
// Code for exponential expansion: rho(t) = rho_0 * exp(-t/tau).
//##############################################################################

#ifdef HYDRO_EXP_EXPANSION

//##############################################################################
// get_nucnet().
//##############################################################################

Libnucnet *
get_nucnet( int argc, char **argv )
{

  Libnucnet * p_nucnet;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--help" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] <<
      " runs a single-zone network calculation for matter expanding" <<
      " exponentially." << std::endl << std::endl;
    std::cout << "For a usage statement, type " << std::endl << std::endl;
    std::cout << argv[0] << " --usage" << std::endl << std::endl;
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
#ifdef HDF5
    std::cout << argv[0] << " ../../data_pub/my_net.xml " <<
      "../../data_pub/zone.xml my_output.h5 \"[z <= 20]\"" <<
      std::endl << std::endl;
#else
    std::cout << argv[0] << " ../../data_pub/my_net.xml " <<
      "../../data_pub/zone.xml my_output.xml \"[z <= 20]\"" <<
      std::endl << std::endl;
#endif
    exit( EXIT_FAILURE );
  }

  if( argc < 4 || argc > 6 || strcmp( argv[1], "--usage" ) == 0 )
  {
    fprintf(
      stderr,
      "\nUsage: %s net_file zone_file out_file xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_file = input network data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_file = input single zone data xml filename\n\n"
    );
#ifdef HDF5
    fprintf(
      stderr, "  out_file = output hdf5 filename\n\n"
    );
#else
    fprintf(
      stderr, "  out_file = output xml filename\n\n"
    );
#endif
    fprintf(
      stderr,
      "  xpath_nuc = nuclear xpath expression (optional--required if xpath_reac specified)\n\n"
    );
    fprintf(
      stderr, "  xpath_reac = reaction xpath expression (optional)\n\n"
    );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Validate input net file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__Net__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet net input!\n" );
      exit( EXIT_FAILURE );
    }
  }

  //============================================================================
  // Validate input zone file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_zone_data_xml( argv[2] ) ) {
      fprintf( stderr, "Not valid libnucnet zone data input!\n" );
      exit( EXIT_FAILURE );
    }
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  p_nucnet = Libnucnet__new();

  if( argc == 4 )
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      NULL,
      NULL
    );
  }
  else if( argc == 5 )
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      argv[4],
      NULL
    );
  }
  else
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      argv[4],
      argv[5]
    );
  }

  Libnucnet__assignZoneDataFromXml( p_nucnet, argv[2], NULL );

  return p_nucnet;

}

//##############################################################################
// initialize_zone().
//##############################################################################

void
initialize_zone( nnt::Zone& zone, char ** argv )
{

  //============================================================================
  // Set output file.
  //============================================================================

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    S_OUTPUT_FILE,
    NULL,
    NULL,
    argv[3]
  );

  //============================================================================
  // Set the initial evolution network.
  //============================================================================

  zone.updateProperty(
    nnt::s_RHO,
    zone.getProperty( nnt::s_RHO_0 )
  );

  zone.updateProperty(
    nnt::s_T9,
    zone.getProperty( nnt::s_T9_0 )
  );

  zone.updateProperty(
    nnt::s_TIME,
    boost::lexical_cast<std::string>( 0. )
  );

}

//##############################################################################
// update_zone_properties().
//##############################################################################

void
update_zone_properties( nnt::Zone& zone )
{

  double d_t, d_tau = 0;

  d_t = boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );

  try
  {
    d_tau = boost::lexical_cast<double>( zone.getProperty( nnt::s_TAU ) );
  }
  catch( const boost::bad_lexical_cast& e )
  {
    if( zone.getProperty( nnt::s_TAU ) == "inf" )
    {
      return;
    }
    else
    {
      std::cerr << "Invalid tau." << std::endl;
    }
  }

  zone.updateProperty(
    nnt::s_RHO,
    boost::lexical_cast<std::string>(
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_RHO_0 ) 
      ) * exp( -d_t / d_tau )
    )
  );

  zone.updateProperty(
    nnt::s_T9,
    boost::lexical_cast<std::string>(
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_T9_0 ) 
      ) * exp( -d_t / ( 3. * d_tau ) )
    )
  );

  if(
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) < 1.e-10
  )
    zone.updateProperty(
      nnt::s_T9,
      boost::lexical_cast<std::string>( 1.e-10 )
    );

  if(
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) < 1.e-30
  )
    zone.updateProperty(
      nnt::s_RHO,
      boost::lexical_cast<std::string>( 1.e-30 )
    );

}

//##############################################################################
// set_zone().
//##############################################################################

int
set_zone( Libnucnet * p_nucnet, nnt::Zone& zone, char ** argv )
{

  if( !argv )
  {
    std::cerr << "Invalid input." << std::endl;
    return 0;
  }

  if( !Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" ) )
  {
    std::cerr << "Zone not found!" << std::endl;
    return 0;
  }

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" )
  );

  return 1;

}

#endif // HYDRO_EXP_EXPANSION

//##############################################################################
// Code for expansion from ascii trajectory file.
//##############################################################################

#ifdef HYDRO_TRAJ

//##############################################################################
// globals.
//##############################################################################

std::vector<double> time, t9, log10_rho;

//##############################################################################
// get_nucnet().
//##############################################################################

Libnucnet *
get_nucnet( int argc, char **argv )
{

  Libnucnet * p_nucnet;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--help" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] <<
      " runs a single-zone network calculation for matter expanding" <<
      std::endl;
    std::cout <<
      "according to data from a trajectory file." << std::endl << std::endl;
    std::cout << "For a usage statement, type " << std::endl << std::endl;
    std::cout << argv[0] << " --usage" << std::endl << std::endl;
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
#ifdef HDF5
    std::cout << argv[0] << " ../../data_pub/my_net.xml " <<
      "../../data_pub/zone.xml traj.txt my_output.h5 \"[z <= 20]\"" <<
      std::endl << std::endl;
#else
    std::cout << argv[0] << " ../../data_pub/my_net.xml " <<
      "../../data_pub/zone.xml traj.txt my_output.xml \"[z <= 20]\"" <<
      std::endl << std::endl;
#endif
    exit( EXIT_FAILURE );
  }

  if ( argc < 5 || argc > 7 || strcmp( argv[1], "--usage" ) == 0 )
  {

    fprintf(
      stderr,
      "\nUsage: %s net_file zone_file traj_file out_file xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_file = input network data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_file = input zone data xml filename\n\n"
    );
    fprintf(
      stderr, "  traj_file = trajectory text file\n\n"
    );
#ifdef HDF5
    fprintf(
      stderr, "  out_file = output hdf5 file\n\n"
    );
#else
    fprintf(
      stderr, "  out_file = output xml file\n\n"
    );
#endif
    fprintf(
      stderr,
      "  xpath_nuc = nuclear xpath expression (optional--required if xpath_reac specified)\n\n"
    );
    fprintf(
      stderr, "  xpath_reac = reaction xpath expression (optional)\n\n"
    );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Validate input net file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__Net__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet input!\n" );
      exit( EXIT_FAILURE );
    }
  }

  //============================================================================
  // Validate input zone file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_zone_data_xml( argv[2] ) ) {
      fprintf( stderr, "Not valid libnucnet zone data input!\n" );
      exit( EXIT_FAILURE );
    }
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  if( argc == 5 ) {
    p_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );
  } else if( argc == 6 ) {
    p_nucnet = Libnucnet__new_from_xml( argv[1], argv[5], NULL, NULL );
  } else {
    p_nucnet = Libnucnet__new_from_xml( argv[1], argv[5], argv[6], NULL );
  }

  //============================================================================
  // Get zone data.
  //============================================================================

  Libnucnet__assignZoneDataFromXml( p_nucnet, argv[2], NULL );

  //============================================================================
  // Get trajectory data.
  //============================================================================

  get_trajectory_data( argv[3] );

  //============================================================================
  // Done.
  //============================================================================

  return p_nucnet;

}

//##############################################################################
// get_trajectory_data().
//##############################################################################

void
get_trajectory_data( char *s_file )
{

  std::ifstream my_file;
  double d_x1, d_x2, d_x3;

  //============================================================================
  // Open file thermodynamics file.
  //============================================================================

  my_file.open( s_file );

  if( my_file.bad() )
  {
    fprintf( stderr, "Couldn't open file.\n" );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Assign vectors.
  //============================================================================

  while( my_file >> d_x1 >> d_x2 >> d_x3 )
  {
    time.push_back( d_x1 );
    t9.push_back( d_x2 );
    log10_rho.push_back( log10( d_x3 ) );
  }

  my_file.close();

}

//##############################################################################
// initialize_zone().
//##############################################################################

void
initialize_zone( nnt::Zone& zone, char ** argv )
{

  //============================================================================
  // Set output file.
  //============================================================================

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    S_OUTPUT_FILE,
    NULL,
    NULL,
    argv[4]
  );

  update_t9_rho_in_zone_by_interpolation(
    zone,
    "spline",
    time,
    t9,
    log10_rho
  );

  return;

}

//##############################################################################
// update_zone_properties().
//##############################################################################

void
update_zone_properties(
  nnt::Zone& zone
)
{

  update_t9_rho_in_zone_by_interpolation(
    zone,
    "spline",
    time,
    t9,
    log10_rho
  );

} 

//##############################################################################
// set_zone().
//##############################################################################

int
set_zone( Libnucnet * p_nucnet, nnt::Zone& zone, char ** argv )
{

  if( !argv )
  {
    std::cerr << "Invalid input." << std::endl;
    return 0;
  }

  if( !Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" ) )
  {
    std::cerr << "Zone not found!" << std::endl;
    return 0;
  }

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" )
  );

  return 1;

}

#endif // HYDRO_TRAJ

//##############################################################################
// Code for expansion by spherical wind.
//##############################################################################

#ifdef HYDRO_WIND

//##############################################################################
// get_nucnet().
//##############################################################################

Libnucnet *
get_nucnet( int argc, char **argv )
{

  Libnucnet * p_nucnet;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--help" ) == 0 )
  {
    std::cout << std::endl;
  }

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
#ifdef HDF5
    std::cout << argv[0] << " ../../data_pub/my_net.xml " <<
      "../../data_pub/zone.xml my_output.h5 \"[z <= 20]\"" <<
      std::endl << std::endl;
#else
    std::cout << argv[0] << " ../../data_pub/my_net.xml " <<
      "../../data_pub/zone.xml my_output.xml \"[z <= 20]\"" <<
      std::endl << std::endl;
#endif
    exit( EXIT_FAILURE );
  }

  if( argc < 4 || argc > 6 )
  {
    std::cout << "Purpose:" << std::endl;
    std::cout << "  " << argv[0] <<
      " runs a single-zone network calculation for matter expanding" <<
      " in a constant velocity wind." << std::endl << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s net_file zone_file out_file xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_file = input network data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_file = input single zone data xml filename\n\n"
    );
#ifdef HDF5
    fprintf(
      stderr, "  out_file = output hdf5 filename\n\n"
    );
#else
    fprintf(
      stderr, "  out_file = output xml filename\n\n"
    );
#endif
    fprintf(
      stderr,
      "  xpath_nuc = nuclear xpath expression (optional--required if xpath_reac specified)\n\n"
    );
    fprintf(
      stderr, "  xpath_reac = reaction xpath expression (optional)\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Validate input net file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__Net__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet net input!\n" );
      exit( EXIT_FAILURE );
    }
  }

  //============================================================================
  // Validate input zone file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_zone_data_xml( argv[2] ) ) {
      fprintf( stderr, "Not valid libnucnet zone data input!\n" );
      exit( EXIT_FAILURE );
    }
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  p_nucnet = Libnucnet__new();

  if( argc == 4 )
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      NULL,
      NULL
    );
  }
  else if( argc == 5 )
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      argv[4],
      NULL
    );
  }
  else
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      argv[4],
      argv[5]
    );
  }

  Libnucnet__assignZoneDataFromXml( p_nucnet, argv[2], NULL );

  return p_nucnet;

}

//##############################################################################
// initialize_zone().
//##############################################################################

void
initialize_zone( nnt::Zone& zone, char ** argv )
{

  //============================================================================
  // Set output file.
  //============================================================================

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    S_OUTPUT_FILE,
    NULL,
    NULL,
    argv[3]
  );

  //============================================================================
  // Set the initial evolution network.
  //============================================================================

  zone.updateProperty(
    nnt::s_RHO,
    zone.getProperty( nnt::s_RHO_0 )
  );

  zone.updateProperty(
    nnt::s_T9,
    zone.getProperty( nnt::s_T9_0 )
  );

  zone.updateProperty(
    nnt::s_TIME,
    boost::lexical_cast<std::string>( 0. )
  );

  zone.updateProperty(
    nnt::s_RADIUS,
    nnt::s_RADIUS_0
  );

}

//##############################################################################
// update_zone_properties().
//##############################################################################

void
update_zone_properties(
  nnt::Zone& zone
)
{

  double d_t = boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );

  zone.updateProperty(
    nnt::s_RHO,
    boost::lexical_cast<std::string>(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO_0 ) ) /
      (
        gsl_pow_2(
          1. +
          (
            d_t /
            boost::lexical_cast<double>( zone.getProperty( nnt::s_TAU ) )
          )
        )
      )
    )
  );

  zone.updateProperty(
    nnt::s_RADIUS,
    boost::lexical_cast<std::string>(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RADIUS_0 ) ) *
      (
        1. +
        (
          d_t /
          boost::lexical_cast<double>( zone.getProperty( nnt::s_TAU ) )
        )
      )
    )
  );

  if( zone.hasProperty( nnt::s_T9_0 ) )
  {
    zone.updateProperty(
      nnt::s_T9,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9_0 ) ) *
        pow(
          boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) /
          boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO_0 ) ),
          1. / 3.
        )
      )
    );
  }

}

//##############################################################################
// set_zone().
//##############################################################################

int
set_zone( Libnucnet * p_nucnet, nnt::Zone& zone, char ** argv )
{

  if( !argv )
  {
    std::cerr << "Invalid input." << std::endl;
    return 0;
  }

  if( !Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" ) )
  {
    std::cerr << "Zone not found!" << std::endl;
    return 0;
  }

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" )
  );

  return 1;

}

#endif // HYDRO_WIND

//##############################################################################
// Code for expansion by shock.
//##############################################################################

#ifdef HYDRO_SHOCK

//##############################################################################
// get_nucnet().
//##############################################################################

Libnucnet *
get_nucnet( int argc, char **argv )
{

  Libnucnet * p_nucnet;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc < 7 || argc > 9 ) {
    fprintf(
      stderr,
      "\nUsage: %s input_xml out_file zone mach t_end steps nuc_xpath reac_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  input_xml = input presn data xml filename\n\n"
    );
#ifdef HDF5
    fprintf(
      stderr, "  out_file = output hdf5 filename\n\n"
    );
#else
    fprintf(
      stderr, "  out_file = output xml filename\n\n"
    );
#endif
    fprintf(
      stderr, "  zone = label of zone\n\n"
    );
    fprintf(
      stderr,
      "  mach = mach number of shock\n\n"
    );
    fprintf(
      stderr,
      "  t_end = duration of calculation (in seconds)\n\n"
    );
    fprintf(
      stderr,
      "  steps = frequency of timestep print out\n\n"
    );
    fprintf(
      stderr,
      "  nuc_xpath = XPath expression for nuclei \
(optional--required if reac_xpath present)\n\n"
    );
    fprintf(
      stderr,
      "  reac_xpath = XPath expression for reactions (optional)\n\n"
    );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Validate input net file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet input!\n" );
      exit( EXIT_FAILURE );
    }
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  if( argc == 7 )
    p_nucnet =
      Libnucnet__new_from_xml(
        argv[1],
        NULL,
        NULL,
        NULL
      );
  else if( argc == 8 )
    p_nucnet =
      Libnucnet__new_from_xml(
        argv[1],
        argv[7],
        NULL,
        NULL
      );
  else
    p_nucnet =
      Libnucnet__new_from_xml(
        argv[1],
        argv[7],
        argv[8],
        NULL
      );

  return p_nucnet;

}

//##############################################################################
// initialize_zone().
//##############################################################################

void
initialize_zone( nnt::Zone& zone, char ** argv )
{

  //============================================================================
  // Set output file.
  //============================================================================

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    S_OUTPUT_FILE,
    NULL,
    NULL,
    argv[2]
  );

  //============================================================================
  // Set the initial evolution network.
  //============================================================================

  if(
    !set_zone_post_shock_conditions( zone, argv[4] )
  )
  {
    std::cerr << "Couldn't set zone conditions." << std::endl;
    exit( EXIT_FAILURE );
  }

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    nnt::s_TEND,
    NULL,
    NULL,
    argv[5]
  );

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    nnt::s_STEPS,
    NULL,
    NULL,
    argv[6]
  );

  zone.normalizeAbundances();

}

//##############################################################################
// update_zone_properties().
//##############################################################################

void
update_zone_properties( nnt::Zone& zone )
{

  double d_t, d_tau;

  d_t = boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );

  d_tau = boost::lexical_cast<double>( zone.getProperty( nnt::s_TAU ) );

  if( d_tau > 0. )
  { 

    zone.updateProperty(
      nnt::s_RHO,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO_0 ) ) *
          exp( -d_t / d_tau )
      )
    );

    zone.updateProperty(
      nnt::s_T9,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9_0 ) ) *
          exp( -d_t / ( 3. * d_tau ) )
      )
    );

  }

  if(
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) < 1.e-10
  )
    zone.updateProperty(
      nnt::s_T9,
      "1.e-10"
    );

  if(
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) < 1.e-30
  )
    zone.updateProperty(
      nnt::s_RHO,
      "1.e-30"
    );

}

//##############################################################################
// set_zone().
//##############################################################################

int
set_zone( Libnucnet * p_nucnet, nnt::Zone& zone, char ** argv )
{

  if( !argv )
  {
    std::cerr << "Invalid input." << std::endl;
    return 0;
  }

  if( !Libnucnet__getZoneByLabels( p_nucnet, argv[3], "0", "0" ) )
  {
    std::cerr << "Zone not found!" << std::endl;
    return 0;
  }

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_nucnet, argv[3], "0", "0" )
  );

  return 1;

}

#endif // HYDRO_SHOCK

} // namespace user

