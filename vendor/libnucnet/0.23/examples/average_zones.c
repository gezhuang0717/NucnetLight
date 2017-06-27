/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Bradley S. Meyer.
//
//     This is free software; you can redistribute it and/or modify it
//     under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     Please see the src/README.txt file in this distribution for more
//     information.
//   </license>
//   <description>
//     <abstract>
//       Example to demonstrate how to use libnucnet routines to parse
//       in zone data from an xml file, average the abundances over
//       those zones, save the average to a new file,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet.h>

int
main( int argc, char **argv ) {

  FILE *p_file;
  double d_mix, d_mix_total = 0.;
  char s_label_1[256], s_label_2[256], s_label_3[256];
  gsl_vector *p_abunds, *p_average;
  Libnucnet *p_libnucnet, *p_libnucnet_average;
  Libnucnet__Zone *p_zone;

  if ( argc != 6 ) {
      fprintf(
        stderr, "\nUsage: %s nuc_filename zone_filename average_filename type new_filename\n\n", argv[0]
      );
      fprintf(
        stderr, "  nuc_filename = input nuclear data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  zone_filename = input zone data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  average_filename = text file with average data.\n\n"
      );
      fprintf(
        stderr,
        "  type = type of averaging to do (by \"number\" or \"mass\")\n\n"
      );
      fprintf(
        stderr, "  new_filename = name of xml file to dump subset to.\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Check averaging.
  //==========================================================================*/

  if(
    strcmp( argv[4], "number" ) != 0 &&
    strcmp( argv[4], "mass" ) != 0
  )
  {
    fprintf( stderr, "\nAveraging type must be \"number\" or \"mass\".\n\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Create a libnucnet structure.
  //==========================================================================*/

  p_libnucnet = Libnucnet__new();

  /*============================================================================
  // Get the nuclei.
  //==========================================================================*/

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc(
      Libnucnet__getNet( p_libnucnet )
    ),
    argv[1],
    NULL
  );

  /*============================================================================
  // Parse in the zone data.
  //==========================================================================*/

  Libnucnet__assignZoneDataFromXml( p_libnucnet, argv[2], NULL );

  /*============================================================================
  // Get average vector.
  //==========================================================================*/

  p_average =
    gsl_vector_calloc( 
      Libnucnet__Nuc__getNumberOfSpecies(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_libnucnet ) )
      )
    );

  /*============================================================================
  // Read in averaging data.
  //==========================================================================*/

  p_file = fopen( argv[3], "r" );

  if( !p_file ) {
    fprintf( stderr, "File not found!\n" );
    return EXIT_FAILURE;
  }

  while( !feof( p_file ) )
  {

    fscanf(
      p_file,
      "%s %s %s %lf\n",
      s_label_1,
      s_label_2,
      s_label_3,
      &d_mix
    );

    p_zone =
      Libnucnet__getZoneByLabels(
        p_libnucnet,
        s_label_1,
        s_label_2,
        s_label_3
      );

    if( strcmp( argv[4], "mass" ) == 0 )
      d_mix *=
        atof( Libnucnet__Zone__getProperty( p_zone, "mass", NULL, NULL ) );
   
    d_mix_total += d_mix;

    p_abunds = Libnucnet__Zone__getAbundances( p_zone );

    gsl_vector_scale( p_abunds, d_mix );

    gsl_vector_add( p_average, p_abunds );

    gsl_vector_free( p_abunds );

  }

  fclose( p_file );

  /*============================================================================
  // Get average.
  //==========================================================================*/

  gsl_vector_scale( p_average, 1. / d_mix_total );

  /*============================================================================
  // Assign average to new structure.
  //==========================================================================*/

  p_libnucnet_average = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc(
      Libnucnet__getNet( p_libnucnet_average )
    ),
    argv[1],
    NULL
  );

  p_zone =
    Libnucnet__Zone__new(
      Libnucnet__getNet( p_libnucnet_average ),
      "average",
      NULL,
      NULL
    );

  Libnucnet__addZone( p_libnucnet_average, p_zone );

  Libnucnet__Zone__updateAbundances( p_zone, p_average );

  /*============================================================================
  // Write the average to a new xml file.
  //==========================================================================*/

  Libnucnet__writeZoneDataToXmlFile( p_libnucnet_average, argv[5] );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  gsl_vector_free( p_average );

  Libnucnet__free( p_libnucnet_average );
  Libnucnet__free( p_libnucnet );
  
  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}
