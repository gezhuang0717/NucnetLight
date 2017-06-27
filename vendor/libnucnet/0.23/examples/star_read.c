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
//       Example to demonstrate how to use libnucnet routines to create
//       a mass fraction xml file from a Rauscher et al. (2002) stellar model.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet.h>

#define I_NUMBER 5000
#define I_WIDTH  12

int main( int argc, char **argv )
{

  Libnucnet *p_nucnet;
  Libnucnet__Zone *p_zone;
  Libnucnet__Species *p_species;

  FILE *p_input;
  char **s_strings;
  char **s_nucname;
  double d_mass_below, d_mass, d_x;
  size_t i, i_strings, i_species;
  char s_zone[10];
  char c = ' ';
  char c_prev = ' ';
  char s_c[2];
  char s_mass_below[256], s_mass[256];

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc != 4 ) {
      fprintf(
        stderr,
        "\nUsage: %s nuc_filename filename new_xml_filename\n\n",
        argv[0]
      );
      fprintf(
        stderr, "  nuc_filename = nuclear data xml filename.\n\n"
      );
      fprintf(
        stderr, "  filename = input stellar model text file.\n\n"
      );
      fprintf(
        stderr, "  new_xml_filename = output xml filename.\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Create libnucnet structure and update nuclear data.
  //==========================================================================*/

  p_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc(
      Libnucnet__getNet( p_nucnet )
    ),
    argv[1],
    NULL
  );

  /*============================================================================
  // Allocate memory for strings.
  //==========================================================================*/

  s_strings = ( char ** ) malloc( I_NUMBER * I_WIDTH * sizeof( char * ) );
  s_nucname = ( char ** ) malloc( I_NUMBER * I_WIDTH * sizeof( char * ) );

  /*============================================================================
  // Open stellar model file.
  //==========================================================================*/

  p_input = fopen( argv[2], "r" );

  /*============================================================================
  // Read strings.
  //==========================================================================*/

  i_strings = 0;

  do
  {

    if( ( c != ' ' ) && ( c_prev == ' ' ) )
    {
      s_strings[i_strings] = ( char * ) malloc( I_WIDTH * sizeof( char * ) );
      sprintf( s_strings[i_strings++], "%c", c );
    }
    else if( ( c != ' ' ) && ( c_prev != ' ' ) )
    {
      sprintf( s_c, "%c", c );
      strcat( s_strings[i_strings-1], s_c );
    }

    c_prev = c;

    fscanf( p_input, "%c", &c );

  } while( c != '\n' );

  /*============================================================================
  // Assign species names.  Rename neutron.
  //==========================================================================*/

  i_species = 0;
  for( i = 0; i < i_strings; i++ )
  {
    if(
        strcmp( s_strings[i], "zone" ) &&
        strcmp( s_strings[i], "mass" ) &&
        strcmp( s_strings[i], "below" )
    )
    {
      s_nucname[i_species] = ( char * ) malloc( I_WIDTH * sizeof( char * ) );
      strcpy( s_nucname[i_species++], s_strings[i] );
    }
  }

  strcpy( s_nucname[0], "n" );

  /*============================================================================
  // Done with s_strings.
  //==========================================================================*/

  for( i = 0; i < i_strings; i++ )
    free( s_strings[i] );

  free( s_strings );

  /*============================================================================
  // Read in mass fraction data.
  //==========================================================================*/

  while( !feof( p_input ) )
  {

    if( 
      fscanf( p_input, "%s %lf %lf", s_zone, &d_mass_below, &d_mass )
      == EOF
    )
      break;

    sprintf( s_mass_below, "%g", d_mass_below );
    sprintf( s_mass, "%g", d_mass );

    p_zone =
      Libnucnet__Zone__new(
        Libnucnet__getNet( p_nucnet ),
        s_zone,
        NULL,
        NULL
      );

    Libnucnet__Zone__updateProperty(
      p_zone,
      "mass below",
      NULL,
      NULL,
      s_mass_below
    );

    Libnucnet__Zone__updateProperty(
      p_zone,
      "mass",
      NULL,
      NULL,
      s_mass
    );

    for( i = 0; i < i_species; i++ )
    {
      fscanf( p_input, "%lf", &d_x );
      if( d_x > 0. )
      {
        p_species =
          Libnucnet__Nuc__getSpeciesByName( 
            p_nucnet->pNet->pNuc, s_nucname[i]
          );
        Libnucnet__Zone__updateSpeciesAbundance(
          p_zone, p_species, d_x / Libnucnet__Species__getA( p_species )
        );
      }
    }

    if( !Libnucnet__addZone( p_nucnet, p_zone ) )
    {
      fprintf( stderr, "Couldn't add zone!\n" );
      return EXIT_FAILURE;
    }

  }

  fclose( p_input );

  /*============================================================================
  // Set format code for zone xml mass fraction output.  Default code
  // is "%g".
  //==========================================================================*/

  Libnucnet__updateZoneXmlMassFractionFormat( p_nucnet, "%.4e" );

  /*============================================================================
  // Write xml file.
  //==========================================================================*/

  Libnucnet__writeZoneDataToXmlFile( p_nucnet, argv[3] );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  for( i = 0; i < i_species; i++ )
    free( s_nucname[i] );

  free( s_nucname );

  Libnucnet__free( p_nucnet );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return EXIT_SUCCESS;

}
