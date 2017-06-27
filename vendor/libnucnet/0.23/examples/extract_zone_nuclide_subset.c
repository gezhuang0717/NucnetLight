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
//       in zone data from an xml file, select a subset of the nuclei within
//       those zones by xpath, save the subset to a new file,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet.h>

/*##############################################################################
// Prototypes.
//############################################################################*/

int
zone_iterator( Libnucnet__Zone *, Libnucnet * );

int
species_iterator( Libnucnet__Species *, void * );

void
add_optional_property(
  const char *, const char *, const char *, const char *, Libnucnet__Zone *
);

/*##############################################################################
// main()
//############################################################################*/

int
main( int argc, char **argv ) {

  Libnucnet *p_libnucnet, *p_libnucnet_subset;

  if ( argc != 6 ) {
      fprintf(
        stderr, "\nUsage: %s nuc_filename zone_filename zone_xpath nuc_xpath new_filename\n\n", argv[0]
      );
      fprintf(
        stderr, "  nuc_filename = input nuclear data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  zone_filename = input zone data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  zone_xpath = xpath expression for zone subset.\n\n"
      );
      fprintf(
        stderr, "  nuc_xpath = xpath expression for nuclide subset.\n\n"
      );
      fprintf(
        stderr, "  new_filename = name of file to dump subset.\n\n"
      );
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

  Libnucnet__assignZoneDataFromXml( p_libnucnet, argv[2], argv[3] );

  /*============================================================================
  // Create a new libnucnet structure.
  //==========================================================================*/

  p_libnucnet_subset = Libnucnet__new();

  Libnucnet__Net__updateFromXml(
    Libnucnet__getNet( p_libnucnet_subset ),
    argv[1],
    argv[4],
    NULL
  );

  /*============================================================================
  // Iterate and assign.
  //==========================================================================*/

  Libnucnet__iterateZones(
    p_libnucnet,
    (Libnucnet__Zone__iterateFunction) zone_iterator,
    p_libnucnet_subset
  );

  /*============================================================================
  // Write the subset to a new xml file.
  //==========================================================================*/

  Libnucnet__writeZoneDataToXmlFile( p_libnucnet_subset, argv[5] );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnucnet__free( p_libnucnet_subset );
  Libnucnet__free( p_libnucnet );
  
  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// zone_iterator()
//############################################################################*/

int
zone_iterator( Libnucnet__Zone *p_old_zone, Libnucnet *p_new )
{

  typedef struct {
    Libnucnet__Zone *pOldZone;
    Libnucnet__Zone *pNewZone;
  } user_data;

  const char *s_label1, *s_label2, *s_label3;
  Libnucnet__Zone *p_new_zone;
  user_data *p_user_data;

  s_label1 = Libnucnet__Zone__getLabel( p_old_zone, 1 );
  s_label2 = Libnucnet__Zone__getLabel( p_old_zone, 2 );
  s_label3 = Libnucnet__Zone__getLabel( p_old_zone, 3 );

  if( !strcmp( s_label1, "0" ) )
    p_new_zone = Libnucnet__Zone__new( p_new->pNet, NULL, NULL, NULL ); 
  else
    if( !strcmp( s_label2, "0" ) )
       p_new_zone = Libnucnet__Zone__new( p_new->pNet, s_label1, NULL, NULL );
    else
       if( !strcmp( s_label3, "0" ) )
         p_new_zone = Libnucnet__Zone__new( p_new->pNet, s_label1, s_label2, NULL );
       else
         p_new_zone =
           Libnucnet__Zone__new( p_new->pNet, s_label1, s_label2, s_label3 );

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  p_user_data->pOldZone = p_old_zone;
  p_user_data->pNewZone = p_new_zone;

  Libnucnet__Zone__iterateOptionalProperties(
    p_old_zone,
    NULL,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
       add_optional_property,
    p_user_data->pNewZone
  );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_new ) ),
    (Libnucnet__Species__iterateFunction) species_iterator,
    p_user_data
  );

  Libnucnet__addZone( p_new, p_new_zone );

  free( p_user_data );

  return 1;

}

/*##############################################################################
// species_iterator()
//############################################################################*/

int
species_iterator( Libnucnet__Species *p_species, void *p_data )
{

  typedef struct {
    Libnucnet__Zone *pOldZone; 
    Libnucnet__Zone *pNewZone; 
  } user_data;

  Libnucnet__Species *p_species_in_old;
  user_data *p_user_data;

  p_user_data = ( user_data * ) p_data;

  p_species_in_old =
    Libnucnet__Nuc__getSpeciesByName(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( p_user_data->pOldZone )
      ),
      Libnucnet__Species__getName( p_species )
    );

  Libnucnet__Zone__updateSpeciesAbundance(
    p_user_data->pNewZone,
    p_species,
    Libnucnet__Zone__getSpeciesAbundance(
      p_user_data->pOldZone, p_species_in_old
    )
  );

  return 1;

}

/*##############################################################################
// add_optional_property()
//############################################################################*/

void
add_optional_property(
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  const char *s_value,
  Libnucnet__Zone *p_zone
)
{

  Libnucnet__Zone__updateProperty(
    p_zone,
    s_name,
    s_tag1,
    s_tag2,
    s_value
  );

}

