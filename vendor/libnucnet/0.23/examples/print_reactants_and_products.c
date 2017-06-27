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
//       a new Libnucnet__Reac structure of nuclear reactions from an
//       input xml file, print out the reactants and products for a reaction
//       chosen by its string, and clear the structure and free the allocated
//       memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>

#define NUCLIDE "nuclide"

typedef struct {
  int iCount;
} user_data;

int
get_number(
  Libnucnet__Reaction__Element *, 
  user_data *
);

int print_element( Libnucnet__Reaction__Element *, const char * );

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;
  Libnucnet__Reaction *p_reaction;
  char s_nuclide[100], s_null[100];
  user_data *p_user_data;
  int i_number_reactants, i_number_nuclide_reactants;
  int i_number_products, i_number_nuclide_products;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if ( argc != 3 ) {
      fprintf(
        stderr, "\nUsage: %s reac_filename reac_string\n\n", argv[0]
      );
      fprintf(
        stderr, "  reac_filename = input reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  reac_string = reaction string\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read reaction file and exit if not present.
  //==========================================================================*/

  p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], NULL );

  if( !p_my_reactions ) {
    fprintf( stderr, "Reaction data not read!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get reaction.
  //==========================================================================*/

  p_reaction = Libnucnet__Reac__getReactionByString( p_my_reactions, argv[2] );

  if( !p_reaction ) {
    fprintf( stderr, "Reaction not found!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get user data and initialize;
  //==========================================================================*/

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data )
  {
    fprintf( stderr, "Could not allocate memory for user data!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get strings.
  //==========================================================================*/

  strcpy( s_nuclide, NUCLIDE );
  strcpy( s_null, " " );

  /*============================================================================
  // Get number of reactant and products.
  //==========================================================================*/

  p_user_data->iCount = 0;
  Libnucnet__Reaction__iterateReactants(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) get_number,
    p_user_data
  );
  i_number_reactants = p_user_data->iCount;

  p_user_data->iCount = 0;
  Libnucnet__Reaction__iterateNuclideReactants(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) get_number,
    p_user_data
  );
  i_number_nuclide_reactants = p_user_data->iCount;

  p_user_data->iCount = 0;
  Libnucnet__Reaction__iterateProducts(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) get_number,
    p_user_data
  );
  i_number_products = p_user_data->iCount;

  p_user_data->iCount = 0;
  Libnucnet__Reaction__iterateNuclideProducts(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) get_number,
    p_user_data
  );
  i_number_nuclide_products = p_user_data->iCount;

  /*============================================================================
  // Print reaction string.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nFor %s:\n",
    Libnucnet__Reaction__getString( p_reaction )
  );

  /*============================================================================
  // Print reactants.
  //==========================================================================*/

  fprintf(
    stdout,
    "\n%d nuclide ",
    i_number_nuclide_reactants
  );

  if( i_number_nuclide_reactants == 0 )
    fprintf( stdout, "reactants\n" );
  else if( i_number_nuclide_reactants == 1 )
    fprintf( stdout, "reactant:\n" );
  else
    fprintf( stdout, "reactants:\n" );

  Libnucnet__Reaction__iterateReactants(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) print_element,
    s_nuclide
  );

  fprintf(
    stdout,
    "\n%d non-nuclide ",
    i_number_reactants - i_number_nuclide_reactants
  );

  if( i_number_reactants - i_number_nuclide_reactants == 0 )
    fprintf( stdout, "reactants\n" );
  else if( i_number_reactants - i_number_nuclide_reactants == 1 )
    fprintf( stdout, "reactant:\n" );
  else
    fprintf( stdout, "reactants:\n" );

  Libnucnet__Reaction__iterateReactants(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) print_element,
    s_null
  );
    
  /*============================================================================
  // Print products.
  //==========================================================================*/

  fprintf(
    stdout,
    "\n%d nuclide ",
    i_number_nuclide_products
  );

  if( i_number_nuclide_products == 0 )
    fprintf( stdout, "products\n" );
  else if( i_number_nuclide_products == 1 )
    fprintf( stdout, "product:\n" );
  else
    fprintf( stdout, "products:\n" );

  Libnucnet__Reaction__iterateProducts(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) print_element,
    s_nuclide
  );

  fprintf(
    stdout,
    "\n%d non-nuclide ",
    i_number_products - i_number_nuclide_products
  );

  if( i_number_products - i_number_nuclide_products == 0 )
    fprintf( stdout, "products\n" );
  else if( i_number_products - i_number_nuclide_products == 1 )
    fprintf( stdout, "product:\n" );
  else
    fprintf( stdout, "products:\n" );

  Libnucnet__Reaction__iterateProducts(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) print_element,
    s_null
  );

  fprintf( stdout, "\n" );
    
  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  free( p_user_data );
  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

}

int
get_number(
  Libnucnet__Reaction__Element *p_element,
  user_data *p_data
)
{

  if( !p_element )
  {
    fprintf( stderr, "Invalid input.\n" );
    exit( EXIT_FAILURE );
  }

  p_data->iCount++;

  return 1;

}

int
print_element(
  Libnucnet__Reaction__Element *p_element,
  const char *s_type
)
{

  if( strcmp( s_type, NUCLIDE ) == 0 )
  {
    if( Libnucnet__Reaction__Element__isNuclide( p_element ) )
      fprintf(
        stdout,
        " %s\n",
        Libnucnet__Reaction__Element__getName( p_element )
      );
  }
  else if( strcmp( s_type, NUCLIDE ) )
  {
    if( !Libnucnet__Reaction__Element__isNuclide( p_element ) )
      fprintf(
        stdout,
        " %s\n",
        Libnucnet__Reaction__Element__getName( p_element )
      );
  }
  else
  {
    fprintf( stderr, "Invalid type!\n" );
    return 0;
  }

  return 1;

}
