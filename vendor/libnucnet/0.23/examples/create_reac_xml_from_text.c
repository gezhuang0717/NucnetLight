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
//       Example to demonstrate how to use libnucnet routines to create a
//       Libnucnet__Reac structure of nuclear reactions, add reactions to the
//       structure, write the data to an xml file,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet__Reac.h>

#define MAX_SIZE 256
#define USER_SUPPLIED "user_supplied_fit"

int main( int argc, char **argv ) {

  Libnucnet__Reac *p_my_reactions;
  Libnucnet__Reaction *p_reaction, *p_copy;
  FILE *p_input_file;
  char s_reaction_type[MAX_SIZE], s_tmp[MAX_SIZE], s_source[MAX_SIZE];
  char s_fit_note[MAX_SIZE], s_key[MAX_SIZE];
  char s_name[MAX_SIZE], s_tag1[MAX_SIZE], s_tag2[MAX_SIZE], s_value[MAX_SIZE];
  char *s_end;
  size_t i, i_number, i_number_points, i_number_fits;
  double d_t9, d_rate, d_sef;
  gsl_vector *p_t9, *p_rate, *p_sef;
  int i_z1, i_a1, i_z2, i_a2, i_sef_flag;
  unsigned int i_number_tags;
  double *a, d_tmp, d_spint, d_spinf, d_TlowHf, d_Tlowfit, d_acc;
  double d_Thighfit = 10.;
  
  if ( argc != 3 ) {
      fprintf(
        stderr, "\nUsage: %s input_file output_file \n", argv[0]
      );
      fprintf(
        stderr, "\n  input_file: name of input text file\n"
      );
      fprintf(
        stderr, "\n  output_file: name of output xml file\n\n"
      );
      return EXIT_FAILURE;
  }

  /*============================================================================
  // Create reaction list.
  //==========================================================================*/

  p_my_reactions = Libnucnet__Reac__new( );

  /*============================================================================
  // Open input file.
  //==========================================================================*/

  if( ( p_input_file = fopen( argv[1], "r" ) ) == NULL ) {
    fprintf( stderr, "Could not open file.\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Read in input.
  //==========================================================================*/

  while( !feof( p_input_file ) ) {

  /*----------------------------------------------------------------------------
  // Create reaction.
  //--------------------------------------------------------------------------*/

    p_reaction = Libnucnet__Reaction__new();

  /*----------------------------------------------------------------------------
  // Read in reaction type.
  //--------------------------------------------------------------------------*/

    fscanf( p_input_file, "%s\n", s_reaction_type );

  /*----------------------------------------------------------------------------
  // If user-supplied rate type, get key and register function (NULL is ok to
  // simply register key).
  //--------------------------------------------------------------------------*/

    if( strcmp( s_reaction_type, USER_SUPPLIED ) == 0 )
    {
      fgets( s_key, MAX_SIZE, p_input_file );

      s_end = strchr( s_key, '\n' );
      if( s_end ) {
        *s_end = '\0';
      }

      Libnucnet__Reaction__setUserRateFunctionKey(
        p_reaction,
        s_key
      );
    }

  /*----------------------------------------------------------------------------
  // Read in reaction source.
  //--------------------------------------------------------------------------*/

    fgets( s_source, MAX_SIZE, p_input_file );

    s_end = strchr( s_source, '\n' );
    if( s_end ) {
      *s_end = '\0';
    }

    Libnucnet__Reaction__updateSource( p_reaction, s_source );

  /*----------------------------------------------------------------------------
  // Read in reactants.
  //--------------------------------------------------------------------------*/

    fscanf( p_input_file, "%lu\n", (unsigned long *) &i_number );

    for( i = 0; i < i_number; i++ )
    {
      fscanf( p_input_file, "%s\n", s_tmp );
      Libnucnet__Reaction__addReactant( p_reaction, s_tmp );
    }

  /*----------------------------------------------------------------------------
  // Read in products.
  //--------------------------------------------------------------------------*/

    fscanf( p_input_file, "%lu\n", (unsigned long *) &i_number );
    
    for( i = 0; i < i_number; i++ )
    {
      fscanf( p_input_file, "%s\n", s_tmp );
      Libnucnet__Reaction__addProduct( p_reaction, s_tmp );
    }

  /*----------------------------------------------------------------------------
  // Read in and store rate data.
  //--------------------------------------------------------------------------*/

    /*..........................................................................
    // Single rate.
    //........................................................................*/

    if( strcmp( s_reaction_type, "single_rate" ) == 0 ) {

        fscanf( p_input_file, "%lf\n", &d_rate );
        Libnucnet__Reaction__updateSingleRate(
          p_reaction, d_rate
        );

    }

    /*..........................................................................
    // Rate table.
    //........................................................................*/

    else if( strcmp( s_reaction_type, "rate_table" ) == 0 ) {

        fscanf( p_input_file, "%lu\n", (unsigned long *) &i_number_points );

        p_t9 = gsl_vector_calloc( i_number_points );
        p_rate = gsl_vector_calloc( i_number_points );
        p_sef = gsl_vector_calloc( i_number_points );

        for( i = 0; i < i_number_points; i++ )
        {

          fscanf(
            p_input_file,
            "%lf%lf%lf\n",
            &d_t9, &d_rate, &d_sef
          );
          
          gsl_vector_set( p_t9, i, d_t9 );
          gsl_vector_set( p_rate, i, d_rate );
          gsl_vector_set( p_sef, i, d_sef );

        }

        Libnucnet__Reaction__updateRateTable(
          p_reaction,
          p_t9,
          p_rate,
          p_sef
        );

        gsl_vector_free( p_t9 );
        gsl_vector_free( p_rate );
        gsl_vector_free( p_sef );

    }

    /*..........................................................................
    // Non Smoker fit.
    //........................................................................*/

    else if( strcmp( s_reaction_type, "non_smoker_fit" ) == 0 ) {

      fscanf(
        p_input_file,
        "%lu\n",
        (unsigned long *) &i_number_fits
      );

      for( i = 0; i < i_number_fits; i++ ) {

        fgets( s_fit_note, MAX_SIZE, p_input_file );

        s_end = strchr( s_fit_note, '\n' );
        if( s_end ) {
          *s_end = '\0';
        }

        fscanf(
          p_input_file,
          "%d%d%d%d%lf%lf%lf%lf%lf%lf%d\n",
          &i_z1, &i_a1, &i_z2, &i_a2, &d_tmp, &d_spint, &d_spinf, &d_TlowHf,
          &d_Tlowfit, &d_acc, &i_sef_flag
        );

        a = ( double * ) malloc( sizeof( double ) * 8 );

        fscanf( p_input_file, "%lf%lf%lf%lf\n", &a[0], &a[1], &a[2], &a[3] );
        fscanf( p_input_file, "%lf%lf%lf%lf\n", &a[4], &a[5], &a[6], &a[7]  );

        Libnucnet__Reaction__addNonSmokerFit(
          p_reaction, 
          s_fit_note,
          a,
          d_spint,
          d_spinf,
          d_TlowHf,
          d_Tlowfit,
          d_Thighfit,
          d_acc
        );

      }

    }

    /*..........................................................................
    // User fit.
    //........................................................................*/

    else if( strcmp( s_reaction_type, USER_SUPPLIED ) == 0 )
    {

      fscanf(
        p_input_file,
        "%lu %u\n",
        (unsigned long *) &i_number_fits,
        &i_number_tags
      );

      for( i = 0; i < i_number_fits; i++ )
      {
        switch( i_number_tags )
        {
          case 0:
            fscanf(
              p_input_file,
              "%s",
              s_name
            );
            fgets( s_value, MAX_SIZE, p_input_file );

            s_end = strchr( s_value, '\n' );
            if( s_end ) {
              *s_end = '\0';
            }

            if(
              !Libnucnet__Reaction__updateUserRateFunctionProperty(
                 p_reaction,
                 s_name,
                 NULL,
                 NULL,
                 s_value
              )
            )
            {
              fprintf( stderr, "Couldn't update property.\n" );
              return EXIT_FAILURE;
            }
            break;
          case 1:
            fscanf(
              p_input_file,
              "%s %s",
              s_name,
              s_tag1
            );
            fgets( s_value, MAX_SIZE, p_input_file );

            s_end = strchr( s_value, '\n' );
            if( s_end ) {
              *s_end = '\0';
            }

            if(
              !Libnucnet__Reaction__updateUserRateFunctionProperty(
                 p_reaction,
                 s_name,
                 s_tag1,
                 NULL,
                 s_value
              )
            )
            {
              fprintf( stderr, "Couldn't update property.\n" );
              return EXIT_FAILURE;
            }
            break;
          case 2:
            fscanf(
              p_input_file,
              "%s %s %s\n",
              s_name,
              s_tag1,
              s_tag2
            );
            fgets( s_value, MAX_SIZE, p_input_file );

            s_end = strchr( s_value, '\n' );
            if( s_end ) {
              *s_end = '\0';
            }

            if(
              !Libnucnet__Reaction__updateUserRateFunctionProperty(
                 p_reaction,
                 s_name,
                 s_tag1,
                 s_tag2,
                 s_value
              )
            )
            {
              fprintf( stderr, "Couldn't update property.\n" );
              return EXIT_FAILURE;
            }
            break;
          default:
            fprintf( stderr, "Invalid number of tags.\n" );
            return EXIT_FAILURE;
        }
            
      }
    }

    else {

        fprintf( stderr, "%s\n", s_reaction_type );
        fprintf( stderr, "No such reaction!\n" );
        return EXIT_FAILURE;

    }

    if( !Libnucnet__Reac__addReaction( p_my_reactions, p_reaction ) )
    {
      fprintf(
        stderr, "Couldn't add reaction %s\n",
        Libnucnet__Reaction__getString( p_reaction )
      );
      return EXIT_FAILURE;
    }

  /*----------------------------------------------------------------------------
  // For illustration, copy reaction and update p_my_reactions with it.
  //--------------------------------------------------------------------------*/

    p_copy = Libnucnet__Reaction__copy( p_reaction );

    if( !Libnucnet__Reac__updateReaction( p_my_reactions, p_copy ) )
    {
      fprintf(
        stderr, "Couldn't update reaction %s\n",
        Libnucnet__Reaction__getString( p_reaction )
      );
      return EXIT_FAILURE;
    }

    fscanf( p_input_file, "\n" );
 
  }
  
  /*============================================================================
  // Close input file.
  //==========================================================================*/

  fclose( p_input_file );

  /*============================================================================
  // Write to XML.
  //==========================================================================*/

  Libnucnet__Reac__writeToXmlFile( p_my_reactions, argv[2] );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnucnet__Reac__free( p_my_reactions );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

