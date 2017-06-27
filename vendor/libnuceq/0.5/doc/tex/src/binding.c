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
//       Compute the two-species minimum mass per nucleon.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet__Nuc.h>

typedef struct {
  Libnucnet__Nuc *pNuc;
  double dBinding;
  double dYe;
  double dY1;
  double dY2;
  Libnucnet__Species *pSpecies;
  Libnucnet__Species *pSpecies1;
  Libnucnet__Species *pSpecies2;
} work;

int max_binding( Libnucnet__Species *, void * );

int pair_binding( Libnucnet__Species *, void * );

int
main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  work *p_work;

  if ( argc != 3 && argc != 4 ) {
      fprintf(
        stderr, "\nUsage: %s filename ye xpath\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  ye = Ye for calculation.\n\n"
      );
      fprintf(
        stderr, "  xpath = xpath expression (optional).\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Create nuclear species collection.
  //==========================================================================*/

  if( argc == 3 ) {
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );
  } else {
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], argv[3] );
  }
    
  /*============================================================================
  // Create.
  //==========================================================================*/

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work )
  {
    fprintf( stderr, "Couldn't allocate memory!\n" );
    return EXIT_FAILURE;
  }

  p_work->pNuc = p_my_nuclei;

  p_work->dBinding = 1.e300;

  p_work->dYe = atof( argv[2] );

  Libnucnet__Nuc__iterateSpecies(
    p_my_nuclei,
    (Libnucnet__Species__iterateFunction) max_binding,
    p_work
  );

  fprintf(
    stdout,
    "%s  %s  %s  %e  %e  %e\n",
    argv[2],
    Libnucnet__Species__getName( p_work->pSpecies1 ),
    Libnucnet__Species__getName( p_work->pSpecies2 ),
    p_work->dY1,
    p_work->dY2,
    p_work->dBinding
  );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  free( p_work );

  Libnucnet__Nuc__free( p_my_nuclei );
  
  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// max_binding().
//############################################################################*/

int
max_binding( Libnucnet__Species *p_species, void *p_data )
{

  work *p_work = ( work * ) p_data;

  p_work->pSpecies = p_species;

  Libnucnet__Nuc__iterateSpecies(
    p_work->pNuc,
    (Libnucnet__Species__iterateFunction) pair_binding,
    p_work
  );    

  return 1; 

}

/*##############################################################################
// pair_binding()
//############################################################################*/

int
pair_binding(
  Libnucnet__Species *p_species,
  void *p_data
)
{

  double d_bind, d_m1, d_m2, d_ye1, d_ye2, d_y1, d_y2;
  int i_denom;
  work * p_work = ( work * ) p_data;

  if(
    Libnucnet__Species__getIndex( p_species ) >
    Libnucnet__Species__getIndex( p_work->pSpecies )
  )
    return 1;

  d_m1 =
    WN_AMU_TO_MEV * Libnucnet__Species__getA( p_species ) +
    Libnucnet__Species__getMassExcess( p_species );

  d_m2 =
    WN_AMU_TO_MEV * Libnucnet__Species__getA( p_work->pSpecies ) +
    Libnucnet__Species__getMassExcess( p_work->pSpecies );

  d_ye1 =
    (double) Libnucnet__Species__getZ( p_species ) /
      (double) Libnucnet__Species__getA( p_species );

  d_ye2 =
    (double) Libnucnet__Species__getZ( p_work->pSpecies ) /
      (double) Libnucnet__Species__getA( p_work->pSpecies );

  if( p_species == p_work->pSpecies || d_ye1 == d_ye2 )
  {
    if( d_ye1 == p_work->dYe )
    {
      d_y1 = 1. / Libnucnet__Species__getA( p_species );
      d_bind = d_m1 * d_y1;
      d_y2 = 0.;
    }
    else return 1;
  }
  else
  {
    i_denom =
      Libnucnet__Species__getZ( p_work->pSpecies ) *
      Libnucnet__Species__getA( p_species )
      -
      Libnucnet__Species__getZ( p_species ) *
      Libnucnet__Species__getA( p_work->pSpecies );

    d_y1 =
      (
        Libnucnet__Species__getZ( p_work->pSpecies ) -
        p_work->dYe * Libnucnet__Species__getA( p_work->pSpecies ) 
      ) /
      (double) i_denom;
    d_y2 =
      (
        p_work->dYe * Libnucnet__Species__getA( p_species ) -
        Libnucnet__Species__getZ( p_species )
      ) /
      (double) i_denom;
    if( d_y1 < 0. || d_y2 < 0. ) return 1;
    if( d_y1 > 1. || d_y2 > 1. ) return 1;
    d_bind = d_m1 * d_y1 + d_m2 * d_y2;
  }

  if( d_bind < p_work->dBinding )
  {
    if( d_y1 >= d_y2 )
    {
      p_work->dY1 = d_y1;
      p_work->dY2 = d_y2;
      p_work->pSpecies1 = p_species;
      p_work->pSpecies2 = p_work->pSpecies;
    }
    else
    {
      p_work->dY1 = d_y2;
      p_work->dY2 = d_y1;
      p_work->pSpecies1 = p_work->pSpecies;
      p_work->pSpecies2 = p_species;
    }
    p_work->dBinding = d_bind;
  }

  return 1;

}
   
