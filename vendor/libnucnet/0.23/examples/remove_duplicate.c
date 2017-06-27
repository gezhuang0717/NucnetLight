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
//       Iterate function to remove duplicate reactions in a
//       Libnucnet__Net structure.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "remove_duplicate.h"

int
get_number_reactants( Libnucnet__Reaction__Element *, int * );

int
remove_duplicate(
  Libnucnet__Reaction *p_reaction, Libnucnet__Net *p_net
)
{

  Libnucnet__Reaction *p_reaction1, *p_reaction2;
  double d_q1, d_q2;
  int i_reactants1, i_reactants2;

  /*============================================================================
  // Get reaction pointers in original network structure.  Return if
  // reaction not valid or parent duplicate has already been removed.
  //==========================================================================*/

  p_reaction1 =
    Libnucnet__Reac__getReactionByString(
      Libnucnet__Net__getReac( p_net ),
      Libnucnet__Reaction__getString( p_reaction )
    );

  p_reaction2 =
    Libnucnet__Reac__getReactionByString(
      Libnucnet__Net__getReac( p_net ),
      Libnucnet__Reaction__getString(
        Libnucnet__Reaction__getParentDuplicate( p_reaction )
      )
    );

  if( !Libnucnet__Net__isValidReaction( p_net, p_reaction1 ) )
    return 1;

  if( !p_reaction2 ) return 1;

  /*============================================================================
  // Check whether reaction has only one nuclide reactant.
  //==========================================================================*/

  i_reactants1 = 0;

  Libnucnet__Reaction__iterateReactants(
    p_reaction1,
    (Libnucnet__Reaction__Element__iterateFunction) get_number_reactants,
    &i_reactants1
  );

  i_reactants2 = 0;

  Libnucnet__Reaction__iterateReactants(
    p_reaction2,
    (Libnucnet__Reaction__Element__iterateFunction) get_number_reactants,
    &i_reactants2
  );

  if( i_reactants1 == 1 && i_reactants2 != 1 )
  {
    fprintf(
      stdout,
      "Removing %s (Data source: %s)\n",
      Libnucnet__Reaction__getString( p_reaction1 ),
      Libnucnet__Reaction__getSource( p_reaction1 )
    );
    Libnucnet__Reac__removeReaction(
      Libnucnet__Net__getReac( p_net ),
      p_reaction1
    );
    return 1;
  }

  if( i_reactants1 != 1 && i_reactants2 == 1 )
  {
    fprintf(
      stdout,
      "Removing %s (Data source: %s)\n",
      Libnucnet__Reaction__getString( p_reaction2 ),
      Libnucnet__Reaction__getSource( p_reaction2 )
    );
    Libnucnet__Reac__removeReaction(
      Libnucnet__Net__getReac( p_net ),
      p_reaction2
    );
    return 1;
  }

  /*============================================================================
  // Calculate the reaction Q values.
  //==========================================================================*/

  d_q1 = Libnucnet__Net__computeReactionQValue( p_net, p_reaction1 );

  d_q2 = Libnucnet__Net__computeReactionQValue( p_net, p_reaction2 );

  /*============================================================================
  // Print the Q values if PRINTOUT (set in remove_duplicate.h) is
  // "yes".
  //==========================================================================*/

  if( strcmp( PRINTOUT, "yes" ) == 0 )
  {

    fprintf(
      stdout,
      "\n%s has Q value %g MeV\n",
      Libnucnet__Reaction__getString( p_reaction1 ),
      d_q1
    ); 

    fprintf(
      stdout,
      "%s has Q value %g MeV\n",
      Libnucnet__Reaction__getString( p_reaction2 ),
      d_q2
    ); 

  }

  /*============================================================================
  // Remove duplicate or its parent depending on which is endothermic.
  // If both endothermic or exothermic, remove the duplicate.
  //==========================================================================*/

  if( d_q1 < 0. && d_q2 > 0. )
  {

    if( strcmp( PRINTOUT, "yes" ) == 0 )
    {
      fprintf(
        stdout,
        "Removing %s (Data source: %s)\n",
        Libnucnet__Reaction__getString( p_reaction1 ),
        Libnucnet__Reaction__getSource( p_reaction1 )
      );
    }

    Libnucnet__Reac__removeReaction(
      Libnucnet__Net__getReac( p_net ),
      p_reaction1
    );

  }
  else if( d_q1 > 0. && d_q2 < 0. )
  {

    if( strcmp( PRINTOUT, "yes" ) == 0 )
    {
      fprintf(
        stdout,
        "Removing %s (Data source: %s)\n",
        Libnucnet__Reaction__getString( p_reaction2 ),
        Libnucnet__Reaction__getSource( p_reaction2 )
      );
    }

    Libnucnet__Reac__removeReaction(
      Libnucnet__Net__getReac( p_net ),
      p_reaction2
    );

  }
  else
  {

    if( strcmp( PRINTOUT, "yes" ) == 0 )
    {
      fprintf(
        stdout,
        "Removing %s (Data source: %s)\n",
        Libnucnet__Reaction__getString( p_reaction1 ),
        Libnucnet__Reaction__getSource( p_reaction1 )
      );
    }

    Libnucnet__Reac__removeReaction(
      Libnucnet__Net__getReac( p_net ),
      p_reaction1
    );

  }

  return 1;

}

int
get_number_reactants(
  Libnucnet__Reaction__Element *p_reactant,
  int *p_number
)
{

  if( Libnucnet__Reaction__Element__isNuclide( p_reactant ) )
    (*p_number)++;

  return 1;

}
