////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2012 Clemson University.
//
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
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief Code for wrapping zones,
//!        reactions, species, reactants and products, and properties.
//!
////////////////////////////////////////////////////////////////////////////////

#include "nnt/wrappers.h"

/**
 * @brief The NucNet Tools namespace.
 */
namespace nnt
{

//##############################################################################
// ReactionElement methods.
//##############################################################################

/**
 * A method that sets the wrapped Libnucnet__Reaction__Element.
 * \param p_reaction_element A pointer to the element to be wrapped.
 */

void
ReactionElement::setNucnetReactionElement(
  Libnucnet__Reaction__Element * p_reaction_element )
{
  pReactionElement = p_reaction_element;
}

/**
 * A method that returns the wrapped Libnucnet__Reaction__Element.
 * \return A pointer to the wrapped Libnucnet__Reaction__Element.
 */

Libnucnet__Reaction__Element * ReactionElement::getNucnetReactionElement()
{
  return pReactionElement;
}

//##############################################################################
// Reaction methods.
//##############################################################################

/**
 * A method that sets the wrapped Libnucnet__Reaction.
 * \param p_reaction A pointer to the Libnucnet__Reaction to be wrapped.
 */

void Reaction::setNucnetReaction( Libnucnet__Reaction * p_reaction )
{
  pReaction = p_reaction;
}

/**
 * A method that returns the wrapped Libnucnet__Reaction.
 * \return A pointer to the wrapped Libnucnet__Reaction.
 */

Libnucnet__Reaction * Reaction::getNucnetReaction()
{
  return pReaction;
}

//##############################################################################
// Zone methods.
//##############################################################################

/**
 * Zone constructor.
 */

Zone::Zone()
{

  pfGuessFunction = NULL;
  pfMatrixFunction = NULL;
  pfNseCorrectionFactorFunction = NULL;

}

/**
 * A method that sets the wrapped Libnucnet__Zone.
 * \param p_zone A pointer to the Libnucnet__Zone to be wrapped.
 */

void Zone::setNucnetZone( Libnucnet__Zone * p_zone )
{
  pZone = p_zone;
}

/**
 * A method that returns the wrapped Libnucnet__Zone.
 * \return A pointer to the wrapped Libnucnet__Zone.
 */

Libnucnet__Zone * Zone::getNucnetZone()
{
  if( !this )
  {
    std::cerr << "Invalid zone." << std::endl;
    exit( EXIT_FAILURE );
  }
  return pZone;
}

/**
 * A method to set the guess function.
 * \param pf_function A pointer to the guess function.
 * \param p_data A pointer to extra data for the guess function;
 */

void
Zone::setGuessFunction( guessFunction pf_function, void * p_data )
{

  pfGuessFunction = pf_function;
  pGuessFunctionData = p_data;

}

/**
 * A method to set the matrix function.
 * \param pf_function A pointer to the matrix function.
 * \param p_data A pointer to extra data for the matrix function;
 */

void
Zone::setMatrixFunction( matrixFunction pf_function, void * p_data )
{

  pfMatrixFunction = pf_function;
  pMatrixFunctionData = p_data;

}

/**
 * A method to update the NSE correction factor function.
 * \param pf_function A pointer to the function.
 * \param p_data A pointer to extra data for the function;
 */

void
Zone::updateNseCorrectionFactorFunction(
  Libnucnet__Species__nseCorrectionFactorFunction pf_function,
  void * p_data
)
{

  pfNseCorrectionFactorFunction = pf_function;
  pNseCorrectionFactorFunctionData = p_data;

}

/**
 * A method to return a Libnucnet__NetView from a zone.  If the view
 *   does not exist or if the underlying network has been updated, the
 *   routine creates, stores, and returns a new view.  To create the
 *   new view, the routine will use the first label as the nuclear
 *   XPath expression and the second as the reaction XPath expression.
 * \param s_label1 The first label for the view.
 * \param s_label2 The second label for the view (optional).
 * \param s_label3 The third label for the view (optional).
 */

Libnucnet__NetView *
Zone::getNetView(
  const char * s_label1,
  const char * s_label2,
  const char * s_label3
)
{

  Libnucnet__NetView * p_view, * p_new_view;

  p_view =
    Libnucnet__Zone__getNetView(
      this->getNucnetZone(),
      s_label1,
      s_label2,
      s_label3
    );

  if( p_view && !Libnucnet__NetView__wasNetUpdated( p_view ) ) return p_view;

  p_new_view =
    Libnucnet__NetView__new(
      Libnucnet__Zone__getNet( this->getNucnetZone() ),
      s_label1,
      s_label2
    );

  Libnucnet__Zone__updateNetView(
    this->getNucnetZone(),
    s_label1,
    s_label2,
    s_label3,
    p_new_view
  );

  return p_new_view; 

}

Libnucnet__NetView *
Zone::getNetView(
  const char * s_label1,
  const char * s_label2
)
{

  return this->getNetView( s_label1, s_label2, NULL );

}

Libnucnet__NetView *
Zone::getNetView(
  const char * s_label1
)
{

  return this->getNetView( s_label1, NULL );

}

//##############################################################################
// Species methods.
//##############################################################################

/**
 * A method that sets the wrapped Libnucnet__Species.
 * \param p_species A pointer to the Libnucnet__Species to be wrapped.
 */

void Species::setNucnetSpecies( Libnucnet__Species * p_species )
{
  pSpecies = p_species;
}

/**
 * A method that returns the wrapped Libnucnet__Species.
 * \return A pointer to the wrapped Libnucnet__Species.
 */

Libnucnet__Species * Species::getNucnetSpecies()
{
  return pSpecies;
}

} //namespace nnt
