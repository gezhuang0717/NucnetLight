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
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file to define classes to wrap zones,
//!        reactions, species, reactants and products, and properties.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_WRAPPERS_H
#define NNT_WRAPPERS_H

#include <iostream>
#include <string>

#include <Libnucnet.h>
#include <Libstatmech.h>
#include <Libnuceq.h>

/**
 * @brief The NucNet Tools namespace.
 */
namespace nnt
{

  //############################################################################
  // Typedefs. 
  //############################################################################

  typedef double ( *quantityFunction )( double, void * );

  typedef std::pair<double,double> ( *guessFunction )( void * );

  typedef void ( *matrixFunction )( WnMatrix *, gsl_vector *, void * );

  //############################################################################
  // ReactionElement.
  //############################################################################

  /**
   * A class that wraps a Libnucnet__Reaction__Element. Member functions are
   * defined in wrappers.cpp.
   */

  class ReactionElement
  {

    public:
      Libnucnet__Reaction__Element * getNucnetReactionElement();
      void setNucnetReactionElement( Libnucnet__Reaction__Element * );

    protected:
      Libnucnet__Reaction__Element * pReactionElement;

  };

  //############################################################################
  // Reaction.
  //############################################################################

  /**
   * A class that wraps a Libnucnet__Reaction.
   */

  class Reaction
  {

    public:
      Libnucnet__Reaction * getNucnetReaction();
      void setNucnetReaction( Libnucnet__Reaction * );

    protected:
      Libnucnet__Reaction * pReaction;

  };

  //############################################################################
  // Zone.
  //############################################################################

  /**
   * A class that wraps a Libnucnet__Zone.
   */

  class Zone
  {

    public:
      Zone();
      Libnucnet__Zone * getNucnetZone();
      void setNucnetZone( Libnucnet__Zone * );
      std::string *
        createOutReactionXPathString( const char *, const char * );
      void normalizeAbundances();
      void printAbundances();
      std::string getProperty( std::string );
      std::string getProperty( std::string, std::string );
      std::string getProperty( std::string, std::string, std::string );
      int updateProperty( std::string, std::string );
      int updateProperty( std::string, std::string, std::string );
      int updateProperty(
          std::string, std::string, std::string, std::string
        );
      int hasProperty( std::string );
      int hasProperty( std::string, std::string );
      int hasProperty( std::string, std::string, std::string );
      matrixFunction getMatrixFunction() const {return pfMatrixFunction; }
      void * getMatrixFunctionData() const {return pMatrixFunctionData;}
      double computeRootFromQuantity( quantityFunction, double, void * );
      double computeRootFromQuantity( quantityFunction, void * );
      void setWeakDetailedBalance();
      void setGuessFunction( guessFunction, void * );
      void setMatrixFunction( matrixFunction, void * );
      void
        updateNseCorrectionFactorFunction(
          Libnucnet__Species__nseCorrectionFactorFunction,
          void *
        );
      Libnucnet__Species__nseCorrectionFactorFunction
        getNseCorrectionFactorFunction()
          const {return pfNseCorrectionFactorFunction;}
      void * getNseCorrectionFactorFunctionData()
          const {return pNseCorrectionFactorFunctionData;}
      Libnucnet__NetView *
        getNetView( const char *, const char *, const char * );
      Libnucnet__NetView * getNetView( const char *, const char * );
      Libnucnet__NetView * getNetView( const char * );
      double getElectronNeutrinoChemicalPotential();

      bool operator==(const Zone &p) const
      {

        std::string s11 = Libnucnet__Zone__getLabel( pZone, 1 );
        std::string s12 = Libnucnet__Zone__getLabel( pZone, 2 );
        std::string s13 = Libnucnet__Zone__getLabel( pZone, 3 );

        std::string s21 = Libnucnet__Zone__getLabel( p.pZone, 1 );
        std::string s22 = Libnucnet__Zone__getLabel( p.pZone, 2 );
        std::string s23 = Libnucnet__Zone__getLabel( p.pZone, 3 );

        return (s11 == s21) && (s12 == s22) && (s13 == s23 );

      }

      bool operator<(const Zone &p) const
      {

        std::string s11 = Libnucnet__Zone__getLabel( pZone, 1 );
        std::string s12 = Libnucnet__Zone__getLabel( pZone, 2 );
        std::string s13 = Libnucnet__Zone__getLabel( pZone, 3 );

        std::string s21 = Libnucnet__Zone__getLabel( p.pZone, 1 );
        std::string s22 = Libnucnet__Zone__getLabel( p.pZone, 2 );
        std::string s23 = Libnucnet__Zone__getLabel( p.pZone, 3 );

        if( s11 != s21 )
          return s11 < s21;

        if( s12 != s22 )
          return s22 < s22;

        return s13 < s23;

      }

    protected:
      Libnucnet__Zone * pZone;
      void setClusterData( Libnuceq * );
      guessFunction pfGuessFunction;
      void * pGuessFunctionData;
      matrixFunction pfMatrixFunction;
      void * pMatrixFunctionData;
      Libnucnet__Species__nseCorrectionFactorFunction
        pfNseCorrectionFactorFunction;
      void * pNseCorrectionFactorFunctionData;

  };

  //############################################################################
  // Species.
  //############################################################################

  /**
   * A class that wraps a Libnucnet__Species.
   */

  class Species
  {

    public:
      Libnucnet__Species * getNucnetSpecies();
      void setNucnetSpecies( Libnucnet__Species * );

    protected:
      Libnucnet__Species * pSpecies;

  };

} // namespace nnt

#endif // NNT_WRAPPERS_H
