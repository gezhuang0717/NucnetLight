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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file for routines to compute weak rates depending
//!        on temperature and mass density of electrons.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_TWO_D_WEAK_RATES_H
#define NNT_TWO_D_WEAK_RATES_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <Libnucnet.h>
#include <Libstatmech.h>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/math.h"

namespace nnt
{

class TwoDWeakQuantity
{

  public:
    TwoDWeakQuantity( Libnucnet__Reaction *, const char * );
    TwoDWeakQuantity( const TwoDWeakQuantity& );
    ~TwoDWeakQuantity();
    std::pair<double,double> computeValue( double, double );
    std::string getReactionString() const { return sReaction; }
    gsl_vector * getT9Vector();
    gsl_vector * getLog10RhoeVector();
    gsl_matrix * getMatrix();

  private:
    Libnucnet__Reaction * pReaction;
    std::string sReaction;
    gsl_vector * pT9Vector;
    gsl_vector * pLog10RhoeVector;
    gsl_matrix * pMatrix;

};

double
ffnIV_I( double, double, double, double );

double
ffnIV_Ie( double, double, double, double, double );

double
ffnIV__fermi_dirac( int, double );

double
ffnIV__compute_Ie(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double,
  double,
  double,
  double
);

} // namespace nnt

#endif /* NNT_TWO_D_WEAK_RATES_H */
