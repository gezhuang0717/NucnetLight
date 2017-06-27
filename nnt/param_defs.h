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
//! \brief A header file to define NucNet Tools parameters.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_PARAM_DEFS_H
#define NNT_PARAM_DEFS_H

/*##############################################################################
// Define some parameters.
//############################################################################*/

namespace nnt
{

  const double d_EQ_DIFF       = 1.e-15;
  const double d_FLOW_MIN      = 1.e-10;
  const double d_EXP_LARGE     = 600.;
  const double d_Y_MIN_PRINT   = 1.e-30;
  const double d_REL_EPS       = 1.e-08;
  const double d_STEP_FRACTION = 1.e-3;

  const size_t i_BUF_SIZE      = 32;

  const double d_ELECTRON_MASS_IN_MEV = \
    GSL_CONST_CGSM_MASS_ELECTRON * \
    GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT * \
    GSL_CONST_NUM_MICRO / GSL_CONST_CGSM_ELECTRON_VOLT;

} // namespace nnt

#endif /* NNT_PARAM_DEFS_H */
