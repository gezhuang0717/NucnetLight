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
//! \brief A header file for computation of approximate weak rates according to
//!        Arcones et al., A&A 522, A25 (2010)
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_AA522A25_H
#define NNT_AA522A25_H

#include <iostream>

#include <Libnucnet.h>

#include "nnt/string_defs.h"
#include "nnt/param_defs.h"
#include "nnt/auxiliary.h"
#include "nnt/two_d_weak_rates.h"

namespace user
{

const double d_AA522A25_WEAK_K    = 6144;         /* Hardy & Towner 2009 */
const size_t i_AA522A25_WORKSPACE = 1000;
const double d_AA522A25_EPSABS    = 1.e-12;
const double d_AA522A25_EPSREL    = 1.e-12;
const size_t i_AA522A25_LIMIT     = 1000;
const double d_AA522A25_EC_CUTOFF = 1.e-100;
const double d_AA522A25_EC_LIMIT  = 1.01;

const char s_AA522A25_B_FACTOR[]            = "aa522a25 B factor";
const char s_AA522A25_BETA_MINUS[]          = "aa522a25 beta minus";
const char s_AA522A25_BETA_PLUS[]           = "aa522a25 beta plus";
const char s_AA522A25_ELECTRON_CAPTURE[]    = "aa522a25 electron capture";
const char s_AA522A25_POSITRON_CAPTURE[]    = "aa522a25 positron capture";

typedef struct
{
  double dEtaF;
  double dMunuekT;
  double dQ;
  double dkT;
  int iNeutrinoExponent;
} aa522a25_function_data;

void
aa522a25__update_net(
  Libnucnet__Net *,
  const char *
);

void
aa522a25__update_net(
  Libnucnet__Net *
);

int
aa522a25__update_reactions(
  Libnucnet__Species *,
  Libnucnet__Net *
);

double
aa522a25__compute_rate(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double,
  double,
  double,
  double
);

double
aa522a25__compute_electron_capture_integral(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double,
  double,
  double,
  double,
  int
);

double
aa522a25__compute_beta_plus_integral(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double,
  double,
  double,
  double,
  int
);

double
aa522a25__compute_positron_capture_integral(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double,
  double,
  double,
  double,
  int
);

double
aa522a25__compute_beta_minus_integral(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double,
  double,
  double,
  double,
  int
);

double
aa522a25__electron_capture_integrand(
  double,
  void *
);

double
aa522a25__beta_plus_integrand(
  double,
  void *
);

double
aa522a25__positron_capture_integrand(
  double,
  void *
);

double
aa522a25__beta_minus_integrand(
  double,
  void *
);

double
aa522a25__compute_reaction_neutrino_energy_loss_rate(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double,
  double,
  double,
  double
);

double
aa522a25__compute_reaction_average_neutrino_energy(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double,
  double,
  double,
  double
);

} // namespace user

#endif // NNT_AA522A25_H
