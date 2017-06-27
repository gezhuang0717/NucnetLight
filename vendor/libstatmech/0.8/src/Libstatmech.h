/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//      Copyright (c) 2008-2012 Clemson University.
//
//      This file is part of the Clemson Webnucleo group's
//      libstatmech module, originally developed by Tianhong Yu
//      and Bradley S. Meyer.  For more information,
//      please see http://www.webnucleo.org.
//
//      This is free software; you can redistribute it and/or modify it
//      under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This software is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this software (please see the "gnu_gpl.txt" file in the doc/
//      directory of this distribution); if not, write to the Free Software
//      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
//      USA
//   </license>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#ifndef LIBSTATMECH_H
#define LIBSTATMECH_H

/*##############################################################################
// Use extern "C" for C++ compilers.
//############################################################################*/

#ifdef __cplusplus
extern "C"
{
#endif

/*##############################################################################
// Standard library #include's.
//############################################################################*/

#include <stdio.h>
#include <math.h>
#include <string.h>

/*##############################################################################
// Libxml #includes.
//############################################################################*/

#include <libxml/tree.h>

/*##############################################################################
// Gsl #includes.
//############################################################################*/

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_min.h>

/*##############################################################################
// Integration definitions.
//############################################################################*/

#define D_EPS_ABSOLUTE      1.e-08    /*  Absolute integration epsilon */
#define D_EPS_RELATIVE      1.e-08    /*  Relative integration epsilon */
#define D_GAMMA_ALPHA       1.e-03    /*  Test for -alpha = gamma */
#define I_LIMIT             5000      /*  Max number of subdivisions */
#define I_WORKSPACE         10000L    /*  Workspace allocation */

/*##############################################################################
// Differentiation definitions.
//############################################################################*/

#define D_STEP_FRACTION     1.e-3     /*  Percentage of argument as the step
                                          value */
/*##############################################################################
// Root-finding definitions.
//############################################################################*/

#define D_EPS_ROOT          1.e-12    /*  Root convergence criterion */
#define I_ITER_MAX          1000      /*  Max number of root-finding iters */ 
#define D_BOSON_ALPHA_LIMIT -1.e-300  /*  Upper limit on boson alpha */

/*##############################################################################
// Macro definitions.
//############################################################################*/

/* To debug, uncomment following line. */
/* #define debug_thermo */

#ifndef debug_thermo

#define LIBSTATMECH__ERROR( s_reason ) \
       do { \
         fprintf( stderr, "ERROR: %s in %s on line %d\n", \
           s_reason, __FILE__, __LINE__ \
         ); \
         exit( EXIT_FAILURE ); \
       } while (0)

#else

#define LIBSTATMECH__ERROR( s_reason ) \
       do { \
         fprintf( stderr, "ERROR: %s in %s on line %d\n", \
           s_reason, __FILE__, __LINE__ \
         ); \
         abort(); \
       } while (0)

#endif

/*##############################################################################
// Definitions.
//############################################################################*/

#define  S_NUMBER_DENSITY           "number density"
#define  S_PRESSURE                 "pressure"
#define  S_ENERGY_DENSITY           "energy density"
#define  S_INTERNAL_ENERGY_DENSITY  "internal energy density"
#define  S_ENTROPY_DENSITY          "entropy density"
#define  S_CHEMICAL_POTENTIAL       "chemical potential"

/*##############################################################################
// Fermions.
//############################################################################*/

/*##############################################################################
// <user_routine name="Libstatmech__Fermion__Integrand()">
//
//   <description>
//     <abstract>
//       User-supplied routine to compute the integrand for the desired
//       thermodynamic function of a fermion gas given the
//       temperature and chemical potential/kT.
//     </abstract>
//     <keywords>
//       Libstatmech, integrand, chemical, potential, user,
//       supplied, function
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Fermion__Integrand(
//   Libstatmech__Fermion *self,
//   double d_x,
//   double d_T,
//   double d_mukT,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="d_x"
//       kind="in,positional,required"
//     >
//       A double giving the value of the integration variable x.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the thermodynamic function.
//     </param>
//
//     <param
//       name="d_mukT"
//       kind="in,positional,required"
//     >
//       A double giving the chemical potential divided by kT at which to
//       compute the thermodynamic function.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return a double giving the integrand for
//       the desired thermodynamic function for the input temperature
//       and chemical potential divided by Boltzmann's constant times the
//       temperature.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef double ( *Libstatmech__Fermion__Integrand ) (
  void *, double, double, double, void *
);

/*##############################################################################
// <user_routine name="Libstatmech__Fermion__Function()">
//
//   <description>
//     <abstract>
//       User-supplied routine to compute the desired
//       thermodynamic function of a fermion gas given the temperature and
//       chemical potential/kT.
//     </abstract>
//     <keywords>
//       Libstatmech, function, chemical, potential, user,
//       supplied, function
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Fermion__Function(
//   Libstatmech__Fermion *self,
//   double d_T,
//   double d_mukT,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the thermodynamic function.
//     </param>
//
//     <param
//       name="d_mukT"
//       kind="in,positional,required"
//     >
//       A double giving the chemical potential divided by kT at which to
//       compute the thermodynamic function.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return a double giving the 
//       desired thermodynamic function for the input temperature
//       and chemical potential divided by Boltzmann's constant times the
//       temperature.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef double ( *Libstatmech__Fermion__Function ) (
  void *, double, double, void *
);

/*##############################################################################
// <class name="Libstatmech__Fermion">
//
//   <description>
//     <abstract>
//       Libstatmech__Fermion is a structure that stores information
//       about a fermion.  The contents of the structure are not made
//       public by the API but rather are accessed by API routines.
//     </abstract>
//     <keywords>
//       statistical, mechanics, fermion, boson
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="tyu" start_date="2008/05/02"/>
//       <author userid="mbradle" start_date="2008/05/02"/>
//     </current>
//     <previous>
//     </previous>
//   </authors>
//   
//   <compatibility>
//     Compiled successfully in gcc (GCC) 4.1.1.
//   </compatibility>
//   
// </class>
//############################################################################*/

typedef struct Libstatmech__Fermion {
  xmlChar *sxName;
  double dRestMass;
  int iMultiplicity;
  double dCharge;
  xmlHashTablePtr pWorkHash;
} Libstatmech__Fermion;

/*##############################################################################
// Bosons. 
//############################################################################*/

/*##############################################################################
// <user_routine name="Libstatmech__Boson__Integrand()">
//
//   <description>
//     <abstract>
//       User-supplied routine to compute the integrand for the desired
//       thermodynamic function of a boson gas given the
//       temperature and chemical potential/kT.
//     </abstract>
//     <keywords>
//       Libstatmech, integrand, chemical, potential, user,
//       supplied, function
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Boson__Integrand(
//   Libstatmech__Boson *self,
//   double d_x,
//   double d_T,
//   double d_mukT,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="d_x"
//       kind="in,positional,required"
//     >
//       A double giving the value of the integration variable x.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the thermodynamic function.
//     </param>
//
//     <param
//       name="d_mukT"
//       kind="in,positional,required"
//     >
//       A double giving the chemical potential divided by kT at which to
//       compute the thermodynamic function.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return a double giving the integrand for
//       the desired thermodynamic function for the input temperature
//       and chemical potential divided by Boltzmann's constant times the
//       temperature.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef double ( *Libstatmech__Boson__Integrand ) (
  void *, double, double, double, void *
);

/*##############################################################################
// <user_routine name="Libstatmech__Boson__Function()">
//
//   <description>
//     <abstract>
//       User-supplied routine to compute the desired
//       thermodynamic function of a boson gas given the temperature and
//       chemical potential/kT.
//     </abstract>
//     <keywords>
//       Libstatmech, function, chemical, potential, user,
//       supplied, function
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Boson__Function(
//   Libstatmech__Boson *self,
//   double d_T,
//   double d_mukT,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the thermodynamic function.
//     </param>
//
//     <param
//       name="d_mukT"
//       kind="in,positional,required"
//     >
//       A double giving the chemical potential divided by kT at which to
//       compute the thermodynamic function.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return a double giving the 
//       desired thermodynamic function for the input temperature
//       and chemical potential divided by Boltzmann's constant times the
//       temperature.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef double ( *Libstatmech__Boson__Function ) (
  void *, double, double, void *
);

/*##############################################################################
// <class name="Libstatmech__Boson">
//
//   <description>
//     <abstract>
//       Libstatmech__Boson is a structure that stores information
//       about a boson.  The contents of the structure are not made public
//       by the API but rather are accessed through API routines.
//     </abstract>
//     <keywords>
//       statistical, mechanics, boson
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="tyu" start_date="2008/05/02"/>
//       <author userid="mbradle" start_date="2008/05/02"/>
//     </current>
//     <previous>
//     </previous>
//   </authors>
//   
//   <compatibility>
//     Compiled successfully in gcc (GCC) 4.1.1.
//   </compatibility>
//   
// </class>
//############################################################################*/

typedef struct Libstatmech__Boson {
  xmlChar *sxName;
  double dRestMass;
  int iMultiplicity;
  double dCharge;
  xmlHashTablePtr pWorkHash;
} Libstatmech__Boson;

/*##############################################################################
// Other type definitiones. 
//############################################################################*/

typedef double ( *Libstatmech__integratorFunction ) (
  void *, double, double, double, void *
);

typedef double ( *Libstatmech__function ) (
  void *, double, double, void *
);

typedef struct Libstatmech__Work {
  double dT;
  double dAlpha;
  Libstatmech__Fermion *pFermion;
  Libstatmech__Boson *pBoson;
  Libstatmech__integratorFunction pfIntegrand;
  Libstatmech__function pfFunction;
  void *pFunctionData;
  void *pIntegrandData;
  double dIntegralLowerLimit;
  double dIntegralUpperLimit;
  double dEpsAbsolute;
  double dEpsRelative;
} Libstatmech__Work;

/*##############################################################################
// API Routines.
//############################################################################*/

/*##############################################################################
// Routines for Fermions.
//############################################################################*/

/*##############################################################################
// <routine name="Libstatmech__Fermion__new()">
//
//   <description>
//     <abstract>
//       Creates a new Libstatmech__Fermion structure.
//     </abstract>
//     <keywords>
//       Libstatmech, fermion
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libstatmech__Fermion
// *Libstatmech__Fermion__new(
//    const char *s_fermion_name,
//    double d_rest_mass,
//    int i_multiplicity,
//    double d_charge
// );
//     </calling_sequence>
//
//     <param
//       name="s_fermion_name"
//       kind="in,positional,required"
//     >
//       A string giving the name of the new fermion.
//     </param>
//
//     <param
//       name="d_rest_mass"
//       kind="in,positional,required"
//     >
//       A double giving the rest mass energy of the fermion in MeV.
//     </param>
//
//     <param
//       name="i_multiplicity"
//       kind="in,positional,required"
//     >
//       An int giving the multiplicity of the fermion.  This is usually
//       twice the spin plus unity (that is, 2J+1, where J is the spin).
//     </param>
//
//     <param
//       name="d_charge"
//       kind="in,positional,required"
//     >
//       A double giving the charge of the fermion.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new Libstatmech__Fermion
//       structure. If the routine cannot allocate sufficient memory,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Create a new fermion, an electron, with mass 0.511 Mev,
//         multiplicity 2, and charge -1.:
//       </synopsis>
//
//       <code>
// p_electron =
//   Libstatmech__Fermion__new(
//     "electron", 0.511, 2, -1.
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libstatmech__Fermion
*Libstatmech__Fermion__new( const char *, double, int, double );

/*##############################################################################
// <routine name="Libstatmech__Fermion__getName()">
//
//   <description>
//     <abstract>
//       Returns the name of a fermion.
//     </abstract>
//     <keywords>
//       Libstatmech, fermion
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// const char *
// Libstatmech__Fermion__getName(
//    const Libstatmech__Fermion *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the name of the fermion.
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the name of the fermion pointed to by p_fermion:
//       </synopsis>
//
//       <code>
// printf(
//   "The name of the fermion is %s\n",
//   Libstatmech__Fermion__getName( p_fermion )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

const char *Libstatmech__Fermion__getName( const Libstatmech__Fermion * );

/*##############################################################################
// <routine name="Libstatmech__Fermion__getRestMass()">
//
//   <description>
//     <abstract>
//       Returns the rest mass of a fermion.
//     </abstract>
//     <keywords>
//       Libstatmech, fermion
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Fermion__getRestMass(
//    const Libstatmech__Fermion *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the rest mass of the fermion.
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the rest mass of the fermion pointed to by p_fermion:
//       </synopsis>
//
//       <code>
// printf(
//   "The rest mass of the fermion is %e\n",
//   Libstatmech__Fermion__getRestMass( p_fermion )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double Libstatmech__Fermion__getRestMass( const Libstatmech__Fermion * );

/*##############################################################################
// <routine name="Libstatmech__Fermion__getMultiplicity()">
//
//   <description>
//     <abstract>
//       Returns the multiplicity of a fermion.
//     </abstract>
//     <keywords>
//       Libstatmech, fermion
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Fermion__getMultiplicity(
//    const Libstatmech__Fermion *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the multiplicity of the fermion.
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the multiplicity of the fermion pointed to by p_fermion:
//       </synopsis>
//
//       <code>
// printf(
//   "The multiplicity of the fermion is %d\n",
//   Libstatmech__Fermion__getMultiplicity( p_fermion )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int Libstatmech__Fermion__getMultiplicity( const Libstatmech__Fermion * );

/*##############################################################################
// <routine name="Libstatmech__Fermion__getCharge()">
//
//   <description>
//     <abstract>
//       Returns the charge of a fermion.
//     </abstract>
//     <keywords>
//       Libstatmech, fermion
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Fermion__getCharge(
//    const Libstatmech__Fermion *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the charge of the fermion.
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the charge of the fermion pointed to by p_fermion:
//       </synopsis>
//
//       <code>
// printf(
//   "The charge of the fermion is %e\n",
//   Libstatmech__Fermion__getCharge( p_fermion )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libstatmech__Fermion__getCharge( const Libstatmech__Fermion * );

/*##############################################################################
// <routine name="Libstatmech__Fermion__free()">
//
//   <description>
//     <abstract>
//      Free the memory allocated for a Libstatmech fermion. 
//     </abstract>
//     <keywords>
//       Libstatmech, fermion, free
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void 
// Libstatmech__Fermion__free(
//    Libstatmech__Fermion *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//
//     <doc kind="post" id="result">
//       Upon successful return, the memory for the fermion has been freed. 
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//        Free the memory for the fermion pointed to by p_fermion: 
//       </synopsis>
//
//       <code>
// Libstatmech__Fermion__free( p_fermion );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void Libstatmech__Fermion__free( Libstatmech__Fermion * );

/*##############################################################################
// <routine name="Libstatmech__Fermion__computeIntegrandValue()">
//
//   <description>
//     <abstract>
//       Routine to compute the thermodynamic integrand values of a fermion
//       gas given the name of thermodynamic quantity, temperature
//       and chemical potential/kT.
//     </abstract>
//     <keywords>
//       Libstatmech, integrand, temperature, chemical, potential
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Fermion__computeIntegrandValue(
//   Libstatmech__Fermion *self,
//   const char *s_integrand_name,
//   double d_x,
//   double d_T,
//   double d_mukT,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="s_integrand_name"
//       kind="in,positional,required"
//     >
//       A string giving the name associated to the integrand.
//     </param>
//
//     <param
//       name="d_x"
//       kind="in,positional,required"
//     >
//       A double giving the value of the integration variable x.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the thermodynamic integrand.
//     </param>
//
//     <param
//       name="d_mukT"
//       kind="in,positional,required"
//     >
//       A double giving the chemical potential divided by kT at which to
//       compute the thermodynamic integrand. 
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a double giving the integrand for the input
//       thermodynamic quantity at the input integration variable (x) value,
//       temperature, and chemical potential divided by Boltzmann's constant
//       times the temperature.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Compute the pressure integrand of a fermion gas with d_T = 1.e7 K 
//         and d_mukT = 1 at integration variable value d_x:
//       </synopsis>
//
//       <code>
// Libstatmech__Fermion__computeIntegrandValue(
//   p_fermion, 
//   S_PRESSURE,
//   d_x,
//   1.e7,
//   1.,
//   NULL
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libstatmech__Fermion__computeIntegrandValue(
  Libstatmech__Fermion *,
  const char *,
  double,
  double,
  double,
  void *
);

/*##############################################################################
// <routine name="Libstatmech__Fermion__computeQuantity()">
//
//   <description>
//     <abstract>
//       Routine to compute the thermodynamic quantity values of a fermion
//       gas given the name of thermodynamic quantity, temperature
//       and chemical potential/kT.
//     </abstract>
//     <keywords>
//       Libstatmech, thermodynamic, quantity, temperature, chemical,
//       potential
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Fermion__computeQuantity(
//   Libstatmech__Fermion *self,
//   const char *s_quantity_name,
//   double d_T,
//   double d_mukT,
//   void *p_function_data,
//   void *p_integrand_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string qiving the quantity name.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the thermodynamic quantity.
//     </param>
//
//     <param
//       name="d_mukT"
//       kind="in,positional,required"
//     >
//       A double giving the chemical potential divided by kT at which to
//       compute the thermodynamic quantity.
//     </param>
//
//     <param
//       name="p_function_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic function.
//     </param>
//
//     <param
//       name="p_integrand_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic integrand.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a double giving the thermodynamic quantity value
//       for the input thermodynamic quantity name,
//       temperature and chemical potential divided by Boltzmann's constant
//       times the temperature.  If input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Compute the pressure of a fermion gas with d_T = 1.e7 K 
//         and d_mukT = 1.: 
//       </synopsis>
//
//       <code>
// Libstatmech__Fermion__computeQuantity(
//   p_fermion, 
//   S_PRESSURE,
//   1.e7,
//   1.,
//   NULL,
//   NULL
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libstatmech__Fermion__computeQuantity(
  Libstatmech__Fermion *,
  const char *,
  double,
  double,
  void *,
  void *
);

/*##############################################################################
// <routine name="Libstatmech__Fermion__computeChemicalPotential()">
//
//   <description>
//     <abstract>
//       Routine to compute the chemical potential divided by kT of a fermion
//       gas given the temperature and number density.
//     </abstract>
//     <keywords>
//       Libstatmech, chemical, potential, temperature, number, density
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Fermion__computeChemicalPotential(
//   Libstatmech__Fermion *self,
//   double d_T,
//   double d_number_density,
//   void *p_function_data,
//   void *p_integrand_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the chemical potential.
//     </param>
//
//     <param
//       name="d_number_density"
//       kind="in,positional,required"
//     >
//       A double giving the number density at which to
//       compute the chemical potential.
//     </param>
//
//     <param
//       name="p_function_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the already integrated part of the thermodynamic quantity.
//     </param>
//
//     <param
//       name="p_integrand_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the integrand part of the thermodynamic quantity.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a double giving the chemical potential divided by kT
//       for the input temperature and number density. 
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the chemical potential divided by kT
//         of a fermion gas with d_T = 1.e7 K 
//         and d_number_density = 1.e25 cm^-3: 
//       </synopsis>
//
//       <code>
// fprintf(
//   stdout,
//   "mu / kT = %e\n",
//   Libstatmech__Fermion__computeChemicalPotential(
//     p_fermion, 
//     1.e7,
//     1.e25,
//     NULL,
//     NULL
//   )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libstatmech__Fermion__computeChemicalPotential(
  Libstatmech__Fermion *, double, double, void *, void *
);

/*##############################################################################
// <routine name="Libstatmech__Fermion__computeTemperatureDerivative()">
//
//   <description>
//     <abstract>
//       Routine to compute the temperature derivatives of input thermodynamic
//       quatities of a fermion gas given the temperature and number density.
//     </abstract>
//     <keywords>
//       Libstatmech, derivatives, temperature, number, density
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Fermion__computeTemperatureDerivative(
//   Libstatmech__Fermion *self,
//   const char *s_function,
//   double d_T,
//   double d_number_density,
//   void *p_function_data,
//   void *p_integrand_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="s_function"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the chemical potential.
//     </param>
//
//     <param
//       name="d_number_density"
//       kind="in,positional,required"
//     >
//       A double giving the number density at which to
//       compute the chemical potential.
//     </param>
//
//     <param
//       name="p_function_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the already integrated part of the thermodynamic quantity.
//     </param>
//
//     <param
//       name="p_integrand_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the integrand part of the thermodynamic quantity.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a double giving the temperature derivative of the 
//       input thermodynamic quantity, at the input temperature and
//       number density.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Compute the temperature derivative of pressure 
//         of a fermion gas with d_T = 1.e7 K 
//         and d_number_density = 1.e25 cm^-3: 
//       </synopsis>
//
//       <code>
// Libstatmech__Fermion__computeTemperatureDerivative(
//   p_fermion, 
//   S_PRESSURE,
//   1.e7,
//   1.e25,
//   NULL,
//   NULL
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libstatmech__Fermion__computeTemperatureDerivative(
  Libstatmech__Fermion *,
  const char *,
  double,
  double,
  void *,
  void *
); 

/*##############################################################################
// <routine name="Libstatmech__Fermion__updateIntegralLowerLimit()">
//
//   <description>
//     <abstract>
//       Routine to update the integral lower limit of a thermodynamic
//       quantity of a fermion gas. 
//     </abstract>
//     <keywords>
//       Libstatmech, update, integral, lower, limit 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Fermion__updateIntegralLowerLimit(
//   Libstatmech__Fermion *self,
//   const char *s_quantity_name,
//   double d_lower_limit 
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_lower_limit"
//       kind="in,positional,required"
//     >
//       A double giving the integral lower limit.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns 1 (true) if the update succeeded and 0 (false) if
//       not. Upon successful return, the integral lower limit has been
//       updated. If the input thermodynamic quantity is not valid,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Update the integral lower limit of number density to d_x0 = 20.: 
//       </synopsis>
//
//       <code>
// d_x0 = 20.;
// Libstatmech__Fermion__updateIntegralLowerLimit(
//   p_fermion, 
//   S_NUMBER_DENSITY,
//   d_x0 
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
Libstatmech__Fermion__updateIntegralLowerLimit(
  Libstatmech__Fermion *,
  const char *,
  double
);

/*##############################################################################
// <routine name="Libstatmech__Fermion__updateQuantityIntegralAccuracy()">
//
//   <description>
//     <abstract>
//       Routine to update the accuracy parameters for the integral
//       of a thermodynamic quantity of a fermion gas. 
//     </abstract>
//     <keywords>
//       Libstatmech, update, integral, accuracy
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Fermion__updateQuantityIntegralAccuracy(
//   Libstatmech__Fermion *self,
//   const char *s_quantity_name,
//   double d_eps_absolute,
//   double d_eps_relative
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_eps_absolute"
//       kind="in,positional,required"
//     >
//       A double giving the absolute accuracy desired for the integral.
//       The default is 1.e-8.
//     </param>
//
//     <param
//       name="d_eps_relative"
//       kind="in,positional,required"
//     >
//       A double giving the relative accuracy desired for the integral.
//       The default is 1.e-8.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns 1 (true) if the update succeeded and 0 (false)
//       if not. Upon successful return, the absolute and relative accuracies
//       for the integral have been updated. 
//       If the input thermodynamic quantity is not valid,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Update the integral absolute accuracy to 1.e-4 and the relative
//         accuracy to 1.e-12 for the number density: 
//       </synopsis>
//
//       <code>
// if(
//    Libstatmech__Fermion__updateQuantityIntegralAccuracy(
//      p_fermion, 
//      S_NUMBER_DENSITY,
//      1.e-4,
//      1.e-12
//    )
// )
//   fprintf(stdout, "Update succeeded.\n" );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
Libstatmech__Fermion__updateQuantityIntegralAccuracy(
  Libstatmech__Fermion *,
  const char *,
  double,
  double
);

/*##############################################################################
// Routines for Bosons.
##############################################################################*/

/*##############################################################################
// <routine name="Libstatmech__Boson__new()">
//
//   <description>
//     <abstract>
//       Creates a new Libstatmech__Boson structure.
//     </abstract>
//     <keywords>
//       Libstatmech, boson 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libstatmech_Boson
// *Libstatmech__Boson__new(
//    const char *s_boson_name,
//    double d_rest_mass,
//    int i_multiplicity,
//    double d_charge
// );
//     </calling_sequence>
//
//     <param
//       name="s_boson_name"
//       kind="in,positional,required"
//     >
//       A string giving the name of the new boson.
//     </param>
//
//     <param
//       name="d_rest_mass"
//       kind="in,positional,required"
//     >
//       A double giving the rest mass energy of the boson in MeV.
//     </param>
//
//     <param
//       name="i_multiplicity"
//       kind="in,positional,required"
//     >
//       An int giving the multiplicity of the boson.  This is usually
//       twice the spin plus unity (that is, 2J+1, where J is the spin).
//     </param>
//
//     <param
//       name="d_charge"
//       kind="in,positional,required"
//     >
//       A double giving the charge of the boson.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new Libstatmech_Boson
//       structure. If the routine cannot allocate sufficient memory,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Create a new boson, a photon, with mass 0. Mev,
//         multiplicity 2, and charge 0.:
//       </synopsis>
//
//       <code>
// p_photon =
//   Libstatmech__Boson__new(
//     "photon", 0., 2, 0.
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libstatmech__Boson
*Libstatmech__Boson__new( const char *, double, int, double );

/*##############################################################################
// <routine name="Libstatmech__Boson__getName()">
//
//   <description>
//     <abstract>
//       Returns the name of a boson.
//     </abstract>
//     <keywords>
//       Libstatmech, boson
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// const char *
// Libstatmech__Boson__getName(
//    const Libstatmech__Boson *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Boson structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the name of the boson.
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the name of the boson pointed to by p_boson:
//       </synopsis>
//
//       <code>
// printf(
//   "The name of the boson is %s\n",
//   Libstatmech__Boson__getName( p_boson )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

const char *Libstatmech__Boson__getName( const Libstatmech__Boson * );

/*##############################################################################
// <routine name="Libstatmech__Boson__getRestMass()">
//
//   <description>
//     <abstract>
//       Returns the rest mass of a boson.
//     </abstract>
//     <keywords>
//       Libstatmech, boson
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Boson__getRestMass(
//    const Libstatmech__Boson *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Boson structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the rest mass of the boson.
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the rest mass of the boson pointed to by p_boson:
//       </synopsis>
//
//       <code>
// printf(
//   "The rest mass of the boson is %e\n",
//   Libstatmech__Boson__getRestMass( p_boson )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double Libstatmech__Boson__getRestMass( const Libstatmech__Boson * );

/*##############################################################################
// <routine name="Libstatmech__Boson__getMultiplicity()">
//
//   <description>
//     <abstract>
//       Returns the multiplicity of a boson.
//     </abstract>
//     <keywords>
//       Libstatmech, boson
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Boson__getMultiplicity(
//    const Libstatmech__Boson *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Boson structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the multiplicity of the boson.
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the multiplicity of the boson pointed to by p_boson:
//       </synopsis>
//
//       <code>
// printf(
//   "The multiplicity of the boson is %d\n",
//   Libstatmech__Boson__getMultiplicity( p_boson )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int Libstatmech__Boson__getMultiplicity( const Libstatmech__Boson * );

/*##############################################################################
// <routine name="Libstatmech__Boson__getCharge()">
//
//   <description>
//     <abstract>
//       Returns the charge of a boson.
//     </abstract>
//     <keywords>
//       Libstatmech, boson
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Boson__getCharge(
//    const Libstatmech__Boson *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Boson structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the charge of the boson.
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the charge of the boson pointed to by p_boson:
//       </synopsis>
//
//       <code>
// printf(
//   "The charge of the boson is %e\n",
//   Libstatmech__Boson__getCharge( p_boson )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double Libstatmech__Boson__getCharge( const Libstatmech__Boson * );

/*##############################################################################
// <routine name="Libstatmech__Boson__free()">
//
//   <description>
//     <abstract>
//      Free the memory allocated for a Libstatmech boson. 
//     </abstract>
//     <keywords>
//       Libstatmech, boson, free
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void 
// Libstatmech__Boson__free(
//    Libstatmech__Boson *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Boson structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//
//     <doc kind="post" id="result">
//       Upon successful return, the memory for the boson has been freed. 
//       If the input pointer is not valid, Libstatmech error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//        Free the memory for the boson pointed to by p_boson: 
//       </synopsis>
//
//       <code>
// Libstatmech__Boson__free( p_boson );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void Libstatmech__Boson__free( Libstatmech__Boson * );

/*##############################################################################
// <routine name="Libstatmech__Boson__computeIntegrandValue()">
//
//   <description>
//     <abstract>
//       Routine to compute the thermodynamic integrand values of a boson
//       gas given the name of thermodynamic quantity, temperature
//       and chemical potential/kT.
//     </abstract>
//     <keywords>
//       Libstatmech, integrand, temperature, chemical, potential
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Boson__computeIntegrandValue(
//   Libstatmech__Boson *self,
//   const char *s_quantity_name, 
//   double d_x,
//   double d_T,
//   double d_mukT,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_x"
//       kind="in,positional,required"
//     >
//       A double giving the value of the integration variable x.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the thermodynamic integrand.
//     </param>
//
//     <param
//       name="d_mukT"
//       kind="in,positional,required"
//     >
//       A double giving the chemical potential divided by kT at which to
//       compute the thermodynamic integrand. 
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a double giving the integrand for the input
//       thermodynamic quantity at the input integration variable (x) value,
//       temperature, and chemical potential divided by Boltzmann's constant
//       times the temperature.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Compute the pressure integrand of a boson gas with d_T = 1.e5 K 
//         and d_mukT = 1 at integration variable value d_x: 
//       </synopsis>
//
//       <code>
// Libstatmech__Boson__computeIntegrandValue(
//   p_boson, 
//   S_PRESSURE,
//   d_x,
//   1.e5,
//   1.,
//   NULL
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libstatmech__Boson__computeIntegrandValue(
  Libstatmech__Boson *,
  const char *,
  double,
  double,
  double,
  void *
);

/*##############################################################################
// <routine name="Libstatmech__Boson__computeQuantity()">
//
//   <description>
//     <abstract>
//       Routine to compute the thermodynamic quantity values of a boson
//       gas given the name of thermodynamic quantity, temperature
//       and chemical potential/kT.
//     </abstract>
//     <keywords>
//       Libstatmech, thermodynamic, quantity, temperature, chemical,
//       potential
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Boson__computeQuantity(
//   Libstatmech__Boson *self,
//   const char *s_quantity_name,
//   double d_T,
//   double d_mukT,
//   void *p_function_data,
//   void *p_integrand_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the thermodynamic quantity.
//     </param>
//
//     <param
//       name="d_mukT"
//       kind="in,positional,required"
//     >
//       A double giving the chemical potential divided by kT at which to
//       compute the thermodynamic quantity.
//     </param>
//
//     <param
//       name="p_function_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic function.
//     </param>
//
//     <param
//       name="p_integrand_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic integrand.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a double giving the thermodynamic quantity value
//       for the input thermodynamic quantity name,
//       temperature and chemical potential divided by Boltzmann's constant
//       times the temperature.  If input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Compute the pressure of a boson gas with d_T = 1.e5 K 
//         and d_mukT = 1.: 
//       </synopsis>
//
//       <code>
// Libstatmech__Boson__computeQuantity(
//   p_boson, 
//   S_PRESSURE,
//   1.e5,
//   1.,
//   NULL,
//   NULL
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libstatmech__Boson__computeQuantity(
  Libstatmech__Boson *,
  const char *,
  double,
  double,
  void *,
  void *
);

/*##############################################################################
// <routine name="Libstatmech__Boson__computeChemicalPotential()">
//
//   <description>
//     <abstract>
//       Routine to compute the chemical potential divided by kT of a boson
//       gas given the temperature and number density.
//     </abstract>
//     <keywords>
//       Libstatmech, chemical, potential, temperature, number, density
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Boson__computeChemicalPotential(
//   Libstatmech__Boson *self,
//   double d_T,
//   double d_number_density,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the chemical potential.
//     </param>
//
//     <param
//       name="d_number_density"
//       kind="in,positional,required"
//     >
//       A double giving the number density at which to
//       compute the chemical potential.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the thermodynamic function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a double giving the chemical potential divided by kT
//       for the input temperature and number density. 
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the chemical potential divided by kT
//         of a boson gas with d_T = 1.e5 K 
//         and d_number_density = 1.e18:
//       </synopsis>
//
//       <code>
// fprintf(
//   stdout,
//   "mu / kT = %e\n",
//   Libstatmech__Boson__computeChemicalPotential(
//     p_boson, 
//     1.e5,
//     1.e18,
//     NULL
//   )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libstatmech__Boson__computeChemicalPotential(
  Libstatmech__Boson *, double, double, void *, void *
);

/*##############################################################################
// <routine name="Libstatmech__Boson__computeTemperatureDerivative()">
//
//   <description>
//     <abstract>
//       Routine to compute the temperature derivatives of thermodynamic
//       quatities of a boson gas given the temperature and number density.
//     </abstract>
//     <keywords>
//       Libstatmech, derivatives, temperature, number, density
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libstatmech__Boson__computeTemperatureDerivative(
//   Libstatmech__Boson *self,
//   const char *s_function,
//   double d_T,
//   double d_number_density,
//   void *p_function_data,
//   void *p_integrand_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="s_function"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_T"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in K at which to
//       compute the chemical potential.
//     </param>
//
//     <param
//       name="d_number_density"
//       kind="in,positional,required"
//     >
//       A double giving the number density at which to
//       compute the chemical potential.
//     </param>
//
//     <param
//       name="p_function_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the already integrated part of the thermodynamic quantity.
//     </param>
//
//     <param
//       name="p_integrand_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined structure containing extra data for
//       the integrand part of the thermodynamic quantity.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a double giving the temperature derivative value of 
//       the thermodynamic quantity name, at the input temperature and
//       number density.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Compute the temperature derivative of pressure
//         of a boson gas with d_T = 1.e7 K 
//         and d_number_density = 1.e25 cm^-3: 
//       </synopsis>
//
//       <code>
// Libstatmech__Boson__computeTemperatureDerivative(
//   p_boson, 
//   S_PRESSURE,
//   1.e7,
//   1.e25,
//   NULL,
//   NULL
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libstatmech__Boson__computeTemperatureDerivative(
  Libstatmech__Boson *,
  const char *,
  double,
  double,
  void *,
  void *
); 

/*##############################################################################
// <routine name="Libstatmech__Boson__updateIntegralLowerLimit()">
//
//   <description>
//     <abstract>
//       Routine to update the integral lower limit of a thermodynamic
//       quantity of a boson gas. 
//     </abstract>
//     <keywords>
//       Libstatmech, update, integral, lower, limit 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Boson__updateIntegralLowerLimit(
//   Libstatmech__Boson *self,
//   const char *s_quantity_name,
//   double d_lower_limit 
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_lower_limit"
//       kind="in,positional,required"
//     >
//       A double giving the integral lower limit.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns 1 (true) if the update succeeded and 0 (false) if
//       not. Upon successful return, the integral lower limit has been
//       updated. If the input thermodynamic quantity is not valid,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Update the integral lower limit of number density to d_x0 = 20.: 
//       </synopsis>
//
//       <code>
// d_x0 = 20.;
// Libstatmech__Boson__updateIntegralLowerLimit(
//   p_boson, 
//   S_NUMBER_DENSITY,
//   d_x0 
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
Libstatmech__Boson__updateIntegralLowerLimit(
  Libstatmech__Boson *,
  const char *,
  double
);

/*##############################################################################
// <routine name="Libstatmech__Fermion__updateQuantity()">
//
//   <description>
//     <abstract>
//       Routine to update a thermodynamic quantity.
//     </abstract>
//     <keywords>
//       Libstatmech, quantity, integrand, function
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Fermion__updateQuantity(
//   Libstatmech__Fermion *self,
//   const char *s_quantity_name,
//   Libstatmech__Fermion__Function pf_function,
//   Libstatmech__Fermion__Integrand pf_integrand
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="pf_function"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Fermion__Function.  If set to
//       NULL, the quantity will be computed only from the integral.
//     </param>
//
//     <param
//       name="pf_integrand"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Fermion__Integrand.  If set to
//       NULL, the quantity will be computed only from the function.
//       If set to DEFAULT_INTEGRAND, the quantity will be computed
//       from the default integral.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns 1 (true) if the quantity update was successful and
//       0 (false) if not.  On successful return, if the quantity previously
//       existed, its function and integrand have been updated, while if the
//       quantity did not previously exist, it has been added.  If any input
//       is not valid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Update the pressure for Libstatmech__Fermion *p_fermion with
//         function my_pressure_function and integrand my_pressure_integrand:
//       </synopsis>
//
//       <code>
// if(
//   Libstatmech__Fermion__updateQuantity(
//     p_fermion, 
//     S_PRESSURE,
//     (Libstatmech__Fermion__Function) my_pressure_function,
//     (Libstatmech__Fermion__Integrand) my_pressure_integrand
//   )
// )
//   fprintf( stdout, "Update succeeded.\n" );
//       </code>
//     </doc>
//
//     <doc kind="example" id="example2">
//       <synopsis>
//         Create a new thermodynamic quantity called "my quantity" for
//         Libstatmech__Fermion *p_fermion with function my_function and integrand
//         my_integrand:
//       </synopsis>
//
//       <code>
// if(
//   Libstatmech__Fermion__updateQuantity(
//     p_fermion, 
//     "my quantity",
//     (Libstatmech__Fermion__Function) my_function,
//     (Libstatmech__Fermion__Integrand) my_integrand
//   )
// )
//   fprintf( stdout, "Update succeeded.\n" );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
Libstatmech__Fermion__updateQuantity(
  Libstatmech__Fermion *,
  const char *,
  Libstatmech__Fermion__Function,
  Libstatmech__Fermion__Integrand
);

/*##############################################################################
// <routine name="Libstatmech__Boson__updateQuantity()">
//
//   <description>
//     <abstract>
//       Routine to update a thermodynamic quantity.
//     </abstract>
//     <keywords>
//       Libstatmech, quantity, integrand, function
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Boson__updateQuantity(
//   Libstatmech__Boson *self,
//   const char *s_quantity_name,
//   Libstatmech__Boson__Function pf_function,
//   Libstatmech__Boson__Integrand pf_integrand
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="pf_function"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Boson__Function.  If set to NULL,
//       the quantity will be computed only from the integral.
//     </param>
//
//     <param
//       name="pf_integrand"
//       kind="in,positional,required"
//     >
//       A pointer to a Libstatmech__Boson__Integrand.  If set to
//       NULL, the quantity will be computed only from the function.
//       If set to DEFAULT_INTEGRAND, the quantity will be computed
//       from the default integral.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns 1 (true) if the quantity update was successful and
//       0 (false) if not.  On successful return, if the quantity previously
//       existed, its function and integrand have been updated, while if the
//       quantity did not previously exist, it has been added.  If any input
//       is not valid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Create a new thermodynamic quantity called "my quantity" for
//         Libstatmech__Boson *p_boson with function my_function and integrand
//         my_integrand:
//       </synopsis>
//
//       <code>
// if(
//   Libstatmech__Boson__updateQuantity(
//     p_boson, 
//     "my quantity",
//     (Libstatmech__Boson__Function) my_function,
//     (Libstatmech__Boson__Integrand) my_integrand
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
Libstatmech__Boson__updateQuantity(
  Libstatmech__Boson *,
  const char *,
  Libstatmech__Boson__Function,
  Libstatmech__Boson__Integrand
);

/*##############################################################################
// <routine name="Libstatmech__Boson__updateQuantityIntegralAccuracy()">
//
//   <description>
//     <abstract>
//       Routine to update the accuracy parameters for the integral
//       of a thermodynamic quantity of a boson gas. 
//     </abstract>
//     <keywords>
//       Libstatmech, update, integral, accuracy
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Boson__updateIntegralLowerLimit(
//   Libstatmech__Boson *self,
//   const char *s_quantity_name,
//   double d_eps_absolute,
//   double d_eps_relative
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_eps_absolute"
//       kind="in,positional,required"
//     >
//       A double giving the absolute accuracy desired for the integral.
//       The default is 1.e-8.
//     </param>
//
//     <param
//       name="d_eps_relative"
//       kind="in,positional,required"
//     >
//       A double giving the relative accuracy desired for the integral.
//       The default is 1.e-8.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns 1 (true) if the update succeeded and 0 (false)
//       if not. Upon successful return, the absolute and relative accuracies
//       for the integral have been updated. 
//       If the input thermodynamic quantity is not valid,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Update the integral absolute accuracy to 1.e-4 and the relative
//         accuracy to 1.e-12 for the number density: 
//       </synopsis>
//
//       <code>
// if(
//    Libstatmech__Boson__updateQuantityIntegralAccuracy(
//      p_boson, 
//      S_NUMBER_DENSITY,
//      1.e-4,
//      1.e-12
//    )
// )
//   fprintf(stdout, "Update succeeded.\n" );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
Libstatmech__Boson__updateQuantityIntegralAccuracy(
  Libstatmech__Boson *,
  const char *,
  double,
  double
);

/*##############################################################################
// <routine name="Libstatmech__Boson__copy()">
//
//   <description>
//     <abstract>
//       Creates a copy of a boson.
//     </abstract>
//     <keywords>
//       Libstatmech, boson 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libstatmech_Boson *
// Libstatmech__Boson__copy(
//   Libstatmech__Boson * self 
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the boson to be copied.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new Libstatmech_Boson
//       structure that has the same rest mass, multiplicity, and charge as
//       the input structure. The new boson has the default quantity
//       functions and integrands.  If the input boson is not valid,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Copy p_old_photon to p_photon:
//       </synopsis>
//
//       <code>
// p_photon =
//   Libstatmech__Boson__copy(
//     p_old_photon
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libstatmech__Boson *
Libstatmech__Boson__copy(
  Libstatmech__Boson *
);

/*##############################################################################
// <routine name="Libstatmech__Fermion__copy()">
//
//   <description>
//     <abstract>
//       Creates a copy of a fermion.
//     </abstract>
//     <keywords>
//       Libstatmech, fermion 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libstatmech_Fermion *
// Libstatmech__Fermion__copy(
//   Libstatmech__Fermion * self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the fermion to be copied.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new Libstatmech_Fermion
//       structure that has the same rest mass, multiplicity, and charge as
//       the input structure. The new fermion has the default quantity
//       functions and integrands.  If the input fermion is not valid,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Copy p_old_fermion to p_fermion:
//       </synopsis>
//
//       <code>
// p_fermion =
//   Libstatmech__Fermion__copy(
//     p_old_fermion
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libstatmech__Fermion *
Libstatmech__Fermion__copy(
  Libstatmech__Fermion *
);

/*##############################################################################
// <routine name="Libstatmech__Fermion__updateIntegralUpperLimit()">
//
//   <description>
//     <abstract>
//       Routine to update the integral upper limit of a thermodynamic
//       quantity of a fermion gas. 
//     </abstract>
//     <keywords>
//       Libstatmech, update, integral, upper, limit 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Fermion__updateIntegralUpperLimit(
//   Libstatmech__Fermion *self,
//   const char *s_quantity_name,
//   double d_upper_limit 
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Fermion structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_upper_limit"
//       kind="in,positional,required"
//     >
//       A double giving the integral upper limit.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns 1 (true) if the update succeeded and 0 (false) if
//       not. Upon successful return, the integral upper limit has been
//       updated. If the input thermodynamic quantity is not valid,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Update the integral upper limit of number density to d_x0 = 20.: 
//       </synopsis>
//
//       <code>
// d_x0 = 20.;
// Libstatmech__Fermion__updateIntegralUpperLimit(
//   p_fermion, 
//   S_NUMBER_DENSITY,
//   d_x0 
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
Libstatmech__Fermion__updateIntegralUpperLimit(
  Libstatmech__Fermion *,
  const char *,
  double
);

/*##############################################################################
// <routine name="Libstatmech__Boson__updateIntegralUpperLimit()">
//
//   <description>
//     <abstract>
//       Routine to update the integral upper limit of a thermodynamic
//       quantity of a boson gas. 
//     </abstract>
//     <keywords>
//       Libstatmech, update, integral, upper, limit 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int 
// Libstatmech__Boson__updateIntegralUpperLimit(
//   Libstatmech__Boson *self,
//   const char *s_quantity_name,
//   double d_upper_limit 
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to the Libstatmech__Boson structure.
//     </param>
//
//     <param
//       name="s_quantity_name"
//       kind="in,positional,required"
//     >
//       A string giving the quantity name.
//     </param>
//
//     <param
//       name="d_upper_limit"
//       kind="in,positional,required"
//     >
//       A double giving the integral upper limit.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns 1 (true) if the update succeeded and 0 (false) if
//       not. Upon successful return, the integral upper limit has been
//       updated. If the input thermodynamic quantity is not valid,
//       Libstatmech error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Update the integral upper limit of number density to d_x0 = 20.: 
//       </synopsis>
//
//       <code>
// d_x0 = 20.;
// Libstatmech__Boson__updateIntegralUpperLimit(
//   p_boson, 
//   S_NUMBER_DENSITY,
//   d_x0 
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
Libstatmech__Boson__updateIntegralUpperLimit(
  Libstatmech__Boson *,
  const char *,
  double
);

/*##############################################################################
// Non-API Routines.
//############################################################################*/

double
Libstatmech__Fermion__defaultPressureIntegrand(
  Libstatmech__Fermion *, double, double, double, void *
);

double
Libstatmech__Fermion__defaultNumberDensityIntegrand(
  Libstatmech__Fermion *, double, double, double, void *
);

double
Libstatmech__Fermion__defaultNumberDensityFunction(
  Libstatmech__Fermion *, double, double, void *
);

double
Libstatmech__Fermion__defaultEnergyDensityIntegrand(
  Libstatmech__Fermion *, double, double, double, void *
);

double
Libstatmech__Fermion__defaultInternalEnergyDensityIntegrand(
  Libstatmech__Fermion *, double, double, double, void *
);

double
Libstatmech__Fermion__defaultEntropyDensityIntegrand(
  Libstatmech__Fermion *, double, double, double, void *
);

double
Libstatmech__Fermion__numberDensityRootFinder(
  double, void *
);

double
Libstatmech__Fermion__differentiate_helper(
  double,
  void *
);

double
Libstatmech__Boson__defaultPressureIntegrand(
  Libstatmech__Boson *, double, double, double, void *
);

double
Libstatmech__Boson__defaultNumberDensityIntegrand(
  Libstatmech__Boson *, double, double, double, void *
);

double
Libstatmech__Boson__defaultEnergyDensityIntegrand(
  Libstatmech__Boson *, double, double, double, void *
);

double
Libstatmech__Boson__defaultInternalEnergyDensityIntegrand(
  Libstatmech__Boson *, double, double, double, void *
);

double
Libstatmech__Boson__defaultEnergyDensitySquareIntegrand(
  Libstatmech__Boson *, double, double, double, void *
);

double
Libstatmech__Boson__defaultEntropyDensityIntegrand(
  Libstatmech__Boson *, double, double, double, void *
);

double
Libstatmech__Boson__numberDensityRootFinder(
  double, void *
);

double
Libstatmech__Boson__differentiate_helper(
  double,
  void *
);

int
Libstatmech__bracket_root_of_function(
  gsl_function, double *, double *
);

double
Libstatmech__integrate_helper(
  double, void *
);

double
Libstatmech__exp(
  double
);

void
Libstatmech__Work__free(
  Libstatmech__Work *, xmlChar *
);

Libstatmech__Work *
Libstatmech__Work__new( void );

double
Libstatmech__Boson__integrator(
  Libstatmech__Work *
);

double
Libstatmech__Fermion__integrator(
  Libstatmech__Work *
);

double
Libstatmech__integrate_1(
  Libstatmech__Work *,
  double,
  double
);

double
Libstatmech__integrate_2(
  Libstatmech__Work *,
  double,
  double
);

int
Libstatmech__is_standard_quantity( const char * );

Libstatmech__Fermion__Integrand
Libstatmech__Fermion__get_default_integrand( const char * );

Libstatmech__Boson__Integrand
Libstatmech__Boson__get_default_integrand( const char * );

double
DEFAULT_INTEGRAND( void *, double, double, double, void * );

int
Libstatmech__value_is_zero( double );

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBSTATMECH_H */
