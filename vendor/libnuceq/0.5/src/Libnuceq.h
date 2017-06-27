/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//      Copyright (c) 2010-2012 Clemson University.
//
//      This file is part of the Clemson Webnucleo group's
//      libnucnet module, originally developed by Bradley S. Meyer
//      and Tianhong Yu.  For more information,
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
//   <description>
//     <abstract>
//        Source code for the Libnucnet__Eq part of the libnucnet module.
//        For documentation, see the header file Libnucnet__Eq.h.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#ifndef LIBNUCEQ_H
#define LIBNUCEQ_H

/*##############################################################################
// Use extern "C" for C++ compilers.
//############################################################################*/

#ifdef __cplusplus
extern "C"
{
#endif

/*##############################################################################
// Includes.
//############################################################################*/

#include <Libnucnet__Nuc.h>
#include <Libstatmech.h>
#include <gsl/gsl_roots.h>
 
/*##############################################################################
// Macro definitions.
//############################################################################*/

/* To debug, uncomment following line. */
/* #define debug_eq */

#ifndef debug_eq

#define LIBNUCEQ__ERROR( s_reason ) \
       do { \
         fprintf( stderr, "ERROR: %s in %s on line %d\n", \
           s_reason, __FILE__, __LINE__ \
         ); \
         exit( EXIT_FAILURE ); \
       } while (0)

#else

#define LIBNUCEQ__ERROR( s_reason ) \
       do { \
         fprintf( stderr, "ERROR: %s in %s on line %d\n", \
           s_reason, __FILE__, __LINE__ \
         ); \
         abort(); \
       } while (0)

#endif

/*##############################################################################
// Some defines.
//############################################################################*/

#define D_EQ_EPS    1.e-20
#define D_EQ_DIFF   1.e-15
#define D_EQ_EXP_MAX   690.
#define I_EQ_ITER_MAX  1000
#define I_EQ_HASH_MAX  1024

/*##############################################################################
// Some string defines.
//############################################################################*/

#define EQ_NEUT  "n"
#define EQ_PROT  "h1"
#define S_NUCEQ_ELECTRON "electron"

/*##############################################################################
// Forward declarations.
//############################################################################*/

typedef struct _Libnuceq Libnuceq;

typedef struct _Libnuceq__Species Libnuceq__Species;

typedef struct _Libnuceq__Cluster Libnuceq__Cluster;

/*##############################################################################
// Functions.
//############################################################################*/

/*##############################################################################
// <user_routine name="Libnuceq__Cluster_constraint_function()">
//
//   <description>
//     <abstract>
//       User-supplied routine to compute the contribution of a species
//       to a constraint on a cluster.
//     </abstract>
//     <keywords>
//       Libnuceq, user, supplied, function, cluster, constraint
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__Cluster_constraint_function(
//   Libnuceq__Species *p_species,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="p_species"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Species.
//     </param>
//
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to user-supplied extra data for the
//       constraint function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return the contribution of the input species
//       to the cluster constraint.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef
  double
  ( *Libnuceq__Cluster__constraint_function )
  ( Libnuceq__Species *, void * );

/*##############################################################################
// <user_routine name="Libnuceq__Cluster_prefactor_function()">
//
//   <description>
//     <abstract>
//       User-supplied routine to compute the contribution of a species
//       to the prefactor of an equilibrium cluster.
//     </abstract>
//     <keywords>
//       Libnuceq, user, supplied, function, prefactor
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__Cluster_prefactor_function(
//   Libnuceq__Cluster *p_cluster,
//   Libnuceq__Species *p_species,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="p_cluster"
//       kind="in,positional,required"
//     >
//       A pointer to an equilibrium cluster.
//     </param>
//
//     <param
//       name="p_species"
//       kind="in,positional,required"
//     >
//       A pointer to a species in the cluster.
//     </param>
//
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to user-supplied extra data for the
//       prefactor function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return the contribution to the cluster 
//       prefactor.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef
  double
  ( *Libnuceq__Cluster__prefactorFunction )
  ( Libnuceq__Cluster *, Libnuceq__Species *, void * );

/*##############################################################################
// <user_routine name="Libnuceq__wseCorrectionFunction()">
//
//   <description>
//     <abstract>
//       User-supplied routine to compute the extra addend to the weak
//       statistical equilibrium relation between neutrons, protons,
//       and electrons.  This function is equal to mu_p/kT + mu_e/kT - mu_n/kT
//       in weak statistical equilibrium.
//     </abstract>
//     <keywords>
//       Libnuceq, user, supplied, function, weak, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__wseCorrectionFunction(
//   Libnuceq *self,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to an equilibrium.
//     </param>
//
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to user-supplied extra data for the
//       weak equilibrium relation function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return the addend to the weak equilibrium
//       relation.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef
  double
  ( *Libnuceq__wseCorrectionFunction )
  ( Libnuceq *, void * );

/*##############################################################################
// Structures.
//############################################################################*/

/*##############################################################################
// <class name="Libnuceq__Species">
//
//   <description>
//     <abstract>
//       Libnuceq__Species is a structure that stores information about a
//       species in an equilibrium calculation.
//     </abstract>
//     <keywords>
//       isotope chain, nuclear, data, nuclide, xml, equilibrium, species
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2010/05/18"/>
//     </current>
//   </authors>
//   
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//   
// </class>
//############################################################################*/

struct _Libnuceq__Species {
  double dNseFactor;
  double dNseCorrectionFactor;
  double dAbundance;
  Libnucnet__Species *pNucSpecies;
};

/*##############################################################################
// <class name="Libnuceq">
//
//   <description>
//     <abstract>
//       Libnuceq is a structure that stores information for an equilibrium
//       calculation.
//     </abstract>
//     <keywords>
//       isotope chain, nuclear, data, nuclide, xml, equilibrium
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2007/12/18"/>
//     </current>
//   </authors>
//   
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//   
// </class>
//############################################################################*/

struct _Libnuceq {
  Libnucnet__Nuc *pNuc;
  xmlHashTablePtr pSpeciesHash;
  xmlHashTablePtr pClusterHash;
  Libnucnet__Species__nseCorrectionFactorFunction pfCorr;
  void *pCorrData;
  Libstatmech__Fermion__Function pfElectronFunction;
  Libstatmech__Fermion__Integrand pfElectronIntegrand;
  Libnuceq__wseCorrectionFunction pfWseCorrectionFunction;
  void *pElectronFunctionData;
  void *pElectronIntegrandData;
  void *pWseCorrectionFunctionData;
  double dT9;
  double dRho;
  double *pYe;
  double dMunkT;
  double dMupkT;
  double *pMuekT;
};

/*##############################################################################
// <class name="Libnuceq__Cluster">
//
//   <description>
//     <abstract>
//       Libnuceq__Cluster is a structure that stores information about an
//       equilibrium cluster, that is, a subset of nuclei in an equilibrium
//       that are in equilibrium with each other.
//     </abstract>
//     <keywords>
//       isotope chain, nuclear, data, nuclide, xml, equilibrium, cluster
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2010/05/18"/>
//     </current>
//   </authors>
//   
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//   
// </class>
//############################################################################*/

struct _Libnuceq__Cluster {
  Libnuceq *pEquil;
  xmlChar *sxXPath;
  double dConstraint, dMukT;
  xmlHashTablePtr pSpeciesHash;
  Libnuceq__Cluster__constraint_function pfConstraintFunction;
  void *pConstraintData;
  Libnuceq__Cluster__prefactorFunction pfPrefactorFunction;
  void *pPrefactorData;
};

/*##############################################################################
// More functions.
//############################################################################*/

/*##############################################################################
// <user_routine name="Libnuceq__Species__iterateFunction()">
//
//   <description>
//     <abstract>
//       User-supplied routine to be applied during an iteration over
//       an equilibrium.
//     </abstract>
//     <keywords>
//       Libnuceq, user, supplied, function, species, iterate
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int
// Libnucnet__Species__iterateFunction(
//   Libnuceq__Species *self,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a libnuceq species.
//     </param>
//
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined data structure containing extra
//       data for the function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return 1 to continue or 0 to stop.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef
int ( *Libnuceq__Species__iterateFunction )( Libnuceq__Species *, void * );

/*##############################################################################
// <user_routine name="Libnuceq__Cluster__iterateFunction()">
//
//   <description>
//     <abstract>
//       User-supplied routine to be applied during an iteration over
//       the clusters in an equilibrium.
//     </abstract>
//     <keywords>
//       Libnuceq, user, supplied, function, cluster, iterate
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int
// Libnucnet__Cluster__iterateFunction(
//   Libnuceq__Cluster *self,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined data structure containing extra
//       data for the function.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return 1 to continue or 0 to stop.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef
void ( *Libnuceq__Cluster__iterateFunction )( Libnuceq__Cluster *, void * );

typedef double ( *Libnuceq__solve_function )( double, void * );

/*##############################################################################
// API Routines.
//############################################################################*/

/*##############################################################################
// <routine name="Libnuceq__new()">
//
//   <description>
//     <abstract>
//       Create a new Libnuceq structure.
//     </abstract>
//     <keywords>
//       Libnucnet, nuclear, nse, statistical, equilibrium, data, xml, hash
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libnuceq *
// Libnuceq__new( Libnucnet__Nuc *p_nuc );
//     </calling_sequence>
//
//     <param
//       name="p_nuc"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnucnet__Nuc structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns a pointer to a new Libnuceq structure.
//       If it is not possible to allocate memory for the new structure,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Create a Libnuceq structure p_my_eq from the
//         Libnucnet__Nuc structure p_my_nuclei:
//       </synopsis>
//
//       <code>
// p_my_eq =
//   Libnuceq__new(
//     p_my_nuclei
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libnuceq *
Libnuceq__new(
  Libnucnet__Nuc *
);

/*##############################################################################
// <routine name="Libnuceq__free()">
//
//   <description>
//     <abstract>
//       Free a Libnuceq structure.
//     </abstract>
//     <keywords>
//       Libnucnet, nuclear, nse, statistical, equilibrium, data, xml, hash
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void Libnucnet__free( Libnuceq *self );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the Libnuceq structure memory has been freed.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Free the memory for Libnuceq *p_my_eq:
//       </synopsis>
//
//       <code>
//   Libnuceq__free( p_my_eq );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void Libnuceq__free( Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__getNuc()">
//
//   <description>
//     <abstract>
//       Return the Libnucnet__Nuc nuclear collection underlying the
//       equilibrium structure.
//     </abstract>
//     <keywords>
//       Libnucnet, nuclear, nse, statistical, equilibrium, data, xml
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libnucnet__Nuc *
// Libnuceq__getNuc( Libnuceq *self );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns a pointer to the underlying nuclear collection.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the pointer to the Libnucnet__Nuc structure
//         underlying Libnuceq *p_my_equil and call it p_my_nuclei:
//       </synopsis>
//
//       <code>
// p_my_nuclei =
//   Libnuceq__getNuc(
//     p_my_equil
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libnucnet__Nuc *
Libnuceq__getNuc( Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__computeEquilibrium()">
//
//   <description>
//     <abstract>
//       Compute the equilibrium for the input temperature and density.
//     </abstract>
//     <keywords>
//       Libnucnet, nuclear, nse, statistical, equilibrium, compute
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__computeEquilibrium(
//   Libnuceq *self
//   double d_t9,
//   double d_rho
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="d_t9"
//       kind="in,positional,required"
//     >
//       A double giving the temperature in billions of Kelvins at which
//       to compute the equilibrium.
//     </param>
//     <param
//       name="d_rho"
//       kind="in,positional,required"
//     >
//       A double giving the density in g/cc at which
//       to compute the equilibrium.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the equilibrium has been computed.  Equilibrium
//       abundances and chemical potentials may be retrieved with other
//       API routines.  If the input is invalid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Compute the equilibrium for the input structure p_my_equil
//         at a temperature of 5 billion K and 1.e7 g/cc:
//       </synopsis>
//
//       <code>
//   Libnuceq__computeEquilibrium(
//     p_my_equil,
//     5.,
//     1.e7
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__computeEquilibrium(
  Libnuceq *,
  double,
  double
);

/*##############################################################################
// <routine name="Libnuceq__getSpeciesByName()">
//
//   <description>
//     <abstract>
//       Retrieve a species by its name from an equilibrium structure.
//     </abstract>
//     <keywords>
//       species, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libnuceq__Species *
// Libnuceq__getSpeciesByName(
//   Libnuceq *self,
//   const char *s_name
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="s_name"
//       kind="in,positional,required"
//     >
//       A string giving the name of the species to be retrieved.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to the species.  If the species
//       does not exist in the equilibrium, routine returns NULL.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the pointer to the Libnuceq__Species *p_eq_species
//         for ni56 from equilibrium p_my_equil:
//       </synopsis>
//
//       <code>
// p_eq_species =
//   Libnuceq__getSpeciesByName(
//     p_my_equil,
//     "ni56"
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libnuceq__Species *
Libnuceq__getSpeciesByName( Libnuceq *, const char * );

/*##############################################################################
// <routine name="Libnuceq__Species__getAbundance()">
//
//   <description>
//     <abstract>
//       Retrieve the equilibrium abundance of a species.
//     </abstract>
//     <keywords>
//       species, nuclear, nse, statistical, equilibrium, abundance
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__Species__getAbundance(
//   Libnuceq__Species *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Species structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a double giving the abundance of the species.
//       If the input species is not valid, Libnuceq error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the abundance of ni56 in the equilibrium p_my_equil:
//       </synopsis>
//
//       <code>
// fprintf(
//   stdout,
//   "Abundance of ni56 is %g.\n",
//   Libnuceq__Species__getAbundance(
//     Libnuceq__getSpeciesByName(
//       p_my_equil,
//       "ni56"
//     )
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__Species__getAbundance( Libnuceq__Species * );

/*##############################################################################
// <routine name="Libnuceq__iterateSpecies()">
//
//   <description>
//     <abstract>
//       Iterate over the species in an equilibrium and apply the
//       user-supplied function.
//     </abstract>
//     <keywords>
//       species, nuclear, nse, statistical, equilibrium, iterate
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__iterateSpecies(
//   const Libnuceq *self,
//   (Libnuceq__Species__iterateFunction) pf_my_function,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="p_my_function"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined function to be applied to the species.
//     </param>
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined data structure for extra data
//       to be applied.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine iterates through the species in the equilibrium
//       and applies the user-supplied routine to each species.
//       If any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Iterate through the species in p_my_equilibrium and apply
//         the function my_iterate_function and the extra data p_user_data:
//       </synopsis>
//
//       <code>
// Libnuceq__iterateSpecies(
//   p_my_equil,
//   (Libnuceq__Species__iterateFunction) my_iterate_function,
//   p_user_data
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__iterateSpecies(
  const Libnuceq *,
  Libnuceq__Species__iterateFunction,
  void *
);

/*##############################################################################
// <routine name="Libnuceq__Species__getNucSpecies()">
//
//   <description>
//     <abstract>
//       Return the Libnucnet__Nuc species corresponding to the given
//       species in the equilibrium structure.
//     </abstract>
//     <keywords>
//       Libnucnet, nuclear, nse, statistical, equilibrium, data, species
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libnucnet__Species *
// Libnuceq__Species__getNucSpecies( Libnuceq__Species *self );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Species structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns a pointer to the underlying Libnucnet__Species
//       structure.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the mass excess of Libnuceq__Species *p_eq_species:
//       </synopsis>
//
//       <code>
// fprintf(
//   stdout,
//   "The mass excess of %s is %g (MeV).\n",
//   Libnucnet__Species__getName(
//     Libnuceq__Species__getNucSpecies( p_eq_species )
//   ),
//   Libnucneq__Species__getMassExcess(
//     Libnuceq__Species__getNucSpecies( p_eq_species )
//   )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libnucnet__Species *
Libnuceq__Species__getNucSpecies( Libnuceq__Species * );

/*##############################################################################
// <routine name="Libnuceq__Cluster__iterateSpecies()">
//
//   <description>
//     <abstract>
//       Iterate over the species in an equilibrium cluster and apply the
//       user-supplied function.
//     </abstract>
//     <keywords>
//       species, nuclear, nse, statistical, equilibrium, iterate, cluster
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__Cluster__iterateSpecies(
//   const Libnuceq__Cluster *self,
//   (Libnuceq__Species__iterateFunction) pf_my_function,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       name="p_my_function"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined function to be applied to the species.
//     </param>
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to a user-defined data structure for extra data
//       to be applied.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine iterates through the species in the equilibrium cluster
//       and applies the user-supplied routine to each species.
//       If any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Iterate through the species in p_my_cluster and apply
//         the function my_iterate_function and the extra data p_user_data:
//       </synopsis>
//
//       <code>
// Libnuceq__Cluster__iterateSpecies(
//   p_my_cluster,
//   (Libnuceq__Species__iterateFunction) my_iterate_function,
//   p_user_data
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__Cluster__iterateSpecies(
  const Libnuceq__Cluster *,
  Libnuceq__Species__iterateFunction,
  void *
);

/*##############################################################################
// <routine name="Libnuceq__getNumberOfClusters()">
//
//   <description>
//     <abstract>
//       Retrieve the current number of clusters in the equilibrium.
//     </abstract>
//     <keywords>
//       species, equilibrium, cluster, number
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t
// Libnuceq__getNumberOfClusters(
//   const Libnuceq__Cluster *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns the number of clusters currently set for
//       the equilibrium.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the number of clusters in the equilibrum p_equil:
//       </synopsis>
//
//       <code>
// printf(
//   "The equilibrium has %lu clusters.\n",
//   (unsigned long) Libnuceq__getNumberOfClusters( p_equil )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t
Libnuceq__getNumberOfClusters(
  const Libnuceq *
);

/*##############################################################################
// <routine name="Libnuceq__setYe()">
//
//   <description>
//     <abstract>
//       Set the Ye constraint for an equilibrium.
//     </abstract>
//     <keywords>
//       Libnucnet, statistical, equilibrium, electron, nucleon, ratio, ye
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__setYe( Libnuceq *self, double d_ye );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="d_ye"
//       kind="in,positional,required"
//     >
//       A double giving the electron-to-nucleon ratio Ye to be set for
//       the equilibrium.  This number must be between zero and one.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the Ye constraint for the equilibrium
//       has been set.  If a previous Ye constraint existed, it has been
//       replaced with the new value.
//       If the input structure is not valid, or the Ye is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Set the Ye constraint on p_my_equil to 0.43:
//       </synopsis>
//
//       <code>
// Libnucnet__setYe(
//   p_my_equil, 0.43
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__setYe(
  Libnuceq *,
  double
);

/*##############################################################################
// <routine name="Libnuceq__clearYe()">
//
//   <description>
//     <abstract>
//       Clear the Ye constraint for an equilibrium.
//     </abstract>
//     <keywords>
//       Libnucnet, statistical, equilibrium, electron, nucleon, ratio, ye
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__clearYe( Libnuceq *self );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the Ye constraint for the equilibrium
//       has been removed.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Clear the Ye constraint on p_my_equil:
//       </synopsis>
//
//       <code>
// Libnucnet__clearYe(
//   p_my_equil
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__clearYe( Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__newCluster()">
//
//   <description>
//     <abstract>
//       Create a new Libnuceq cluster.
//     </abstract>
//     <keywords>
//       Libnucnet, nuclear, nse, statistical, equilibrium, data, cluster
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libnuceq__Cluster *
// Libnuceq__newCluster(
//   Libnuceq *self,
//   const char *s_xpath
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="s_xpath"
//       kind="in,positional,required"
//     >
//       A string giving the XPath expression defining the cluster.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns a pointer to a new Libnuceq__Cluster structure.
//       If it is not possible to allocate memory for the new structure,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Create a new equilibrium cluster within the equilibrium p_my_equil
//         for all nuclei with charge
//         greater than or equal to 6:
//       </synopsis>
//
//       <code>
// p_my_cluster =
//   Libnuceq__newCluster(
//     p_my_equil,
//     "[z >= 6]"
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libnuceq__Cluster *
Libnuceq__newCluster(
  Libnuceq *,
  const char *
);

/*##############################################################################
// <routine name="Libnuceq__removeCluster()">
//
//   <description>
//     <abstract>
//       Remove a cluster from an equilibrium.
//     </abstract>
//     <keywords>
//       Libnucnet, nuclear, nse, statistical, equilibrium, remove, cluster
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int
// Libnuceq__removeCluster(
//   Libnuceq *self,
//   Libnuceq__Cluster *p_cluster
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="p_cluster"
//       kind="in,positional,required"
//     >
//       A pointer to the Libnuceq__Cluster to be removed.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns 1 (true) if the removal succeeded and 0 (false)
//       if not.  If the input is invalid, Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Remove the cluster containing all nuclear with atomic number
//         Z = 12 from the equilibrium p_my_equil:
//       </synopsis>
//
//       <code>
// p_cluster =
//   Libnuceq__getCluster(
//     p_my_equil,
//     "[z = 12]"
//   );
// if( !Libnuceq__removeCluster( p_my_equil, p_cluster ) )
//   fprintf(
//     stderr,
//     "Couldn't remove cluster!\n"
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
Libnuceq__removeCluster(
  Libnuceq *,
  Libnuceq__Cluster *
);

/*##############################################################################
// <routine name="Libnuceq__Cluster__free()">
//
//   <description>
//     <abstract>
//       Free a Libnuceq cluster.
//     </abstract>
//     <keywords>
//       Libnucnet, nuclear, nse, statistical, equilibrium, data, cluster
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__Cluster__free(
//   Libnuceq__Cluster *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the cluster has been freed.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Free the cluster Libnuceq__Cluster *p_my_cluster:
//       </synopsis>
//
//       <code>
// Libnuceq__Cluster__free(
//     p_my_cluster
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__Cluster__free( Libnuceq__Cluster * );

/*##############################################################################
// <routine name="Libnuceq__getMunkT()">
//
//   <description>
//     <abstract>
//       Retrieve the neutron chemical potential from an equilibrium.
//     </abstract>
//     <keywords>
//       Libnucnet, neutron, statistical, equilibrium, chemical, potential
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__getMunkT(
//   const Libnuceq *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns the neutron chemical potential divided by kT for
//       the equilibrium.
//       If the input equilibrium is invalid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the neutron chemical potential (divided by kT) for
//         p_my_equil:
//       </synopsis>
//
//       <code>
//  fprintf(
//    stdout,
//    "The mu_n / kT for the equilibrium is %g.\n",
//    Libnuceq__getMunkT( p_my_equil )
//  );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__getMunkT( const Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__getMupkT()">
//
//   <description>
//     <abstract>
//       Retrieve the proton chemical potential from an equilibrium.
//     </abstract>
//     <keywords>
//       Libnucnet, proton, statistical, equilibrium, chemical, potential
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__getMupkT(
//   const Libnuceq *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns the proton chemical potential divided by kT for
//       the equilibrium.
//       If the input equilibrium is invalid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the proton chemical potential (divided by kT) for
//         p_my_equil:
//       </synopsis>
//
//       <code>
//  fprintf(
//    stdout,
//    "The mu_p / kT for the equilibrium is %g.\n",
//    Libnuceq__getMupkT( p_my_equil )
//  );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__getMupkT( const Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__getMuekT()">
//
//   <description>
//     <abstract>
//       Retrieve the electron chemical potential from an equilibrium.
//     </abstract>
//     <keywords>
//       Libnucnet, electron, statistical, equilibrium, chemical, potential
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__getMuekT(
//   const Libnuceq *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns the electron chemical potential divided by kT for
//       the equilibrium.
//       If the input equilibrium is invalid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the electron chemical potential (divided by kT) for
//         p_my_equil:
//       </synopsis>
//
//       <code>
//  fprintf(
//    stdout,
//    "The mu_e / kT for the equilibrium is %g.\n",
//    Libnuceq__getMuekT( p_my_equil )
//  );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__getMuekT( const Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__Cluster__getMukT()">
//
//   <description>
//     <abstract>
//       Retrieve the chemical potential of an equilibrium cluster.
//     </abstract>
//     <keywords>
//       Libnucnet, cluster, statistical, equilibrium, chemical, potential
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__Cluster__getMukT(
//   const Libnuceq__Cluster *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns the chemical potential divided by kT for
//       the equilibrium cluster.
//       If the input equilibrium cluster is invalid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the chemical potential (divided by kT) for
//         the cluster p_my_cluster:
//       </synopsis>
//
//       <code>
//  fprintf(
//    stdout,
//    "The mu / kT for the equilibrium cluster is %g.\n",
//    Libnuceq__Cluster__getMukT( p_my_equil )
//  );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__Cluster__getMukT( const Libnuceq__Cluster * );

/*##############################################################################
// <routine name="Libnuceq__getCluster()">
//
//   <description>
//     <abstract>
//       Retrieve a cluster from an equilibrium structure.
//     </abstract>
//     <keywords>
//       cluster, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// Libnucnet__Cluster *
// Libnuceq__getCluster(
//   Libnuceq *self,
//   const char *s_xpath
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="s_xpath"
//       kind="in,positional,required"
//     >
//       A string giving the XPath expression used to define the cluster.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns a pointer to the cluster.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//       If the cluster is not found, routine returns NULL.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the pointer to the Libnuceq__Cluster *p_my_cluster
//         from Libnuceq *p_my_equil that was defined by the XPath expression
//         "[z >= 6]";
//       </synopsis>
//
//       <code>
// p_my_cluster
//   Libnuceq__getCluster(
//     p_my_equil,
//     "[z >= 6]"
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

Libnuceq__Cluster *
Libnuceq__getCluster( const Libnuceq *, const char * );

/*##############################################################################
// <routine name="Libnuceq__getT9()">
//
//   <description>
//     <abstract>
//       Retrieve the temperature (in billions of K) at which an equilibrium
//       was computed.
//     </abstract>
//     <keywords>
//       temperature, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__getT9(
//   const Libnuceq *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns a double giving the temperature at which the
//       equilibrium was computed.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//       If the equilibrium has not yet been computed, the routine returns
//       zero.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the temperature at which the equilibrium p_my_equil
//         was computed:
//       </synopsis>
//
//       <code>
// d_t9 = Libnuceq__getT9( p_my_equil );
// if( d_t9 )
//   printf(
//     "Equilibrium computed at T9 = %g\n",
//     d_t9
//   );
// else
//   printf( "Equilibrium not yet computed.\n" );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__getT9( const Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__getRho()">
//
//   <description>
//     <abstract>
//       Retrieve the density (in g/cc) at which an equilibrium
//       was computed.
//     </abstract>
//     <keywords>
//       density, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__getRho(
//   const Libnuceq *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns a double giving the density at which the
//       equilibrium was computed.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//       If the equilibrium has not yet been computed, the routine returns
//       zero.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the density at which the equilibrium p_my_equil
//         was computed:
//       </synopsis>
//
//       <code>
// d_rho = Libnuceq__getRho( p_my_equil );
// if( d_rho )
//   printf(
//     "Equilibrium computed at rho = %g g/cc\n",
//     d_rho
//   );
// else
//   printf( "Equilibrium not yet computed.\n" );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__getRho( const Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__Cluster__getXPathString()">
//
//   <description>
//     <abstract>
//       Retrieve the XPath string defining a cluster.
//     </abstract>
//     <keywords>
//       XPath, cluster, string, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// const char *
// Libnuceq__Cluster__getXPathString(
//   const Libnuceq__Cluster *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns the string defining the input cluster.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the XPath string for p_my_cluster:
//       </synopsis>
//
//       <code>
// printf(
//   "Cluster XPath is %s\n",
//   Libnuceq__Cluster__getXPathString( p_my_cluster )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

const char *
Libnuceq__Cluster__getXPathString( const Libnuceq__Cluster * );

/*##############################################################################
// <routine name="Libnuceq__Cluster__getConstraint()">
//
//   <description>
//     <abstract>
//       Retrieve the abundance constraint for a cluster.
//     </abstract>
//     <keywords>
//       XPath, cluster, constraint, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__Cluster__getConstraint(
//   const Libnuceq__Cluster *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns a double giving the abundance constraint on a
//       Cluster. If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the constraint for p_my_cluster:
//       </synopsis>
//
//       <code>
// printf(
//   "Cluster constraint is %g\n",
//   Libnuceq__Cluster__getConstraint( p_my_cluster )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__Cluster__getConstraint(
  const Libnuceq__Cluster *
);

/*##############################################################################
// <routine name="Libnuceq__Cluster__updateConstraint()">
//
//   <description>
//     <abstract>
//       Update the abundance constraint for a cluster.
//     </abstract>
//     <keywords>
//       cluster, constraint, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__Cluster__updateConstraint(
//   Libnuceq__Cluster *self,
//   double d_constraint
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       name="d_constraint"
//       kind="in,positional,required"
//     >
//       A double giving the new abundance constraint for the cluster.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the abundance constraint on a cluster has
//       been updated.  If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Update the constraint for p_my_cluster to 0.01:
//       </synopsis>
//
//       <code>
// Libnuceq__Cluster__updateConstraint(
//   p_my_cluster,
//   0.01
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__Cluster__updateConstraint(
  Libnuceq__Cluster *,
  double
);

/*##############################################################################
// <routine name="Libnuceq__Cluster__updateConstraintFunction()">
//
//   <description>
//     <abstract>
//       Update the constraint function for a cluster and its associated
//       data.
//     </abstract>
//     <keywords>
//       cluster, constraint, function, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__Cluster__updateConstraintFunction(
//   Libnuceq__Cluster *self,
//   Libnuceq__Cluster__constraint_function pf_func,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       name="pf_func"
//       kind="in,positional,required"
//     >
//       The name of the constraint function.
//     </param>
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to a data structure giving extra data for the constraint
//       function.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the cluster's constraint function and
//       associated data have been updated.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Update Libnucnet__Cluster *p_my_cluster with the constraint
//         function pf_my_func and associated data *p_data:
//       </synopsis>
//
//       <code>
// Libnuceq__Cluster__updateConstraintFunction(
//   p_my_cluster,
//   (Libnuceq__Cluster__constraint_function) pf_my_func,
//   p_data
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__Cluster__updateConstraintFunction(
  Libnuceq__Cluster *,
  Libnuceq__Cluster__constraint_function,
  void *
);

/*##############################################################################
// <routine name="Libnuceq__Cluster__updatePrefactorFunction()">
//
//   <description>
//     <abstract>
//       Update the prefactor function for a cluster and its associated
//       data.
//     </abstract>
//     <keywords>
//       cluster, prefactor, function, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__Cluster__updatePrefactorFunction(
//   Libnuceq__Cluster *self,
//   Libnuceq__Cluster__prefactorFunction pf_func,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       name="pf_func"
//       kind="in,positional,required"
//     >
//       The name of the prefactor function.
//     </param>
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to a data structure giving extra data for the prefactor
//       function.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the cluster's prefactor function and
//       associated data have been updated.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Update Libnucnet__Cluster *p_my_cluster with the prefactor
//         function pf_my_func and associated data *p_data:
//       </synopsis>
//
//       <code>
// Libnuceq__Cluster__updatePrefactorFunction(
//   p_my_cluster,
//   (Libnuceq__Cluster__prefactorFunction) pf_my_func,
//   p_data
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__Cluster__updatePrefactorFunction(
  Libnuceq__Cluster *,
  Libnuceq__Cluster__prefactorFunction,
  void *
);

/*##############################################################################
// <routine name="Libnuceq__getAbundances()">
//
//   <description>
//     <abstract>
//       Retrieve the abundances of the species in an equilibrium.
//     </abstract>
//     <keywords>
//       abundances, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// Libnuceq__getAbundances(
//   const Libnuceq__Cluster *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a new gsl_vector containing the abundances in the
//       equilibrium.  The species are sorted according to the sorting
//       set for the underlying Libnucnet__Nuc structure.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the abundances from the equilibrium p_my_equil and store
//         them in the new gsl_vector p_my_abundances:
//       </synopsis>
//
//       <code>
// p_my_abundances =
//   Libnuceq__getAbundances( p_my_equil );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
Libnuceq__getAbundances( const Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__iterateClusters()">
//
//   <description>
//     <abstract>
//       Iterate the clusters in an equilibrium and apply the user-defined
//       function and associated data.
//     </abstract>
//     <keywords>
//       iterate, cluster, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__iterateClusters(
//   Libnuceq *self,
//   Libnuceq__Cluster__iterateFunction pf_func,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="pf_func"
//       kind="in,positional,required"
//     >
//       The name of the function to be applied to the clusters.
//     </param>
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to the extra data associated with the cluster iterate
//       function.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine iterates over the clusters in the equilibrium and applies
//       the user-supplied function and data.
//       If any input is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Iterate over the clusters in the equilibrium p_my_equil and
//         apply the function p_my_cluster_function and the data structure
//         p_my_data.
//       </synopsis>
//
//       <code>
// Libnuceq__iterateClusters(
//   p_my_equil,
//   (Libnuceq__Cluster__iterateFunction) pf_my_cluster_function,
//   p_data
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__iterateClusters(
  Libnuceq *,
  Libnuceq__Cluster__iterateFunction,
  void *
);

/*##############################################################################
// <routine name="Libnuceq__setNseCorrectionFactorFunction()">
//
//   <description>
//     <abstract>
//       Set the NSE correction factor function and the associated data
//       for the equilibrium.
//     </abstract>
//     <keywords>
//       nuclear, nse, correction, factor, equilibrium, function
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__setNseCorrectionFactorFunction(
//   Libnuceq *self,
//   Libnucnet__Species__nseCorrectionFactorFunction pf_func,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="pf_func"
//       kind="in,positional,required"
//     >
//       The name of the correction factor function to be applied to the
//       species.
//     </param>
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to the extra data associated with the correction factor
//       function.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the correction factor function and the
//       associated data have been set to the input values.
//       If any input is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         For equilibrium p_my_equil, set the correction factor function 
//         to my_correction_function and set the associated data to point to the
//         data structure my_data:
//       </synopsis>
//
//       <code>
// Libnuceq__setNseCorrectionFactorFunction(
//   p_my_equil,
//   (Libnucnet__Species__nseCorrectionFactorFunction) my_correction_function,
//   p_data
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__setNseCorrectionFactorFunction(
  Libnuceq *,
  Libnucnet__Species__nseCorrectionFactorFunction,
  void *
);

/*##############################################################################
// <routine name="Libnuceq__clearNseCorrectionFactorFunction()">
//
//   <description>
//     <abstract>
//       Clear the NSE correction factor function and the associated data
//       for the equilibrium.
//     </abstract>
//     <keywords>
//       nuclear, nse, correction, factor, equilibrium, function
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__clearNseCorrectionFactorFunction(
//   Libnuceq *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the correction factor function has been
//       reset to the default (no correction).
//       If the input equilibrium is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         For equilibrium p_my_equil, clear the correction factor function: 
//       </synopsis>
//
//       <code>
// Libnuceq__clearNseCorrectionFactorFunction(
//   p_my_equil
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__clearNseCorrectionFactorFunction(
  Libnuceq *
);

/*##############################################################################
// <routine name="Libnuceq__computeZMoment()">
//
//   <description>
//     <abstract>
//       Compute the moment of the abundances of the equilibrium about the
//       charge of each species raised to a power.
//     </abstract>
//     <keywords>
//       nuclear, nse, equilibrium, Z, moment
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__computeZMoment(
//   const Libnuceq *self,
//   unsigned int i
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="i"
//       kind="in,positional,required"
//     >
//       An unsigned int giving the power to which to raise the atomic
//       number of each species in computing the moment.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the sum of the abundance of each species times the
//       species charge raised to the power i.
//       If the input equilibrium is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print Ye for the equilibrium p_my_equil:
//       </synopsis>
//
//       <code>
// fprintf(
//   stdout,
//   "Ye = %e\n",
//   Libnuceq__computeZMoment( p_my_equil, 1 )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__computeZMoment(
  const Libnuceq *,
  unsigned int
);

/*##############################################################################
// <routine name="Libnuceq__computeAMoment()">
//
//   <description>
//     <abstract>
//       Compute the moment of the abundances of the equilibrium about the
//       mass number of each species raised to a power.
//     </abstract>
//     <keywords>
//       nuclear, nse, equilibrium, A, moment
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double
// Libnuceq__computeAMoment(
//   const Libnuceq *self,
//   unsigned int i
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="i"
//       kind="in,positional,required"
//     >
//       An unsigned int giving the power to which to raise the mass
//       number of each species in computing the moment.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the sum of the abundance of each species times the
//       species mass number raised to the power i.
//       If the input equilibrium is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Print the sum of abundances for the equilibrium p_my_equil:
//       </synopsis>
//
//       <code>
// fprintf(
//   stdout,
//   "Sum of abundances = %e\n",
//   Libnuceq__computeAMoment( p_my_equil, 0 )
// );
//       </code>
//     </doc>
//
//     <doc kind="example" id="example2">
//       <synopsis>
//         Print the sum of mass fractions for the equilibrium p_my_equil:
//       </synopsis>
//
//       <code>
// fprintf(
//   stdout,
//   "Sum of mass fractions = %e\n",
//   Libnuceq__computeAMoment( p_my_equil, 1 )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double
Libnuceq__computeAMoment(
  const Libnuceq *, unsigned int
);

/*##############################################################################
// <routine name="Libnuceq__getNumberOfSpecies()">
//
//   <description>
//     <abstract>
//       Retrieve the number of species in an equilibrium.
//     </abstract>
//     <keywords>
//       number, species, nuclear, nse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t
// Libnuceq__getNumberOfSpecies(
//   const Libnuceq *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the number of species in the equilibrium.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the number of species in equilibrium p_my_equil:
//       </synopsis>
//
//       <code>
// printf(
//   "The number of species in the equilibrium is %lu:\n"
//   (unsigned long) Libnuceq__getNumberOfSpecies( p_my_equil )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t
Libnuceq__getNumberOfSpecies( const Libnuceq * );

/*##############################################################################
// <routine name="Libnuceq__Cluster__getNumberOfSpecies()">
//
//   <description>
//     <abstract>
//       Retrieve the number of species in an equilibrium cluster.
//     </abstract>
//     <keywords>
//       number, species, nuclear, nse, statistical, equilibrium, cluster
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t
// Libnuceq__Cluster__getNumberOfSpecies(
//   const Libnuceq__Cluster *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq__Cluster structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the number of species in the equilibrium cluster.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the number of species in the cluster defined by nuclei with
//         Z >= 6 in equilibrium p_my_equil:
//       </synopsis>
//
//       <code>
// p_cluster = Libnuceq__getCluster( p_my_equil, "[z >= 6]" );
// printf(
//   "The number of species in the equilibrium cluster is %lu:\n"
//   (unsigned long) Libnuceq__Cluster__getNumberOfSpecies( p_cluster )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t
Libnuceq__Cluster__getNumberOfSpecies( const Libnuceq__Cluster * );

/*##############################################################################
// <routine name="Libnuceq__updateUserElectronNumberDensity()">
//
//   <description>
//     <abstract>
//       Set the function and integrand and associated data
//       for computing the electron
//       number density from its corresponding chemical potential in an
//       equilibrium to the user-supplied values.
//     </abstract>
//     <keywords>
//       number, species, nuclear, nse, statistical, equilibrium, electron
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__updateUserElectronNumberDensity(
//   Libnuceq *self,
//   Libstatmech__Fermion__Function pf_func,
//   Libstatmech__Fermion__Integrand pf_integrand,
//   void *p_function_data,
//   void *p_integrand_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="pf_func"
//       kind="in,positional,required"
//     >
//       The name of the user-supplied Libstatmech__Fermion__Function
//       to compute the number density
//       from the temperature and fermion chemical potential.
//     </param>
//     <param
//       name="pf_integrand"
//       kind="in,positional,required"
//     >
//       The name of the user-supplied Libstatmech__Fermion__Integrand
//       to compute the number density
//       from the temperature and fermion chemical potential.
//     </param>
//     <param
//       name="p_function_data"
//       kind="in,positional,required"
//     >
//       A pointer to the data structure for extra data to the function.
//     </param>
//     <param
//       name="p_integrand_data"
//       kind="in,positional,required"
//     >
//       A pointer to the data structure for extra data to the integrand.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return the number density function and integrand for
//       the electrons has been set to the input values along with
//       the extra data.
//       If both the function and integrand are NULL,
//       the default (non-interacting relativistic electrons) is used.
//       If the either function or integrand is not NULL,
//       it is used in place of the default.  If the input structure is
//       not valid, Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Set the number density function for electrons to the user-supplied
//         function my_function that uses extra data in my_data and do not
//         use an integrand:
//       </synopsis>
//
//       <code>
// Libnuceq__updateUserElectronNumberDensityFunction(
//   p_my_equil,
//   (Libstatmech__Fermion__function) my_function,
//   NULL,
//   &my_data,
//   NULL
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__updateUserElectronNumberDensity(
  Libnuceq *,
  Libstatmech__Fermion__Function,
  Libstatmech__Fermion__Integrand,
  void *,
  void *
);

/*##############################################################################
// <routine name="Libnuceq__updateWseCorrectionFunction()">
//
//   <description>
//     <abstract>
//       Update the WSE relation correction function and its associated
//       data.
//     </abstract>
//     <keywords>
//       cluster, function, nuclear, wse, statistical, equilibrium
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__updateWseCorrectionFunction(
//   Libnuceq *self,
//   Libnuceq__wseCorrectionFunction pf_func,
//   void *p_data
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//     >
//       A pointer to a Libnuceq structure.
//     </param>
//     <param
//       name="pf_func"
//       kind="in,positional,required"
//     >
//       The name of the correction function.
//     </param>
//     <param
//       name="p_data"
//       kind="in,positional,required"
//     >
//       A pointer to a data structure giving extra data for the correction
//       function.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the wse correction function and
//       associated data have been updated.
//       If the input structure is not valid,
//       Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Update Libnucnet *p_my_equilibrium with the wse correction
//         function pf_my_func and associated data *p_data:
//       </synopsis>
//
//       <code>
// Libnuceq__updateWseCorrectionFunction(
//   p_my_equilibrium,
//   (Libnuceq__wseCorrectionFunction) pf_my_func,
//   p_data
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__updateWseCorrectionFunction(
  Libnuceq *,
  Libnuceq__wseCorrectionFunction,
  void *
);

/*##############################################################################
// <routine name="Libnuceq__copy_clusters()">
//
//   <description>
//     <abstract>
//       Copy the clusters from one equilibrium to another.
//     </abstract>
//     <keywords>
//       number, species, nuclear, nse, statistical, equilibrium, cluster
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// Libnuceq__copy_clusters(
//   Libnuceq * destination, Libnuceq * source
// );
//     </calling_sequence>
//
//     <param
//       name="destination"
//       kind="in,positional,required"
//     >
//       A pointer to the destination Libnuceq equilibrim.
//     </param>
//     <param
//       name="source"
//       kind="in,positional,required"
//     >
//       A pointer to the source Libnuceq equilibrim.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the clusters in the source equilibrium
//       have been copied to the destination.  Before copying, the routine
//       first clears any clusters in the destination equilibrium.
//       If the source or destination structure is not valid,
//       or if the source and destination structures do not share the same
//       parent nuclear collection, Libnuceq error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Copy the equilibrium clusters from equilibrium p_source to
//         p_destination.
//       </synopsis>
//
//       <code>
// Libnuceq__copy_clusters( p_destination, p_source );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
Libnuceq__copy_clusters(
  Libnuceq *,
  Libnuceq *
);

/*##############################################################################
// Non-API Routines.
//############################################################################*/

Libnuceq__Species *
Libnuceq__Species__new(
  Libnucnet__Species *
);

void
Libnuceq__Species__free(
  Libnuceq__Species *,
  xmlChar *
);

void
Libnuceq__setSpeciesNseFactors( Libnuceq * );

void
Libnuceq__setSpeciesNseFactorsCallback(
  Libnuceq__Species *,
  Libnuceq *,
  xmlChar *
);

void
Libnuceq__setSpeciesNseCorrectionFactors( Libnuceq * );

void
Libnuceq__setSpeciesNseCorrectionFactorsCallback(
  Libnuceq__Species *,
  Libnuceq *,
  xmlChar *
);

void
Libnuceq__Species__zeroNseCorrectionFactor(
  Libnuceq__Species *,
  void *,
  xmlChar *
);

double
Libnuceq__solveEquilibrium( Libnuceq__solve_function, void * );

double
Libnuceq__A_function( 
  double,
  void *
);

double
Libnuceq__WSE_function( 
  double,
  void *
);

double
Libnuceq__NSE_function( 
  double,
  void *
);

void
Libnuceq__assign_species(
  Libnucnet__Species *,
  Libnuceq *,
  xmlChar *
);

double
Libnuceq__computeSpeciesBaseLogAbundance(
  Libnuceq *,
  Libnuceq__Species *
);

void
Libnuceq__A_function_callback(
  Libnuceq__Species *,
  void *,
  xmlChar *
);

void
Libnuceq__function2_callback(
  Libnuceq__Species *,
  void *,
  xmlChar *
);

Libnuceq__Species **
Libnuceq__createSpeciesArray( const Libnuceq * );

Libnuceq__Species **
Libnuceq__Cluster__createSpeciesArray ( const Libnuceq__Cluster * );

void
Libnuceq__create_species_array_callback(
  Libnuceq__Species *,
  void *,
  const xmlChar *
);

int
Libnuceq__sorter(
  const void *,
  const void *
);

int
Libnuceq__bracket_root_of_function(
  gsl_function, double *, double *, void *
);

void
Libnuceq__Cluster__hash_free( Libnuceq__Cluster *, xmlChar * );

void
Libnuceq__assignSpeciesToCluster(
  Libnuceq *,
  Libnuceq__Cluster *
);

void
Libnuceq__assignSpeciesToClusterCallback(
  Libnucnet__Species *,
  void *,
  xmlChar *
);

void
Libnuceq__computeAbundances( Libnuceq * );

void
Libnuceq__computeAbundancesCallback(
  Libnuceq__Species *,
  Libnuceq *,
  xmlChar *
);

void
Libnuceq__set_cluster_abundances_to_zero(
  Libnuceq__Species *,
  void *,
  xmlChar *
);

void
Libnuceq__compute_cluster_abundances(
  Libnuceq__Cluster *,
  void *,
  xmlChar *
);

double
Libnuceq__cluster_abundance_function(
  double,
  Libnuceq__Cluster *
);

void
Libnuceq__cluster_abundance_function_callback(
  Libnuceq__Species *,
  void *,
  xmlChar *
);

void
Libnuceq__computeZMomentCallback(
  Libnuceq__Species *,
  void *,
  xmlChar *
);

void
Libnuceq__computeAMomentCallback(
  Libnuceq__Species *,
  void *,
  xmlChar *
);

void
Libnuceq__getAbundancesCallback(
  Libnuceq__Species *,
  gsl_vector *,
  xmlChar *
);

void
Libnuceq__iterateClustersHelper(
  Libnuceq__Cluster *,
  void *,
  xmlChar *
);

double
Libnuceq__computeElectronNumberDensity(
  Libnuceq *
);

void
Libnuceq__initializeMuekT(
  Libnuceq *
);

int
Libnuceq__copy_species_to_new_cluster(
  Libnuceq__Species *,
  Libnuceq__Cluster *
);

int
Libnuceq__copy_cluster(
  Libnuceq__Cluster *,
  Libnuceq *
);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBNUCEQ_H */

