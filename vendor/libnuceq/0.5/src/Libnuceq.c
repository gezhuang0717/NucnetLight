/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//      Copyright (c) 2010-2012 Clemson University.
//
//      This file is part of the Clemson Webnucleo group's
//      libnuceq module, originally developed by Bradley S. Meyer
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
//        Source code for the Libnuceq module.
//        For documentation, see the header file Libnuceq.h.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Include header.
//############################################################################*/

#include "Libnuceq.h"

/*##############################################################################
// Libnuceq__new()
//############################################################################*/

Libnuceq *Libnuceq__new(
  Libnucnet__Nuc *p_nuc
)
{

  /*---------------------------------------------------------------------------
  // Types.
  //--------------------------------------------------------------------------*/

  Libnuceq *self;

  /*---------------------------------------------------------------------------
  // Check input.
  //--------------------------------------------------------------------------*/

  if( !p_nuc  ) {
      LIBNUCEQ__ERROR( "Invalid input" );
  }
 
  /*---------------------------------------------------------------------------
  // Allocate Equil structure.
  //--------------------------------------------------------------------------*/

  self = ( Libnuceq * ) malloc( sizeof( Libnuceq ) );

  if( !self ) {
    LIBNUCEQ__ERROR( "Couldn't allocate memory for structure" );
  }

  /*---------------------------------------------------------------------------
  // Assign structure elements.
  //--------------------------------------------------------------------------*/

  self->pNuc = p_nuc;
  self->dT9 = 0.;
  self->dRho = 0.;
  self->pYe = NULL;
  self->pMuekT = NULL;

  self->pfCorr = NULL;
  self->pCorrData = NULL;

  self->pfElectronFunction = NULL;
  self->pElectronFunctionData = NULL;

  self->pfElectronIntegrand = NULL;
  self->pElectronIntegrandData = NULL;

  self->pfWseCorrectionFunction = NULL;
  self->pWseCorrectionFunctionData = NULL;

  /*---------------------------------------------------------------------------
  // Create and assign hash.
  //--------------------------------------------------------------------------*/

  self->pSpeciesHash = xmlHashCreate( 0 );

  xmlHashScan(
    self->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnuceq__assign_species,
    self
  );

  self->pClusterHash = xmlHashCreate( I_EQ_HASH_MAX );

  return self;

}

/*##############################################################################
// Libnuceq__assign_species().
//############################################################################*/

void
Libnuceq__assign_species(
  Libnucnet__Species *p_nuc_species,
  Libnuceq *p_equil,
  xmlChar *sx_name
)
{

  Libnuceq__Species *p_eq_species;

  p_eq_species = Libnuceq__Species__new( p_nuc_species );

  xmlHashAddEntry(
    p_equil->pSpeciesHash,
    sx_name,
    p_eq_species
  );

}

/*##############################################################################
// Libnuceq__free()
//############################################################################*/

void
Libnuceq__free( Libnuceq *self )
{

  xmlHashFree(
    self->pSpeciesHash,
    (xmlHashDeallocator) Libnuceq__Species__free
  );

  if( self->pYe ) free( self->pYe );

  if( self->pMuekT ) free( self->pMuekT );

  xmlHashFree(
    self->pClusterHash,
    (xmlHashDeallocator) Libnuceq__Cluster__hash_free
  );

  free( self );

}

/*##############################################################################
// Libnuceq__getNuc()
//############################################################################*/

Libnucnet__Nuc *
Libnuceq__getNuc( Libnuceq *self ) {

  return self->pNuc;

}

/*##############################################################################
// Libnuceq__computeEquilibrium().
//############################################################################*/
  
void
Libnuceq__computeEquilibrium(
  Libnuceq *self,
  double d_t9,
  double d_rho
)
{
  
  /*---------------------------------------------------------------------------
  // Set t9 and rho.
  //--------------------------------------------------------------------------*/

  self->dT9 = d_t9;
  self->dRho = d_rho;

  /*---------------------------------------------------------------------------
  // Initialize electron chemical potential.
  //--------------------------------------------------------------------------*/
  
  self->pMuekT = NULL; 

  /*---------------------------------------------------------------------------
  // Set NSE factors.
  //--------------------------------------------------------------------------*/

  Libnuceq__setSpeciesNseFactors( self );

  /*---------------------------------------------------------------------------
  // Solve.
  //--------------------------------------------------------------------------*/

  Libnuceq__solveEquilibrium(
    ( Libnuceq__solve_function ) Libnuceq__A_function,
    self
  );

}

/*##############################################################################
// Libnuceq__solveEquilibrium().
//############################################################################*/
  
double
Libnuceq__solveEquilibrium(
  Libnuceq__solve_function pf_func,
  void *p_params
)
{

  /*---------------------------------------------------------------------------
  // Solvers and related parameters.
  //--------------------------------------------------------------------------*/

  const gsl_root_fsolver_type *p_solver_type; 
  int i_iter = 0, i_status;
  gsl_root_fsolver *p_solver;
  gsl_function fn_equil;
  double d_x_lo = -10., d_x_hi = 10., d_result;

  /*---------------------------------------------------------------------------
  // Initialize function.
  //--------------------------------------------------------------------------*/

  fn_equil.function = pf_func;
  fn_equil.params = p_params;

  /*---------------------------------------------------------------------------
  // Bracket root.
  //--------------------------------------------------------------------------*/

  if(
      !Libnuceq__bracket_root_of_function(
         fn_equil, &d_x_lo, &d_x_hi, fn_equil.params
      )
  )
     LIBNUCEQ__ERROR( "Couldn't bracket root" );

  /*---------------------------------------------------------------------------
  // Set up solver.
  //--------------------------------------------------------------------------*/

  p_solver_type = gsl_root_fsolver_brent;
  p_solver = gsl_root_fsolver_alloc( p_solver_type );

  gsl_root_fsolver_set( p_solver, &fn_equil, d_x_lo, d_x_hi );

  /*---------------------------------------------------------------------------
  // Loop on solution.
  //--------------------------------------------------------------------------*/

  do
    {
      i_iter++;
      i_status = gsl_root_fsolver_iterate( p_solver );
      d_x_lo = gsl_root_fsolver_x_lower( p_solver );
      d_x_hi = gsl_root_fsolver_x_upper( p_solver );
      d_result = gsl_root_fsolver_root( p_solver );
      i_status =
        gsl_root_test_interval( d_x_lo, d_x_hi, D_EQ_EPS, D_EQ_EPS );

    }
    while( i_status != GSL_SUCCESS && i_iter < I_EQ_ITER_MAX );
           
  /*---------------------------------------------------------------------------
  // Clean up and return.
  //--------------------------------------------------------------------------*/

  gsl_root_fsolver_free( p_solver );

  return d_result;

}

/*##############################################################################
// Libnuceq__setSpeciesNseFactors()
//############################################################################*/

void
Libnuceq__setSpeciesNseFactors( Libnuceq *self )
{

  xmlHashScan(
    self->pSpeciesHash,
    (xmlHashScanner) Libnuceq__setSpeciesNseFactorsCallback,
    self
  );

}

/*##############################################################################
// Libnuceq__setSpeciesNseFactorsCallback()
//############################################################################*/

void
Libnuceq__setSpeciesNseFactorsCallback(
  Libnuceq__Species *p_eq_species,
  Libnuceq *p_equil,
  xmlChar *sx_name
)
{

  if( !sx_name ) LIBNUCEQ__ERROR( "Invalid species" );

  p_eq_species->dNseFactor =
    Libnucnet__Nuc__computeSpeciesNseFactor(
      p_equil->pNuc,
      p_eq_species->pNucSpecies,
      p_equil->dT9,
      p_equil->dRho
    );

}

/*##############################################################################
// Libnuceq__setSpeciesNseCorrectionFactors()
//############################################################################*/

void
Libnuceq__setSpeciesNseCorrectionFactors( Libnuceq *self )
{

  xmlHashScan(
    self->pSpeciesHash,
    (xmlHashScanner) Libnuceq__setSpeciesNseCorrectionFactorsCallback,
    self
  );

}

/*##############################################################################
// Libnuceq__setSpeciesNseCorrectionFactorsCallback()
//############################################################################*/

void
Libnuceq__setSpeciesNseCorrectionFactorsCallback(
  Libnuceq__Species *p_eq_species,
  Libnuceq *p_equil,
  xmlChar *sx_name
)
{

  if( !sx_name ) LIBNUCEQ__ERROR( "Invalid species" );

  p_eq_species->dNseCorrectionFactor =
    p_equil->pfCorr(
      p_eq_species->pNucSpecies,
      p_equil->dT9,
      p_equil->dRho,
      Libnuceq__computeZMoment( p_equil, 1 ),
      p_equil->pCorrData
    );

}

/*##############################################################################
// Libnuceq__Species__zeroNseCorrectionFactor()
//############################################################################*/

void
Libnuceq__Species__zeroNseCorrectionFactor(
  Libnuceq__Species *p_eq_species,
  void *p_data,
  xmlChar *sx_name
)
{

  if( p_data || !sx_name ) LIBNUCEQ__ERROR( "Invalid input" );

  p_eq_species->dNseCorrectionFactor = 0;

}

/*##############################################################################
// Libnuceq__Species__new()
//############################################################################*/

Libnuceq__Species *
Libnuceq__Species__new(
  Libnucnet__Species *p_nuc_species
)
{

  Libnuceq__Species *self =
    ( Libnuceq__Species * ) malloc( sizeof( Libnuceq__Species ) );

  self->pNucSpecies = p_nuc_species;
  self->dAbundance = 0.;
  self->dNseFactor = 0.;
  self->dNseCorrectionFactor = 0.;

  return self;

}

/*##############################################################################
// Libnuceq__Species__free()
//############################################################################*/

void
Libnuceq__Species__free(
  Libnuceq__Species *p_eq_species,
  xmlChar *sx_name
)
{

  if( !sx_name ) LIBNUCEQ__ERROR( "Invalid species" );

  free( p_eq_species );

}

/*##############################################################################
// Libnuceq__getAbundances().
//############################################################################*/

gsl_vector *
Libnuceq__getAbundances(
  const Libnuceq *self
)
{

  gsl_vector *p_vector;

  if( !self ) LIBNUCEQ__ERROR( "Invalid species" );

  p_vector =
    gsl_vector_calloc( Libnuceq__getNumberOfSpecies( self ) );

  xmlHashScan(
    self->pSpeciesHash,
    (xmlHashScanner) Libnuceq__getAbundancesCallback,
    p_vector
  );

  return p_vector;

}

/*##############################################################################
// Libnuceq__getAbundancesCallback().
//############################################################################*/

void
Libnuceq__getAbundancesCallback(
  Libnuceq__Species *p_eq_species,
  gsl_vector *p_vector,
  xmlChar *sx_name
)
{

  if( !sx_name ) LIBNUCEQ__ERROR( "No such species" );

  gsl_vector_set(
    p_vector,
    Libnucnet__Species__getIndex(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ),
    p_eq_species->dAbundance
  );

}

/*##############################################################################
// Libnuceq__Species__getAbundance().
//############################################################################*/

double
Libnuceq__Species__getAbundance(
  Libnuceq__Species *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid species" );

  return self->dAbundance;

}

/*##############################################################################
// Libnuceq__Species__getNucSpecies().
//############################################################################*/

Libnucnet__Species *
Libnuceq__Species__getNucSpecies( Libnuceq__Species *self )
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid species" );

  return self->pNucSpecies;

}

/*##############################################################################
// Libnuceq__computeAbundances()
//############################################################################*/

void
Libnuceq__computeAbundances(
  Libnuceq *self
)
{

  xmlHashScan(
    self->pSpeciesHash,
    (xmlHashScanner) Libnuceq__computeAbundancesCallback,
    self
  );

  xmlHashScan(
    self->pClusterHash,
    (xmlHashScanner) Libnuceq__compute_cluster_abundances,
    NULL
  );

} 

/*##############################################################################
// Libnuceq__computeAbundancesCallback()
//############################################################################*/

void
Libnuceq__computeAbundancesCallback(
  Libnuceq__Species *p_species,
  Libnuceq *p_equil,
  xmlChar *sx_name
)
{

  if( !sx_name ) LIBNUCEQ__ERROR( "Invalid species" );

  p_species->dAbundance =
    exp(
      GSL_MIN(
        D_EQ_EXP_MAX,
        Libnuceq__computeSpeciesBaseLogAbundance( p_equil, p_species )
      )
    );

}

/*##############################################################################
// Libnuceq__computeSpeciesBaseLogAbundance()
//############################################################################*/

double
Libnuceq__computeSpeciesBaseLogAbundance(
  Libnuceq *self,
  Libnuceq__Species *p_eq_species
) {

  double d_exp;

  d_exp =
    Libnucnet__Species__getZ(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ) * ( self->dMupkT - self->dMunkT )
    +
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ) * self->dMunkT;

  d_exp += p_eq_species->dNseFactor + p_eq_species->dNseCorrectionFactor;

  return d_exp;

}

/*##############################################################################
// Libnuceq__compute_cluster_abundances()
//############################################################################*/

void
Libnuceq__compute_cluster_abundances(
  Libnuceq__Cluster *p_cluster,
  void *p_data,
  xmlChar *sx_xpath
)
{

  if( !sx_xpath || p_data ) LIBNUCEQ__ERROR( "Invalid input" );

  if( WnMatrix__value_is_zero( p_cluster->dConstraint ) )
  {
    p_cluster->dMukT = GSL_NEGINF;
    xmlHashScan(
      p_cluster->pSpeciesHash,
      (xmlHashScanner) Libnuceq__set_cluster_abundances_to_zero,
      NULL
    );
  }
  else
    p_cluster->dMukT =
      Libnuceq__solveEquilibrium(
        (Libnuceq__solve_function) Libnuceq__cluster_abundance_function,
        p_cluster
      );
  
}

/*##############################################################################
// Libnuceq__set_cluster_abundances_to_zero()
//############################################################################*/

void
Libnuceq__set_cluster_abundances_to_zero(
  Libnuceq__Species *p_species,
  void *p_data,
  xmlChar *sx_name
)
{

  if( p_data || !sx_name ) LIBNUCEQ__ERROR( "Invalid input" );

  p_species->dAbundance = 0.;

}

/*##############################################################################
// Libnuceq__cluster_abundance_function()
//############################################################################*/

double
Libnuceq__cluster_abundance_function(
  double d_x,
  Libnuceq__Cluster *p_cluster
)
{

  typedef struct {
    Libnuceq__Cluster *pCluster;
    double dResult;
  } work;

  double d_result;
  work *p_work;

  p_cluster->dMukT = d_x;

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_work->pCluster = p_cluster;

  p_work->dResult = -p_cluster->dConstraint;

  xmlHashScan(
    p_cluster->pSpeciesHash,
    (xmlHashScanner) Libnuceq__cluster_abundance_function_callback,
    p_work
  );

  d_result = p_work->dResult;

  free( p_work );

  return d_result;

}

/*##############################################################################
// Libnuceq__cluster_abundance_function_callback()
//############################################################################*/

void
Libnuceq__cluster_abundance_function_callback(
  Libnuceq__Species *p_species,
  void *p_data,
  xmlChar *sx_name
)
{

  typedef struct {
    Libnuceq__Cluster *pCluster;
    double dResult;
  } work;

  double d_exp;
  work *p_work = ( work * ) p_data;

  if(
    xmlStrcmp( (const xmlChar *) EQ_NEUT, sx_name ) == 0 ||
    xmlStrcmp( (const xmlChar *) EQ_PROT, sx_name ) == 0
  )
    LIBNUCEQ__ERROR( "Cluster must not include neutron or proton" );

  d_exp =
    Libnuceq__computeSpeciesBaseLogAbundance(
      p_work->pCluster->pEquil,
      p_species
    );
  
  if( p_work->pCluster->pfPrefactorFunction )
    d_exp +=
      p_work->pCluster->pfPrefactorFunction(
        p_work->pCluster,
        p_species,
        p_work->pCluster->pPrefactorData
      );
  else
    d_exp += p_work->pCluster->dMukT;

  p_species->dAbundance = exp( GSL_MIN( D_EQ_EXP_MAX, d_exp ) );

  if( p_work->pCluster->pfConstraintFunction )
    p_work->dResult +=
      p_work->pCluster->pfConstraintFunction(
        p_species,
        p_work->pCluster->pConstraintData
      );
  else
    p_work->dResult += Libnuceq__Species__getAbundance( p_species );

}

/*##############################################################################
// Libnuceq__A_function()
//############################################################################*/

double
Libnuceq__A_function(
  double d_x,
  void *p_params
)
{

  typedef struct {
    double dResult1;
    double dResult2;
    Libnuceq *pEquil;
  } work;

  double d_result;
  work *p_work;

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work ) LIBNUCEQ__ERROR( "Couldn't allocate work structure" );

  p_work->pEquil = ( Libnuceq * ) p_params;
  p_work->dResult1 = -1.;

  p_work->pEquil->dMunkT = d_x;

  if( p_work->pEquil->pfCorr )
    Libnuceq__setSpeciesNseCorrectionFactors( p_work->pEquil );

  if( p_work->pEquil->pYe )
    p_work->pEquil->dMupkT =
      Libnuceq__solveEquilibrium(
        (Libnuceq__solve_function) Libnuceq__NSE_function,
        p_work
      );
  else
  {
    Libnuceq__initializeMuekT( p_work->pEquil );
    p_work->pEquil->dMupkT =
      Libnuceq__solveEquilibrium(
        (Libnuceq__solve_function) Libnuceq__WSE_function,
        p_work
      );
  }

  xmlHashScan(
    p_work->pEquil->pSpeciesHash,
    (xmlHashScanner) Libnuceq__A_function_callback,
    p_work
  );

  if( !gsl_finite( p_work->dResult1 ) ) {
      return GSL_FAILURE;
  }

  d_result = p_work->dResult1;

  free( p_work );

  return d_result;

}

/*##############################################################################
// Libnuceq__A_function_callback()
//############################################################################*/

void
Libnuceq__A_function_callback(
  Libnuceq__Species *p_eq_species,
  void *p_data,
  xmlChar *sx_name
)
{

  typedef struct {
    double dResult1;
    double dResult2;
    Libnuceq *pEquil;
  } work;

  work *p_work = ( work * ) p_data;    

  if( !sx_name ) LIBNUCEQ__ERROR( "Invalid species" );

  p_work->dResult1 +=
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ) * Libnuceq__Species__getAbundance( p_eq_species );
    
}

/*##############################################################################
// Libnuceq__WSE_function()
//############################################################################*/

double
Libnuceq__WSE_function(
  double d_x,
  void *p_params
)
{

  typedef struct {
    double dResult1;
    double dResult2;
    Libnuceq *pEquil;
  } work;

  work *p_work;

  p_work = ( work * ) p_params;

  p_work->dResult2 = 0.;

  p_work->pEquil->dMupkT = d_x;

  Libnuceq__computeAbundances( p_work->pEquil );

  xmlHashScan(
    p_work->pEquil->pSpeciesHash,
    (xmlHashScanner) Libnuceq__function2_callback,
    p_work
  );

  if( !gsl_finite( p_work->dResult2 ) ) {
      return GSL_FAILURE;
  }

  *(p_work->pEquil->pMuekT) =
    p_work->pEquil->dMunkT -
    p_work->pEquil->dMupkT +
    (
      Libnucnet__Species__getMassExcess(
        Libnucnet__Nuc__getSpeciesByName(
          Libnuceq__getNuc( p_work->pEquil ),
          EQ_NEUT
        )
      ) -
      Libnucnet__Species__getMassExcess(
        Libnucnet__Nuc__getSpeciesByName(
          Libnuceq__getNuc( p_work->pEquil ),
          EQ_PROT
        )
      )
    ) /
    (
      GSL_CONST_CGSM_BOLTZMANN *
      Libnuceq__getT9( p_work->pEquil ) *
      GSL_CONST_NUM_GIGA *
      WN_ERGS_TO_MEV
    );

  if( p_work->pEquil->pfWseCorrectionFunction )
    *(p_work->pEquil->pMuekT) +=
      p_work->pEquil->pfWseCorrectionFunction(
         p_work->pEquil,
         p_work->pEquil->pWseCorrectionFunctionData
      );

  if( *p_work->pEquil->pMuekT > -3000 )
    p_work->dResult2 -=
      Libnuceq__computeElectronNumberDensity( p_work->pEquil ) /
      ( p_work->pEquil->dRho * GSL_CONST_NUM_AVOGADRO );
        
  return p_work->dResult2;

}

/*##############################################################################
// Libnuceq__computeElectronNumberDensity().
//############################################################################*/

double
Libnuceq__computeElectronNumberDensity(
  Libnuceq *p_equil
)
{

  Libstatmech__Fermion *p_electron;
  double d_result;

  p_electron =
    Libstatmech__Fermion__new(
      S_NUCEQ_ELECTRON,
      WN_MASS_ELECTRON_MEV,
      2,
      -1.
  );

  if( p_equil->pfElectronFunction || p_equil->pfElectronIntegrand )
    Libstatmech__Fermion__updateQuantity(
      p_electron,
      S_NUMBER_DENSITY,
      (Libstatmech__Fermion__Function) p_equil->pfElectronFunction,
      (Libstatmech__Fermion__Integrand) p_equil->pfElectronIntegrand
    );

  d_result =
    Libstatmech__Fermion__computeQuantity(
      p_electron,
      S_NUMBER_DENSITY,
      Libnuceq__getT9( p_equil ) * GSL_CONST_NUM_GIGA,
      Libnuceq__getMuekT( p_equil ),
      p_equil->pElectronFunctionData,
      p_equil->pElectronIntegrandData
    );

  Libstatmech__Fermion__free( p_electron );

  return d_result;

}

/*##############################################################################
// Libnuceq__NSE_function()
//############################################################################*/

double
Libnuceq__NSE_function(
  double d_x,
  void *p_params
)
{

  typedef struct {
    double dResult1;
    double dResult2;
    Libnuceq *pEquil;
  } work;

  work *p_work;

  p_work = ( work * ) p_params;

  p_work->dResult2 = -(*(p_work->pEquil->pYe));

  p_work->pEquil->dMupkT = d_x;

  Libnuceq__computeAbundances( p_work->pEquil );

  xmlHashScan(
    p_work->pEquil->pSpeciesHash,
    (xmlHashScanner) Libnuceq__function2_callback,
    p_work
  );

  if( !gsl_finite( p_work->dResult2 ) ) {
      return GSL_FAILURE;
  }

  return p_work->dResult2;

}

/*##############################################################################
// Libnuceq__function2_callback()
//############################################################################*/

void
Libnuceq__function2_callback(
  Libnuceq__Species *p_eq_species,
  void *p_data,
  xmlChar *sx_name
)
{

  typedef struct {
    double dResult1;
    double dResult2;
    Libnuceq *pEquil;
  } work;

  work *p_work = ( work * ) p_data;    

  if( !sx_name ) LIBNUCEQ__ERROR( "Invalid species" );

  p_work->dResult2 +=
    Libnucnet__Species__getZ(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ) * Libnuceq__Species__getAbundance( p_eq_species );

}

/*##############################################################################
// Libnuceq__iterateClusters().
//############################################################################*/

void
Libnuceq__iterateClusters(
  Libnuceq *self,
  Libnuceq__Cluster__iterateFunction pf_func,
  void *p_data
)
{

  typedef struct {
    Libnuceq__Cluster__iterateFunction pfFunc;
    void *pData;
  } work;

  work *p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_work->pfFunc = pf_func;
  p_work->pData = p_data;

  xmlHashScan(
    self->pClusterHash,
    (xmlHashScanner) Libnuceq__iterateClustersHelper,
    p_work
  );

  free( p_work );

}

/*##############################################################################
// Libnuceq__iterateClustersHelper().
//############################################################################*/

void
Libnuceq__iterateClustersHelper(
  Libnuceq__Cluster *p_cluster,
  void *p_data,
  xmlChar *sx_cluster
)
{

  typedef struct {
    Libnuceq__Cluster__iterateFunction pfFunc;
    void *pData;
  } work;

  work *p_work = ( work * ) p_data;  

  if( !sx_cluster ) LIBNUCEQ__ERROR( "Invalid cluster" );

  p_work->pfFunc( p_cluster, p_work->pData );

}

/*##############################################################################
// Libnuceq__iterateSpecies().
//############################################################################*/

void
Libnuceq__iterateSpecies(
  const Libnuceq *self,
  Libnuceq__Species__iterateFunction pf_func,
  void *p_data
)
{

  Libnuceq__Species **p_species_array;
  size_t i, i_species;

  p_species_array =
    Libnuceq__createSpeciesArray( self );

  i_species = Libnuceq__getNumberOfSpecies( self );

  for( i = 0; i < i_species; i++ )
    if( pf_func( p_species_array[i], p_data ) == 0 ) break;

  free( p_species_array );

} 

/*##############################################################################
// Libnuceq__createSpeciesArray().
//############################################################################*/

Libnuceq__Species **
Libnuceq__createSpeciesArray(
  const Libnuceq *self
)
{

  typedef struct {
    size_t iIndex;
    Libnuceq__Species **pSpeciesArray;
  } work;

  Libnuceq__Species **p_species_array;

  work *p_work;

  if( !self ) LIBNUCEQ__ERROR( "Invalid input" );

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_work->pSpeciesArray =
    ( Libnuceq__Species ** )
    malloc(
      sizeof( Libnuceq__Species * ) *
      Libnuceq__getNumberOfSpecies( self )
    );

  if( !p_work->pSpeciesArray )
    LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_work->iIndex = 0;

  xmlHashScan(
    self->pSpeciesHash,
    (xmlHashScanner) Libnuceq__create_species_array_callback,
    p_work
  );

  p_species_array = p_work->pSpeciesArray;

  free( p_work );

  qsort(
    p_species_array,
    Libnuceq__getNumberOfSpecies( self ),
    sizeof( Libnuceq__Species * ),
    Libnuceq__sorter
  );

  return p_species_array;

}
  
/*##############################################################################
// Libnuceq__create_species_array_callback().
//############################################################################*/

void
Libnuceq__create_species_array_callback(
  Libnuceq__Species *p_species,
  void *p_data,
  const xmlChar *sx_species
)
{

  typedef struct {
    size_t iIndex;
    Libnuceq__Species **pSpeciesArray;
  } work;

  work *p_work = ( work * ) p_data;

  if( !sx_species )
    LIBNUCNET__NUC__ERROR( "No such species" );

  p_work->pSpeciesArray[p_work->iIndex++] = p_species;

}

/*##############################################################################
// Libnuceq__sorter().
//############################################################################*/

int
Libnuceq__sorter(
  const void *p_1,
  const void *p_2
)
{

  Libnuceq__Species *p_species_1, *p_species_2;

  p_species_1 = *(Libnuceq__Species * const *) p_1;
  p_species_2 = *(Libnuceq__Species * const *) p_2;

  if( p_species_1->pNucSpecies->iIndex < p_species_2->pNucSpecies->iIndex )
    return -1;
  else
    return 1;

}

/*##############################################################################
// Libnuceq__Cluster__iterateSpecies().
//############################################################################*/

void
Libnuceq__Cluster__iterateSpecies(
  const Libnuceq__Cluster *self,
  Libnuceq__Species__iterateFunction pf_func,
  void *p_data
)
{

  Libnuceq__Species **p_species_array;
  size_t i, i_species;

  p_species_array =
    Libnuceq__Cluster__createSpeciesArray( self );

  i_species = Libnuceq__Cluster__getNumberOfSpecies( self );

  for( i = 0; i < i_species; i++ )
    if( pf_func( p_species_array[i], p_data ) == 0 ) break;

  free( p_species_array );

}

/*##############################################################################
// Libnuceq__Cluster__createSpeciesArray().
//############################################################################*/

Libnuceq__Species **
Libnuceq__Cluster__createSpeciesArray(
  const Libnuceq__Cluster *self
)
{

  typedef struct {
    size_t iIndex;
    Libnuceq__Species **pSpeciesArray;
  } work;

  Libnuceq__Species **p_species_array;

  work *p_work;

  if( !self ) LIBNUCEQ__ERROR( "Invalid input" );

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_work->pSpeciesArray =
    ( Libnuceq__Species ** )
    malloc(
      sizeof( Libnuceq__Species * ) *
      Libnuceq__Cluster__getNumberOfSpecies( self )
    );

  if( !p_work->pSpeciesArray )
    LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_work->iIndex = 0;

  xmlHashScan(
    self->pSpeciesHash,
    (xmlHashScanner) Libnuceq__create_species_array_callback,
    p_work
  );

  p_species_array = p_work->pSpeciesArray;

  free( p_work );

  qsort(
    p_species_array,
    Libnuceq__Cluster__getNumberOfSpecies( self ),
    sizeof( Libnuceq__Species * ),
    Libnuceq__sorter
  );

  return p_species_array;

}
  
/*##############################################################################
// Libnuceq__getNumberOfSpecies().
//############################################################################*/

size_t
Libnuceq__getNumberOfSpecies(
  const Libnuceq *self
)
{

  return
    (size_t) xmlHashSize( self->pSpeciesHash );

}

/*##############################################################################
// Libnuceq__Cluster__getNumberOfSpecies().
//############################################################################*/

size_t
Libnuceq__Cluster__getNumberOfSpecies(
  const Libnuceq__Cluster *self
)
{

  return
    (size_t) xmlHashSize( self->pSpeciesHash );

}

/*##############################################################################
// Libnuceq__getSpeciesByName().
//############################################################################*/

Libnuceq__Species *
Libnuceq__getSpeciesByName(
  Libnuceq *self,
  const char *s_name
)
{

  return
    ( Libnuceq__Species * )
    xmlHashLookup( self->pSpeciesHash, ( const xmlChar * ) s_name );

}

/*##############################################################################
// Libnuceq__bracket_of_root_function()
//############################################################################*/

int
Libnuceq__bracket_root_of_function(
  gsl_function F, double *p_x1, double *p_x2, void *params
)
{

  double d_f1, d_f2, d_factor = 1.6;
  int i_iter = 0, i_max_iter = 1000;

  d_f1 = F.function( *p_x1, params );
  d_f2 = F.function( *p_x2, params );

  while( i_iter < i_max_iter )
  {
    i_iter++;
    if( d_f1 * d_f2 < 0. ) return 1;
    if( fabs( d_f1 ) < fabs( d_f2 ) )
    {
      *p_x1 += d_factor * ( *p_x1 - *p_x2 );
      d_f1 = F.function( *p_x1, params );
    }
    else if( gsl_fcmp( d_f1, d_f2, D_EQ_DIFF ) == 0 )
    {
      *p_x1 += d_factor * ( *p_x1 - *p_x2 );
      *p_x2 += d_factor * ( *p_x2 - *p_x1 );
      d_f1 = F.function( *p_x1, params );
      d_f2 = F.function( *p_x2, params );
    }
    else
    {
      *p_x2 += d_factor * ( *p_x2 - *p_x1 );
      d_f2 = F.function( *p_x2, params );
    }
  }

  return 0;

}

/*##############################################################################
// Libnuceq__setYe().
//############################################################################*/

void
Libnuceq__setYe(
  Libnuceq *self,
  double d_ye
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  if( d_ye < 0. || d_ye > 1. ) LIBNUCEQ__ERROR( "Invalid Ye" );

  if( self->pYe ) free( self->pYe );

  self->pYe = ( double * ) malloc( sizeof( double ) );

  if( !self->pYe ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  *self->pYe = d_ye;

}

/*##############################################################################
// Libnuceq__initializeMuekT().
//############################################################################*/

void
Libnuceq__initializeMuekT(
  Libnuceq *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  if( self->pMuekT ) free( self->pMuekT );

  self->pMuekT = ( double * ) malloc( sizeof( double ) );

  if( !self->pMuekT ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

}

/*##############################################################################
// Libnuceq__clearYe().
//############################################################################*/

void
Libnuceq__clearYe(
  Libnuceq *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  free( self->pYe );
  self->pYe = NULL;

}

/*##############################################################################
// Libnuceq__newCluster().
//############################################################################*/

Libnuceq__Cluster *
Libnuceq__newCluster(
  Libnuceq *self,
  const char *s_xpath
)
{

  Libnuceq__Cluster *p_cluster;

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  p_cluster = ( Libnuceq__Cluster * ) malloc( sizeof( Libnuceq__Cluster ) );

  p_cluster->pEquil = self;

  if( !p_cluster ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_cluster->sxXPath = xmlCharStrdup( s_xpath );

  if( !p_cluster->sxXPath ) LIBNUCEQ__ERROR( "Couldn't allocate memory " );

  p_cluster->pSpeciesHash = xmlHashCreate( 0 );

  Libnuceq__assignSpeciesToCluster( self, p_cluster );

  p_cluster->dConstraint = 0;

  p_cluster->pfConstraintFunction = NULL;

  p_cluster->pConstraintData = NULL;

  p_cluster->pfPrefactorFunction = NULL;

  p_cluster->pPrefactorData = NULL;

  if(
    xmlHashAddEntry(
      self->pClusterHash,
      p_cluster->sxXPath,
      p_cluster
    ) == -1
  ) LIBNUCEQ__ERROR( "Couldn't add cluster" );

  return p_cluster;

}

/*##############################################################################
// Libnuceq__Cluster__updateConstraint().
//############################################################################*/

void
Libnuceq__Cluster__updateConstraint(
  Libnuceq__Cluster *self,
  double d_constraint
)
{

  self->dConstraint = d_constraint;

}

/*##############################################################################
// Libnuceq__removeCluster().
//############################################################################*/

int
Libnuceq__removeCluster(
  Libnuceq *self,
  Libnuceq__Cluster *p_cluster
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  return
    xmlHashRemoveEntry(
      self->pClusterHash,
      p_cluster->sxXPath,
      (xmlHashDeallocator) Libnuceq__Cluster__hash_free
    ) + 1;

}

/*##############################################################################
// Libnuceq__assignSpeciesToCluster().
//############################################################################*/

void
Libnuceq__assignSpeciesToCluster(
  Libnuceq *self,
  Libnuceq__Cluster *p_cluster
)
{

  typedef struct {
    Libnuceq *pEquil;
    Libnuceq__Cluster *pCluster;
  } work;

  Libnucnet__NucView *p_cluster_nuc_view;

  work *p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_work->pEquil = self;
  p_work->pCluster = p_cluster;

  p_cluster_nuc_view =
    Libnucnet__NucView__new(
      self->pNuc,
      (const char *) p_cluster->sxXPath
    );

  xmlHashScan(
    p_cluster_nuc_view->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnuceq__assignSpeciesToClusterCallback,
    p_work
  );

  Libnucnet__NucView__free( p_cluster_nuc_view );

  free( p_work );

}

/*##############################################################################
// Libnuceq__assignSpeciesToClusterCallback().
//############################################################################*/

void
Libnuceq__assignSpeciesToClusterCallback(
  Libnucnet__Species *p_nuc_species,
  void *p_data,
  xmlChar *sx_name
)
{

  typedef struct {
    Libnuceq *pEquil;
    Libnuceq__Cluster *pCluster;
  } work;

  Libnuceq__Species *p_species;
  work *p_work = ( work * ) p_data;

  if( !p_nuc_species ) LIBNUCEQ__ERROR( "Invalid species" );

  p_species =
    (Libnuceq__Species *)
    xmlHashLookup(
      p_work->pEquil->pSpeciesHash,
      sx_name
    );

  if( !p_species ) LIBNUCEQ__ERROR( "Couldn't look up species" );

  if(
    xmlHashAddEntry(
      p_work->pCluster->pSpeciesHash,
      sx_name,
      p_species
    ) == -1
  )
    LIBNUCEQ__ERROR( "Couldn't add species" );

}

/*##############################################################################
// Libnuceq__Cluster__free().
//############################################################################*/

void
Libnuceq__Cluster__free( Libnuceq__Cluster *self )
{

  Libnuceq__Cluster__hash_free( self, self->sxXPath );

}
  
/*##############################################################################
// Libnuceq__Cluster__hash_free().
//############################################################################*/

void
Libnuceq__Cluster__hash_free( Libnuceq__Cluster *self, xmlChar *sx_xpath )
{

  if( !sx_xpath ) LIBNUCEQ__ERROR( "No such xpath" );

  xmlFree( self->sxXPath );
  
  xmlHashFree(
    self->pSpeciesHash,
    NULL
  );

  free( self );

} 

/*##############################################################################
// Libnuceq__getMunkT()
//############################################################################*/

double
Libnuceq__getMunkT(
  const Libnuceq *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  return self->dMunkT;

}

/*##############################################################################
// Libnuceq__getMupkT()
//############################################################################*/

double
Libnuceq__getMupkT(
  const Libnuceq *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  return self->dMupkT;

}

/*##############################################################################
// Libnuceq__getMuekT()
//############################################################################*/

double
Libnuceq__getMuekT(
  const Libnuceq *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  if( !self->pMuekT )
    LIBNUCEQ__ERROR( "Electron chemical potential not available" );

  return *(self->pMuekT);

}

/*##############################################################################
// Libnuceq__Cluster__getXPathString()
//############################################################################*/

const char *
Libnuceq__Cluster__getXPathString(
  const Libnuceq__Cluster *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid cluster" );

  return (const char *) self->sxXPath;

}

/*##############################################################################
// Libnuceq__Cluster__getConstraint()
//############################################################################*/

double
Libnuceq__Cluster__getConstraint(
  const Libnuceq__Cluster *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid cluster" );

  return self->dConstraint;

}

/*##############################################################################
// Libnuceq__Cluster__getMukT()
//############################################################################*/

double
Libnuceq__Cluster__getMukT(
  const Libnuceq__Cluster *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid cluster" );

  return self->dMukT;

}

/*##############################################################################
// Libnuceq__getT9()
//############################################################################*/

double
Libnuceq__getT9(
  const Libnuceq *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  return self->dT9;

}

/*##############################################################################
// Libnuceq__getRho()
//############################################################################*/

double
Libnuceq__getRho(
  const Libnuceq *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  return self->dRho;

}

/*##############################################################################
// Libnuceq__getCluster()
//############################################################################*/

Libnuceq__Cluster *
Libnuceq__getCluster(
  const Libnuceq *self,
  const char *s_xpath
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  return
    ( Libnuceq__Cluster * )
    xmlHashLookup(
      self->pClusterHash,
      ( const xmlChar * ) s_xpath
    );
  
}

/*##############################################################################
// Libnuceq__computeZMoment()
//############################################################################*/

double
Libnuceq__computeZMoment(
  const Libnuceq *self,
  unsigned int i_index
)
{

  typedef struct {
    double dResult;
    unsigned int iIndex;
  } work;

  double d_result;
  work *p_work;

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  p_work = ( work * ) malloc( sizeof( work ) );
  if( !p_work ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_work->dResult = 0.;
  p_work->iIndex = i_index;

  xmlHashScan(
    self->pSpeciesHash,
    (xmlHashScanner) Libnuceq__computeZMomentCallback,
    p_work
  );

  d_result = p_work->dResult;

  free( p_work );

  return d_result;

}

/*##############################################################################
// Libnuceq__computeZMomentCallback()
//############################################################################*/

void
Libnuceq__computeZMomentCallback(
  Libnuceq__Species *p_species,
  void *p_data,
  xmlChar *sx_name
)
{

  typedef struct {
    double dResult;
    unsigned int iIndex;
  } work;

  work *p_work = ( work * ) p_data;

  if( !sx_name ) LIBNUCEQ__ERROR( "No such species" );

  p_work->dResult +=
    Libnuceq__Species__getAbundance( p_species ) *
    pow(
      (double)
      Libnucnet__Species__getZ(
        Libnuceq__Species__getNucSpecies( p_species )
      ),
      (double) p_work->iIndex
    );

}

/*##############################################################################
// Libnuceq__computeAMoment()
//############################################################################*/

double
Libnuceq__computeAMoment(
  const Libnuceq *self,
  unsigned int i_index
)
{

  typedef struct {
    double dResult;
    unsigned int iIndex;
  } work;

  double d_result;
  work *p_work;

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  p_work = ( work * ) malloc( sizeof( work ) );
  if( !p_work ) LIBNUCEQ__ERROR( "Couldn't allocate memory" );

  p_work->dResult = 0.;
  p_work->iIndex = i_index;

  xmlHashScan(
    self->pSpeciesHash,
    (xmlHashScanner) Libnuceq__computeAMomentCallback,
    p_work
  );

  d_result = p_work->dResult;

  free( p_work );

  return d_result;

}

/*##############################################################################
// Libnuceq__computeAMomentCallback()
//############################################################################*/

void
Libnuceq__computeAMomentCallback(
  Libnuceq__Species *p_species,
  void *p_data,
  xmlChar *sx_name
)
{

  typedef struct {
    double dResult;
    unsigned int iIndex;
  } work;

  work *p_work = ( work * ) p_data;

  if( !sx_name ) LIBNUCEQ__ERROR( "No such species" );

  p_work->dResult +=
    Libnuceq__Species__getAbundance( p_species ) *
    pow(
      (double)
      Libnucnet__Species__getA(
        Libnuceq__Species__getNucSpecies( p_species )
      ),
      (double) p_work->iIndex
    );

}

/*##############################################################################
// Libnuceq__Cluster__updateConstraintFunction()
//############################################################################*/

void
Libnuceq__Cluster__updateConstraintFunction(
  Libnuceq__Cluster *self,
  Libnuceq__Cluster__constraint_function pf_func,
  void *p_data
)
{

  self->pfConstraintFunction = pf_func;
  self->pConstraintData = p_data;

} 

/*##############################################################################
// Libnuceq__Cluster__updatePrefactorFunction()
//############################################################################*/

void
Libnuceq__Cluster__updatePrefactorFunction(
  Libnuceq__Cluster *self,
  Libnuceq__Cluster__prefactorFunction pf_func,
  void *p_data
)
{

  self->pfPrefactorFunction = pf_func;
  self->pPrefactorData = p_data;

} 

/*##############################################################################
// Libnuceq__updateWseCorrectionFunction()
//############################################################################*/

void
Libnuceq__updateWseCorrectionFunction(
  Libnuceq *self,
  Libnuceq__wseCorrectionFunction pf_func,
  void *p_data
)
{

  self->pfWseCorrectionFunction = pf_func;
  self->pWseCorrectionFunctionData = p_data;

} 

/*##############################################################################
// Libnuceq__setNseCorrectionFactorFunction().
//############################################################################*/

void
Libnuceq__setNseCorrectionFactorFunction(
  Libnuceq *self,
  Libnucnet__Species__nseCorrectionFactorFunction pf_func,
  void *p_data
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  self->pfCorr = pf_func;
  self->pCorrData = p_data;

} 

/*##############################################################################
// Libnuceq__clearNseCorrectionFactorFunction().
//############################################################################*/

void
Libnuceq__clearNseCorrectionFactorFunction(
  Libnuceq *self
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid equilibrium" );

  self->pfCorr = NULL;
  self->pCorrData = NULL;

  xmlHashScan(
    self->pSpeciesHash,
    (xmlHashScanner) Libnuceq__Species__zeroNseCorrectionFactor,
    NULL
  );

} 

/*##############################################################################
// Libnuceq__updateUserElectronNumberDensity().
//############################################################################*/

void
Libnuceq__updateUserElectronNumberDensity(
  Libnuceq *self,
  Libstatmech__Fermion__Function pf_func,
  Libstatmech__Fermion__Integrand pf_integrand,
  void *p_function_data,
  void *p_integrand_data
)
{

  if( !self ) LIBNUCEQ__ERROR( "Invalid input" );

  self->pfElectronFunction = pf_func;
  self->pElectronFunctionData = p_function_data;
  self->pfElectronIntegrand = pf_integrand;
  self->pElectronIntegrandData = p_integrand_data;

}

/*##############################################################################
// Libnuceq__getNumberOfClusters().
//############################################################################*/

size_t
Libnuceq__getNumberOfClusters( const Libnuceq * self )
{

  return (size_t) xmlHashSize( self->pClusterHash );

}

/*##############################################################################
// Libnuceq__copy_species_to_new_cluster()
//############################################################################*/

int
Libnuceq__copy_species_to_new_cluster(
  Libnuceq__Species * p_eq_species,
  Libnuceq__Cluster * p_cluster
)
{

  Libnuceq__Species * p_new_eq_species;

  p_new_eq_species =
    Libnuceq__getSpeciesByName(
      p_cluster->pEquil,
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      )
    );

  if( !p_new_eq_species ) LIBNUCEQ__ERROR( "Couldn't find species" );

  if(
    xmlHashAddEntry(
      p_cluster->pSpeciesHash,
      (const xmlChar *)
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      ),
      p_new_eq_species
    ) == -1
  )
    LIBNUCEQ__ERROR( "Couldn't add species" );

  return 1;

}

/*##############################################################################
// Libnuceq__copy_cluster()
//############################################################################*/

int
Libnuceq__copy_cluster(
  Libnuceq__Cluster * p_cluster,
  Libnuceq * p_equil
)
{

  Libnuceq__Cluster * p_new_cluster =
    (Libnuceq__Cluster *) malloc( sizeof( Libnuceq__Cluster ) );

  p_new_cluster->pEquil = p_equil;
  p_new_cluster->sxXPath = xmlStrdup( p_cluster->sxXPath );

  p_new_cluster->dConstraint = 0;
  p_new_cluster->dMukT = 0;

  p_new_cluster->pfConstraintFunction = NULL;
  p_new_cluster->pConstraintData = NULL;
  p_new_cluster->pfPrefactorFunction = NULL;
  p_new_cluster->pPrefactorData = NULL;

  p_new_cluster->pSpeciesHash = xmlHashCreate( 0 );

  Libnuceq__Cluster__iterateSpecies(
    p_cluster,
    (Libnuceq__Species__iterateFunction) Libnuceq__copy_species_to_new_cluster,
    p_new_cluster
  );

  xmlHashAddEntry(
    p_equil->pClusterHash,
    p_new_cluster->sxXPath,
    p_new_cluster
  );

  return 1;

}
  
/*##############################################################################
// Libnuceq__copy_clusters()
//############################################################################*/

void
Libnuceq__copy_clusters( Libnuceq * p_destination, Libnuceq * p_source )
{

  if( !p_destination || !p_source ) LIBNUCEQ__ERROR( "Invalid input" );

  if( p_destination->pNuc != p_source->pNuc )
    LIBNUCEQ__ERROR(
      "Destination nuclear collection does not equal that of source"
    );

  xmlHashFree(
    p_destination->pClusterHash,
    (xmlHashDeallocator) Libnuceq__Cluster__hash_free
  );

  p_destination->pClusterHash = xmlHashCreate( I_EQ_HASH_MAX );

  Libnuceq__iterateClusters(
    p_source,
    (Libnuceq__Cluster__iterateFunction) Libnuceq__copy_cluster,
    p_destination
  );
  
}

