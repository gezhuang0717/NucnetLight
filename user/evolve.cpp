//////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Code for evolving a network.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes. 
//##############################################################################

#include "user/evolve.h"

/**
 * @brief A namespace for user-defined functions.
 */
namespace user
{

//##############################################################################
// evolve()
//##############################################################################

int
evolve( 
  nnt::Zone& zone
) {

  WnMatrix *p_matrix; 
  size_t i_iter;
  gsl_vector *p_y_old, *p_rhs, *p_sol, *p_work;
  double d_dt;
  std::pair<double,double> check;
  
  //==========================================================================
  // Evolve NSE + weak rates, if appropriate.
  //==========================================================================

  if( zone.hasProperty( nnt::s_USE_HI_T_EQUIL ) )
  {

    if( zone.getProperty( nnt::s_USE_HI_T_EQUIL ) == "yes" )
    {

      if(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) >
        boost::lexical_cast<double>( zone.getProperty( nnt::s_HI_T9_EQUIL ) )
      )
      {
        evolve_nse_plus_weak_rates( zone );
        return 1;
      }

    }

  }

  //============================================================================
  // Get timestep. 
  //============================================================================

  d_dt = boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  //============================================================================
  // Save the old abundances.
  //============================================================================

  p_y_old = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  //============================================================================
  // Newton-Raphson Iterations.
  //============================================================================

  for( i_iter = 1; i_iter <= I_ITMAX; i_iter++ ) {

    //--------------------------------------------------------------------------
    // Get matrix and rhs vector.
    //--------------------------------------------------------------------------

    boost::tie( p_matrix, p_rhs ) = get_evolution_matrix_and_vector( zone );

    //--------------------------------------------------------------------------
    // Add 1/dt to diagonal.
    //--------------------------------------------------------------------------

    WnMatrix__addValueToDiagonals(
      p_matrix,
      1.0 / boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
    );

    //--------------------------------------------------------------------------
    // Correct vector for iteration.
    //--------------------------------------------------------------------------

    p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );
    gsl_vector_sub( p_work, p_y_old );
    gsl_vector_scale( p_work, 1. / d_dt );
    gsl_vector_sub( p_rhs, p_work );

    gsl_vector_free( p_work );

    //--------------------------------------------------------------------------
    // Solve matrix equation.
    //--------------------------------------------------------------------------

    p_sol = solve_matrix_for_zone( zone, p_matrix, p_rhs );

    //--------------------------------------------------------------------------
    // Check solution.
    //--------------------------------------------------------------------------

    check = check_matrix_solution( zone, p_sol );

    //--------------------------------------------------------------------------
    // Update abundances.
    //--------------------------------------------------------------------------

    p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

    gsl_vector_add( p_work, p_sol );

    Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_work );

    gsl_vector_free( p_work );

    //--------------------------------------------------------------------------
    // Free matrix, p_rhs, and p_sol. Remember
    // Libnucnet__Zone__computeJacobianMatrix returns a new matrix and
    // Libnucnet__computeFlowVector and WnMatrix__solve return new gsl_vectors
    // each time they are called.
    //--------------------------------------------------------------------------

    WnMatrix__free( p_matrix ); 
    gsl_vector_free( p_rhs );
    gsl_vector_free( p_sol );

    //--------------------------------------------------------------------------
    // Exit iterations if converged.
    //--------------------------------------------------------------------------

    if( zone.hasProperty( nnt::s_NEWTON_RAPHSON_CONVERGE ) )
    {
      if(
        check.first <
        boost::lexical_cast<double>(
          zone.getProperty( nnt::s_NEWTON_RAPHSON_CONVERGE )
        )
      )
        break;
    }
    else
    {
      if( check.first < D_MIN ) break;
    }

    //--------------------------------------------------------------------------
    // Return with negative value if large negative abundances.
    //--------------------------------------------------------------------------

    if( zone.hasProperty( nnt::s_LARGE_NEG_ABUND_THRESHOLD ) )
    {
      if( !is_nonneg_abunds( zone ) ) return -1;
    }
      
  }

  //==========================================================================
  // Update abundance changes.
  //==========================================================================

  p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  gsl_vector_sub( p_work, p_y_old );

  Libnucnet__Zone__updateAbundanceChanges( zone.getNucnetZone(), p_work );

  gsl_vector_free( p_work );
  
  //==========================================================================
  // Free allocated memory and return.
  //==========================================================================

  gsl_vector_free( p_y_old );

  return (int) i_iter;

}

//##############################################################################
// safe_evolve()
//##############################################################################

void
safe_evolve(
  nnt::Zone& zone,
  double d_dt,
  double d_tend
)
{

  double d_xsum;
  gsl_vector *p_y, *p_y_old;
  double d_t;

  Libnucnet__Zone * p_zone = zone.getNucnetZone();

  p_y_old = Libnucnet__Zone__getAbundances( p_zone );

  zone.updateProperty(
    nnt::s_DTIME,
    boost::lexical_cast<std::string>( d_dt )
  );

  evolve( zone );

  d_xsum = 1. - Libnucnet__Zone__computeAMoment( p_zone, 1 );

  if( fabs( d_xsum ) > D_X_EPS || !is_nonneg_abunds( zone ) )
  {

    while(
      ( fabs( d_xsum ) > D_X_EPS || !is_nonneg_abunds( zone ) ) &&
      d_dt > 1.e-20
    )
    {
      Libnucnet__Zone__updateAbundances( p_zone, p_y_old );
      d_dt /= 10;
      zone.updateProperty(
        nnt::s_DTIME,
        boost::lexical_cast<std::string>( d_dt )
      );
      evolve( zone );
      d_xsum = 1. - Libnucnet__Zone__computeAMoment( p_zone, 1 );
    }

  }

  d_t = d_dt;

  while( d_t < d_tend )
  {
    evolve( zone );
    d_t += d_dt;
    d_dt *= 1.15;
    if( d_dt + d_t > d_tend ) d_dt = d_tend - d_t + 1.e-30;
    zone.updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( d_dt )
    );
  }

  p_y = Libnucnet__Zone__getAbundances( p_zone );

  gsl_vector_sub( p_y, p_y_old );

  Libnucnet__Zone__updateAbundanceChanges( p_zone, p_y );

  zone.updateProperty(
    nnt::s_DTIME,
    boost::lexical_cast<std::string>( d_tend )
  );

  gsl_vector_free( p_y );
  gsl_vector_free( p_y_old );

}

//##############################################################################
// is_nonneg_abunds().
//##############################################################################

bool
is_nonneg_abunds( nnt::Zone& zone )
{

  double d_abund_min = 0.;

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      )
    );

  if( zone.hasProperty( nnt::s_LARGE_NEG_ABUND_THRESHOLD ) )
  {
    d_abund_min =
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_LARGE_NEG_ABUND_THRESHOLD )
      );
  }

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    double d_abund =
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      );

    if( d_abund < 0 )
    {
      if( fabs( d_abund ) < d_abund_min )
      {
        Libnucnet__Zone__updateSpeciesAbundance(
          zone.getNucnetZone(),
          species.getNucnetSpecies(),
          0.
        );
      }
      else
      {
        return false;
      }
    }

  }

  return true;

}

//##############################################################################
// evolve_zone().
//##############################################################################

int
evolve_zone( nnt::Zone& zone, double d_t_end )
{

  Libnucnet__Zone * p_zone;
  Libnuceq *p_equil;
  gsl_vector *p_abundances;
  double d_t = 0., d_dt, d_t9, d_rho; 
  int i_steps = 0;

  //============================================================================
  // Initializations.
  //============================================================================

  p_zone = zone.getNucnetZone();

  fprintf(
    stdout,
    "Start Zone %s with T9 = %s\n",
    Libnucnet__Zone__getLabel( p_zone, 1 ),
    Libnucnet__Zone__getProperty( p_zone, nnt::s_T9_0, NULL, NULL )
  );

  fflush( stdout );
  
  zone.updateProperty(
    nnt::s_TIME,
    boost::lexical_cast<std::string>( d_t )
  );

  d_t9 = boost::lexical_cast<double>( zone.getProperty( nnt::s_T9_0 ) );

  d_rho = boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO_0 ) );


  //============================================================================
  // Normalize abundances so that Xsum = 1.
  //============================================================================

  zone.normalizeAbundances( );

  //============================================================================
  // Set to NSE if t9 > 7.
  //============================================================================

  if( d_t9 > 7. )
  {

    d_t9 = 7.;

    d_rho =
      gsl_pow_3(
        d_t9 /
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9_0 ) )
      ) *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO_0 ) );

    d_t =
      boost::lexical_cast<double>( zone.getProperty( nnt::s_TAU ) )
      *
      log(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO_0 ) )
        / d_rho
      );
        
    p_equil =
      Libnuceq__new(
        Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( p_zone ) )
      );

    Libnuceq__setYe(
      p_equil, 
      Libnucnet__Zone__computeZMoment( p_zone, 1 )
    );

    Libnuceq__computeEquilibrium(
      p_equil, 
      d_t9,
      d_rho
    );

    p_abundances = Libnuceq__getAbundances( p_equil );

    Libnucnet__Zone__updateAbundances( p_zone, p_abundances );

    gsl_vector_free( p_abundances );

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );

    zone.updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( 1.e-6 )
    );


  }

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  while( d_t < d_t_end )
  {

    i_steps++;

    d_dt = boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );
    d_t += d_dt;

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );

    d_t9 =
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9_0 ) ) *
      exp(
        -d_t /
        (
          3. * boost::lexical_cast<double>( zone.getProperty( nnt::s_TAU ) )
        )
      );

    if( d_t9 < 1.e-6 ) d_t9 = 1.e-6;

    zone.updateProperty(
      nnt::s_T9,
      boost::lexical_cast<std::string>( d_t9 )
    );

    d_rho =
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO_0 ) ) *
      exp(
        -d_t /
        boost::lexical_cast<double>( zone.getProperty( nnt::s_TAU ) )
      );

    if( d_rho < 1.e-18 ) d_rho = 1.e-18;

    zone.updateProperty(
      nnt::s_RHO,
      boost::lexical_cast<std::string>( d_rho )
    );

  //============================================================================
  // Evolve system.
  //============================================================================

    if( !evolve( zone ) )
    {
      fprintf( stderr, "Problem converging.\n" );
      exit( EXIT_FAILURE );
    }

  //============================================================================
  // Update timestep.
  //============================================================================

    Libnucnet__Zone__updateTimeStep(
      p_zone,
      &d_dt,
      0.15,
      0.15,
      1.e-10
    );

    if ( d_t + d_dt > d_t_end ) {

      d_dt = d_t_end - d_t;

    }

    zone.updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( d_dt )
    );

    zone.updateProperty(
      nnt::s_STEPS,
      boost::lexical_cast<std::string>( i_steps )
    );

  //============================================================================
  // Limit the network.
  //============================================================================

    limit_evolution_network( zone );

  }  

  fprintf(
    stdout,
    "                          End Zone %s, %d steps, 1 - xsum = %g\n",
    Libnucnet__Zone__getLabel( p_zone, 1 ),
    i_steps,
    1. - Libnucnet__Zone__computeAMoment( p_zone, 1 )
  );

  return 0;

}

//##############################################################################
// evolve_nse_plus_weak_rates().
//##############################################################################

void
evolve_nse_plus_weak_rates( nnt::Zone zone )
{

  struct {
    nnt::Zone * pZone;
    double dYeOld;
  } work;

  int i_status;
  int i_iter = 0, i_max_iter = 100;
  const gsl_root_fsolver_type * p_T;
  gsl_root_fsolver * p_s;
  gsl_function F;
  double d_x_lo, d_x_hi;

  gsl_vector * p_old_abundances, * p_abundance_changes;

  p_old_abundances = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  work.pZone = &zone;
  work.dYeOld = Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );
     
  F.function = &nse_plus_weak_rates_function;
  F.params = &work;

  d_x_lo = work.dYeOld * 0.9;
  d_x_hi = work.dYeOld * 1.1;

  if(
      !nnt::bracket_root_of_function(
         F, &d_x_lo, &d_x_hi
      )
  )
  {
    fprintf( stderr, "Couldn't bracket root of function.\n" );
    exit( EXIT_FAILURE );
  }

  p_T = gsl_root_fsolver_brent;
  p_s = gsl_root_fsolver_alloc( p_T );
  gsl_root_fsolver_set( p_s, &F, d_x_lo, d_x_hi );
     
  do
  {
    i_iter++;
    i_status = gsl_root_fsolver_iterate( p_s );
    gsl_root_fsolver_root( p_s );
    d_x_lo = gsl_root_fsolver_x_lower( p_s );
    d_x_hi = gsl_root_fsolver_x_upper( p_s );
    i_status = gsl_root_test_interval( d_x_lo, d_x_hi, 0, 0.001 );
  } while( i_status == GSL_CONTINUE && i_iter < i_max_iter );
     
  gsl_root_fsolver_free( p_s );

  p_abundance_changes =
    Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  gsl_vector_sub( p_abundance_changes, p_old_abundances );

  Libnucnet__Zone__updateAbundanceChanges(
    zone.getNucnetZone(),
    p_abundance_changes
  );

  gsl_vector_free( p_old_abundances );
  gsl_vector_free( p_abundance_changes );
  
}

//##############################################################################
// nse_plus_weak_rates_function().
//##############################################################################

double
nse_plus_weak_rates_function(
  double d_ye,
  void * p_data
)
{

  typedef struct {
    nnt::Zone * pZone;
    double dYeOld;
  } work;

  boost::tuple<double, double, double, double, double> T;
  double d_result;

  work * p_work = ( work * ) p_data;

  p_work->pZone->updateProperty(
    nnt::s_YE,
    boost::lexical_cast<std::string>( d_ye )
  ); 

  nnt::set_zone_abundances_to_equilibrium( *p_work->pZone );

  T = compute_all_yedot( d_ye, p_work->pZone );

  d_result =
    ( d_ye - p_work->dYeOld ) /
      boost::lexical_cast<double>( p_work->pZone->getProperty( nnt::s_DTIME ) )
    -
    T.get<4>();

  return d_result;

}

//##############################################################################
// set_zone_for_evolution().
//##############################################################################

void
set_zone_for_evolution( nnt::Zone& zone )
{

  screening_data my_screening_data;
  coul_corr_data my_coul_corr_data;

  if(
    zone.hasProperty( nnt::s_USE_SCREENING ) &&
    zone.getProperty( nnt::s_USE_SCREENING )  == "yes"
  )
  {
  
    my_screening_data = get_screening_data( zone );

    Libnucnet__Zone__setScreeningFunction(
      zone.getNucnetZone(),
      (Libnucnet__Net__screening_function) my_screening_function,
      &my_screening_data
    );
    
    my_coul_corr_data = get_coulomb_corr_data();

    Libnucnet__Zone__setNseCorrectionFactorFunction(
      zone.getNucnetZone(),
      (Libnucnet__Species__nseCorrectionFactorFunction)
        my_coulomb_correction,
      &my_coul_corr_data
    );

  }

  //--------------------------------------------------------------------------
  // Update data for reactions.
  //--------------------------------------------------------------------------

  update_my_rate_functions_data( zone );

  //--------------------------------------------------------------------------
  // Compute rates.
  //--------------------------------------------------------------------------

  Libnucnet__Zone__computeRates(
    zone.getNucnetZone(),
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
  ); 

  //--------------------------------------------------------------------------
  // Set weak detailed balance.
  //--------------------------------------------------------------------------

  zone.setWeakDetailedBalance();

  //--------------------------------------------------------------------------
  // Modify rates.
  //--------------------------------------------------------------------------

  modify_rates( zone );

  //--------------------------------------------------------------------------
  // Zero out small rates.
  //--------------------------------------------------------------------------

  if( zone.hasProperty( nnt::s_SMALL_RATES_THRESHOLD ) )
  {
    zero_out_small_rates(
      zone,
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_SMALL_RATES_THRESHOLD )
      )
    );
  }

}

//##############################################################################
// get_evolution_matrix_and_vector().
//##############################################################################

std::pair< WnMatrix *, gsl_vector *>
get_evolution_matrix_and_vector( nnt::Zone& zone )
{

  set_zone_for_evolution( zone );

  //--------------------------------------------------------------------------
  // Return pair.
  //--------------------------------------------------------------------------

  return
    std::make_pair(
      Libnucnet__Zone__computeJacobianMatrix( zone.getNucnetZone() ),
      Libnucnet__Zone__computeFlowVector( zone.getNucnetZone() )
    );

}

//##############################################################################
// get_evolution_matrix().
//##############################################################################

WnMatrix *
get_evolution_matrix( nnt::Zone& zone )
{

  set_zone_for_evolution( zone );

  //--------------------------------------------------------------------------
  // Return matrix.
  //--------------------------------------------------------------------------

  return Libnucnet__Zone__computeJacobianMatrix( zone.getNucnetZone() );

}

//##############################################################################
// check_matrix_solution().
//##############################################################################

std::pair<double,double>
check_matrix_solution( nnt::Zone& zone, gsl_vector * p_sol )
{

  double d_checkT = 0, d_check = 0, d_abund, d_total = 0;
  size_t i_index;

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( zone.getNucnetZone() ) )
    );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    d_abund =
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      );

    i_index = Libnucnet__Species__getIndex( species.getNucnetSpecies() );

    if( zone.hasProperty( nnt::s_NEWTON_RAPHSON_ABUNDANCE ) )
    {
      if(
         d_abund >
         boost::lexical_cast<double>(
           zone.getProperty( nnt::s_NEWTON_RAPHSON_ABUNDANCE )
         )
      )
      {
        d_checkT = fabs( gsl_vector_get( p_sol, i_index ) / d_abund );
        if( d_checkT > d_check ) d_check = d_checkT;
      }
    }
    else if( d_abund > D_Y_MIN )
    {
      d_checkT = fabs( gsl_vector_get( p_sol, i_index ) / d_abund );
      if( d_checkT > d_check ) d_check = d_checkT;
    }

    d_total +=
      gsl_pow_2(
        gsl_vector_get( p_sol, i_index ) *
        Libnucnet__Species__getA( species.getNucnetSpecies() )
      );
    
  }

  return std::make_pair( d_check, sqrt( d_total ) );

}
       
//##############################################################################
// my_log10_t9_from_entropy_root(). 
//##############################################################################

double
my_log10_t9_from_entropy_root(
  double d_x,
  nnt::Zone * p_nnt_zone
)
{

  double d_result;

  gsl_vector * p_abundances, * p_abundance_changes;

  p_nnt_zone->updateProperty(
    nnt::s_T9,
    boost::lexical_cast<std::string>( pow( 10., d_x ) )
  );

  p_abundances =
    Libnucnet__Zone__getAbundances(
      p_nnt_zone->getNucnetZone()
    );

  p_abundance_changes =
    Libnucnet__Zone__getAbundanceChanges(
      p_nnt_zone->getNucnetZone()
    );

  safe_evolve(
    *p_nnt_zone,
    boost::lexical_cast<double>( p_nnt_zone->getProperty( nnt::s_DTIME ) ),
    boost::lexical_cast<double>( p_nnt_zone->getProperty( nnt::s_DTIME ) )
 );

  d_result =
    compute_thermo_quantity(
      *p_nnt_zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    )
    -      
    boost::lexical_cast<double>(
      p_nnt_zone->getProperty( nnt::s_ENTROPY_PER_NUCLEON )
    );

  Libnucnet__Zone__updateAbundances(
    p_nnt_zone->getNucnetZone(),
    p_abundances
  );

  Libnucnet__Zone__updateAbundanceChanges(
    p_nnt_zone->getNucnetZone(),
    p_abundance_changes
  );

  gsl_vector_free( p_abundances );
  gsl_vector_free( p_abundance_changes );

   std::cout << pow( 10., d_x ) << "  " << d_result << std::endl;

  return d_result;

}
    
//##############################################################################
// my_density_from_entropy_root(). 
//##############################################################################

double
my_density_from_entropy_root(
  double d_x,
  nnt::Zone * p_nnt_zone
)
{

  p_nnt_zone->updateProperty(
    nnt::s_RHO,
    boost::lexical_cast<std::string>( d_x )
  );

  return
    compute_thermo_quantity(
      *p_nnt_zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    )
    -      
    boost::lexical_cast<double>(
      p_nnt_zone->getProperty( nnt::s_ENTROPY_PER_NUCLEON )
    );

}

//##############################################################################
// my_log10_t9_guess().
//##############################################################################

std::pair<double,double>
my_log10_t9_guess(
  nnt::Zone * p_nnt_zone
)
{

  return
    std::make_pair<double,double>(
      log10( 
	0.999 *
	  boost::lexical_cast<double>(
	    p_nnt_zone->getProperty( nnt::s_T9 )
	  )
      ),
      log10(
	1.001 *
	  boost::lexical_cast<double>(
	    p_nnt_zone->getProperty( nnt::s_T9 )
	  )
      )
    );

}

} // namespace user
