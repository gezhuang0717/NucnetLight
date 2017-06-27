#include <boost/function.hpp>
#include <boost/bind.hpp>

namespace user
{

//##############################################################################
// exp_multi_zone().
//##############################################################################

template<class Function1, class Function2>
int
exp_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 phi_fn,
  Function2 const_vec_fn,
  boost::unordered_map<std::string,std::string>& param_map,
  double d_t
)
{

  WnSparseSolve__Phi * p_phi;
  boost::unordered_map<std::string,std::string>::iterator it;
  gsl_vector * p_sol, * p_old, * p_constant;
  std::vector<WnMatrix *> mix_matrices;
  WnMatrix * p_jacobian_matrix;
  size_t i, i_species, i_offset, i_delta;

  //============================================================================
  // Create solver and set parameters.
  //============================================================================

  p_phi = WnSparseSolve__Phi__new();

  it = param_map.find( nnt::s_ITER_SOLVER_DEBUG );

  if( it != param_map.end() )
    WnSparseSolve__Phi__setDebug( p_phi ); 

  it = param_map.find( nnt::s_ITER_SOLVER_REL_TOL );

  if( it != param_map.end() )
    WnSparseSolve__Phi__updateTolerance(
      p_phi,
      boost::lexical_cast<double>( it->second )
    );

  it = param_map.find( nnt::s_ITER_SOLVER_MAX_ITERATIONS );

  if( it != param_map.end() )
    WnSparseSolve__Phi__updateMaximumIterations(
      p_phi,
      boost::lexical_cast<int>( it->second )
    );

  //============================================================================
  // Get number of species and abundance vectors.
  //============================================================================

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  i_offset = i_species;

  if( multi_mass_calculation_check( zones ) )
    i_delta = 1;
  else
    i_delta = 0;

  i_offset += i_delta;

  p_old = create_full_vector( zones );

  try{
    p_constant = const_vec_fn();
  }
  catch( boost::bad_function_call& ex )
  {
    p_constant = gsl_vector_calloc( i_offset * zones.size() );
  }

  //============================================================================
  // Exponential solution.
  //============================================================================

  mix_matrices = phi_fn();

  p_jacobian_matrix = get_evolution_matrix( zones[0] );
  WnMatrix__scaleMatrix( p_jacobian_matrix, -1. );

  #pragma omp parallel for schedule( dynamic, 1 )
    for( i = 0; i < mix_matrices.size(); i++ )
    {
      WnMatrix__insertMatrix(
        mix_matrices[i],
        p_jacobian_matrix,
        (i * i_offset) + i_delta + 1,
        (i * i_offset ) + i_delta + 1
      );

    }

  WnMatrix__free( p_jacobian_matrix );

  //============================================================================
  // Exponential solution.
  //============================================================================

  p_sol =
    phi__solve__parallel(
      p_phi,
      mix_matrices,
      p_old,
      p_constant,
      d_t
    );

  for( size_t i = 0; i < mix_matrices.size(); i++ )
    WnMatrix__free( mix_matrices[i] );

  WnSparseSolve__Phi__free( p_phi );
  gsl_vector_free( p_constant );

  if( !p_sol ) 
  {
    update_from_full_vector(
      zones,
      p_old,
      "abundances"
    );
    gsl_vector_free( p_old );
    return 0;
  }

  update_from_full_vector(
    zones,
    p_sol,
    "abundances"
  );

  gsl_vector_free( p_sol );
  gsl_vector_free( p_old );

  return 1;

}

template<class Function>
int
exp_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function phi_fn,
  boost::unordered_map<std::string,std::string>& param_map,
  double d_t
)
{

  boost::function< gsl_vector*( ) > f;
  f = 0;

  return
    exp_multi_zone(
      zones,
      phi_fn,
      boost::bind( f ),
      param_map,
      d_t
    );

}

//##############################################################################
// safe_exp_multi_zone().
//##############################################################################

template<class Function1, class Function2>
int
safe_exp_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 phi_fn,
  Function2 const_vec_fn,
  boost::unordered_map<std::string,std::string>& param_map,
  double d_dt
)
{

  double d_dt_e, d_dt_cum = 0.;

  gsl_vector * p_old, * p_new;

  p_old = create_full_vector( zones );

  if( !exp_multi_zone( zones, phi_fn, const_vec_fn, param_map, d_dt ) )
  {

    d_dt_e = d_dt / D_DT_E_STEP;

    while( !exp_multi_zone( zones, phi_fn, const_vec_fn, param_map, d_dt_e ) )
    {
      d_dt_e /= D_DT_E_STEP;
    }

    d_dt_cum = d_dt_e;

    while( gsl_fcmp( d_dt_cum, d_dt, D_DT_E_CMP ) )
    { 

      std::cout << d_dt_cum << "  " << d_dt << std::endl;

      if( exp_multi_zone( zones, phi_fn, const_vec_fn, param_map, d_dt_e ) )
      {
	d_dt_cum += d_dt_e;
	d_dt_e *= D_DT_E_REG;
	if( d_dt_cum + d_dt_e > d_dt ) d_dt_e = d_dt - d_dt_cum;
      }
      else
      {
	d_dt_e /= D_DT_E_REG;
        if( d_dt_e < 1.e-100 ) return 0;
      }

    }

  }

  p_new = create_full_vector( zones );

  gsl_vector_sub( p_new, p_old );

  update_from_full_vector(
    zones,
    p_new,
    "abundance changes"
  );

  gsl_vector_free( p_old );
  gsl_vector_free( p_new );

  return 1;

}

template<class Function>
int
safe_exp_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function phi_fn,
  boost::unordered_map<std::string,std::string>& param_map,
  double d_t
)
{

  boost::function< gsl_vector*( ) > f;
  f = 0;

  return
    safe_exp_multi_zone(
      zones,
      phi_fn,
      boost::bind( f ),
      param_map,
      d_t
    );

}

//##############################################################################
// evolve_multi_zone().
//##############################################################################

template<class Function1, class Function2>
int
evolve_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 fn1,
  Function2 fn2,
  boost::unordered_map<std::string,std::string>& param_map,
  double d_dt
)
{

  WnMatrix * p_matrix;
  gsl_vector * p_vector, * p_sol, * p_old, * p_current, * p_work;
  gsl_vector_view view;
  std::vector<WnMatrix *> jacobians;
  std::vector<gsl_vector *> rhs_vectors;
  size_t i, i_species, i_iter, i_offset, i_delta;
  int i_fn2;
  double d_check;

  //============================================================================
  // Get number of species and abundance vectors.
  //============================================================================

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  i_offset = i_species;

  if( multi_mass_calculation_check( zones ) )
    i_delta = 1;
  else
    i_delta = 0;

  i_offset += i_delta;

  p_old = create_full_vector( zones );

  p_current = gsl_vector_alloc( WnMatrix__get_gsl_vector_size( p_old ) );

  gsl_vector_memcpy( p_current, p_old );

  //============================================================================
  // Newton-Raphson loop.
  //============================================================================

  for( i_iter = 0; i_iter < IT_MAX; i_iter++ )
  {

    boost::tie( p_matrix, p_vector ) = fn1( p_current );

    boost::tie( jacobians, rhs_vectors ) =
      get_zone_jacobians_and_rhs_vectors( zones );

    i = 0;
    BOOST_FOREACH( WnMatrix * p_sub_matrix, jacobians )
    {
      
      WnMatrix__insertMatrix(
	p_matrix,
	p_sub_matrix,
	(i * i_offset) + i_delta + 1,
	(i * i_offset ) + i_delta + 1
      );

      WnMatrix__free( p_sub_matrix );

      i++;

    }

    WnMatrix__addValueToDiagonals( p_matrix, 1. / d_dt );

    i = 0;
    BOOST_FOREACH( gsl_vector * p_sub_vector, rhs_vectors )
    {
      
      view =
        gsl_vector_subvector( p_vector, ( i * i_offset ) + i_delta, i_species );
      gsl_vector_add( &view.vector, p_sub_vector );
      gsl_vector_free( p_sub_vector );
      i++;

    }

    p_work = gsl_vector_alloc( WnMatrix__get_gsl_vector_size( p_current ) );
    gsl_vector_memcpy( p_work, p_current );
    gsl_vector_sub( p_work, p_old );
    gsl_vector_scale( p_work, 1. / d_dt );
    gsl_vector_sub( p_vector, p_work );
    gsl_vector_free( p_work );

    p_sol =
      solve_sparse_matrix_with_ilu_preconditioner(
        p_matrix,
        p_vector,
        param_map
      ); 

    WnMatrix__free( p_matrix );
    gsl_vector_free( p_vector );

    if( !p_sol ) 
    {
      update_from_full_vector(
        zones,
        p_old,
        "abundances"
      );
      gsl_vector_free( p_old );
      gsl_vector_free( p_current );
      return 0;
    }

    gsl_vector_add( p_current, p_sol );

    update_from_full_vector(
      zones,
      p_current,
      "abundances"
    );

    d_check =
      sqrt( gsl_blas_dnrm2( p_sol ) ) / sqrt( gsl_blas_dnrm2( p_current ) );

    gsl_vector_free( p_sol );

    if( d_check < 1.e-4 ) break;

  }

  try
  {
    i_fn2 = fn2();
  }
  catch( boost::bad_function_call& ex )
  {
    i_fn2 = 1;
  }

  if( i_fn2 == 0 )
  {
    update_from_full_vector(
      zones,
      p_old,
      "abundances"
    );
    gsl_vector_free( p_old );
    gsl_vector_free( p_current );
    return 0;
  }

  gsl_vector_free( p_old );
  gsl_vector_free( p_current );

  return 1;

}

template<class Function1>
int
evolve_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 fn,
  boost::unordered_map<std::string,std::string>& param_map,
  double d_dt
)
{

  boost::function< int() > f;
  f = 0;

  return
    evolve_multi_zone(
      zones,
      fn,
      boost::bind( f ),
      param_map,
      d_dt
    );

}

//##############################################################################
// safe_evolve_multi_zone().
//##############################################################################

template<class Function1, class Function2>
void
safe_evolve_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 fn1,
  Function2 fn2,
  boost::unordered_map<std::string,std::string>& param_map,
  double d_dt
)
{

  double d_dt_e, d_dt_cum = 0.;

  gsl_vector * p_old, * p_new;

  p_old = create_full_vector( zones );

  if( !evolve_multi_zone( zones, fn1, fn2, param_map, d_dt ) )
  {

    d_dt_e = d_dt / D_DT_E_STEP;

    while( !evolve_multi_zone( zones, fn1, fn2, param_map, d_dt_e ) )
    {
      d_dt_e /= D_DT_E_STEP;
    }

    d_dt_cum = d_dt_e;

    while( gsl_fcmp( d_dt_cum, d_dt, D_DT_E_CMP ) )
    { 

      std::cout << d_dt_cum << "  " << d_dt << std::endl;

      if( evolve_multi_zone( zones, fn1, param_map, d_dt_e ) )
      {
	d_dt_cum += d_dt_e;
	d_dt_e *= D_DT_E_REG;
	if( d_dt_cum + d_dt_e > d_dt ) d_dt_e = d_dt - d_dt_cum;
      }
      else
      {
	d_dt_e /= D_DT_E_REG;
      }

    }

  }

  p_new = create_full_vector( zones );

  gsl_vector_sub( p_new, p_old );

  update_from_full_vector(
    zones,
    p_new,
    "abundance changes"
  );

  gsl_vector_free( p_old );
  gsl_vector_free( p_new );

}


} // namespace user
