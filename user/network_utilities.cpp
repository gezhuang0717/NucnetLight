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
//! \brief Code for useful network routines.
////////////////////////////////////////////////////////////////////////////////

#include "user/network_utilities.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

//##############################################################################
// set_zone_post_shock_conditions()
//##############################################################################

int
set_zone_post_shock_conditions( nnt::Zone zone, std::string s_mach )
{

  double d_mach, d_t9, d_rho;

  d_mach = boost::lexical_cast<double>( s_mach );

  d_rho = boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO1 ) );

  d_rho *=
    ( D_GAMMA + 1. ) * gsl_pow_2( d_mach ) /
    (
      ( D_GAMMA - 1. ) * gsl_pow_2( d_mach ) + 2.
    );

  zone.updateProperty(
    nnt::s_RHO_0,
    boost::lexical_cast<std::string>( d_rho )
  );
  
  d_t9 =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T1 ) ) /
    GSL_CONST_NUM_GIGA;

  d_t9 *=
    (
      2 * D_GAMMA * gsl_pow_2( d_mach ) - ( D_GAMMA - 1. )
    ) *
    (
      ( D_GAMMA - 1. ) * gsl_pow_2( d_mach ) + 2.
    ) /
    gsl_pow_2(
      ( D_GAMMA + 1. ) * d_mach
    );
     
  zone.updateProperty(
    nnt::s_T9_0,
    boost::lexical_cast<std::string>( d_t9 )
  );
  
  zone.updateProperty(
    nnt::s_T9,
    boost::lexical_cast<std::string>( d_t9 )
  );
  
  zone.updateProperty(
    nnt::s_MACH,
    s_mach
  );
  
  zone.updateProperty(
    nnt::s_TAU,
    boost::lexical_cast<std::string>( 446. / sqrt( d_rho ) )
  );
  
  return 1;
  
}

//##############################################################################
// compute_luminosity().
//##############################################################################

double
compute_luminosity(
  nnt::Zone zone
)
{
  double d_result;

  d_result =  
    4. *
    M_PI * 
    GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT *
    gsl_pow_2( 
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RADIUS ) ) 
    ) * 
    gsl_pow_4( 
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
      GSL_CONST_NUM_GIGA
    );

  return d_result;

}

//##############################################################################
// compute_photon_entropy_loss_rate().
//##############################################################################

double
compute_photon_entropy_loss_rate(
  nnt::Zone zone
)
{

  double d_r;
  double d_opacity_factor;

  d_r = 
    1. /
    (
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
      Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ) *
      GSL_CONST_NUM_AVOGADRO
    ) /
    GSL_CONST_CGSM_THOMSON_CROSS_SECTION;

  d_opacity_factor =  
    1. -
    gsl_pow_3( 
      1. -
      d_r / 
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RADIUS ) )
    );

  return
    compute_luminosity( zone ) /
    (
      4. / 3. *
      M_PI *
      gsl_pow_3(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RADIUS ) ) 
      ) *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
      GSL_CONST_NUM_AVOGADRO
    ) /
    ( 
      GSL_CONST_CGSM_BOLTZMANN *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
      GSL_CONST_NUM_GIGA
    ) *
    d_opacity_factor;

} 

//##############################################################################
// compute_v_Thermal().
//##############################################################################

double
compute_v_Thermal(
  Libnucnet__Species * p_species,
  double d_t9
)
{

  if( !p_species )
  {
    std::cerr << "No such species." << std::endl;
    exit( EXIT_FAILURE );
  }

  return
    sqrt(
      GSL_CONST_CGSM_BOLTZMANN * d_t9 *GSL_CONST_NUM_GIGA /
      ( 
        GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS *
          Libnucnet__Species__getA( p_species )
        +
        Libnucnet__Species__getMassExcess( p_species ) *
          GSL_CONST_CGSM_ELECTRON_VOLT *
          GSL_CONST_NUM_MEGA /
          gsl_pow_2( GSL_CONST_CGSM_SPEED_OF_LIGHT )
      )
    );

}

//##############################################################################
// update_exposures().
//##############################################################################

void
update_exposures( nnt::Zone& zone )
{

  Libnucnet__Zone__iterateOptionalProperties(
    zone.getNucnetZone(),
      nnt::s_EXPOSURE,
      NULL,
      NULL,
      (Libnucnet__Zone__optional_property_iterate_function)
         update_exposures_callback,
      &zone
  );

}
    
//##############################################################################
// update_exposures_callback().
//##############################################################################

void
update_exposures_callback(
  const char * s_property,
  const char * s_species,
  const char * s_tag2,
  const char * s_value,
  nnt::Zone& zone
)
{

  Libnucnet__Species * p_species =
    Libnucnet__Nuc__getSpeciesByName(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      s_species
    );

  if( !p_species ) return;

  zone.updateProperty(
    nnt::s_EXPOSURE,
    std::string( s_species ),
    boost::lexical_cast<std::string>(
      boost::lexical_cast<double>( s_value )
      +
      compute_v_Thermal(
        p_species,
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) )
      ) *
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        p_species
      ) *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
      GSL_CONST_NUM_AVOGADRO *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) ) *
      GSL_CONST_CGSM_BARN *
      GSL_CONST_NUM_MILLI
    )
  );
}

//##############################################################################
// compute_Ysum().
//##############################################################################

double
compute_Ysum(
  nnt::Zone& zone,
  nnt::species_list_t& species_list
)
{

  double d_Ysum = 0.;

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    d_Ysum +=
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      );

  }

  return d_Ysum;

}

//##############################################################################
// update_t9_rho_in_zone_by_interpolation().
//##############################################################################

void
update_t9_rho_in_zone_by_interpolation(
  nnt::Zone& zone,
  std::string s_interpolation_type,
  const std::vector<double> &time,
  const std::vector<double> &t9,
  const std::vector<double> &log10_rho
)
{

  gsl_vector * p_time, * p_t9, * p_log10_rho;
  double d_t9 = 0, d_rho = 0, d_change_t9, d_change_rho;
  double d_time, d_dt;

  //==========================================================================
  // Get gsl vectors.
  //==========================================================================

  p_time = nnt::get_new_gsl_vector_from_std_vector( time );
  p_t9 = nnt::get_new_gsl_vector_from_std_vector( t9 );
  p_log10_rho = nnt::get_new_gsl_vector_from_std_vector( log10_rho );

  //==========================================================================
  // Get time and dt.
  //==========================================================================

  d_time =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );
  d_dt =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  //==========================================================================
  // Interpolate.  If change in t9 or rho too large, decrease time step by
  // factor of two and interpolate again.
  //==========================================================================

  while( d_dt > 0. )
  {

     if( s_interpolation_type == "spline" )
     {
       d_t9 = nnt::spline_interpolation( p_time, p_t9, d_time );
       d_rho =
         pow( 10., nnt::spline_interpolation( p_time, p_log10_rho, d_time ) );
     }
     else if( s_interpolation_type == "linear" )
     {
       d_t9 = nnt::linear_interpolation( p_time, p_t9, d_time );
       d_rho =
         pow( 10., nnt::linear_interpolation( p_time, p_log10_rho, d_time ) );
     }
     else
     {
       std::cerr << "Invalid interpolation type." << std::endl;
       exit( EXIT_FAILURE );
     }

     if( !zone.hasProperty( "interpolation flag" ) )
     {
       zone.updateProperty(
         nnt::s_OLD_T9, 
         boost::lexical_cast<std::string>( d_t9 )
       );
       zone.updateProperty(
         nnt::s_OLD_RHO, 
         boost::lexical_cast<std::string>( d_rho )
       );
       d_change_t9 = 0.;
       d_change_rho = 0.;
       zone.updateProperty( "interpolation flag", "set" );
     }
     else
     {

       d_change_t9 =
         fabs(
           d_t9 -
           boost::lexical_cast<double>( zone.getProperty( nnt::s_OLD_T9 ) )
         ) /
         boost::lexical_cast<double>( zone.getProperty( nnt::s_OLD_T9 ) );

       d_change_rho =
         fabs(
           d_rho -
           boost::lexical_cast<double>( zone.getProperty( nnt::s_OLD_RHO ) )
         ) /
         boost::lexical_cast<double>( zone.getProperty( nnt::s_OLD_RHO ) );
     }

     if( d_change_t9 < D_EPS && d_change_rho < D_EPS )
       break;

     d_time -= d_dt;
     d_dt /= 2;
     d_time += d_dt;

     if( d_dt < D_DT_MIN )
       break;
        
  }
  
  zone.updateProperty(
    nnt::s_T9,
    boost::lexical_cast<std::string>( d_t9 )
  );

  zone.updateProperty(
    nnt::s_RHO,
    boost::lexical_cast<std::string>( d_rho )
  );

  zone.updateProperty(
    nnt::s_TIME,
    boost::lexical_cast<std::string>( d_time )
  );

  zone.updateProperty(
    nnt::s_DTIME,
    boost::lexical_cast<std::string>( d_dt )
  );

  zone.updateProperty(
    nnt::s_OLD_T9,
    boost::lexical_cast<std::string>( d_t9 )
  );

  zone.updateProperty(
    nnt::s_OLD_RHO,
    boost::lexical_cast<std::string>( d_rho )
  );

  gsl_vector_free( p_time );
  gsl_vector_free( p_t9 );
  gsl_vector_free( p_log10_rho );

}

//############################################################################
// copy_zone_abundances_as_properties().
//##########################################################################//

void
copy_zone_abundances_as_properties(
  nnt::Zone& zone1,
  nnt::Zone& zone2,
  const char * s_name
)
{

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zone1.getNucnetZone() )
      )
    );

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    zone2.updateProperty(
      s_name,
      Libnucnet__Species__getName( species.getNucnetSpecies() ),
      boost::lexical_cast<std::string>(
        Libnucnet__Zone__getSpeciesAbundance(
          zone1.getNucnetZone(),
          species.getNucnetSpecies()
        )
      )
    );
  }

}

//##############################################################################
// compute_cluster_abundance_moment().
//##############################################################################

double
compute_cluster_abundance_moment(
  nnt::Zone& zone,
  const char * s_nuc_xpath,
  std::string s_nuc,
  int i_moment
)
{

  Libnucnet__NetView * p_view;
  double d_yc = 0;

  assert( s_nuc == "z" || s_nuc == "n" || s_nuc == "a" );

  p_view = zone.getNetView( s_nuc_xpath, "" );

  nnt::species_list_t species_list =
    nnt::make_species_list( 
      Libnucnet__Net__getNuc(
        Libnucnet__NetView__getNet( p_view )
      )
    );

  if( s_nuc == "z" )
  {
    BOOST_FOREACH( nnt::Species species, species_list )
    {
      d_yc +=
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          species.getNucnetSpecies()
        ) *
        pow(
          (double) Libnucnet__Species__getZ( species.getNucnetSpecies() ),
          (double) i_moment
        );   
    }
  }
  else if( s_nuc == "n" )
  {
    BOOST_FOREACH( nnt::Species species, species_list )
    {
      d_yc +=
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          species.getNucnetSpecies()
        ) *
        pow(
          (double)
            Libnucnet__Species__getA( species.getNucnetSpecies() ) -
            Libnucnet__Species__getZ( species.getNucnetSpecies() ),
          (double) i_moment
        );   
    }
  }
  else
  {
    BOOST_FOREACH( nnt::Species species, species_list )
    {
      d_yc +=
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          species.getNucnetSpecies()
        ) *
        pow(
          (double) Libnucnet__Species__getA( species.getNucnetSpecies() ),
          (double) i_moment
        );   
    }
  }

  return d_yc;

}

//##############################################################################
// set_specific_species().
//##############################################################################

void
set_specific_species(
  WnMatrix * p_matrix,
  gsl_vector * p_rhs,
  void * p_data
)
{

  nnt::Zone zone = *(nnt::Zone *) p_data;

  double d_dt =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  Libnucnet__Species * p_fixed_species =
    Libnucnet__Nuc__getSpeciesByName(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      zone.getProperty( nnt::s_SPECIFIC_SPECIES ).c_str()
    );

  size_t i_index = Libnucnet__Species__getIndex( p_fixed_species );

  WnMatrix__removeRow( p_matrix, i_index + 1 );

  gsl_vector * p_row =
    gsl_vector_calloc( WnMatrix__getNumberOfColumns( p_matrix ) );

  gsl_vector_set(
    p_row,
    i_index,
    1. / d_dt
  );

  WnMatrix__insertRow( p_matrix, i_index + 1, p_row );

  gsl_vector_free( p_row );

  gsl_vector_set(
    p_rhs,
    i_index,
    boost::lexical_cast<double>(
      zone.getProperty( nnt::s_SPECIFIC_ABUNDANCE )
    ) / d_dt -
    Libnucnet__Zone__getSpeciesAbundance(
      zone.getNucnetZone(),
      p_fixed_species
    ) / d_dt
  );

}

} // namespace user
