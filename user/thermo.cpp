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
//! \brief Code for computing relevant thermodynamics quantities.
//!
////////////////////////////////////////////////////////////////////////////////

#include "user/thermo.h"

namespace user
{

//##############################################################################
// get_thermo_species_list().
//##############################################################################

nnt::species_list_t
get_thermo_species_list( nnt::Zone& zone )
{

  if( zone.hasProperty( nnt::s_THERMO_NUC_VIEW ) )
  {
    return
      nnt::make_species_list(
        Libnucnet__Net__getNuc(
          Libnucnet__NetView__getNet(
            zone.getNetView(
              zone.getProperty( nnt::s_THERMO_NUC_VIEW ).c_str()
            )  
          )
        )
      );
  }
  else
  {
    return
      nnt::make_species_list(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      );
  }

}

//############################################################################
// compute_electron_chemical_potential_kT().
//############################################################################

/**
  \brief Compute the electron chemical potential (less the rest mass)
	 divided by kT.

  \param zone A NucNet Tools zonel
  \return A double giving the chemical potential / kT.
*/

double
compute_electron_chemical_potential_kT(
  nnt::Zone& zone
)
{

  Libstatmech__Fermion *p_electron;
  double d_result;

  double d_t9 = boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) );

  double d_rho = boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );

  double d_ye = Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );

  p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  d_result =
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      d_t9 * GSL_CONST_NUM_GIGA,
      d_rho * d_ye * GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return d_result;

}

//##############################################################################
// compute_baryon_internal_energy_density(). 
//##############################################################################

double
compute_baryon_internal_energy_density( nnt::Zone& zone )
{

  double d_eB = 0;
  double d_n, d_T;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  d_n =
    GSL_CONST_NUM_AVOGADRO *
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    d_eB +=
      ( 3. / 2. ) *
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      ) *
      d_n *
      GSL_CONST_CGSM_BOLTZMANN *
      d_T;
  }

  return d_eB;

}

//##############################################################################
// compute_electron_internal_energy_density().
//##############################################################################

double
compute_electron_internal_energy_density( nnt::Zone& zone )
{

  Libstatmech__Fermion * p_electron;
  double d_T, d_eE = 0;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    nnt::s_INTERNAL_ENERGY_DENSITY,
    0.,
    nnt::d_REL_EPS
  );

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;

  d_eE =
    Libstatmech__Fermion__computeQuantity(
      p_electron,
      nnt::s_INTERNAL_ENERGY_DENSITY,
      d_T,
      compute_electron_chemical_potential_kT( zone ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return d_eE;

}
  
//##############################################################################
// compute_photon_internal_energy_density().
//##############################################################################

double
compute_photon_internal_energy_density( nnt::Zone& zone )
{

  double d_eP = 0, d_T;
  Libstatmech__Boson * p_photon;

  p_photon = Libstatmech__Boson__new( "photon", 0., 2, 0. );

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;
  
  d_eP =
    Libstatmech__Boson__computeQuantity(
      p_photon,
      nnt::s_INTERNAL_ENERGY_DENSITY,
      d_T,
      0.,
      NULL,
      NULL
    );

  Libstatmech__Boson__free( p_photon );

  return d_eP;

}

//##############################################################################
// compute_internal_energy_density().
//##############################################################################

double
compute_internal_energy_density( nnt::Zone& zone )
{

  return
    compute_baryon_internal_energy_density( zone ) +
    compute_electron_internal_energy_density( zone ) +
    compute_photon_internal_energy_density( zone );

}
  
//##############################################################################
// compute_rest_mass_energy_density().
//##############################################################################

/**
 * \brief Compute the total rest mass energy density in a zone.
 * \return The rest mass energy density (cgs).
 */


double
compute_rest_mass_energy_density( nnt::Zone& zone )
{

  double d_n, d_result = 0;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  d_n =
    GSL_CONST_NUM_AVOGADRO *
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    d_result +=
      d_n *
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      ) *
      (
        Libnucnet__Species__getA( species.getNucnetSpecies() ) *
        GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS *
        gsl_pow_2( GSL_CONST_CGSM_SPEED_OF_LIGHT ) +
        Libnucnet__Species__getMassExcess( species.getNucnetSpecies() ) *
        GSL_CONST_CGSM_ELECTRON_VOLT *
        GSL_CONST_NUM_MEGA
      );
       
  }

  return d_result;

}

//##############################################################################
// compute_baryon_pressure().
//##############################################################################

double
compute_baryon_pressure( nnt::Zone& zone )
{

  double d_pB = 0;
  double d_n, d_T;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  d_n =
    GSL_CONST_NUM_AVOGADRO *
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    d_pB +=
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      ) *
      d_n *
      GSL_CONST_CGSM_BOLTZMANN *
      d_T;
  }

  return d_pB;

}

//##############################################################################
// compute_electron_pressure().
//##############################################################################

double
compute_electron_pressure( nnt::Zone& zone )
{

  Libstatmech__Fermion * p_electron;
  double d_T, d_pE = 0;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    nnt::s_PRESSURE,
    0.,
    nnt::d_REL_EPS
  );

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;

  d_pE =
    Libstatmech__Fermion__computeQuantity(
      p_electron,
      nnt::s_PRESSURE,
      d_T,
      compute_electron_chemical_potential_kT( zone ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return d_pE;

}
  
//##############################################################################
// compute_photon_pressure().
//##############################################################################

double
compute_photon_pressure( nnt::Zone& zone )
{

  double d_pP = 0, d_T;
  Libstatmech__Boson * p_photon;

  p_photon = Libstatmech__Boson__new( "photon", 0., 2, 0. );

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;
  
  d_pP =
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_PRESSURE,
      d_T,
      0.,
      NULL,
      NULL
    );

  Libstatmech__Boson__free( p_photon );

  return d_pP;

}

//##############################################################################
// compute_pressure().
//##############################################################################

double
compute_pressure( nnt::Zone& zone )
{

  return
    compute_baryon_pressure( zone ) +
    compute_electron_pressure( zone ) +
    compute_photon_pressure( zone );

}
  
//##############################################################################
// compute_baryon_dPdT().
//##############################################################################

double
compute_baryon_dPdT( nnt::Zone& zone )
{

  double d_B_dPdT = 0;
  double d_n;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  d_n =
    GSL_CONST_NUM_AVOGADRO *
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    d_B_dPdT +=
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      ) *
      d_n *
      GSL_CONST_CGSM_BOLTZMANN;
  }

  return d_B_dPdT;

}

//##############################################################################
// compute_electron_dPdT().
//##############################################################################

double
compute_electron_dPdT( nnt::Zone& zone )
{

  Libstatmech__Fermion * p_electron;
  double d_E_dPdT = 0;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    S_PRESSURE,
    0.,
    nnt::d_REL_EPS
  );

  d_E_dPdT =
    Libstatmech__Fermion__computeTemperatureDerivative(
      p_electron,
      S_PRESSURE,
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
        GSL_CONST_NUM_GIGA,
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
        GSL_CONST_NUM_AVOGADRO *
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return d_E_dPdT;

}
  
//##############################################################################
// compute_photon_dPdT().
//##############################################################################

double
compute_photon_dPdT( nnt::Zone& zone )
{

  double d_T;

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;
  
  return
    4. / 3. *
    (
      4. *
      GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT /
      GSL_CONST_CGSM_SPEED_OF_LIGHT
    ) *
    gsl_pow_3( d_T );

}

//##############################################################################
// compute_dPdT().
//##############################################################################

double
compute_dPdT( nnt::Zone& zone )
{

  return
    compute_baryon_dPdT( zone ) +
    compute_electron_dPdT( zone ) +
    compute_photon_dPdT( zone );

}
  
//##############################################################################
// compute_baryon_entropy_per_nucleon().
//##############################################################################

double
compute_baryon_entropy_per_nucleon( nnt::Zone& zone )
{

  double d_sB = 0;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    if( 
	Libnucnet__Zone__getSpeciesAbundance(
	  zone.getNucnetZone(),
	  species.getNucnetSpecies()
	) > 1.e-30
    )
    {

      d_sB +=
	(
	  2.5 *
	  Libnucnet__Zone__getSpeciesAbundance(
	    zone.getNucnetZone(),
	    species.getNucnetSpecies()
	  )
	)
	-
	(
	  Libnucnet__Zone__getSpeciesAbundance(
	    zone.getNucnetZone(),
	    species.getNucnetSpecies()
	  )
	  *
	  log(
	    Libnucnet__Zone__getSpeciesAbundance(
	      zone.getNucnetZone(),
	      species.getNucnetSpecies()
	    ) /
	    Libnucnet__Species__computeQuantumAbundance(
	      species.getNucnetSpecies(),
	      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
	      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
	    )
	  )
	);

    }

  }

  return d_sB;

}

//##############################################################################
// compute_electron_entropy_per_nucleon().
//##############################################################################

double
compute_electron_entropy_per_nucleon( nnt::Zone& zone )
{

  Libstatmech__Fermion * p_electron;
  double d_T, d_sE = 0;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    S_ENTROPY_DENSITY,
    0.,
    nnt::d_REL_EPS
  );

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;

  d_sE =
    Libstatmech__Fermion__computeQuantity(
      p_electron,
      S_ENTROPY_DENSITY,
      d_T,
      compute_electron_chemical_potential_kT( zone ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return
    d_sE /
    (
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
      GSL_CONST_NUM_AVOGADRO *
      GSL_CONST_CGSM_BOLTZMANN
    );

}
  
//##############################################################################
// compute_photon_entropy_per_nucleon().
//##############################################################################

double
compute_photon_entropy_per_nucleon( nnt::Zone& zone )
{

  double d_sP = 0, d_T;
  Libstatmech__Boson * p_photon;

  p_photon = Libstatmech__Boson__new( "photon", 0., 2, 0. );

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;
  
  d_sP =
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_ENTROPY_DENSITY,
      d_T,
      0.,
      NULL,
      NULL
    );

  Libstatmech__Boson__free( p_photon );

  return
    d_sP /
    (
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
      GSL_CONST_NUM_AVOGADRO *
      GSL_CONST_CGSM_BOLTZMANN
    );


}

//##############################################################################
// compute_entropy_per_nucleon().
//##############################################################################

double
compute_entropy_per_nucleon( nnt::Zone& zone )
{

  return
    compute_baryon_entropy_per_nucleon( zone ) +
    compute_electron_entropy_per_nucleon( zone ) +
    compute_photon_entropy_per_nucleon( zone );

}
  
//##############################################################################
// compute_baryon_specific_heat_per_nucleon().
//##############################################################################

double
compute_baryon_specific_heat_per_nucleon( nnt::Zone& zone )
{

  double d_cvb = 0;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    d_cvb +=
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      )
      *
      (
        1.5
        +
        (
          boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
          GSL_CONST_NUM_GIGA *
          compute_dlnG_dT(
            species.getNucnetSpecies(), 
            boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
              GSL_CONST_NUM_GIGA
          )
        )
      );

  }

  return d_cvb;

}

//##############################################################################
// compute_photon_specific_heat_per_nucleon().
//##############################################################################

double
compute_photon_specific_heat_per_nucleon( nnt::Zone& zone )
{

  return
    4.
    *
    (
      4. * GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT /
      GSL_CONST_CGSM_SPEED_OF_LIGHT
    )
    *
    gsl_pow_3(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
      GSL_CONST_NUM_GIGA
    )
    /
    (
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
        GSL_CONST_NUM_AVOGADRO *
        GSL_CONST_CGSM_BOLTZMANN
    );

}

//##############################################################################
// compute_electron_specific_heat_per_nucleon().
//##############################################################################

double
compute_electron_specific_heat_per_nucleon( nnt::Zone& zone )
{
  
  Libstatmech__Fermion * p_electron;
  double d_T, d_cve = 0, d_ne;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    S_ENTROPY_DENSITY,
    0.,
    nnt::d_REL_EPS
  );

  d_T =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
    GSL_CONST_NUM_GIGA;

  d_ne =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
    GSL_CONST_NUM_AVOGADRO *
    Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );

  d_cve =
    d_T *
    Libstatmech__Fermion__computeTemperatureDerivative(
      p_electron,
      S_ENTROPY_DENSITY,
      d_T,
      d_ne,
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return
    d_cve /
    (
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
      GSL_CONST_NUM_AVOGADRO *
      GSL_CONST_CGSM_BOLTZMANN
    );

}
  
//##############################################################################
// compute_specific_heat_per_nucleon().
//##############################################################################

double
compute_specific_heat_per_nucleon( nnt::Zone& zone )
{

  return
    compute_baryon_specific_heat_per_nucleon( zone ) +
    compute_electron_specific_heat_per_nucleon( zone ) +
    compute_photon_specific_heat_per_nucleon( zone );

}
  
//##############################################################################
// compute_thermo_quantity().
//##############################################################################

/**
 * \brief Compute a thermodynamic quantity from a zone.
 * \param zone A Nucnet-Tools zone.
 * \param s_quantity A string giving the quantity to compute.
 * \param s_particle The particle ("baryon", "electron", or "photon")
 *     to compute, or the total ("total").
 * \return The computed quantity.
 */

double
compute_thermo_quantity(
  nnt::Zone& zone,
  std::string s_quantity,
  std::string s_particle
)
{

  if( s_quantity == nnt::s_PRESSURE )
  {
    if( s_particle == nnt::s_BARYON )
      return compute_baryon_pressure( zone );
    else if( s_particle == nnt::s_ELECTRON )
      return compute_electron_pressure( zone );
    else if( s_particle == nnt::s_PHOTON )
      return compute_photon_pressure( zone );
    else if( s_particle == nnt::s_TOTAL )
      return compute_pressure( zone );
    else
    {
      std::cerr <<
        std::endl << s_particle << " not a valid particle type." << std::endl;
      exit( EXIT_FAILURE );
    }
  }
  else if( s_quantity == nnt::s_ENTROPY_PER_NUCLEON )
  {
    if( s_particle == nnt::s_BARYON )
      return compute_baryon_entropy_per_nucleon( zone );
    else if( s_particle == nnt::s_ELECTRON )
      return compute_electron_entropy_per_nucleon( zone );
    else if( s_particle == nnt::s_PHOTON )
      return compute_photon_entropy_per_nucleon( zone );
    else if( s_particle == nnt::s_TOTAL )
      return compute_entropy_per_nucleon( zone );
    else
    {
      std::cerr <<
        std::endl << s_particle << " not a valid particle type." << std::endl;
      exit( EXIT_FAILURE );
    }
  }
  else if( s_quantity == nnt::s_DPDT )
  {
    if( s_particle == nnt::s_BARYON )
      return compute_baryon_dPdT( zone );
    else if( s_particle == nnt::s_ELECTRON )
      return compute_electron_dPdT( zone );
    else if( s_particle == nnt::s_PHOTON )
      return compute_photon_dPdT( zone );
    else if( s_particle == nnt::s_TOTAL )
      return compute_dPdT( zone );
    else
    {
      std::cerr <<
        std::endl << s_particle << " not a valid particle type." << std::endl;
      exit( EXIT_FAILURE );
    }
  }
  else if( s_quantity == nnt::s_SPECIFIC_HEAT_PER_NUCLEON )
  {
    if( s_particle == nnt::s_BARYON )
      return compute_baryon_specific_heat_per_nucleon( zone );
    else if( s_particle == nnt::s_ELECTRON )
      return compute_electron_specific_heat_per_nucleon( zone );
    else if( s_particle == nnt::s_PHOTON )
      return compute_photon_specific_heat_per_nucleon( zone );
    else if( s_particle == nnt::s_TOTAL )
      return compute_specific_heat_per_nucleon( zone );
    else
    {
      std::cerr <<
        std::endl << s_particle << " not a valid particle type." << std::endl;
      exit( EXIT_FAILURE );
    }
  }
  else if( s_quantity == nnt::s_INTERNAL_ENERGY_DENSITY )
  {
    if( s_particle == nnt::s_BARYON )
      return compute_baryon_internal_energy_density( zone );
    else if( s_particle == nnt::s_ELECTRON )
      return compute_electron_internal_energy_density( zone );
    else if( s_particle == nnt::s_PHOTON )
      return compute_photon_internal_energy_density( zone );
    else if( s_particle == nnt::s_TOTAL )
      return compute_internal_energy_density( zone );
    else
    {
      std::cerr <<
        std::endl << s_particle << " not a valid particle type." << std::endl;
      exit( EXIT_FAILURE );
    }
  }

  else
  {
    std::cerr <<
      std::endl << s_quantity << " not a valid quantity type." << std::endl;
    exit( EXIT_FAILURE );
  }

}
  
//##############################################################################
// compute_log10_t9_entropy_root_with_equilibrium(). 
//##############################################################################

double
compute_log10_t9_entropy_root_with_equilibrium(
  double d_log10_t9,
  nnt::Zone * p_zone
)
{

  p_zone->updateProperty(
    nnt::s_T9,
    boost::lexical_cast<std::string>( pow( 10., d_log10_t9 ) )
  );

  nnt::set_zone_abundances_to_equilibrium( *p_zone );

  return
    compute_thermo_quantity(
      *p_zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    ) -
    boost::lexical_cast<double>(
      p_zone->getProperty( nnt::s_ENTROPY_PER_NUCLEON )
    );

}
    
//##############################################################################
// compute_log10_t9_entropy_root().
//##############################################################################

double
compute_log10_t9_entropy_root(
  double d_log10_t9,
  nnt::Zone * p_zone
)
{

  p_zone->updateProperty(
    nnt::s_T9,
    boost::lexical_cast<std::string>( pow( 10., d_log10_t9 ) )
  );

  return
    compute_thermo_quantity(
      *p_zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    ) -
    boost::lexical_cast<double>(
      p_zone->getProperty( nnt::s_ENTROPY_PER_NUCLEON )
    );

}

//##############################################################################
// compute_log10_density_entropy_root_with_equilibrium(). 
//##############################################################################

double
compute_log10_density_entropy_root_with_equilibrium(
  double d_log10_rho,
  nnt::Zone * p_zone 
)
{

  p_zone->updateProperty(
    nnt::s_RHO,
    boost::lexical_cast<std::string>( pow( 10., d_log10_rho ) )
  );

  nnt::set_zone_abundances_to_equilibrium( *p_zone );

  return
    compute_thermo_quantity(
      *p_zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    )
    -
    boost::lexical_cast<double>(
      p_zone->getProperty( nnt::s_ENTROPY_PER_NUCLEON )
    );

}

//##############################################################################
// compute_specific_heat_per_nucleon_at_constant_pressure().
//##############################################################################

double
compute_specific_heat_per_nucleon_at_constant_pressure( nnt::Zone& zone )
{

  zone.updateProperty(
    nnt::s_PRESSURE,
    boost::lexical_cast<std::string>(
      compute_thermo_quantity(
        zone,
        nnt::s_PRESSURE,
        nnt::s_TOTAL
      )
    )
  );

  std::string s_t9 = zone.getProperty( nnt::s_T9 );
  
  double d_cp =
    nnt::compute_temperature_derivative(
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_T9 )
      ) * GSL_CONST_NUM_GIGA,
      (nnt::quantityFunction) dS_dT_at_constant_pressure,
      &zone
    );

  zone.updateProperty( nnt::s_T9, s_t9 );

  return boost::lexical_cast<double>( s_t9 ) * GSL_CONST_NUM_GIGA * d_cp;

}
        
//##############################################################################
// dS_dT_at_constant_pressure().
//##############################################################################

double dS_dT_at_constant_pressure(
  double d_T,
  void * p_data
)
{

  nnt::Zone zone = *(nnt::Zone *) p_data;

  std::string s_rho = zone.getProperty( nnt::s_RHO );

  zone.updateProperty(
    nnt::s_T9, 
    boost::lexical_cast<std::string>( d_T / GSL_CONST_NUM_GIGA ) 
  );

  zone.setGuessFunction(
    (nnt::guessFunction) cp_rho_guess_function,
    &zone
  );

  zone.updateProperty(
    nnt::s_RHO, 
    boost::lexical_cast<std::string>(
      pow(
        10.,
        zone.computeRootFromQuantity(
          (nnt::quantityFunction) compute_log10_rho_pressure_root,
          &zone
        )
      )
    )
  );

  return compute_entropy_per_nucleon( zone );

}

//##############################################################################
// cp_rho_guess_function().
//##############################################################################

std::pair<double, double>
cp_rho_guess_function(
  nnt::Zone *p_zone
)
{

  double d_rho =
    boost::lexical_cast<double>( p_zone->getProperty( nnt::s_RHO ) );

  return
    std::make_pair<double,double>(
      log10( 0.9 * d_rho ),
      log10( 1.1 * d_rho )
    );

}

//##############################################################################
// compute_log10_rho_pressure_root().
//##############################################################################

double
compute_log10_rho_pressure_root(
  double d_log10_rho,
  void *p_data
)
{

  nnt::Zone zone = *(nnt::Zone *) p_data;

  zone.updateProperty(
    nnt::s_RHO,
    boost::lexical_cast<std::string>( pow( 10., d_log10_rho ) )
  );

  double d_result =
    compute_pressure( zone )
    -
    boost::lexical_cast<double>( zone.getProperty( nnt::s_PRESSURE ) );

  return d_result;

}

//##############################################################################
// compute_sound_speed().
//##############################################################################

double
compute_sound_speed( nnt::Zone& zone )
{

  double d_cs;

  zone.updateProperty(
    nnt::s_ENTROPY_PER_NUCLEON,
    boost::lexical_cast<std::string>(
      compute_thermo_quantity(
        zone,
        nnt::s_ENTROPY_PER_NUCLEON,
        nnt::s_TOTAL
      )
    )
  );

  std::string s_rho = zone.getProperty( nnt::s_RHO );
  
  d_cs =
    nnt::compute_density_derivative(
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_RHO )
      ),
      (nnt::quantityFunction) dP_drho,
      &zone
    );

  zone.updateProperty( nnt::s_RHO, s_rho );

  return sqrt( d_cs );

}
        
//##############################################################################
// dP_drho().
//##############################################################################

double dP_drho(
  double d_rho,
  void * p_data
)
{

  double d_t9;

  nnt::Zone * p_zone = (nnt::Zone *) p_data;

  std::string s_t9 = p_zone->getProperty( nnt::s_T9 );

  p_zone->updateProperty(
    nnt::s_RHO, 
    boost::lexical_cast<std::string>( d_rho ) 
  );

  p_zone->setGuessFunction(
    (nnt::guessFunction) sound_speed_t9_guess_function,
    p_zone
  );

  d_t9 =
    pow(
      10.,
      p_zone->computeRootFromQuantity(
        (nnt::quantityFunction) compute_log10_t9_entropy_root,
        p_data
      )
    );

  p_zone->updateProperty(
    nnt::s_T9, 
    boost::lexical_cast<std::string>( d_t9 ) 
  );

  double d_result = compute_pressure( *p_zone );

  p_zone->updateProperty( nnt::s_T9, s_t9 );

  return d_result;

}

//##############################################################################
// sound_speed_t9_guess_function().
//##############################################################################

std::pair<double, double>
sound_speed_t9_guess_function(
  nnt::Zone *p_zone
)
{

  double d_t9 = boost::lexical_cast<double>( p_zone->getProperty( nnt::s_T9 ) );

  return
    std::make_pair<double,double>(
      log10( 0.9 * d_t9 ),
      log10( 1.1 * d_t9 )
    );

}

//##############################################################################
// compute_dlnG_dT()
//##############################################################################

double
compute_dlnG_dT(
  Libnucnet__Species * p_species,
  double d_T
)
{

  gsl_function F;
  double d_result, d_abserr;

  F.function = &compute_lnG;
  F.params = p_species;

  gsl_deriv_central( &F, d_T, 1.e-6, &d_result, &d_abserr );

  return d_result;

}

//##############################################################################
// compute_lnG().
//##############################################################################

double
compute_lnG( double d_T, void * p_param )
{

  Libnucnet__Species * p_species = ( Libnucnet__Species * ) p_param;

  return
    log(
      Libnucnet__Species__computePartitionFunction(
        p_species,
        d_T / GSL_CONST_NUM_GIGA
      )
    );

}

} // namespace user
