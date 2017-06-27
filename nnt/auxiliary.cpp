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
//! \brief Code file for a variety of useful auxiliary utilities.
//!
////////////////////////////////////////////////////////////////////////////////

#include "nnt/auxiliary.h"

namespace nnt
{

//############################################################################
// Zone::hasProperty().
//############################################################################

/**
 * \brief Determines whether a property is present in a zone.
 *
 * \param s_name A string giving the name of the property.
 * \param s_tag1 A string giving an 
 *                 extra tag associated with the property (optional).
 * \param s_tag2 A string giving a second extra tag associated with the
 *		      property (optional).
 * \return 1 (true) if the property is present and 0 (false) if not.
 */

int
Zone::hasProperty(
  std::string s_name,
  std::string s_tag1,
  std::string s_tag2
)
{

  if(
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      s_tag1.c_str(),
      s_tag2.c_str()
    )
  )
    return 1;
  else
    return 0;

} 

int
Zone::hasProperty(
  std::string s_name,
  std::string s_tag1
)
{

  if(
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      s_tag1.c_str(),
      NULL
    )
  )
    return 1;
  else
    return 0;

} 

int
Zone::hasProperty(
  std::string s_name
)
{

  if(
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      NULL,
      NULL
    )
  )
    return 1;
  else
    return 0;

} 

//############################################################################
// Zone::updateProperty().
//############################################################################

/**
 * \brief Updates a property in a zone from a double.
 *
 * \param s_name A string giving the name of the property.
 * \param s_tag1 A string giving an 
 *                 extra tag associated with the property (optional).
 * \param s_tag2 A string giving a second extra tag associated with the
 *		      property (optional).
 * \param s_value A string giving the value to be updated.
 * \return 1 (true) if update succeeded or 0 (false) if not.
 */

int
Zone::updateProperty(
  std::string s_name,
  std::string s_tag1,
  std::string s_tag2,
  std::string s_value
)
{

  return
    Libnucnet__Zone__updateProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      s_tag1.c_str(),
      s_tag2.c_str(),
      s_value.c_str()
    );

} 

int
Zone::updateProperty(
  std::string s_name,
  std::string s_tag1,
  std::string s_value
)
{

  return
    Libnucnet__Zone__updateProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      s_tag1.c_str(),
      NULL,
      s_value.c_str()
    );

} 

int
Zone::updateProperty(
  std::string s_name,
  std::string s_value
)
{

  return
    Libnucnet__Zone__updateProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      NULL,
      NULL,
      s_value.c_str()
    );

} 

//############################################################################
// Zone::getProperty().
//############################################################################

/**
  \brief Returns the current value of a property in a zone.

  \param s_name A string giving the name of the property.
  \param s_tag1 A string giving an
		  extra tag associated with the property (optional).
  \param s_tag2 A string giving a second extra tag associated with the
		    property (optional).
  \return The value of the property.
*/

std::string
Zone::getProperty(
  std::string s_name,
  std::string s_tag1,
  std::string s_tag2
)
{

  const char * s_value =
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(), s_name.c_str(), s_tag1.c_str(), s_tag2.c_str()
    );

  if( s_value )
    return boost::lexical_cast<std::string>( s_value );
  else
  {

    std::cerr << "Property " <<
      s_name << " " << s_tag1 <<
      " not present." << std::endl;

    exit( EXIT_FAILURE );

  }

}

std::string
Zone::getProperty(
  std::string s_name,
  std::string s_tag1
)
{

  const char * s_value =
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(), s_name.c_str(), s_tag1.c_str(), NULL
    );

  if( s_value )
    return boost::lexical_cast<std::string>( s_value );
  else
  {

    std::cerr << "Property " <<
      s_name << " " << s_tag1 <<
      " not present." << std::endl;

    exit( EXIT_FAILURE );

  }

}

std::string
Zone::getProperty(
  std::string s_name
)
{

  const char * s_value =
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(), s_name.c_str(), NULL, NULL
    );

  if( s_value )
    return boost::lexical_cast<std::string>( s_value );
  else
  {

    std::cerr << "Property " <<
      s_name <<
      " not present." << std::endl;

    exit( EXIT_FAILURE );

  }

}

//############################################################################
// print_zone_abunds()
//############################################################################

/**
  \brief Output the abundances in a zone.

*/

void
Zone::printAbundances()
{

  double d_t, d_dt, d_t9, d_rho;

  d_t = boost::lexical_cast<double>( this->getProperty( s_TIME ) );

  d_dt = boost::lexical_cast<double>( this->getProperty( s_DTIME ) );

  d_t9 = boost::lexical_cast<double>( this->getProperty( s_T9 ) );

  d_rho = boost::lexical_cast<double>( this->getProperty( s_RHO ) );

  printf(
    "t = %10.4e, dt = %10.4e, t9 = %10.4e, rho (g/cc) = %10.4e\n\n",
    d_t, d_dt, d_t9, d_rho
  );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__Zone__getNet( this->getNucnetZone() )
    ),
    (Libnucnet__Species__iterateFunction) print_abundance,
    this->getNucnetZone()
  );

  fprintf(
    stdout,
    "1 - xsum = %e\n\n",
    1. - Libnucnet__Zone__computeAMoment( this->getNucnetZone(), 1 )
  );

  fprintf(
    stdout,
    "Ye = %f\n\n",
    Libnucnet__Zone__computeZMoment( this->getNucnetZone(), 1 )
  );

}

//############################################################################
// print_abundance()
//############################################################################

int
print_abundance(
  Libnucnet__Species *p_species,
  Libnucnet__Zone *p_zone
)
{

  double d_abund;

  if( 
      ( d_abund =
	 Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species )
      ) > d_Y_MIN_PRINT
  )
  {
     printf( "%5lu%5lu%16.6e%16.6e\n",
       (unsigned long) Libnucnet__Species__getZ( p_species ),
       (unsigned long) Libnucnet__Species__getA( p_species ),
       d_abund,
       Libnucnet__Zone__getSpeciesAbundanceChange(
	 p_zone, p_species
     )
    );
  }

  return 1;

}

//############################################################################
// compare_reactions_by_string().
//############################################################################

/**
  \brief Standard reaction comparison function for network calculations
	  for reaction strings sorted alphabetically.

  \param p_reaction1 The first reaction.
  \param p_reaction2 The second reaction.
  \return 1 if reaction 1 is after reaction 2, -1 if reaction 1 is before
	  reaction 2, and 0 if the two reactions are the same.
*/

int
compare_reactions_by_string(
  const Libnucnet__Reaction *p_reaction1,
  const Libnucnet__Reaction *p_reaction2
)
{

  return
    strcmp(
      Libnucnet__Reaction__getString( p_reaction1 ),
      Libnucnet__Reaction__getString( p_reaction2 )
  );
 
}

//############################################################################
// compute_fermi_dirac_factor().
//############################################################################

/**
  \brief Compute the Fermi-Dirac factor.

  \param d_E_kT The energy divided by kT.
  \param d_mu_kT The chemical potential divided by kT.
  \return The Fermi-Dirac factor 1. / ( exp( d_E_kT - d_mu_kT ) + 1 ).
*/

double
compute_fermi_dirac_factor(
  double d_E_kT,
  double d_mu_kT
)
{

  double d_x = d_E_kT - d_mu_kT;

  if( gsl_isinf( d_x ) == -1 )
    return 1.;
  else if( gsl_isinf( d_x ) == 1 )
    return 0.;
  else
  {
    if( d_x > d_EXP_LARGE )
      return 0.;
    else if( d_x > 0 && d_x <= d_EXP_LARGE )
      return exp( -d_x ) / ( 1. + exp( -d_x ) );
    else
      return 1. / ( exp( d_x ) + 1. );
  }

}

//############################################################################
// compute_one_minus_fermi_dirac_factor().
//############################################################################

/**
  \brief Compute one minus the Fermi-Dirac factor.

  \param d_E_kT The energy divided by kT.
  \param d_mu_kT The chemical potential divided by kT.
  \return The Fermi-Dirac factor 1. - 1. / ( exp( d_E_kT - d_mu_kT ) + 1 ).
*/

double
compute_one_minus_fermi_dirac_factor(
  double d_E_kT,
  double d_mu_kT
)
{

  double d_x = d_E_kT - d_mu_kT;

  if( gsl_isinf( d_x ) == -1 )
    return 0.;
  else if( gsl_isinf( d_x ) == 1 )
    return 1.;
  else
  {
    if( d_x > d_EXP_LARGE )
      return 1.;
    else
      return exp( d_x ) / ( 1. + exp( d_x ) );
  }

}

//############################################################################
// compute_kT_in_MeV().
//############################################################################

/**
  \brief Compute kT, Boltzmann's constant k times the temperature T, in
	 MeV.

  \param d_t9 The temperature in billions of Kelvins.
  \return A double giving kT in MeV.
*/

double
compute_kT_in_MeV(
  double d_t9
)
{

  return
    d_t9 *
      GSL_CONST_CGSM_BOLTZMANN *
      GSL_CONST_NUM_KILO /
      GSL_CONST_CGSM_ELECTRON_VOLT;

}

//############################################################################
// compute_electron_chemical_potential_kT().
//############################################################################

/**
  \brief Compute the electron chemical potential (less the rest mass)
	 divided by kT.

  \param d_t9 The temperature in billions of Kelvins.
  \param d_rho The mass density in g/cc.
  \param d_ye The electron-to-nucleon ratio.
  \return A double giving the chemical potential / kT.
*/

double
compute_electron_chemical_potential_kT(
  double d_t9,
  double d_rho,
  double d_ye
)
{

  Libstatmech__Fermion *p_electron;
  double d_result;

  p_electron =
    Libstatmech__Fermion__new(
      s_ELECTRON,
      d_ELECTRON_MASS_IN_MEV,
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

//############################################################################
// zone_compare_by_first_label()
//############################################################################

/**
  \brief Standard zone comparison function which sorts the zones by
	  the numerical value associated with their first label.

  \param p_zone1 The first zone.
  \param p_zone2 The second zone.
  \return 1 if zone 1 is after zone 2, -1 if zone 1 is before
	  zone 2, and 0 if the two zones are the same.
*/

int
zone_compare_by_first_label(
  const Libnucnet__Zone *p_zone1,
  const Libnucnet__Zone *p_zone2
)
{

  if( 
      atof( Libnucnet__Zone__getLabel( p_zone1, 1 ) ) <
      atof( Libnucnet__Zone__getLabel( p_zone2, 1 ) )
  )
     return -1;
  else
     return 1;
	  
}

//############################################################################
// species_sort_by_z_then_a()
//############################################################################

/**
  \brief Standard species comparison function (sorted by ascending Z and,
         within Z, by A).

  \param p_species1 The first species.
  \param p_species2 The second species.
  \return 1 if species 1 is after species 2, -1 if species 1 is before
	  species 2, and 0 if the two species are the same.
*/

int
species_sort_by_z_then_a(
  const Libnucnet__Species *p_species1,
  const Libnucnet__Species *p_species2
)
{

  int i;

  if( 
      Libnucnet__Species__getZ( p_species1 ) <
      Libnucnet__Species__getZ( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getZ( p_species1 ) >
      Libnucnet__Species__getZ( p_species2 )
  )
    return 1;

  if( 
      Libnucnet__Species__getA( p_species1 ) <
      Libnucnet__Species__getA( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getA( p_species1 ) >
      Libnucnet__Species__getA( p_species2 )
  )
    return 1;

  i =
    strcmp(
      Libnucnet__Species__getName( p_species1 ),
      Libnucnet__Species__getName( p_species2 )
    );

  if( i == 0 ) {
    return 0;
  } else {
    return GSL_SIGN( i );
  }

}

//############################################################################
// species_sort_function()
//############################################################################

/**
  \brief Standard species comparison function for network calculations
	 with arrow solver (sorted by ascending Z and A but with
	 n, p, and he4 at the end).

  \param p_species1 The first species.
  \param p_species2 The second species.
  \return 1 if species 1 is after species 2, -1 if species 1 is before
	  species 2, and 0 if the two species are the same.
*/

int
species_sort_function(
  const Libnucnet__Species *p_species1,
  const Libnucnet__Species *p_species2
)
{

  int i;

  if( !strcmp( Libnucnet__Species__getName( p_species1 ), "he4" ) )
    return 1;
  else if( !strcmp( Libnucnet__Species__getName( p_species2 ), "he4" ) )
    return -1;

  if( !strcmp( Libnucnet__Species__getName( p_species1 ), "h1" ) )
    return 1;
  else if( !strcmp( Libnucnet__Species__getName( p_species2 ), "h1" ) )
    return -1;

  if( !strcmp( Libnucnet__Species__getName( p_species1 ), "n" ) )
    return 1;
  else if( !strcmp( Libnucnet__Species__getName( p_species2 ), "n" ) )
    return -1;

  if( 
      Libnucnet__Species__getZ( p_species1 ) <
      Libnucnet__Species__getZ( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getZ( p_species1 ) >
      Libnucnet__Species__getZ( p_species2 )
  )
    return 1;

  if( 
      Libnucnet__Species__getA( p_species1 ) <
      Libnucnet__Species__getA( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getA( p_species1 ) >
      Libnucnet__Species__getA( p_species2 )
  )
    return 1;

  i =
    strcmp(
      Libnucnet__Species__getName( p_species1 ),
      Libnucnet__Species__getName( p_species2 )
    );

  if( i == 0 ) {
    return 0;
  } else {
    return GSL_SIGN( i );
  }

}

//##############################################################################
// set_zone_abundances_to_equilibrium().
//##############################################################################

void
set_zone_abundances_to_equilibrium( nnt::Zone& zone )
{

  Libnuceq * p_my_equil;
  gsl_vector * p_abundances;

  p_my_equil = get_zone_equilibrium( zone );

  Libnuceq__computeEquilibrium(
    p_my_equil, 
    boost::lexical_cast<double>( zone.getProperty( s_T9 ) ),
    boost::lexical_cast<double>( zone.getProperty( s_RHO ) )
  );

  set_zone_cluster_data( zone, p_my_equil );

  p_abundances = Libnuceq__getAbundances( p_my_equil );

  Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_abundances );

  zone.updateProperty(
    nnt::s_MUPKT,
    boost::lexical_cast<std::string>( Libnuceq__getMupkT( p_my_equil ) )
  );

  zone.updateProperty(
    nnt::s_MUNKT,
    boost::lexical_cast<std::string>( Libnuceq__getMunkT( p_my_equil ) )
  );

  gsl_vector_free( p_abundances );

  Libnuceq__free( p_my_equil );

}

//##############################################################################
// get_zone_equilibrium(). 
//##############################################################################

Libnuceq *
get_zone_equilibrium( nnt::Zone& zone )
{

  Libnuceq * p_equil;

  p_equil =
    Libnuceq__new(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      )
    );

  update_equilibrium_with_zone_data( p_equil, zone );

  return p_equil;

}

//##############################################################################
// get_clusters().
//##############################################################################

int
get_clusters(
  Libnuceq__Cluster * p_cluster,
  std::vector<Libnuceq__Cluster *> * p_clusters
)
{

  p_clusters->push_back( p_cluster );

  return 1;

}

//##############################################################################
// update_zone_equilibrium(). 
//##############################################################################

/**
 * \brief Updates an equilibrium with data from a zone.
 *
 * \param p_equil An equilibrium.
 * \param zone A Nucnet Tools zone with the data.
 * \return On successful return, the equilibrium has been updated with electro
 *         fraction and cluster data from the zone.
 */

void
update_equilibrium_with_zone_data( Libnuceq * p_equil, nnt::Zone& zone )
{

  Libnuceq__Cluster *p_cluster;
  int i;

  if( zone.getNseCorrectionFactorFunction() )
    Libnuceq__setNseCorrectionFactorFunction(
      p_equil,
      (Libnucnet__Species__nseCorrectionFactorFunction)
         zone.getNseCorrectionFactorFunction(),
      zone.getNseCorrectionFactorFunctionData()
    );

  if( zone.hasProperty( s_YE ) )
    Libnuceq__setYe(
      p_equil, 
      boost::lexical_cast<double>( zone.getProperty( s_YE ) )
    );

  if( Libnuceq__getNumberOfClusters( p_equil ) != 0 )
  {

    std::vector<Libnuceq__Cluster *> clusters;

    Libnuceq__iterateClusters(
      p_equil,
      (Libnuceq__Cluster__iterateFunction) get_clusters,
      &clusters
    );

    for( size_t i = 0; i < clusters.size(); i++ )
    {
      Libnuceq__removeCluster( p_equil, clusters[i] );
    }

  }

  if(
    zone.hasProperty( s_NUMBER_CLUSTERS )
  )
  {
    for(
      i = 0;
      i <
        boost::lexical_cast<int>( zone.getProperty( s_NUMBER_CLUSTERS ) );
      i++
    )
    {
      p_cluster =
        Libnuceq__newCluster(
          p_equil,
          zone.getProperty(
            s_CLUSTER_XPATH,
            boost::lexical_cast<std::string>( i )
          ).c_str()
        );
      Libnuceq__Cluster__updateConstraint(
        p_cluster,
        boost::lexical_cast<double>(
          zone.getProperty(
            s_CLUSTER_CONSTRAINT,
            boost::lexical_cast<std::string>( i )
          )
        )
      );
    }
  }

}

//##############################################################################
// set_zone_cluster_data().
//##############################################################################

void
set_zone_cluster_data( nnt::Zone& zone, Libnuceq * p_equil )
{

  int i;

  if(
      !zone.hasProperty( s_NUMBER_CLUSTERS )
  )
    return;

  for(
    i = 0;
    i < boost::lexical_cast<int>( zone.getProperty( s_NUMBER_CLUSTERS ) );
    i++
  )
  {
    zone.updateProperty(
      s_CLUSTER_CHEMICAL_POTENTIAL,
      boost::lexical_cast<std::string>( i ),
      boost::lexical_cast<std::string>(
        Libnuceq__Cluster__getMukT(
          Libnuceq__getCluster(
            p_equil,
            zone.getProperty(
              s_CLUSTER_XPATH,
              boost::lexical_cast<std::string>( i )
            ).c_str()
          )
        )
      )
    );
  }

}

//############################################################################
// compute_reaction_nuclear_Qvalue().
//############################################################################

double
compute_reaction_nuclear_Qvalue(
  Libnucnet__Net *p_net,
  Libnucnet__Reaction *p_reaction,
  double d_mec2
)
{

  double d_result;

  d_result = Libnucnet__Net__computeReactionQValue( p_net, p_reaction );

  if( is_beta_minus_reaction( p_reaction ) )
    d_result += d_mec2;
  else if( is_electron_capture_reaction( p_reaction ) )
    d_result -= d_mec2;
  else if( is_beta_plus_reaction( p_reaction ) )
    d_result += d_mec2;
  else if( is_positron_capture_reaction( p_reaction ) )
    d_result -= d_mec2;
    
  return d_result;

}

//############################################################################
// is_electron_capture_reaction().
//############################################################################

int
is_electron_capture_reaction(
  Libnucnet__Reaction *p_reaction
)
{

  if(
    strstr( Libnucnet__Reaction__getString( p_reaction ), "electron" ) &&
    strstr( Libnucnet__Reaction__getString( p_reaction ), "neutrino_e" )
  )
  {
    if(
      strstr( Libnucnet__Reaction__getString( p_reaction ), "anti-neutrino_e" )
    )
      return 0;
    else
      return 1;
  }

  return 0;

}

//############################################################################
// is_beta_minus_reaction().
//############################################################################

int
is_beta_minus_reaction(
  Libnucnet__Reaction *p_reaction
)
{

  if(
    strstr( Libnucnet__Reaction__getString( p_reaction ), "electron" ) &&
    strstr( Libnucnet__Reaction__getString( p_reaction ), "anti-neutrino_e" )
  )
    return 1;
  else
    return 0;

}

//############################################################################
// is_positron_capture_reaction().
//############################################################################

int
is_positron_capture_reaction(
  Libnucnet__Reaction *p_reaction
)
{

  if(
    strstr( Libnucnet__Reaction__getString( p_reaction ), "positron" ) &&
    strstr( Libnucnet__Reaction__getString( p_reaction ), "anti-neutrino_e" )
  )
    return 1;
  else
    return 0;

}

//############################################################################
// is_beta_plus_reaction().
//############################################################################

int
is_beta_plus_reaction(
  Libnucnet__Reaction *p_reaction
)
{

  if(
    strstr( Libnucnet__Reaction__getString( p_reaction ), "positron" ) &&
    strstr( Libnucnet__Reaction__getString( p_reaction ), "neutrino_e" )
  )
  {
    if(
      strstr( Libnucnet__Reaction__getString( p_reaction ), "anti-neutrino_e" )
    )
      return 0;
    else
      return 1;
  }

  return 0;

}

//##############################################################################
// Zone::createOutReactionXPathString().
//##############################################################################

/**
 * A method that creates an xpath string of out reactions.
 * \return A pointer to a std::string giving the XPath expression.
 */

std::string *
Zone::createOutReactionXPathString(
  const char * s_nuc_xpath,
  const char * s_reac_xpath
)
{

  Libnucnet__NetView * p_view1, * p_view2;
  std::string * s_out_reaction_xpath, * s_reaction;
  int i_c;

  p_view1 = this->getNetView( "", s_reac_xpath );

  p_view2 = this->getNetView( s_nuc_xpath, s_reac_xpath );

  reaction_list_t reaction_list =
    make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_view1 ) )
    );

  boost::ptr_list<std::string> out_reaction_list;

  for(
    boost::ptr_list<nnt::Reaction>::iterator it = reaction_list.begin();
    it != reaction_list.end();
    it++
  )
  {

    i_c = 0;

    reaction_element_list_t reactant_list =
      make_reaction_nuclide_reactant_list( it->getNucnetReaction() );

    for(
      boost::ptr_list<nnt::ReactionElement>::iterator
        iter = reactant_list.begin();
      iter != reactant_list.end();
      iter++
    )
    {

      if(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_view2 ) ),
          Libnucnet__Reaction__Element__getName(
            iter->getNucnetReactionElement()
          )
        )
      )
        i_c++;

    }

    reaction_element_list_t product_list =
      make_reaction_nuclide_product_list( it->getNucnetReaction() );

    for(
      boost::ptr_list<nnt::ReactionElement>::iterator
        iter = product_list.begin();
      iter != product_list.end();
      iter++
    )
    {

      if(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_view2 ) ),
          Libnucnet__Reaction__Element__getName(
            iter->getNucnetReactionElement()
          )
        )
      )
        i_c--;

    }

    if( i_c != 0 )
    {
      s_reaction = new std::string;
      *s_reaction = Libnucnet__Reaction__getString( it->getNucnetReaction() );
      out_reaction_list.push_back( s_reaction );
    }

  }

  s_out_reaction_xpath =
    create_reac_xpath_from_list(
      Libnucnet__Zone__getNet( this->getNucnetZone() ),
      out_reaction_list
    );

  return s_out_reaction_xpath;

}

//##############################################################################
// create_nuc_xpath_from_set().
//##############################################################################

std::string *
create_nuc_xpath_from_set(
  Libnucnet__Net * p_net,
  boost::unordered_set<std::string> * p_nuc_set
)
{

  std::string * s_nuc_xpath = new std::string;
  std::stringstream s_tmp;

  *s_nuc_xpath = "[";

  for(
    boost::unordered_set<std::string>::iterator it = p_nuc_set->begin();
    it != p_nuc_set->end();
  )
  {
   
    s_tmp <<
      "(z=" <<
      Libnucnet__Species__getZ(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( p_net ),
          (*it).c_str()
        )
      ) <<
      " and a=" <<
      Libnucnet__Species__getA(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( p_net ),
          (*it).c_str()
        )
      ) << ")";
    
    if( ++it != p_nuc_set->end() ) s_tmp << " or ";

  }

  *s_nuc_xpath += s_tmp.str();

  *s_nuc_xpath += "]";

  return s_nuc_xpath;

}

//##############################################################################
// create_reac_xpath_from_list().
//##############################################################################

std::string *
create_reac_xpath_from_list(
  Libnucnet__Net * p_net,
  boost::ptr_list<std::string> reac_list
)
{

  std::string * s_reac_xpath = new std::string;
  std::stringstream s_tmp;

  *s_reac_xpath = "[";

  for(
    boost::ptr_list<std::string>::iterator it = reac_list.begin();
    it != reac_list.end();
  )
  {

    reaction_element_list_t reactant_list =
      make_reaction_reactant_list(
        Libnucnet__Reac__getReactionByString(
          Libnucnet__Net__getReac( p_net ),
          (*it).c_str()
        )
      );

    s_tmp << "(";

    for(
      boost::ptr_list<ReactionElement>::iterator iter =
        reactant_list.begin();
      iter != reactant_list.end();
      iter++
    )
    {   
   
      s_tmp <<
	"reactant = " << "'" <<
        Libnucnet__Reaction__Element__getName(
          iter->getNucnetReactionElement()
        ) <<
        "' and ";

    }

    reaction_element_list_t product_list =
      make_reaction_product_list(
        Libnucnet__Reac__getReactionByString(
          Libnucnet__Net__getReac( p_net ),
          (*it).c_str()
        )
      );

    for(
       reaction_element_list_t::iterator iter = product_list.begin();
       iter != product_list.end();
    )
    {   
     
       s_tmp <<
       "product = " << "'" <<
       Libnucnet__Reaction__Element__getName(
         iter->getNucnetReactionElement()
       ) <<
       "'";

       if( ++iter != product_list.end() ) s_tmp << " and ";

     }

     s_tmp << ")";

     if( ++it != reac_list.end() ) s_tmp << " or ";

  }

  *s_reac_xpath += s_tmp.str();

  *s_reac_xpath += "]";

  return s_reac_xpath;

}

//##############################################################################
// Zone::normalizeAbundances()
//##############################################################################

void
Zone::normalizeAbundances( )
{

  gsl_vector *p_abundances;
  double d_xm;

  p_abundances =
    Libnucnet__Zone__getAbundances( this->getNucnetZone() );

  d_xm = Libnucnet__Zone__computeAMoment( this->getNucnetZone(), 1 );

  if( d_xm <= 0. )
  {
    std::cerr << "Invalid abundances in zone." << std::endl;
    exit( EXIT_FAILURE );
  }

  gsl_vector_scale(
    p_abundances,
    1. / d_xm
  );

  Libnucnet__Zone__updateAbundances( this->getNucnetZone(), p_abundances );

  gsl_vector_free( p_abundances );

}

//##############################################################################
// Zone::getElectronNeutrinoChemicalPotential()
//##############################################################################

double
Zone::getElectronNeutrinoChemicalPotential()
{

/**
 * \brief Returns the neutrino chemical potential for a zone.
 *
 * \return The neutrino chemical potential divided by kT.  If this value is
 *          -infinity, it is returns as GSL_NEGINF.
 */


  double d_mu_nue_kT;

  try{
    d_mu_nue_kT =
      boost::lexical_cast<double>( this->getProperty( s_MU_NUE_KT ) );
  }
  catch( const boost::bad_lexical_cast& e )
  {
    if( this->getProperty( s_MU_NUE_KT ) == "-inf" )
    {
      d_mu_nue_kT = GSL_NEGINF;
    }
    else
    {
      std::cerr << "Neutrino chemical potential " <<
        this->getProperty( s_MU_NUE_KT ) <<
        " not valid." << std::endl;
      throw e;
    }
  }

  return d_mu_nue_kT;

}

//##############################################################################
// get_property_array_size()
//##############################################################################

void
get_property_array_size(
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  const char *s_value,
  size_t *p_count
)
{

  if(
    ( s_tag1 && !s_name ) ||
    ( s_tag2 && !s_name ) ||
    ( !s_value || !s_name ))
  {
    fprintf( stderr, "Invalid input.\n" );
    exit( EXIT_FAILURE );
  }

  (*p_count)++;

}

//##############################################################################
// get_new_gsl_vector_from_std_vector();
//##############################################################################

gsl_vector *
get_new_gsl_vector_from_std_vector(
  const std::vector<double> &old_vector
)
{

  size_t i = 0;
  gsl_vector * p_new_vector = gsl_vector_alloc( old_vector.size() );

  BOOST_FOREACH( double element, old_vector )
  {
    gsl_vector_set( p_new_vector, i++, element );
  }

  return p_new_vector;

}

//##############################################################################
// get_property_array()
//##############################################################################

void
get_property_array(
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  const char *s_value,
  gsl_vector *p_vector
)
{

  if( s_tag2 || !s_name )
  {
    fprintf( stderr, "Invalid input.\n" );
    exit( EXIT_FAILURE );
  }

  gsl_vector_set(
    p_vector,
    boost::lexical_cast<size_t>( s_tag1 ),
    boost::lexical_cast<double>( s_value )
  );

}

//##############################################################################
// get_property_matrix()
//##############################################################################

void
get_property_matrix(
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  const char *s_value,
  gsl_matrix *p_matrix
)
{

  if( !s_name )
  {
    fprintf( stderr, "Invalid input.\n" );
    exit( EXIT_FAILURE );
  }

  gsl_matrix_set(
    p_matrix,
    boost::lexical_cast<size_t>( s_tag1 ),
    boost::lexical_cast<size_t>( s_tag2 ),
    atof( s_value )
  );

}

//##############################################################################
// get_rate_function_property_gsl_vector()
//##############################################################################

gsl_vector *
get_rate_function_property_gsl_vector(
  Libnucnet__Reaction *p_reaction,
  const char *s_name
)
{

  size_t i_size = 0;
  gsl_vector *p_vector;

  Libnucnet__Reaction__iterateUserRateFunctionProperties(
    p_reaction,
    s_name,
    NULL,
    NULL,
    (Libnucnet__Reaction__user_rate_property_iterate_function)
       get_property_array_size,
    &i_size
  );

  p_vector = gsl_vector_alloc( i_size );

  Libnucnet__Reaction__iterateUserRateFunctionProperties(
    p_reaction,
    s_name,
    NULL,
    NULL,
    (Libnucnet__Reaction__user_rate_property_iterate_function)
       get_property_array,
    p_vector
  );

  return p_vector;

}

//##############################################################################
// reaction_element_count().
//##############################################################################

/**
 * \brief Counts the number of occurrences of a reaction element in a reaction.
 *
 * \param reaction A NucNet Tools reaction.
 * \param s_type A string determining whether the element is a "reactant"
 *               or "product".
 * \param s_element The name of the element.
 * \return The number of occurrences of the element in the reaction.
 */

size_t
reaction_element_count(
  Reaction& reaction,
  std::string s_type,
  std::string s_element
)
{

  size_t i_count = 0;
  nnt::reaction_element_list_t element_list;

  if( s_type == "reactant" )
  {
    element_list =
      nnt::make_reaction_reactant_list( reaction.getNucnetReaction() );
  }
  else if( s_type == "product" )
  {
    element_list =
      nnt::make_reaction_product_list( reaction.getNucnetReaction() );
  }
  else
  {
    std::cerr << "No such reaction element type." << std::endl;
    exit( EXIT_FAILURE );
  }

  BOOST_FOREACH( nnt::ReactionElement element, element_list )
  {
    if(
      std::string(
        Libnucnet__Reaction__Element__getName(
          element.getNucnetReactionElement()
        )
      ) == s_element
    )
    {
      ++i_count;
    }
  }

  return i_count;

}

//##############################################################################
// get_nuc_anchors().
//##############################################################################

std::vector<Libnucnet__Species *>
get_nuc_anchors( Libnucnet__Nuc * p_nuc )
{

  unsigned int iz1, in2, iz3, in4;
  unsigned int izmin, izmax, inmin, inmax;
  Libnucnet__Species *p1, *p2, *p3, *p4;

  std::vector<Libnucnet__Species *> my_vector;

  p1 = p2 = p3 = p4 = NULL;

  izmax = Libnucnet__Nuc__getLargestNucleonNumber( p_nuc, "z" );
  inmax = Libnucnet__Nuc__getLargestNucleonNumber( p_nuc, "n" );

  izmin = izmax;
  inmin = inmax;

  species_list_t species_list = make_species_list( p_nuc );

  BOOST_FOREACH( Species sp, species_list )
  {

    unsigned int i_z = Libnucnet__Species__getZ( sp.getNucnetSpecies() );
    unsigned int i_n =
      Libnucnet__Species__getA( sp.getNucnetSpecies() ) - i_z;

    if( i_z < izmin ) izmin = i_z; 
    if( i_n < inmin ) inmin = i_n; 

  }

  iz1 = izmax;
  in2 = inmin;
  iz3 = izmin;
  in4 = inmax;

  BOOST_FOREACH( Species sp, species_list )
  {

    unsigned int i_z = Libnucnet__Species__getZ( sp.getNucnetSpecies() );
    unsigned int i_n =
      Libnucnet__Species__getA( sp.getNucnetSpecies() ) - i_z;

    if( i_z <= iz1 && i_n == inmin )
    {
      iz1 = i_z;
      p1 = sp.getNucnetSpecies();
    }

    if( i_z == izmax && i_n > in2 )
    {
      in2 = i_n;
      p2 = sp.getNucnetSpecies();
    }

    if( i_z >= iz3 && i_n == inmax )
    {
      iz3 = i_z;
      p3 = sp.getNucnetSpecies();
    }

    if( i_z == izmin && i_n < in4 )
    {
      in4 = i_n;
      p4 = sp.getNucnetSpecies();
    }

  }

  my_vector.push_back( p1 );
  my_vector.push_back( p2 );
  my_vector.push_back( p3 );
  my_vector.push_back( p4 );

  return my_vector;

}

//##############################################################################
// get_stable_species().
//##############################################################################

/**
 * \brief Get a vector of the naturally-occurring species.
 *
 * \return A vector with the stable species.
 */

std::vector<std::string>
get_stable_species()
{

  const char *s_stables[] =
    {
      "h1","h2",
      "he3","he4",
      "li6","li7",
      "be9",
      "b10","b11",
      "c12","c13",
      "n14","n15",
      "o16","o17","o18",
      "f19",
      "ne20","ne21","ne22",
      "na23",
      "mg24","mg25","mg26",
      "al27",
      "si28","si29","si30",
      "p31",
      "s32","s33","s34","s36",
      "cl35","cl37",
      "ar36","ar38","ar40",
      "k39","k40","k41",
      "ca40","ca42","ca43","ca44","ca46","ca48",
      "sc45",
      "ti46","ti47","ti48","ti49","ti50",
      "v50","v51",
      "cr50","cr52","cr53","cr54",
      "mn55",
      "fe54","fe56","fe57","fe58",
      "co59",
      "ni58","ni60","ni61","ni62","ni64",
      "cu63","cu65",
      "zn64","zn66","zn67","zn68","zn70",
      "ga69","ga71",
      "ge70","ge72","ge73","ge74","ge76",
      "as75",
      "se74","se76","se77","se78","se80","se82",
      "br79","br81",
      "kr78","kr80","kr82","kr83","kr84","kr86",
      "rb85","rb87",
      "sr84","sr86","sr87","sr88",
      "y89",
      "zr90","zr91","zr92","zr94","zr96",
      "nb93",
      "mo92","mo94","mo95","mo96","mo97","mo98","mo100",
      "ru96","ru98","ru99","ru100","ru101","ru102","ru104",
      "rh103",
      "pd102","pd104","pd105","pd106","pd108","pd110",
      "ag107","ag109",
      "cd106","cd108","cd110","cd111","cd112","cd113","cd114","cd116",
      "in113","in115",
      "sn112","sn114","sn115","sn116","sn117","sn118","sn119","sn120","sn122",
         "sn124",
      "sb121","sb123",
      "te120","te122","te123","te124","te125","te126","te128","te130",
      "i127",
      "xe124","xe126","xe128","xe129","xe130","xe131","xe132","xe134","xe136",
      "cs133",
      "ba130","ba132","ba134","ba135","ba136","ba137","ba138",
      "la138","la139",
      "ce136","ce138","ce140","ce142",
      "pr141",
      "nd142","nd143","nd144","nd145","nd146","nd148","nd150",
      "sm144","sm147","sm148","sm149","sm150","sm152","sm154",
      "eu151","eu153",
      "gd152","gd154","gd155","gd156","gd157","gd158","gd160",
      "tb159",
      "dy156","dy158","dy160","dy161","dy162","dy163","dy164",
      "ho165",
      "er162","er164","er166","er167","er168","er170",
      "tm169",
      "yb168","yb170","yb171","yb172","yb173","yb174","yb176",
      "lu175","lu176",
      "hf174","hf176","hf177","hf178","hf179","hf180",
      "ta180","ta181",
      "w180","w182","w183","w184","w186",
      "re185","re187",
      "os184","os186","os187","os188","os189","os190","os192",
      "ir191","ir193",
      "pt190","pt192","pt194","pt195","pt196","pt198",
      "au197",
      "hg196","hg198","hg200","hg201","hg202","hg204",
      "tl203","tl205",
      "pb204","pb206","pb207","pb208",
      "bi209",
      "th232",
      "u235","u238"
    };

  size_t i, i_species;
  std::vector<std::string> stable_species;

  i_species = sizeof( s_stables ) / sizeof( s_stables[0] );

  for( i = 0; i < i_species; i++ )
    stable_species.push_back( s_stables[i] );

  return stable_species;

} 

} // namespace nnt
