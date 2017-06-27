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
//!
//! \file
//! \brief Code for modifying selected rates.
//!
////////////////////////////////////////////////////////////////////////////////

#include "rate_modifiers.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

//##############################################################################
// modify_rates_in_view().
//##############################################################################

void
modify_rates_in_view(
  nnt::Zone& zone,
  const char * s_nuc_xpath,
  const char * s_reac_xpath,
  double d_factor
)
{

  Libnucnet__NetView * p_view;
  Libnucnet__Reaction * p_reaction;
  double d_forward, d_reverse;

  p_view = zone.getNetView( s_nuc_xpath, s_reac_xpath );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_view ) )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    //==========================================================================
    // Check that the reaction in is the evolution view.
    //==========================================================================

    p_reaction =
      Libnucnet__Reac__getReactionByString(
        Libnucnet__Net__getReac(
          Libnucnet__NetView__getNet(
            zone.getNetView(
              EVOLUTION_NETWORK,
              NULL
            )
          )
        ),
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
      );

    if( p_reaction )
    {

      Libnucnet__Zone__getRatesForReaction(
        zone.getNucnetZone(),
        reaction.getNucnetReaction(),
        &d_forward,
        &d_reverse
      );

      Libnucnet__Zone__updateRatesForReaction(
        zone.getNucnetZone(),
        reaction.getNucnetReaction(),
        d_forward * d_factor,
        d_reverse * d_factor
      );

    }

  }

}

//##############################################################################
// modify_rates().
//##############################################################################

void
modify_rates( nnt::Zone& zone )
{

  view_multi views; 

  Libnucnet__Zone__iterateOptionalProperties(
    zone.getNucnetZone(),
    NULL,
    nnt::s_RATE_MODIFICATION_VIEW,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) modify_rates_views,
    &views
  );

  BOOST_FOREACH( view view, views )
  {

    modify_rates_in_view(
      zone,
      view.nuc_xpath.c_str(),
      view.reac_xpath.c_str(),
      view.factor
    );

  }

}
    
//##############################################################################
// print_modified_reactions().
//##############################################################################

void
print_modified_reactions( nnt::Zone& zone )
{

  view_multi views; 

  Libnucnet__Zone__iterateOptionalProperties(
    zone.getNucnetZone(),
    NULL,
    nnt::s_RATE_MODIFICATION_VIEW,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) modify_rates_views,
    &views
  );

  BOOST_FOREACH( view my_view, views )
  {

    std::cout << std::endl <<
                 "Rate modification view "<< my_view.id << ":" <<
                 std::endl << std::endl;

    fprintf( stdout, "\t\t\tReaction\t\t\t Modification Factor \n" );

    Libnucnet__Reac__setReactionCompareFunction(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet(
          zone.getNetView(
            my_view.nuc_xpath.c_str(),
            my_view.reac_xpath.c_str()
          )
        )
      ),
      (Libnucnet__Reaction__compare_function)
         nnt::compare_reactions_by_string
    );

    nnt::reaction_list_t reaction_list =
      nnt::make_reaction_list(
        Libnucnet__Net__getReac(
          Libnucnet__NetView__getNet(
            zone.getNetView(
              my_view.nuc_xpath.c_str(),
              my_view.reac_xpath.c_str()
            )
          )
        )
      );

    BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
    {

      fprintf(
        stdout,
        "%-55s%12.4e\n",
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
        my_view.factor
      );

    }

  }

}
    
//##############################################################################
// modify_rates_views().
//##############################################################################

int
modify_rates_views(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  view_multi * p_views
)
{

  std::string s_tmp, s_empty = "";

  if( strcmp( s_tag1, nnt::s_RATE_MODIFICATION_VIEW ) )
  {
    std::cerr << "Invalid input for reaction rate modification view." <<
      std::endl;
    exit( EXIT_FAILURE );
  }

  view_multi::iterator it = p_views->find( s_tag2 ); 

  if( it == p_views->end() )
  {
    s_tmp = s_tag2;
    p_views->insert( view( s_tmp, s_empty, s_empty, 1. ) );
  }

  view_multi::iterator it_work = p_views->find( s_tag2 ); 

  if( strcmp( s_name, nnt::s_RATE_MODIFICATION_NUC_XPATH ) == 0 )
  {
    view change = *it_work;
    s_tmp = s_value;
    p_views->modify( it_work, change_nuc_xpath( s_tmp ) );
  }
  else if( strcmp( s_name, nnt::s_RATE_MODIFICATION_REAC_XPATH ) == 0 )
  {
    view change = *it_work;
    s_tmp = s_value;
    p_views->modify( it_work, change_reac_xpath( s_tmp ) );
  }
  else if( strcmp( s_name, nnt::s_RATE_MODIFICATION_FACTOR ) == 0 )
  {
    view change = *it_work;
    p_views->modify( it_work, change_factor( atof( s_value ) ) );
  }
  else
  {
    std::cerr << "Invalid entry for reaction view modification." << std::endl;
    exit( EXIT_FAILURE );
  }

  return 1;

}
  
//##############################################################################
// modify_rates_for_reaction().
//##############################################################################

void
modify_rates_for_reaction(
  nnt::Zone& zone,
  Libnucnet__Reaction * p_reaction,
  double& d_forward,
  double& d_reverse
)
{

  view_multi views;

  Libnucnet__Zone__iterateOptionalProperties(
    zone.getNucnetZone(),
    NULL,
    nnt::s_RATE_MODIFICATION_VIEW,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) modify_rates_views,
    &views
  );

  BOOST_FOREACH( view my_view, views )
  {

    //==========================================================================
    // Check that the reaction in is the view.
    //==========================================================================

    if(
      Libnucnet__Reac__getReactionByString(
        Libnucnet__Net__getReac(
          Libnucnet__NetView__getNet(
            zone.getNetView(
              my_view.nuc_xpath.c_str(),
              my_view.reac_xpath.c_str()
             )
          )
        ),
        Libnucnet__Reaction__getString( p_reaction )
      )
    )
    {

      d_forward *= my_view.factor;
      d_reverse *= my_view.factor;

    }

  }

}

} // namespace user
