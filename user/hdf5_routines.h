////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2013 Clemson University.
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
//! \brief A header file for routines to work with hdf5.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef HDF5_ROUTINES_H
#define HDF5_ROUTINES_H

//##############################################################################
// Includes.
//##############################################################################

#include <map>

#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/format.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include <H5Cpp.h>

#include "nnt/auxiliary.h"

#include "user/containers.h"

namespace user
{

namespace hdf5
{

//##############################################################################
// Defines.
//##############################################################################

#define I_HDF5_BUF   256

#define S_A                "A"
#define S_INDEX            "index"
#define S_LABEL_1          "Label 1"
#define S_LABEL_2          "Label 2"
#define S_LABEL_3          "Label 3"
#define S_MASS_EXCESS      "Mass Excess"
#define S_MASS_FRACTIONS   "Mass Fractions"
#define S_NAME             "Name"
#define S_NUCLIDE_DATA     "Nuclide Data"
#define S_STEP_GROUP       "Step"
#define S_TAG1             "Tag 1"
#define S_TAG2             "Tag 2"
#define S_VALUE            "Value"
#define S_Z                "Z"
#define S_ZONE_LABELS      "Zone Labels"
#define S_ZONE_PROPERTIES  "Zone Properties"

//##############################################################################
// Classes.
//##############################################################################

class zone_labels_tuple
{

  public:
    boost::tuple<std::string,std::string,std::string> tuple;
    zone_labels_tuple(
      const boost::tuple<std::string,std::string,std::string> t
    ) :
      tuple( t ) {}

    bool operator==(const zone_labels_tuple &p) const
    {
      return tuple == p.tuple;
    }

    friend std::size_t hash_value(zone_labels_tuple const& p)
    {
        std::size_t seed = 0;

        boost::hash_combine( seed, p.tuple.get<0>() );
        boost::hash_combine( seed, p.tuple.get<1>() );
        boost::hash_combine( seed, p.tuple.get<2>() );

        return seed;
    }

};

//##############################################################################
// Typedefs.
//##############################################################################

typedef boost::multi_array<double, 2> mass_fraction_array_t;

typedef struct nuclide
{
  char sName[I_HDF5_BUF];
  hsize_t iIndex;
  unsigned int iZ;
  unsigned int iA;
  double dMassExcess;
} nuclide;

typedef struct zone_labels_t
{
  char sLabel1[I_HDF5_BUF];
  char sLabel2[I_HDF5_BUF];
  char sLabel3[I_HDF5_BUF];
} zone_labels_t;

typedef struct zone_properties_t
{
  char sName[I_HDF5_BUF];
  char sTag1[I_HDF5_BUF];
  char sTag2[I_HDF5_BUF];
  char sValue[I_HDF5_BUF];
} zone_properties_t;

typedef std::map<std::string,nuclide> nuclide_map;

typedef std::pair<std::string,nuclide> nuclide_map_entry;

typedef
hashed_unique_bimap<
  zone_labels_tuple,
  size_t
>::type zone_labels_bimap;

typedef
struct zone_properties_hash_entry
{

  std::string sName;
  std::string sTag1;
  std::string sTag2;
  std::string sValue;
  zone_properties_hash_entry(
    const std::string& s1,
    const std::string& s2,
    const std::string& s3,
    const std::string& s4
  ) : sName( s1 ), sTag1( s2 ), sTag2( s3 ), sValue( s4 ) {}

} zone_properties_hash_entry;
    
typedef boost::multi_index::multi_index_container<
  zone_properties_hash_entry,
  boost::multi_index::indexed_by<
    boost::multi_index::hashed_unique< 
      boost::multi_index::composite_key<
        zone_properties_hash_entry,
        boost::multi_index::member<
          zone_properties_hash_entry,
          std::string,
          &zone_properties_hash_entry::sName
        >,
        boost::multi_index::member<
          zone_properties_hash_entry,
          std::string,
          &zone_properties_hash_entry::sTag1
        >,
        boost::multi_index::member<
          zone_properties_hash_entry,
          std::string,
          &zone_properties_hash_entry::sTag2
        >
      >
    >
  >
> zone_properties_hash;

//##############################################################################
// Prototypes.
//############################################################################*/

void create_output( const char *, Libnucnet * );

void append_zones( const char *, Libnucnet * );

void
populate_zone_properties(
  const char *,
  const char *,
  const char *,
  const char *,
  void *
);

void write_nuclide_data( H5::H5File&, Libnucnet * );

extern "C" 
{
herr_t
dataset_counter( hid_t, const char *, const H5L_info_t *, void * );
}

extern "C" 
{
herr_t
group_counter( hid_t, const char *, const H5L_info_t *, void * );
}

nuclide_map get_network( const char * );

boost::multi_array<double, 2>
get_step_mass_fractions( const char *, const char * );

zone_labels_bimap
get_zone_labels( const char *, const char * );

zone_properties_hash
get_zone_properties(
  const char *,
  const char *,
  const char *,
  const char *,
  const char *
);

size_t
count_groups_in_file(
  const char *
);

size_t
count_datasets_in_group(
  const char *
);

std::string get_step_label( size_t );

} // namespace hdf5
 
} // namespace user
 
#endif // HDF5_ROUTINES_H
