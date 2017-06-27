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
//
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A file containing routines to work with hdf5.
//!
////////////////////////////////////////////////////////////////////////////////

#include "hdf5_routines.h"

namespace user
{

namespace hdf5
{

//##############################################################################
// create_output()
//##############################################################################

void
create_output(
  const char * s_file,
  Libnucnet * p_nucnet
)
{

  H5::H5File file( s_file, H5F_ACC_TRUNC );

  write_nuclide_data( file, p_nucnet );

}

//##############################################################################
// append_zones()
//##############################################################################

void
append_zones(
  const char * s_file,
  Libnucnet * p_nucnet
)
{

  H5::DataSpace dataspace;
  H5::DataSet dataset;
  H5::CompType my_type;
  hsize_t * dims;

  //===========================================================================
  // Get string type for use in output.
  //===========================================================================

  H5::StrType string_type( H5::PredType::C_S1, I_HDF5_BUF );

  //===========================================================================
  // Open file and count groups.
  //===========================================================================

  H5::H5File file( s_file, H5F_ACC_RDWR );

  std::string s_g = get_step_label( count_groups_in_file( s_file ) );

  H5::Group group = file.createGroup( s_g.c_str() );

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) )
    );

  Libnucnet__setZoneCompareFunction(
    p_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_nucnet );

  //===========================================================================
  // Mass Fractions.
  //===========================================================================

  std::vector<double> mass_fractions;

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
    gsl_vector * p_abunds =
      Libnucnet__Zone__getAbundances( zone.getNucnetZone() );
    BOOST_FOREACH( nnt::Species species, species_list )
    {
      mass_fractions.push_back(
        gsl_vector_get(
          p_abunds,
          Libnucnet__Species__getIndex( species.getNucnetSpecies() )
        ) * Libnucnet__Species__getA( species.getNucnetSpecies() )
      );
    }
    gsl_vector_free( p_abunds );
  }

  dims = new hsize_t[2];
  dims[0] = zone_list.size();
  dims[1] = species_list.size();
  dataspace = H5::DataSpace( 2, dims );
  delete[] dims;

  dataset =
    group.createDataSet(
      S_MASS_FRACTIONS,
      H5::PredType::NATIVE_DOUBLE,
      dataspace
  );

  dataset.write( mass_fractions.data(), H5::PredType::NATIVE_DOUBLE );

  //===========================================================================
  // Zone labels.
  //===========================================================================

  std::vector<zone_labels_t> zone_labels;

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    zone_labels_t my_zone_labels = {};

    strcpy(
      my_zone_labels.sLabel1,
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 )
    );

    strcpy(
      my_zone_labels.sLabel2,
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 2 )
    );

    strcpy(
      my_zone_labels.sLabel3,
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 3 )
    );

    zone_labels.push_back( my_zone_labels );

  }

  dims = new hsize_t[1];
  dims[0] = zone_list.size();
  dataspace = H5::DataSpace( 1, dims );
  delete[] dims;

  my_type = H5::CompType( sizeof( zone_labels_t ) );

  my_type.insertMember(
    S_LABEL_1, HOFFSET( zone_labels_t, sLabel1 ), string_type
  );

  my_type.insertMember(
    S_LABEL_2, HOFFSET( zone_labels_t, sLabel2 ), string_type
  );

  my_type.insertMember(
    S_LABEL_3, HOFFSET( zone_labels_t, sLabel3 ), string_type
  );

  dataset = group.createDataSet( S_ZONE_LABELS, my_type, dataspace );

  dataset.write( zone_labels.data(), my_type );

  //===========================================================================
  // Zone properties.
  //===========================================================================

  H5::Group subgroup = group.createGroup( S_ZONE_PROPERTIES );

  size_t i_zone = 0;
  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    std::vector<zone_properties_t> zone_properties;

    Libnucnet__Zone__iterateOptionalProperties(
      zone.getNucnetZone(),
      NULL,
      NULL,
      NULL,
      (Libnucnet__Zone__optional_property_iterate_function)
        populate_zone_properties,
      &zone_properties
    );

    dims = new hsize_t[1];
    dims[0] = zone_properties.size();
    dataspace = H5::DataSpace( 1, dims );
    delete[] dims;

    my_type = H5::CompType( sizeof( zone_properties_t ) );

    my_type.insertMember(
      S_NAME, HOFFSET( zone_properties_t, sName ), string_type
    );

    my_type.insertMember(
      S_TAG1, HOFFSET( zone_properties_t, sTag1 ), string_type
    );

    my_type.insertMember(
      S_TAG2, HOFFSET( zone_properties_t, sTag2 ), string_type
    );

    my_type.insertMember(
      S_VALUE, HOFFSET( zone_properties_t, sValue ), string_type
    );

    std::string s_d = boost::lexical_cast<std::string>( i_zone++ );

    dataset = subgroup.createDataSet( s_d, my_type, dataspace );

    dataset.write( zone_properties.data(), my_type );

  }

}

//##############################################################################
// populate_zone_properties().
//##############################################################################

void
populate_zone_properties(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  void * p_data
)
{

  std::vector<zone_properties_t> * p_zone_properties =
    ( std::vector<zone_properties_t> * ) p_data;
    
  zone_properties_t my_zone_properties = {};

  strcpy( my_zone_properties.sName, s_name );

  if( s_tag1 )
    strcpy( my_zone_properties.sTag1, s_tag1 );
  else
    strcpy( my_zone_properties.sTag1, "0" );
 
  if( s_tag2 )
    strcpy( my_zone_properties.sTag2, s_tag2 );
  else
    strcpy( my_zone_properties.sTag2, "0" );
 
  strcpy( my_zone_properties.sValue, s_value );

  p_zone_properties->push_back( my_zone_properties );

} 
  
//##############################################################################
// write_nuclide_data()
//##############################################################################

void
write_nuclide_data( H5::H5File& file, Libnucnet * p_nucnet )
{

  H5::DataSet dataset;
  std::vector<nuclide> nuclides;

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) )
    ); 

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    nuclide my_nuclide = {};

    strcpy(
      my_nuclide.sName,
      Libnucnet__Species__getName( species.getNucnetSpecies() )
    );
    my_nuclide.iIndex = nuclides.size();
    my_nuclide.iZ = Libnucnet__Species__getZ( species.getNucnetSpecies() );
    my_nuclide.iA = Libnucnet__Species__getA( species.getNucnetSpecies() );
    my_nuclide.dMassExcess =
      Libnucnet__Species__getMassExcess( species.getNucnetSpecies() );

    nuclides.push_back( my_nuclide );
        
  }

  //===========================================================================
  // Create dataspace.
  //===========================================================================

  hsize_t dims[1] = {species_list.size()};
  H5::DataSpace dataspace( 1, dims );

  //===========================================================================
  // Create stringtype.
  //===========================================================================

  H5::StrType string_type( H5::PredType::C_S1, I_HDF5_BUF );

  //===========================================================================
  // Create memory datatype and write to file.
  //===========================================================================

  H5::CompType my_type( sizeof( nuclide ) );

  my_type.insertMember( S_NAME, HOFFSET( nuclide, sName ), string_type );
  my_type.insertMember(
    S_INDEX, HOFFSET(nuclide, iIndex), H5::PredType::NATIVE_HSIZE
  );
  my_type.insertMember( S_Z, HOFFSET(nuclide, iZ), H5::PredType::NATIVE_UINT);
  my_type.insertMember( S_A, HOFFSET(nuclide, iA), H5::PredType::NATIVE_UINT);
  my_type.insertMember(
    S_MASS_FRACTIONS,
    HOFFSET(nuclide, dMassExcess),
    H5::PredType::NATIVE_DOUBLE
  );

  dataset = file.createDataSet( S_NUCLIDE_DATA, my_type, dataspace );

  dataset.write( nuclides.data(), my_type );

}

//##############################################################################
// group_counter()
//##############################################################################

herr_t
group_counter(
  hid_t loc_id,
  const char *name,
  const H5L_info_t *info,
  void * p_data
)
{

  H5O_info_t      infobuf;

  int * p_count = ( int * ) p_data;

  H5Oget_info_by_name( loc_id, name, &infobuf, H5P_DEFAULT );

  if( infobuf.type == H5O_TYPE_GROUP )
  {
    (*p_count)++;
  }

  return 0;

}

//##############################################################################
// dataset_counter()
//##############################################################################

herr_t
dataset_counter(
  hid_t loc_id,
  const char *name,
  const H5L_info_t *info,
  void * p_data
)
{

  H5O_info_t      infobuf;

  int * p_count = ( int * ) p_data;

  H5Oget_info_by_name( loc_id, name, &infobuf, H5P_DEFAULT );

  if( infobuf.type == H5O_TYPE_DATASET )
  {
    (*p_count)++;
  }

  return 0;

}

//##############################################################################
// get_network().
//##############################################################################

nuclide_map
get_network(
  const char * s_file
)
{

  nuclide_map my_nuclides;

  hsize_t * dims;
  H5::DataSet dataset;
  H5::DataSpace dataspace;

  H5::H5File my_file( s_file, H5F_ACC_RDONLY );

  dataset = my_file.openDataSet( S_NUCLIDE_DATA );

  //===========================================================================
  // Create stringtype.
  //===========================================================================

  H5::StrType string_type( H5::PredType::C_S1, I_HDF5_BUF );

  //===========================================================================
  // Create memory datatype.
  //===========================================================================

  H5::CompType my_type( sizeof( nuclide ) );

  my_type.insertMember( S_NAME, HOFFSET( nuclide, sName ), string_type );
  my_type.insertMember(
    S_INDEX, HOFFSET(nuclide, iIndex), H5::PredType::NATIVE_HSIZE
  );
  my_type.insertMember( S_Z, HOFFSET(nuclide, iZ), H5::PredType::NATIVE_UINT);
  my_type.insertMember( S_A, HOFFSET(nuclide, iA), H5::PredType::NATIVE_UINT);
  my_type.insertMember(
    S_MASS_FRACTIONS,
    HOFFSET(nuclide, dMassExcess),
    H5::PredType::NATIVE_DOUBLE
  );

  //===========================================================================
  // Get data.
  //===========================================================================

  dataspace = dataset.getSpace();

  dims = new hsize_t[dataspace.getSimpleExtentNdims()];
  dataspace.getSimpleExtentDims( dims );
  nuclide * nuclides = new nuclide[ dims[0] ];

  dataset.read( nuclides, my_type );

  for( size_t i = 0; i < dims[0]; i++ )
  {
    my_nuclides[nuclides[i].sName] = nuclides[i];
  }

  delete[] dims;
  delete[] nuclides;

  return my_nuclides;

}

//##############################################################################
// get_step_mass_fractions().
//##############################################################################

mass_fraction_array_t
get_step_mass_fractions(
  const char * s_file,
  const char * s_step
)
{

  mass_fraction_array_t my_mass_fractions;
  double * mass_fractions;

  hsize_t * dims;
  H5::Group group;
  H5::DataSet dataset;
  H5::DataSpace dataspace;

  H5::H5File my_file( s_file, H5F_ACC_RDONLY );

  group = my_file.openGroup( s_step );

  dataset = group.openDataSet( S_MASS_FRACTIONS );

  dataspace = dataset.getSpace();

  dims = new hsize_t[dataspace.getSimpleExtentNdims()];
  dataspace.getSimpleExtentDims( dims );
  mass_fractions = new double[ dims[0] * dims[1] ];

  dataset.read( mass_fractions, H5::PredType::NATIVE_DOUBLE );

  my_mass_fractions.resize( boost::extents[ (int) dims[0] ][ (int) dims[1] ] );

  for( size_t i = 0; i < dims[0]; i++ )
  {
    for( size_t j = 0; j < dims[1]; j++ )
    {
      my_mass_fractions[i][j] = mass_fractions[i * dims[1] + j];
    }
  }

  delete[] dims;
  delete[] mass_fractions;

  return my_mass_fractions;

}

//##############################################################################
// get_zone_labels().
//##############################################################################

zone_labels_bimap
get_zone_labels(
  const char * s_file,
  const char * s_step
)
{

  zone_labels_bimap my_zone_labels;

  hsize_t * dims;
  H5::Group group;
  H5::DataSet dataset;
  H5::DataSpace dataspace;
  H5::CompType my_type;

  H5::H5File my_file( s_file, H5F_ACC_RDONLY );

  group = my_file.openGroup( s_step );

  dataset = group.openDataSet( S_ZONE_LABELS );

  dataspace = dataset.getSpace();

  dims = new hsize_t[dataspace.getSimpleExtentNdims()];
  dataspace.getSimpleExtentDims( dims );

  //===========================================================================
  // Create stringtype.
  //===========================================================================

  H5::StrType string_type( H5::PredType::C_S1, I_HDF5_BUF );

  //===========================================================================
  // Create memory datatype and read.
  //===========================================================================

  my_type = H5::CompType( sizeof( zone_labels_t ) );

  my_type.insertMember(
    S_LABEL_1, HOFFSET( zone_labels_t, sLabel1 ), string_type
  );

  my_type.insertMember(
    S_LABEL_2, HOFFSET( zone_labels_t, sLabel2 ), string_type
  );

  my_type.insertMember(
    S_LABEL_3, HOFFSET( zone_labels_t, sLabel3 ), string_type
  );

  zone_labels_t * zone_labels = new zone_labels_t[ dims[0] ];

  dataset.read( zone_labels, my_type );

  for( size_t i = 0; i < dims[0]; i++ )
  {
    my_zone_labels.insert(
      zone_labels_bimap::value_type(
        zone_labels_tuple(
          boost::make_tuple(
            zone_labels[i].sLabel1,
            zone_labels[i].sLabel2,
            zone_labels[i].sLabel3
          )
        ),
        i
      )
    );
  }
    
  delete[] dims;
  delete[] zone_labels;

  return my_zone_labels;

}

//##############################################################################
// get_zone_properties().
//##############################################################################

zone_properties_hash
get_zone_properties(
  const char * s_file,
  const char * s_step,
  const char * s_label1,
  const char * s_label2,
  const char * s_label3
)
{

  zone_labels_bimap my_zone_labels;
  zone_properties_hash my_properties_hash;

  hsize_t * dims;
  H5::Group group, subgroup;
  H5::DataSet dataset;
  H5::DataSpace dataspace;
  H5::CompType my_type;

  H5::H5File my_file( s_file, H5F_ACC_RDONLY );

  group = my_file.openGroup( s_step );

  subgroup = group.openGroup( S_ZONE_PROPERTIES );

  my_zone_labels = get_zone_labels( s_file, s_step );

  zone_labels_bimap::index<from>::type::iterator it =
    my_zone_labels.get<from>().find(
      zone_labels_tuple(
        boost::make_tuple( s_label1, s_label2, s_label3 )
      )
    ); 

  if( it == my_zone_labels.get<from>().end() )
  {
    std::cerr << "Zone not found." << std::endl;
    exit( EXIT_FAILURE );
  }

  std::string s_d = boost::lexical_cast<std::string>( it->second );

  dataset = subgroup.openDataSet( s_d.c_str() );

  dataspace = dataset.getSpace();

  dims = new hsize_t[dataspace.getSimpleExtentNdims()];
  dataspace.getSimpleExtentDims( dims );

  //===========================================================================
  // Create stringtype.
  //===========================================================================

  H5::StrType string_type( H5::PredType::C_S1, I_HDF5_BUF );

  //===========================================================================
  // Create memory datatype and read.
  //===========================================================================

  my_type = H5::CompType( sizeof( zone_properties_t ) );

  my_type.insertMember(
    S_NAME, HOFFSET( zone_properties_t, sName ), string_type
  );

  my_type.insertMember(
    S_TAG1, HOFFSET( zone_properties_t, sTag1 ), string_type
  );

  my_type.insertMember(
    S_TAG2, HOFFSET( zone_properties_t, sTag2 ), string_type
  );

  my_type.insertMember(
    S_VALUE, HOFFSET( zone_properties_t, sValue ), string_type
  );

  zone_properties_t * zone_properties = new zone_properties_t[ dims[0] ];

  dataset.read( zone_properties, my_type );

  for( size_t i = 0; i < dims[0]; i++ )
  {
    my_properties_hash.insert(
      zone_properties_hash_entry(
        zone_properties[i].sName,
        zone_properties[i].sTag1,
        zone_properties[i].sTag2,
        zone_properties[i].sValue
      )
    );
  }

  delete[] dims;
  delete[] zone_properties;

  return my_properties_hash;

}

//##############################################################################
// count_groups_in_file().
//##############################################################################

size_t
count_groups_in_file(
  const char * s_file
)
{

  H5::H5File file( s_file, H5F_ACC_RDONLY );

  size_t i_groups = 0;

  H5Literate(
    file.getId(),
    H5_INDEX_NAME,
    H5_ITER_INC,
    NULL,
    group_counter,
    &i_groups
  );

  return i_groups;

}

//##############################################################################
// count_datasets_in_group().
//##############################################################################

size_t
count_datasets_in_group(
  const char * s_file,
  const char * s_group
)
{

  H5::H5File file( s_file, H5F_ACC_RDONLY );

  size_t i_datasets = 0;

  H5::Group group = file.createGroup( s_group );

  H5Literate(
    group.getId(),
    H5_INDEX_NAME,
    H5_ITER_INC,
    NULL,
    dataset_counter,
    &i_datasets
  );

  return i_datasets;

}

//##############################################################################
// get_step_label()
//##############################################################################

std::string
get_step_label( size_t i )
{

  boost::format fmt( "%05d" );
  fmt % i;
  return std::string( S_STEP_GROUP ) + " " + fmt.str();

}

}  // namespace hdf5

}  // namespace user
