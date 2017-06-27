////////////////////////////////////////////////////////////////////////////////
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
//! \brief Header file for useful containers.
////////////////////////////////////////////////////////////////////////////////

#ifndef CONTAINERS_H
#define CONTAINERS_H

//##############################################################################
// Includes.
//##############################################################################

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>

namespace user
{

struct from{};
struct to{};

template<typename FromType,typename ToType>
class hashed_unique_bimap
{

  public:
    struct value_type
    {
      value_type(
	const FromType& first_,
	const ToType& second_
      ):
	first(first_), second(second_)
      {}

      FromType first;
      ToType   second;

    };

    typedef boost::multi_index::multi_index_container<
      value_type,
      boost::multi_index::indexed_by<
	boost::multi_index::hashed_unique<
          boost::multi_index::tag<from>,
	  boost::multi_index::member<value_type, FromType, &value_type::first>
	>,
	boost::multi_index::hashed_unique<
          boost::multi_index::tag<to>,
	  boost::multi_index::member<value_type, ToType, &value_type::second>
	>
      >
    > type;

};

} // namespace user

#endif // CONTAINERS_H
