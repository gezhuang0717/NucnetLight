#//////////////////////////////////////////////////////////////////////////////
# This file was originally written by Bradley S. Meyer.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#//////////////////////////////////////////////////////////////////////////////
 
#///////////////////////////////////////////////////////////////////////////////
#//! \file
#//! \brief Bash include script to generate examples.
#///////////////////////////////////////////////////////////////////////////////

cd ../..
svn up
cd $MY_HOME
if make $TARGET
then
  echo Make succeeded.
else
  make clean
  if make $TARGET
  then
    echo Make succeeded.
  else
    make clean_all
    if make $TARGET
    then
      echo Make succeeded.
    else
      echo Make failed.  Contact project maintainer.
    fi
  fi
fi
