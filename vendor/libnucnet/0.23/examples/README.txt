////////////////////////////////////////////////////////////////////////////////
// <file type = "public">
//   <license>
//      This is free software; you can redistribute it and/or modify it
//      under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This software is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      Please see the src/README.txt file in this distribution for more
//      information.
//   </license>
//
//   <description>
//     <abstract>
//       README file for the examples/ directory.
//     </abstract>
//     <keywords>
//       README, libnucnet, examples, code, source
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="dcadams" start_date="2006/08/09" />
//       <author userid="mbradle" start_date="2006/08/09" />
//     </current>
//     <previous>
//     </previous>
//   </authors>
//
// </file>
////////////////////////////////////////////////////////////////////////////////

################################################################################
# Overview of directory contents.
################################################################################

./README.txt

  Overview of directory contents, compatibility notes, and compilation
  procedure.

*.c, *.h

  These are the C example codes in the distribution.  See the individual files
  for explanations of these example codes.

coul_corr.c coul_corr.h

  These are example codes for computing and applying coulomb corrections
  to NSE factors.

remove_duplicate.c remove_duplicate.h

  These are example codes for removing duplicate reactions from a network.

screen.c screen.h

  These are example codes for computing and applying screening correction
  factors to reactions.

user_rate_functions.c user_rate_functions.h

  These are example codes for computing and applying user-supplied reaction
  rate functions.

Makefile

  This is the C Makefile for the example codes in the distribution.

################################################################################
# Compatibility notes for the general user.
################################################################################

We've tested with gnu gcc 4.3.2 and Cygwin 1.5.24-2.  We have also tested with
the Intel compilers icc and icpc.  For these latter compilers, please use the
warning suppression flags -wd9 -wd981 -wd1419 -wd10148 -wd10156.

################################################################################
# Example dependencies.
################################################################################

The data for the libnucnet example codes are most easily downloaded using
wget. To check whether wget is installed on a particular unix or linux system,
type:

wget --version

To install wget, please see http://www.gnu.org/software/wget/

################################################################################
# Download the data.
################################################################################

If you have wget installed, type:

make data

Alternatively, carry out the following steps:

1)  Create a data_pub directory:

cd ..
mkdir data_pub

2)  Download the data tarball from

http://libnucnet.sf.net/data_pub/2011-05-28/data.tar.gz

and place it in the data_pub directory.

3)  Uncompress the data:

cd data_pub
gunzip data.tar.gz
tar xvf data.tar 

4)  Return to the examples directory:

cd ../examples

You can remove the data by typing:

make clean_data

################################################################################
# Compilation procedure.
################################################################################

Steps to compile examples:

1) Edit Makefile so that MATRIXSRCDIR points to the directory that contains
your local installation of the Webnucleo module wn_matrix.  To obtain that
module, please visit http://www.webnucleo.org.  Other parameters in the
Makefile (GC, LIBNUCNETDIR, MATRIXSRCDIR, VALGRIND, and PROFILE) should not need
changing, although they can be if desired.

2) Run make:

make create_nuc_collection

3) Run the example:

./create_nuc_collection

The other example codes are compiled similarly.  Just replace
create_nuc_collection in the above instructions with the example code you
wish to compile.  Or you can run make all.  See http://webnucleo.org for more
information on compiling and running the examples.
