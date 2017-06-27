////////////////////////////////////////////////////////////////////////////////
// <file type="public">
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
//       README file for the src/examples/ directory.
//     </abstract>
//     <keywords>
//       README, libstatmech, examples, code, source
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="tyu" start_date="2008/05/20" />
//       <author userid="mbradle" start_date="2008/05/20" />
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

*.c

  These are the C example codes in the distribution.  See the individual files
  for explanations of these example codes.

Makefile

  This is the C Makefile for the example codes in the distribution.

################################################################################
# Compatibility notes for the general user.
################################################################################

We've tested with gnu gcc 3.4.4 and Cygwin 1.5.24-2.  We have also tested with
the Intel compilers icc and icpc.  For these latter compilers, please use the
warning suppression flags -wd9 -wd981 -wd10148 -wd10156 -wd1419.

################################################################################
# Compilation procedure.
################################################################################

Steps to compile example 1:

1) Parameters in the Makefile (GC, STATMECHSRCDIR, and VALGRIND)
should not need changing, although they can be if desired.

2) Run make:

>make create_fermion

3) Run the example:

>create_fermion electron 0.511 2 -1

The other example codes are compiled similarly.  Just replace example1 in the
above instructions with the example code you wish to compile.  Or you can run
make all.  See http://webnucleo.org for more information on compiling and
running the examples.
