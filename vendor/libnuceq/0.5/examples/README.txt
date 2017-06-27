////////////////////////////////////////////////////////////////////////////////
// <file type="public">
//
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
//       README, libnuceq, examples, code, source
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2010/12/12" />
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

We've tested with gnu gcc 4.3.2 and Cygwin 1.5.24-2.  We have also tested with
the Intel compilers icc and icpc.  For these latter compilers, please use the
warning suppression flags -wd9 -wd981 -wd1419 -wd10148 -wd10156.

################################################################################
# Compilation procedure.
################################################################################

1) Edit Makefile so that MATRIXSRCDIR, LIBSTATMECHDIR, and
LIBNUCNETDIR point to the directory that contains your local installation
of these Webnucleo modules.  To obtain those modules, please visit
http://www.webnucleo.org.  Other parameters in the
Makefile (GC, VALGRIND, and PROFILE) should not need changing, although
they can be if desired.

2) Run make:

>make compute_wse

3) Run the example:

>./compute_wse ../data_pub/example_nuc.xml 4. 1.e9 0.

The other example codes are compiled similarly.  Just replace compute_wse in the
above instructions with the example code you wish to compile.  Or you can run
make all.  See http://webnucleo.org for more information on compiling and
running the examples.
