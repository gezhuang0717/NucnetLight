////////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//      This directory contains the example source code
//      for the Clemson Webnucleo group's
//      wn_matrix module, originally developed by David C. Adams and 
//      Bradley S. Meyer.
//
//      This is free software; you can redistribute it and/or modify it
//      under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This software is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//   </license>
//
//   <description>
//     <abstract>
//       README file for the examples/ directory.
//     </abstract>
//     <keywords>
//       README, wn_matrix, examples, code, source
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

*.c

  These are the C example codes in the distribution.  See the individual files
  for explanations of these example codes.

examples.xml

  An xml file for use with regression tests.

Makefile

  This is the C Makefile for the example codes in the distribution.

regression_test.sh

  This is a script that generates the regression tester.

################################################################################
# Compatibility notes for the general user.
################################################################################

We've tested with gnu gcc 3.4.4 and gnu gcc 4.3.2 and Cygwin 1.5.24-2.
We have also tested with the Intel compilers icc and icpc.  For these latter
compilers, please use the warning suppression flags -wd9 -wd981 -wd10148
-wd10156 -wd1419.

################################################################################
# Compilation procedure.
################################################################################

Steps to compile example 1:

1) Edit Makefile so that MATRIXSRCDIR points to the directory that contains
your local installation of the Webnucleo module wn_matrix.  To obtain that
module, please visit http://www.webnucleo.org.  Note that if you have
not altered the directory structure of the distribution, you should not need
to perform this step since MATRIXSRCDIR should already point to the correct
directory.

2) Run make:

>make create_matrix

3) Run the example:

>./create_matrix

The other example codes are compiled similarly.  Just replace create_matrix in
the above instructions with the example code you wish to compile.  Or you can
run make all.  See http://webnucleo.org for more information on compiling and
running the examples.
