/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//      Copyright (c) 2006-2012 Clemson University.
//
//      This distribution directory contains the source code
//      for the Clemson Webnucleo group's
//      wn_matrix module, originally developed by David Adams and 
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
//
//      You should have received a copy of the GNU General Public License
//      along with this software (please see the "gnu_gpl.txt" file in the doc/
//      directory of this distribution); if not, write to the Free Software
//      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
//      USA
//
//      All wn_matrix documentation is free documentation; permission is
//      granted to copy, distribute and/or modify the documentation under the
//      terms of the GNU Free Documentation License, Version 1.2 or any later
//      version published by the Free Software Foundation.  A copy of the
//      license is included in the file "gnu_fdl.txt" in the doc/ directory
//      in this distribution.
//   </license>
//
//   <description>
//     <abstract>
//       README file for the wn_matrix src/ directory.  The src/ directory
//       contains WnMatrix.h, the wn_matrix header file, and WnMatrix.c, the
//       collection of wn_matrix codes.
//     </abstract>
//     <keywords>
//       README, wn_matrix, code, source
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="dcadams" start_date="2006/09/20" />
//       <author userid="mbradle" start_date="2006/09/20" />
//     </current>
//     <previous>
//     </previous>
//   </authors>
//
// </file>
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/

################################################################################
# Overview of directory contents.
################################################################################

./README.txt

  Overview of directory contents and compatibility notes.

WnMatrix.c

  The wn_matrix C code in the distribution.

WnMatrix.h

  The header file for WnMatrix.c

################################################################################
# Dependencies.
################################################################################

wn_matrix codes require libxml2 version 2.6.18 or later.  Earlier versions
of libxml2 may also work, but we haven't tested them.  To check whether
libxml2 is installed on a particular unix or linux system, type:

xml2-config --version

To install libxml2, please see http://www.libxml.org.

wn_matrix codes also require gsl version 1.9 or later.  To check whether
gsl is installed on a particular unix or linux system, type:

gsl-config --version

To install gsl, please see http://www.gnu.org/software/gsl/.

################################################################################
# Compatibility notes for the general user.
################################################################################

We've tested with gnu gcc 4.3.4 and gnu gcc 4.3.2 and Cygwin 1.7.5.
We have also tested with the Intel compilers icc and icpc.  For these latter
compilers, please use the warning suppression flags -wd9 -wd981 -wd10148
-wd10156 -wd1419.  For example, in the Makefile, use the assignment:

GC=icc -wd9 -wd981 -wd10148 -wd10156 -wd1419

or

GC=icpc -wd9 -wd981 -wd10148 -wd10156 -wd1419.


