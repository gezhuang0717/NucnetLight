/*//////////////////////////////////////////////////////////////////////////////
// <file type = "public">
//
//   <license>
//      Copyright (c) 2010-2015 Clemson University.
//
//      This distribution contains the source code
//      for the Clemson Webnucleo group's
//      libnuceq module, originally developed by Bradley S. Meyer.
//      and Tianhong Yu.
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
//      All libnuceq documentation is free documentation; permission is
//      granted to copy, distribute and/or modify the documentation under the
//      terms of the GNU Free Documentation License, Version 1.2 or any later
//      version published by the Free Software Foundation; with the Invariant
//      Sections, the Front-Cover Texts, and the Back-Cover Texts identified
//      in the documentation itself.  A copy of the license is included in the
//      file "gnu_fdl_v1.2.txt" in the doc/ directory in this distribution.
//   </license>
//
//   <description>
//     <abstract>
//       README file for the libnuceq src/ directory.
//     </abstract>
//     <keywords>
//       README, libnuceq, code, source
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2010/05/18" />
//       <author userid="tyu" start_date="2010/05/18" />
//     </current>
//     <previous>
//     </previous>
//   </authors>
//
// </file>
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/

################################################################################
# libnuceq source.
################################################################################

This distribution contains release 0.3.  For further details, see
http://www.webnucleo.org.

################################################################################
# Overview of directory contents.
################################################################################

./README.txt

  Overview of directory contents and compatibility notes.

Libnuceq.c

  The libnuceq C code in the distribution for handling equilibria.

Libnuceq.h

  The header file for Libnuceq.c

################################################################################
# Dependencies.
################################################################################

libnuceq codes require libxml2 version 2.6.18 or later.  Earlier versions
of libxml2 may also work, but we haven't tested them.  To check whether
libxml2 is installed on a particular unix or linux system, type:

xml2-config --version

To install libxml2, please see http://www.libxml.org.

libnuceq codes require gsl version 1.9 or later.  Earlier versions
of gsl may also work, but we haven't tested them.  To check whether
gsl is installed on a particular unix or linux system, type:

gsl-config --version

To install gsl, please see http://www.gnu.org/software/gsl/

This libnuceq distribution requires wn_matrix version 0.9 or later.  If you
have not installed wn_matrix, please see

http://www.webnucleo.org/home/modules/wn_matrix/

This libnuceq distribution requires libnucnet version 0.9 or later.  If you
have not installed libnucnet, please see

http://www.webnucleo.org/home/modules/libnucnet/

This libnuceq distribution requires libstatmech version 0.5 or later.  If you
have not installed libstatmech, please see

http://www.webnucleo.org/home/modules/libstatmech/

################################################################################
# Compatibility notes for the general user.
################################################################################

We've tested with gnu gcc 3.4.4 and Cygwin 1.5.24-2.

