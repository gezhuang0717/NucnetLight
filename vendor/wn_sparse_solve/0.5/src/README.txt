/*//////////////////////////////////////////////////////////////////////////////
// <file type = "public">
//
//   <license>
//      Copyright (c) 2006-2013 Clemson University.
//
//      This distribution contains the source code
//      for the Clemson Webnucleo group's
//      wn_sparse_solve module, originally developed by Bradley S. Meyer.
//      and Lih-Sin The.
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
//      All wn_sparse_solve documentation is free documentation; permission is
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
//       README file for the wn_sparse_solve src/ directory.
//     </abstract>
//     <keywords>
//       README, wn_sparse_solve, code, source
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2010/05/18" />
//       <author userid="tlihsin" start_date="2010/05/18" />
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

  Overview of directory contents and compatibility notes.

WnSparseSolve.c

  The wn_sparse_solve C code in the distribution.

WnSparseSolve.h

  The header file for WnSparseSolve.c

blas1.f, exppro.f, itaux.f, phipro.f

  Unsupported SPARSKIT2 routines modified for use with wn_sparse_solve.

################################################################################
# Dependencies.
################################################################################

wn_sparse_solve codes require wn_matrix version 0.14 or later.  Please see
http://www.webnucleo.org for downloading information.  wn_sparse_solve
routines also require SPARSKIT2.  To download, see
http://www-users.cs.umn.edu/~saad/software/SPARSKIT/sparskit.html.

################################################################################
# Compatibility notes for the general user.
################################################################################

We've tested with gnu gcc 4.5 and Cygwin 1.5.24-2.

