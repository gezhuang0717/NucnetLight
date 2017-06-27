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
//       README file for the wn_sparse_solve examples/ directory.
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

  Overview of directory contents, compatibility notes, and compilation
  procedure.

*.c/h

  These are the C example codes in the distribution.  See the individual files
  for explanations of these example codes.

Makefile

  This is the Makefile for the example codes in the distribution.

Makefile.inc

  This is the Makefile for handling sparskit.

my_convergence_tester.c, my_convergence_tester.h

  Code to demonstrate how to supply a user-defined convergence tester.

examples.xml

  An XML file giving information about the example codes.

ilu_solvers.c, ilu_solvers.h

  Code to demonstrate how to interface wn_sparse_solve code with SPARSIT2's
  incomplete lu decomposition preconditioner set and solves.

print_out.c

  Routine to print out matrix and vector data for the examples.

regression_test.sh

  A regression test script.

solution_check.c

  Routine to check iterative solutions.

################################################################################
# Compatibility notes for the general user.
################################################################################

We've tested with gnu gcc 4.5 and Cygwin 1.5.24-2.

################################################################################
# Compilation procedure.
################################################################################

Steps to compile example1:

1) Obtain SPARSKIT2.  If you have wget install, simply type:

>make sparskit

If you don't have wget, visit http://www-users.cs.umn.edu/~saad/software/ and
install SPARSKIT2 in the ..  directory.

2) Edit Makefile so that MATRIXDIR points to the directory that contains your
local installation of the Webnucleo module wn_matrix.  To obtain that module,
please visit http://www.webnucleo.org.
Note that if you have installed wn_matrix, wn_sparse_solve, and SPARSKIT2
according to the download tutorials at http://www.webnucleo.org,
you should not need to modify the Makefile since MATRIXDIR and SPARSKITDIR
will already point to the correct directories.

3) Run make:

>make basic_solve_sparse

To compile all the examples, simply type:

>make all

4) Run the example:

>./basic_solve_sparse ../data_pub/matrix.xml ../data_pub/rhs.xml gmres 10 1.e-4 1.e-4

You can simply type

>./basic_solve_sparse

for a usage statement.  The other example codes are compiled similarly.  Just
replace basic_solve_sparse in the above instructions with the example code you
wish to compile.  Or you can run make all.  See http://webnucleo.org for more
information on compiling and running the examples.
