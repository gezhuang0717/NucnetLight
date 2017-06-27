/*//////////////////////////////////////////////////////////////////////////////
// <file type = "public">
//   <license>
//      Copyright (c) 2006-2011 Clemson University.
//
//      This file contains the source code header file
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
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#ifndef WN_MATRIX_H
#define WN_MATRIX_H

/*##############################################################################
// Use extern "C" for C++ compilers.
//############################################################################*/

#ifdef __cplusplus
extern "C"
{
#endif

/*##############################################################################
// Includes.
//############################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>

#include <libxml/hash.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxml/xmlschemas.h>
#include <libxml/xmlschemastypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

/*##############################################################################
// Parameter defines.
//############################################################################*/

#define INDEX_MAX ULONG_MAX

/*##############################################################################
// String defines.
//############################################################################*/

#define WNMATRIX_VERSION "0.13"
#define XSD_VERSION "2011-05-09"
#define XML_VERSION "1.0"

#define MATRIX "matrix"

#define COO_MATRIX "coordinate_" MATRIX
#define CSR_MATRIX "csr_" MATRIX
#define YALE_MATRIX "yale_" MATRIX
#define VECTOR "vector"

#define COO_PREFIX "coo"
#define CSR_PREFIX "csr"
#define YALE_PREFIX "yale"
#define VECTOR_PREFIX "vec"

#define COO_MATRIX_WITH_PREFIX COO_PREFIX ":" MATRIX
#define CSR_MATRIX_WITH_PREFIX CSR_PREFIX ":" MATRIX
#define YALE_MATRIX_WITH_PREFIX YALE_PREFIX ":" MATRIX
#define VECTOR_WITH_PREFIX VECTOR_PREFIX ":" VECTOR

#define VALUE_FORMAT "%g"

#define WNMATRIX_COLUMN "column"
#define WNMATRIX_COO_DEC "xmlns:" COO_PREFIX
#define WNMATRIX_CSR_DEC "xmlns:" CSR_PREFIX
#define WNMATRIX_ELEMENT "element"
#define WNMATRIX_ELEMENTS "elements"
#define WNMATRIX_IJA "ija"
#define WNMATRIX_ROW "row"
#define WNMATRIX_ROW_POINTER "row_pointer"
#define WNMATRIX_ROW_POINTERS "row_pointers"
#define WNMATRIX_SA "sa"
#define WNMATRIX_VALUE "value"
#define WNMATRIX_XPATH "//element"
#define WNMATRIX_VECTOR_DEC "xmlns:" VECTOR_PREFIX
#define WNMATRIX_VECTOR_ELEMENT "element"
#define WNMATRIX_VECTOR_XPATH "//element"
#define WNMATRIX_YALE_DEC "xmlns:" YALE_PREFIX

#define XSI_DEC "xmlns:xsi"
#define XSI_LOCATION "xsi:schemaLocation"

/*##############################################################################
// WnMatrix defines.
//############################################################################*/

#define WNMATRIX_BUF_SIZE 32

#define W3C__NAMESPACE "http://www.w3.org/2001/XMLSchema-instance"

#define WNMATRIX__PREFIX \
  "http://wnmatrix.sf.net/xsd_pub/" XSD_VERSION "/"

#define COO__NAMESPACE COO_MATRIX "/"
#define COO__SCHEMA COO_MATRIX ".xsd"

#define WNMATRIX__COO__SCHEMA WNMATRIX__PREFIX COO__SCHEMA
#define WNMATRIX__COO__NAMESPACE WNMATRIX__PREFIX COO__NAMESPACE
#define WNMATRIX__COO__SCHEMALOCATION1 WNMATRIX__COO__NAMESPACE " "
#define WNMATRIX__COO__SCHEMALOCATION \
  WNMATRIX__COO__SCHEMALOCATION1 WNMATRIX__COO__SCHEMA

#define CSR__NAMESPACE CSR_MATRIX "/"
#define CSR__SCHEMA CSR_MATRIX ".xsd"

#define WNMATRIX__CSR__SCHEMA WNMATRIX__PREFIX CSR__SCHEMA
#define WNMATRIX__CSR__NAMESPACE WNMATRIX__PREFIX CSR__NAMESPACE
#define WNMATRIX__CSR__SCHEMALOCATION1 WNMATRIX__CSR__NAMESPACE " "
#define WNMATRIX__CSR__SCHEMALOCATION \
  WNMATRIX__CSR__SCHEMALOCATION1 WNMATRIX__CSR__SCHEMA

#define YALE__NAMESPACE YALE_MATRIX "/"
#define YALE__SCHEMA YALE_MATRIX ".xsd"

#define WNMATRIX__YALE__SCHEMA WNMATRIX__PREFIX YALE__SCHEMA
#define WNMATRIX__YALE__NAMESPACE WNMATRIX__PREFIX YALE__NAMESPACE
#define WNMATRIX__YALE__SCHEMALOCATION1 WNMATRIX__YALE__NAMESPACE " "
#define WNMATRIX__YALE__SCHEMALOCATION \
  WNMATRIX__YALE__SCHEMALOCATION1 WNMATRIX__YALE__SCHEMA

#define VECTOR__NAMESPACE VECTOR "/"
#define VECTOR__SCHEMA VECTOR ".xsd"

#define WNMATRIX__VECTOR__SCHEMA WNMATRIX__PREFIX VECTOR__SCHEMA
#define WNMATRIX__VECTOR__NAMESPACE WNMATRIX__PREFIX VECTOR__NAMESPACE
#define WNMATRIX__VECTOR__SCHEMALOCATION1 WNMATRIX__VECTOR__NAMESPACE " "
#define WNMATRIX__VECTOR__SCHEMALOCATION \
  WNMATRIX__VECTOR__SCHEMALOCATION1 WNMATRIX__VECTOR__SCHEMA

/*##############################################################################
// WnMatrix defines.
//############################################################################*/

/* Uncomment following line to debug */ 
/* #define DEBUG */

#ifndef DEBUG

#define WNMATRIX__ERROR( s_reason ) \
       do { \
         fprintf( stderr, "ERROR: %s in %s on line %d\n", \
           s_reason, __FILE__, __LINE__ \
         ); \
         exit( EXIT_FAILURE ); \
       } while (0)

#else

#define WNMATRIX__ERROR( s_reason ) \
       do { \
         fprintf( stderr, "ERROR: %s in %s on line %d\n", \
           s_reason, __FILE__, __LINE__ \
         ); \
         abort(); \
       } while (0)

#endif

/*##############################################################################
// WnMatrix structures.
//############################################################################*/

typedef struct WnMatrix__Element {
  size_t iRow;
  size_t iCol;
  double dValue;
} WnMatrix__Element;

typedef struct WnMatrix__Elements {
  size_t iCount;
  WnMatrix__Element **pData;
} WnMatrix__Elements;

/*##############################################################################
// <class name="WnMatrix__Line">
//
//   <description>
//     <abstract>
//       WnMatrix__Line is a structure for storing and managing row or column
//       data in a WnMatrix.
//     </abstract>
//     <keywords>
//       sparse, matrix, webnucleo
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2007/06/17"/>
//     </current>
//   </authors>
//
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//
// </class>
//############################################################################*/

typedef struct WnMatrix__Line {
  size_t iCount;
  size_t *a1;
  double       *a2;
} WnMatrix__Line;

/*##############################################################################
// <class name="WnMatrix">
//
//   <description>
//     <abstract>
//       WnMatrix is a structure for storing and managing sparse matrix data.
//       It provides functionality for easily adding and removing elements.
//     </abstract>
//     <keywords>
//       sparse, matrix, webnucleo
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="dcadams" start_date="2005/10/17"/>
//       <author userid="mbradle" start_date="2005/10/17"/>
//     </current>
//   </authors>
//
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//
// </class>
//############################################################################*/

typedef struct WnMatrix {
  xmlHashTablePtr pMatrixHash;
  size_t iRows;
  size_t iCols;
#ifdef WN_USE_MATRIX_LOOKUP
  xmlChar * pLookup;
  size_t iMax;
#endif
} WnMatrix;

/*##############################################################################
// <class name="WnMatrix__Arrow">
//
//   <description>
//     <abstract>
//       WnMatrix__Arrow is a structure for storing and managing sparse matrix
//       data in arrow matrix format, that is, a matrix with non-zero
//       entries in a central band and in the last several columns and
//       last several rows.
//     </abstract>
//     <keywords>
//       sparse, matrix, webnucleo
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2009/02/26"/>
//     </current>
//   </authors>
//
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//
// </class>
//############################################################################*/

typedef struct WnMatrix__Arrow {
  size_t iBand;
  size_t iBandWidth;
  size_t iRows; 
  size_t iWingWidth;
  size_t iCurrentRow;
  double **a;
  double **b;
  double **c;
  double **d;
} WnMatrix__Arrow;
    
/*##############################################################################
// <class name="WnMatrix__Coo">
//
//   <description>
//     <abstract>
//       WnMatrix__Coo is a structure for storing and managing sparse matrix
//       data in coordinate matrix format.
//     </abstract>
//     <keywords>
//       sparse, matrix, webnucleo
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2007/09/03"/>
//     </current>
//   </authors>
//
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//
// </class>
//############################################################################*/

typedef struct WnMatrix__Coo {
  size_t iCount;
  size_t *aRow;
  size_t *aCol;
  double *aVal;
} WnMatrix__Coo;

/*##############################################################################
// <class name="WnMatrix__Csr">
//
//   <description>
//     <abstract>
//       WnMatrix__Csr is a structure for storing and managing sparse matrix
//       data in compressed sparse row matrix format.
//     </abstract>
//     <keywords>
//       sparse, matrix, webnucleo
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2007/09/03"/>
//     </current>
//   </authors>
//
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//
// </class>
//############################################################################*/

typedef struct WnMatrix__Csr {
  size_t iCount;
  size_t iRowCount;
  size_t *aRowptr;
  size_t *aCol;
  double *aVal;
} WnMatrix__Csr;

/*##############################################################################
// <class name="WnMatrix__Yale">
//
//   <description>
//     <abstract>
//       WnMatrix__Yale is a structure for storing and managing sparse matrix
//       data in Yale sparse matrix format.
//     </abstract>
//     <keywords>
//       sparse, matrix, webnucleo, Yale
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2007/09/28"/>
//     </current>
//   </authors>
//
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//
// </class>
//############################################################################*/

typedef struct WnMatrix__Yale {
  size_t *aIja;
  double *aVal;
} WnMatrix__Yale;

/*##############################################################################
// API routines.
//############################################################################*/

/*##############################################################################
// <routine name="WnMatrix__new()">
//
//   <description>
//     <abstract>
//       Initializes a WnMatrix structure.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix *WnMatrix__new(
//   size_t i_rows,
//   size_t i_columns
// );
//     </calling_sequence>
//
//     <param
//       name="i_rows"
//       kind="in,positional,required"
//       doc="row"
//     >
//       A size_t integer containing the number of rows for the matrix.
//     </param>
//     <param
//       name="i_columns"
//       kind="in,positional,required"
//       doc="col"
//     >
//       A size_t integer containing the number of columns for the
//       matrix.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="pre" id="row">
//       i_rows > 0.
//     </doc>
//     <doc kind="pre" id="col">
//       i_columns > 0.
//     </doc>
//
//     <doc kind="post" id="result">
//       The routine returns a pointer to a struct WnMatrix containing the 
//       initialized matrix. The struct is initialized to contain the number
//       of rows and columns specified by the input parameters.  If the
//       routine is unable to allocate memory for the structure, the routine
//       returns NULL.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Initialize a matrix with 50 rows and 100 columns and return
//         a reference p_matrix to it:
//       </synopsis>
//       <code>
// WnMatrix *p_matrix;
// p_matrix = WnMatrix__new( 50, 100 );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix *WnMatrix__new( size_t, size_t );

/*##############################################################################
// <routine name="WnMatrix__assignElement()">
//
//   <description>
//     <abstract>
//       Assigns an element to a WnMatrix structure.  If an
//       element already exists at that (row,col), the new value is added to
//       the existing value.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, assign, element 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void WnMatrix__assignElement(
//   WnMatrix *self, size_t i_row,
//   size_t i_col, double d_val
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param
//       name="i_row"
//       kind="in,positional,required"
//       doc="row"
//     >
//       A size_t integer containing the row index for the assigned 
//       element.
//     </param>
//     <param
//       name="i_col"
//       kind="in,positional,required"
//       doc="col"
//     >
//       A size_t integer containing the column index for the
//       assigned element.
//     </param>
//     <param 
//       name="d_val"
//       kind="in,positional,required"
//       doc="val"
//     >
//       A double containing the value for the assigned element.
//     </param>
//     <param 
//       kind="return" 
//       doc="matrix_out" 
//     />
//
//     <doc kind="pre" id="row">
//       i_row > 0 and i_row <= number of rows in self.
//     </doc>
//     <doc kind="pre" id="col">
//       i_col > 0 and i_col <= number of columns in self.
//     </doc>
//
//     <doc kind="post" id="matrix_out">
//       If a value did not previously exist at i_row and i_col in 
//       self, the value there is now d_val.
//       Otherwise, the value there is now the sum of the previous value and 
//       d_val.  If the input matrix pointer or the row or column is
//       invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Assign the value 3.5 to row 1, column 2 of the matrix stored in
//         WnMatrix * p_matrix:
//       </synopsis>
//       <code>
// size_t i_row = 1;
// size_t i_col = 2;
// double d_val = 3.5;
// WnMatrix__assignElement( p_matrix, i_row, i_col, d_val );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int WnMatrix__assignElement(
  WnMatrix *, size_t, size_t, double
);

/*##############################################################################
// <routine name="WnMatrix__updateElement()">
//
//   <description>
//     <abstract>
//       Updates an element in a WnMatrix structure.  If an element already
//       exists at that (row,col), it is replaced by the new value.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, assign, element 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int
// WnMatrix__updateElement(
//   WnMatrix *self,
//   size_t i_row,
//   size_t i_col,
//   double d_val
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param
//       name="i_row"
//       kind="in,positional,required"
//       doc="row"
//     >
//       A size_t integer containing the row index for the assigned 
//       element.
//     </param>
//     <param
//       name="i_col"
//       kind="in,positional,required"
//       doc="col"
//     >
//       A size_t integer containing the column index for the
//       assigned element.
//     </param>
//     <param 
//       name="d_val"
//       kind="in,positional,required"
//       doc="val"
//     >
//       A double containing the value for the assigned element.
//     </param>
//     <param 
//       kind="return" 
//       doc="matrix_out" 
//     />
//
//     <doc kind="pre" id="row">
//       i_row > 0 and i_row <= number of rows in self.
//     </doc>
//     <doc kind="pre" id="col">
//       i_col > 0 and i_col <= number of columns in self.
//     </doc>
//
//     <doc kind="post" id="matrix_out">
//       Routine returns 0 if the value at i_row and i_col has been
//       successful updated to the input value and -1 if not.
//       If the input matrix pointer or the row or column is
//       invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Update the value at row 1, column 2 of WnMatrix * p_matrix to 3.5:
//       </synopsis>
//       <code>
// WnMatrix__updateElement( p_matrix, 1L, 2L, 3.5 );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
WnMatrix__updateElement(
  WnMatrix *, size_t, size_t, double
);

/*##############################################################################
// <routine name="WnMatrix__getElement()">
//
//   <description>
//     <abstract>
//        Retrieves the value of the specified matrix element.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, element, value 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double WnMatrix__getElement(
//   WnMatrix *self, size_t i_row, size_t i_col
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="i_row" 
//       kind="in,positional,required" 
//       doc="row" 
//     >
//       A size_t integer containing the row index of the
//       element to be
//       retrieved.
//     </param>
//     <param 
//       name="i_col" 
//       kind="in,positional,required" 
//       doc="col" 
//     >
//       A size_t integer containing the column index of the
//       element to be retrieved.
//     </param>
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="pre" id="row">
//       i_row > 0 and i_row <= number of rows in self.
//     </doc>
//     <doc kind="pre" id="col">
//       i_col > 0 and i_col <= number of columns in self.
//     </doc>
//
//     <doc kind="post" id="result">
//       Routine returns a double that is the value of the retrieved element. 
//       If the input matrix pointer or the row or column is
//       invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug">
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the value val of the element at row 1, column 2 from the
//         WnMatrix * p_my_matrix:
//       </synopsis>
//       <code>
// i_row = 1;
// i_col = 2;
// double val;
// val = WnMatrix__getElement( p_my_matrix, i_row, i_col );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double WnMatrix__getElement( WnMatrix *, size_t, size_t );

/*##############################################################################
// <routine name="WnMatrix__removeElement()">
//
//   <description>
//     <abstract>
//        Removes the specified matrix element from the matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, element, remove 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int
// WnMatrix__removeElement(
//   WnMatrix *self, size_t i_row, size_t i_col
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="i_row" 
//       kind="in,positional,required" 
//       doc="row" 
//     >
//       A size_t integer containing the row index of the
//       element to be removed.
//     </param>
//     <param 
//       name="i_col" 
//       kind="in,positional,required" 
//       doc="col" 
//     >
//       A size_t integer containing the column index of the
//       element to be removed.
//     </param>
//     <param 
//       name="self"
//       kind="out,positional,required" 
//       doc="matrix_out"
//     />
//
//     <doc kind="pre" id="row">
//       i_row > 0 and i_row <= number of rows in self.
//     </doc>
//     <doc kind="pre" id="col">
//       i_col > 0 and i_col <= number of rows in self.
//     </doc>
//
//     <doc kind="post" id="matrix_out">
//       On return matrix longer contains the specified element.  If removal
//       is successful, routine returns 0.  If removal is not successful,
//       or if entry not found, routine returns -1.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Remove the element at row 1, column 2 from the WnMatrix *p_matrix:
//       </synopsis>
//
//       <code>
// WnMatrix__removeElement( p_matrix, 1, 2 );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int WnMatrix__removeElement( WnMatrix *, size_t, size_t );

/*##############################################################################
// <routine name="WnMatrix__free()">
//
//   <description>
//     <abstract>
//       Frees the memory allocated for a  WnMatrix structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, remove, free, clear 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void WnMatrix__free( WnMatrix *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="self"
//       kind="out,positional,required" 
//       doc="matrix_out"
//     />
//
//     <doc kind="post" id="matrix_out">
//       On successful return, the matrix and its elements have all
//       been cleared and freed.  The matrix pointer self is no longer
//       available for
//       use--it must be reallocated using WnMatrix__new().
//        
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Delete WnMatrix *p_matrix and free the allocated memory:
//       </synopsis>
//
//       <code>
// WnMatrix__free( p_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void WnMatrix__free( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__clear()">
//
//   <description>
//     <abstract>
//       Zeroes all the entries in a WnMatrix structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, remove, delete, clear 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void WnMatrix__clear( WnMatrix *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="self"
//       kind="out,positional,required" 
//       doc="matrix_out"
//     />
//
//     <doc kind="post" id="matrix_out">
//       On successful return, all elements of the matrix are zero.  The
//       matrix is still available for use.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Zero out all the entries in *p_matrix:
//       </synopsis>
//
//       <code>
// WnMatrix__clear( p_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void WnMatrix__clear( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__writeMatrixToAsciiFile()">
//
//   <description>
//     <abstract>
//       Writes out the elements of a matrix stored in a WnMatrix
//       structure larger than a cutoff value to a file.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, write, file, elements
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int WnMatrix__writeMatrixToAsciiFile(
//   WnMatrix *self, char * s_ascii_filename, double d_cutoff
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param
//       name="s_ascii_filename"
//       kind="in,positional,required"
//     >
//       A string giving the name of the output ascii file.
//     </param>
//     <param 
//       name="d_cutoff"
//       kind="in,positional,required" 
//     >
//       A double giving the threshold for not writing out an element to
//       the file.
//     </param>
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       If the routine is successful, it returns 1 itself and the matrix in 
//       coordinate form in the file with the name s_ascii_filename.  If
//       the routine is unsuccessful, it returns 0.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Write out all matrix elements in p_matrix to the file my_matrix.dat:
//       </synopsis>
//
//       <code>
// if(
//   WnMatrix__writeMatrixToAsciiFile(
//     p_matrix, "my_matrix.dat", 0.
//   ) == 1
// ){
//    printf("Successfully wrote matrix.\n");
// }
//       </code>
//     </doc>
//
//     <doc kind="example" id="example2">
//       <synopsis>
//         Write out all matrix elements in p_matrix with absolute value
//         greater than 1.e-6 to the file my_matrix.dat:
//       </synopsis>
//
//       <code>
// if(
//   WnMatrix__writeMatrixToAsciiFile(
//     p_matrix, "my_matrix.dat", 1.e-6
//   )
// ){
//    printf("Successfully wrote matrix.\n");
// }
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int WnMatrix__writeMatrixToAsciiFile(
  WnMatrix *, const char *, double
);

/*##############################################################################
// <routine name="WnMatrix__getNumberOfElements()">
//
//   <description>
//     <abstract>
//       Returns the number of non-zero elements in the matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, elements, number
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t WnMatrix__getNumberOfElements( WnMatrix *self );
//     </calling_sequence>
//
//     <param
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the number of elements in self.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the number of elements in WnMatrix *p_matrix and store it
//         in i_elements:
//       </synopsis>
//
//       <code>
// size_t i_elements;
// i_elements = WnMatrix__getNumberOfElements( p_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t WnMatrix__getNumberOfElements( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__getNumberOfRows()">
//
//   <description>
//     <abstract>
//       Returns the number of rows in the matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, elements, number, rows
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t WnMatrix__getNumberOfRows( WnMatrix *self );
//     </calling_sequence>
//
//     <param
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the number of rows in self.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the number of rows from WnMatrix *p_matrix and stores it in
//         l_rows:
//       </synopsis>
//
//       <code>
// size_t i_rows;
// i_rows = WnMatrix__getNumberOfRows( p_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t WnMatrix__getNumberOfRows( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__getNumberOfColumns()">
//
//   <description>
//     <abstract>
//       Returns the number of columns in the matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, elements, number, column
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t WnMatrix__getNumberOfColumns( WnMatrix *self );
//     </calling_sequence>
//
//     <param
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the number of columns in self.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the number of columns from WnMatrix *p_matrix and store it
//         in l_columns:
//       </synopsis>
//
//       <code>
// size_t i_columns;
// i_columns = WnMatrix__getNumberOfColumns( p_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t WnMatrix__getNumberOfColumns( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__computeMatrixTimesVector()">
//
//   <description>
//     <abstract>
//       Computes the vector y from the equation y = Ax, given matrix A
//       and vector x.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, multiply, vector
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnMatrix__computeMatrixTimesVector(
//   WnMatrix *self,
//   gsl_vector *p_input_vector
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="p_input_vector" 
//       kind="in,positional,required" 
//       doc="x" 
//     >
//       A pointer to a gsl_vector structure containing the vector to
//       be multiplied by the matrix.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new gsl_vector structure containing
//       the product of the matrix and input vector.  It is the caller's
//       responsibility to free the output vector with gsl_vector_free
//       when done with it.  If the number of columns
//       in the input matrix does not equal the length of the input vector, or
//       if any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Multiply the vector p_in by the WnMatrix * p_matrix
//         and return the result as p_out:
//       </synopsis>
//
//       <code>
// p_out =
//   WnMatrix__computeMatrixTimesVector(
//     p_matrix, p_in
//   );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
WnMatrix__computeMatrixTimesVector(
  WnMatrix *, gsl_vector *
);

/*##############################################################################
// <routine name="WnMatrix__computeTransposeMatrixTimesVector()">
//
//   <description>
//     <abstract>
//       Computes the matrix equation Y = A(transpose)x given A and x.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, transpose, multiply, vector
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnMatrix__computeTransposeMatrixTimesVector(
//   WnMatrix *self,
//   gsl_vector *p_input_vector
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="p_input_vector" 
//       kind="in,positional,required" 
//       doc="x" 
//     >
//       A pointer to a gsl_vector structure containing the vector to
//       be multiplied by the matrix.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new gsl_vector structure containing
//       the product of the transpose of the
//       matrix and input vector.  It is the caller's responsibility to
//       free the output vector with gsl_vector_free when done with it.
//       If the number of rows
//       in the input matrix (that is, the number of columns in the transpose
//       matrix) does not equal the length of the input vector, or
//       if any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Multiply the vector p_in by the transpose of WnMatrix * p_matrix
//         and return the result as p_out:
//       </synopsis>
//
//       <code>
// p_out =
//   WnMatrix__computeTransposeMatrixTimesVector(
//     p_matrix, p_in
//   );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
WnMatrix__computeTransposeMatrixTimesVector(
  WnMatrix *, gsl_vector *
);

/*##############################################################################
// <routine name="WnMatrix__insertMatrix()">
//
//   <description>
//     <abstract>
//       Inserts a smaller WnMatrix into a larger one.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, insert
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__insertMatrix(
//   WnMatrix *self, WnMatrix *p_inserted_matrix,
//   size_t i_row, size_t i_col
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param
//       name="p_inserted_matrix"
//       kind="in,positional,required"
//     >
//       A pointer to the WnMatrix that will be inserted into self.
//     </param>
//     <param 
//       name="i_row"
//       kind="in,positional,required" 
//       doc="row"
//     >
//       A size_t integer representing the first row where the matrix will be
//       inserted.
//     </param>
//     <param
//       name="i_col"
//       kind="in,positional,required"
//       doc="col"
//     >
//       A size_t integer representing the first column where the
//       matrix will be inserted.
//     </param>
//     <param 
//       kind="return" 
//       doc="result"
//     />
//     <param
//       name="self"
//       kind="out,positional,required"
//       doc="matrix_out"
//     />
//
//     <doc kind="pre" id="row">
//       i_row > 0 and i_row <= WnMatrix__getNumberOfRows( self )
//     </doc>
//     <doc kind="pre" id="col">
//       i_col > 0 and i_col <= WnMatrix__getNumberOfColumns( self )
//     </doc>
//
//     <doc kind="post" id="matrix_out">
//       self contains p_inserted_matrix beginning at row i_row and column 
//       i_col.
//     </doc>
//     <doc kind="post" id="result">
//       Routine returns with the smaller matrix inserted into self.
//       If a matrix element in the larger matrix already exists prior to
//       insertion of the smaller matrix at a location where the smaller
//       matrix will be added, the final element at that location will
//       be the sum of the prior element and the element of the smaller matrix.
//       If either matrix pointer is invalid, error handling is invoked.
//       If the row or column of the insertion point is invalid, or if
//       the inserted matrix would extend beyond the bounds of the final
//       matrix, error handling is invoked.  The routine does not free
//       the inserted matrix.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Insert WnMatrix *p_inserted_matrix into WnMatrix *self at row 3,
//         column 5:
//       </synopsis>
//
//       <code>
// WnMatrix__insertMatrix( self, p_inserted_matrix, 3, 5 );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void WnMatrix__insertMatrix(
  WnMatrix *, WnMatrix *, size_t, size_t
);

/*##############################################################################
// <routine name="WnMatrix__getDiagonalElements()">
//
//   <description>
//     <abstract>
//        Returns a gsl_vector containing the diagonal elements.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, diagonals, value 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnMatrix__getDiagonalElements(
//   WnMatrix *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix" 
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="diag_out"
//     />
//
//     <doc kind="post" id="diag_out">
//       On successful return, the routine returns a new gsl_vector
//       of length equal to the number of rows in the matrix and containing
//       the diagonal elements (including zeros) in order.  The caller must
//       free the memory for the returned vector.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the diagonal elements of WnMatrix *p_matrix:
//       </synopsis>
//
//       <code>
// gsl_vector *p_diagonals;
// p_diagonals = WnMatrix__getDiagonalElements( p_matrix );
// for(
//   i = 1;
//   i <= WnMatrix__get_gsl_vector_size( p_diagonals );
//   i++
// )
// {
//   printf(
//     "row = %d  diag element = %e\n", 
//     i,
//     gsl_vector_get( p_diagonals, i - 1 );
//   );
// }
// gsl_vector_free( p_diagonals );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
WnMatrix__getDiagonalElements( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__getTransferMatrix()">
//
//   <description>
//     <abstract>
//       Returns the F matrix ( Fij = aij/ajj, but no diagonal elements ).
//     </abstract>
//     <keywords>
//       sparse, webnucleo, F matrix, transfer
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix * WnMatrix__getTransferMatrix( WnMatrix *self );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, routine returns a new WnMatrix * that contains the
//       transfer F matrix.  Routine does not free the input matrix.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the transfer matrix from WnMatrix *p_matrix and store
//         it in WnMatrix *p_transfer_matrix:
//       </synopsis>
//       <code>
// p_transfer_matrix = WnMatrix__getTransferMatrix( p_matrix );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix *WnMatrix__getTransferMatrix( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__getCoo()">
//
//   <description>
//     <abstract>
//        Retrieves the coordinate matrix format of a matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, compressed sparse row, COO
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix__Coo *WnMatrix__getCoo( WnMatrix *self );
//     </calling_sequence>
//
//     <param 
//       name="self"
//       kind="in,positional,required"
//       doc="matrix" 
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param
//       doc="result"
//       kind="return"
//     />
//
//     <doc id="result">
//       Routine returns a pointer to a new WnMatrix__Coo structure that
//       contains the matrix in coordinate format.
//       If the input is invalid or if the memory for the coordinate
//       matrix cannot be allocated, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Store the elements of WnMatrix * p_matrix in Coordinate
//         format in WnMatrix__Coo * p_coo_matrix:
//       </synopsis>
//
//       <code>
// p_coo_matrix = WnMatrix__getCoo( p_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix__Coo *WnMatrix__getCoo( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__Coo__getRowVector()">
//
//   <description>
//     <abstract>
//       Retrieves the row vector for a coordinate matrix structure.
//     </abstract>
//     <keywords>
//        row, webnucleo, coordinate, matrix
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t *
// WnMatrix__Coo__getRowVector(
//   WnMatrix__Coo *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Coo structure.
//     </param>
//     <param 
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the row vector for the coordinate matrix, that is,
//       the array of row numbers of the non-zone elements. The user does
//       not allocate memory for the returned array.  The array
//       memory is freed with WnMatrix__Coo__free().
//       If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the row vector
//         of the matrix stored in coordinate format in p_coo_matrix:
//       </synopsis>
//
//       <code>
// size_t *a_row;
// a_row = WnMatrix__Coo__getRowVector( p_coo_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t *WnMatrix__Coo__getRowVector( WnMatrix__Coo * );

/*##############################################################################
// <routine name="WnMatrix__Coo__getColumnVector()">
//
//   <description>
//     <abstract>
//       Retrieves the column vector for a coordinate
//       matrix structure.
//     </abstract>
//     <keywords>
//        column, webnucleo, coordinate, matrix
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t *
// WnMatrix__Coo__getColumnVector(
//   WnMatrix__Coo *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Coo structure.
//     </param>
//     <param 
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the column vector for the coordinate matrix, that
//       is, an array containing the column numbers of the non-zero elements
//       of the matrix.  The user does not allocate memory for
//       the returned array.  The memory for the array is freed with
//       WnMatrix__Coo_free().
//       If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the column vector
//         of the matrix stored in coordinate format in p_coo_matrix:
//       </synopsis>
//
//       <code>
// a_column = WnMatrix__Coo__getColumnVector( p_coo_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t *WnMatrix__Coo__getColumnVector( WnMatrix__Coo * );

/*##############################################################################
// <routine name="WnMatrix__Coo__getValueVector()">
//
//   <description>
//     <abstract>
//       Retrieves the value vector for a coordinate
//       matrix structure.
//     </abstract>
//     <keywords>
//        value, webnucleo, coordinate, matrix
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double *
// WnMatrix__Coo__getValueVector(
//   WnMatrix__Coo *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Coo structure.
//     </param>
//     <param 
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the value vector for the coordinate matrix.  The
//       user does not allocate memory for the array.
//       If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the value vector
//         of the matrix stored in coordinate format in p_coo_matrix:
//       </synopsis>
//
//       <code>
// a_value = WnMatrix__Coo__getValueVector( p_coo_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double *WnMatrix__Coo__getValueVector( WnMatrix__Coo * );

/*##############################################################################
// <routine name="WnMatrix__Coo__free()">
//
//   <description>
//     <abstract>
//       Frees the memory allocated for a WnMatrix__Coo structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, remove, free, clear 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__Coo__free( WnMatrix__Coo *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Coo structure.
//     </param>
//     <param 
//       name="self"
//       kind="out,positional,required" 
//       doc="matrix_out"
//     />
//
//     <doc kind="post" id="matrix_out">
//       On successful return, the Coo matrix and its elements have all
//       been cleared and freed.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Delete WnMatrix__Coo *p_coo_matrix and free the allocated memory:
//       </synopsis>
//
//       <code>
// WnMatrix__Coo__free( p_coo_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void WnMatrix__Coo__free( WnMatrix__Coo * );

/*##############################################################################
// <routine name="WnMatrix__getCsr()">
//
//   <description>
//     <abstract>
//        Retrieves the  Compressed Sparse Row format of a matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, compressed, row, CSR 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
//  WnMatrix__Csr * WnMatrix__getCsr( WnMatrix *self );
//     </calling_sequence>
//
//     <param 
//       name="self"
//       kind="in,positional,required"
//       doc="matrix" 
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new WnMatrix__Csr structure that
//       contains the matrix in compressed sparse row format.
//       If the input is invalid or if the memory for the compressed
//       sparse row matrix cannot be allocated, error handling is
//       invoked.
//     </doc>
//      
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Store the elements of WnMatrix * p_matrix in Compressed 
//         Sparse Row format in the structure pointed to by p_csr_matrix:
//       </synopsis>
//
//       <code>
// p_csr_matrix = WnMatrix__getCsr( p_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix__Csr *WnMatrix__getCsr( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__Csr__getRowPointerVector()">
//
//   <description>
//     <abstract>
//       Retrieves the row pointer vector for a matrix stored in
//       compressed sparse row matrix format.
//     </abstract>
//     <keywords>
//        row, webnucleo, compressed, sparse, matrix
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t *
// WnMatrix__Csr__getRowPointerVector(
//   WnMatrix__Csr *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Csr structure.
//     </param>
//     <param 
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the row pointer array for a compressed sparse
//       row matrix.  The row pointer array gives the element number of
//       the first non-zero element of the given row.  The returned array
//       has N + 1 elements, where N is the number of rows in the matrix.
//       The user does not allocate memory for the array, and the memory
//       for the array is freed when the user frees the Csr matrix with
//       WnMatrix__Csr__free().
//       If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the row pointer array
//         of the matrix stored in compressed sparse row format in p_csr_matrix:
//       </synopsis>
//
//       <code>
// a_rowptr = WnMatrix__Csr__getRowPointerVector( p_csr_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t *WnMatrix__Csr__getRowPointerVector( WnMatrix__Csr * );

/*##############################################################################
// <routine name="WnMatrix__Csr__getColumnVector()">
//
//   <description>
//     <abstract>
//       Retrieves the column number vector for a matrix in a compressed
//       sparse row matrix structure.
//     </abstract>
//     <keywords>
//        column, webnucleo, compressed, sparse, row, matrix
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t *
// WnMatrix__Csr__getColumnVector(
//   WnMatrix__Csr *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Csr structure.
//     </param>
//     <param 
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the column number vector of the matrix, that is,
//       the array containing the column numbers of the non-zero elements
//       in the matrix.  The user does not allocate memory for the returned
//       array, and the array's memory is freed when the user calls
//       WnMatrix__Csr__free().
//       If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the column number vector
//         of the matrix stored in compressed sparse row format in p_csr_matrix:
//       </synopsis>
//
//       <code>
// a_col = WnMatrix__Csr__getColumnVector( p_csr_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t *WnMatrix__Csr__getColumnVector( WnMatrix__Csr * );

/*##############################################################################
// <routine name="WnMatrix__Csr__getValueVector()">
//
//   <description>
//     <abstract>
//       Retrieves the value vector of a matrix in a compressed sparse
//       row matrix structure.
//     </abstract>
//     <keywords>
//        value, webnucleo, compressed, sparse, row, matrix
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double *
// WnMatrix__Csr__getValue(
//   WnMatrix__Csr *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Csr structure.
//     </param>
//
//     <param 
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the value vector of the element stored in the
//       compressed sparse row format.  This array contains the value
//       of the non-zero matrix elements.  The user does not allocate
//       memory for the returned array, and the memory is freed when the
//       user calls WnMatrix__Csr__free().
//       If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the value vector
//         of the matrix stored in compressed sparse row format in
//         p_csr_matrix:
//       </synopsis>
//
//       <code>
// a_value = WnMatrix__Csr__getValue( p_csr_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double *WnMatrix__Csr__getValueVector( WnMatrix__Csr * );

/*##############################################################################
// <routine name="WnMatrix__Csr__free()">
//
//   <description>
//     <abstract>
//       Frees the memory allocated for a WnMatrix__Csr structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, remove, free, clear 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void WnMatrix__Csr__free( WnMatrix *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Csr structure.
//     </param>
//     <param 
//       name="self"
//       kind="out,positional,required" 
//       doc="matrix_out"
//     />
//
//     <doc kind="post" id="matrix_out">
//       On successful return, the memory for the Csr matrix has
//       been freed.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Delete WnMatrix__Csr *p_csr_matrix and free the allocated memory:
//       </synopsis>
//
//       <code>
// WnMatrix__Csr__free( p_csr_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void WnMatrix__Csr__free( WnMatrix__Csr * );

/*##############################################################################
// <routine name="WnMatrix__Yale__free()">
//
//   <description>
//     <abstract>
//       Frees the memory allocated for a WnMatrix__Yale structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, remove, free, clear, Yale 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void WnMatrix__Yale__free( WnMatrix__Yale *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Yale structure.
//     </param>
//     <param 
//       name="self"
//       kind="out,positional,required" 
//       doc="matrix_out"
//     />
//
//     <doc kind="post" id="matrix_out">
//       On successful return, the memory for the Yale matrix
//       has been freed.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Free WnMatrix__Yale *p_yale_matrix:
//       </synopsis>
//
//       <code>
// WnMatrix__Yale__free( p_yale_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void WnMatrix__Yale__free( WnMatrix__Yale * );

/*##############################################################################
// <routine name="WnMatrix__getYale()">
//
//   <description>
//     <abstract>
//        Retrieves the Yale sparse form of the matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, yale, YSM 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix__Yale *
// WnMatrix__getYaleSparseMatrix(
//   WnMatrix *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix" 
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new WnMatrix__Yale structure that
//       contains the matrix in Yale sparse matrix format.
//       If the input is invalid or if the memory for the Yale
//       sparse matrix cannot be allocated, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Store the WnMatrix *p_matrix in Yale Sparse format
//         in WnMatrix__Yale *p_yale:
//       </synopsis>
//
//       <code>
// p_yale = WnMatrix__getYale( p_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix__Yale *WnMatrix__getYale( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__Yale__getPointerVector()">
//
//   <description>
//     <abstract>
//       Retrieves the element pointer vector for a matrix stored in
//       Yale sparse matrix format.
//     </abstract>
//     <keywords>
//        row, webnucleo, Yale, sparse, matrix
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t *
// WnMatrix__Yale__getPointerVector(
//   WnMatrix__Yale *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Yale structure.
//     </param>
//     <param 
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the element pointer array for a Yale
//       sparse matrix.  The returned array has N + M + 1 elements, where
//       N is the number of rows and M is the number of non-zero matrix
//       elements. For a matrix with N rows, the first N elements
//       of the array point to the index of the same array that contain the
//       first non-zero element of the given row. 
//       The user does not allocate memory for the array, and the memory
//       for the array is freed when the user frees the Yale matrix with
//       WnMatrix__Yale__free().
//       If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the element pointer array
//         of the matrix stored in Yale sparse format in p_yale_matrix:
//       </synopsis>
//
//       <code>
// a_ptr = WnMatrix__Yale__getPointerVector( p_yale_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t *WnMatrix__Yale__getPointerVector( WnMatrix__Yale * );

/*##############################################################################
// <routine name="WnMatrix__Yale__getValueVector()">
//
//   <description>
//     <abstract>
//       Retrieves the value vector for a matrix stored in
//       Yale sparse matrix format.
//     </abstract>
//     <keywords>
//        value, webnucleo, Yale, sparse, matrix
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double *
// WnMatrix__Yale__getValueVector(
//   WnMatrix__Yale *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Yale structure.
//     </param>
//     <param 
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the value array for a Yale
//       sparse matrix.  The returned array has N + M + 1 elements, where
//       N is the number of rows and M is the number of non-zero matrix
//       elements. For a matrix with N rows, the first N elements
//       of the array contain the diagonal matrix elements.
//       The elements for index >= N + 2, where N = number of rows, contain
//       the off-diagonal elements ordered by rows and, within each row,
//       ordered by column.
//       The user does not allocate memory for the array, and the memory
//       for the array is freed when the user frees the Yale matrix with
//       WnMatrix__Yale__free().
//       If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the value array
//         of the matrix stored in Yale sparse format in p_yale_matrix:
//       </synopsis>
//
//       <code>
// a_val = WnMatrix__Yale__getValueVector( p_yale_matrix );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double *WnMatrix__Yale__getValueVector( WnMatrix__Yale * );

/*##############################################################################
// <routine name="WnMatrix__addValueToDiagonals()">
//
//   <description>
//     <abstract>
//       Adds a value to the diagonal elements of the matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, add, diagonals, element 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__addValueToDiagonals(
//   WnMatrix *self, double d_val
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param
//       name="d_val"
//       kind="in,positional,required"
//     >
//       A double containing the value to add to the diagonals.
//     </param>
//     <param 
//       name="self"
//       kind="out,positional,required" 
//       doc="matrix_out"
//     />
//
//     <doc kind="post" id="matrix_out">
//       All diagonal elements are increased by d_val. 
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Add 0.5 to the diagonals of WnMatrix *p_matrix:
//       </synopsis>
//
//       <code>
// WnMatrix__addValueToDiagonals( p_matrix, 0.5 );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void WnMatrix__addValueToDiagonals( WnMatrix *, double );

/*##############################################################################
// <routine name="WnMatrix__getRow()">
//
//   <description>
//     <abstract>
//       Gets a row from a WnMatrix structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, row, get
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix__Line *
// WnMatrix__getRow(
//   WnMatrix *self, size_t i_row
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="i_row" 
//       kind="in,positional,required" 
//       doc="i_row"
//     >
//       A size_t integer giving the row to get.
//     </param>
//     <param 
//       name="return"
//       kind="out,positional,required" 
//       doc="p_row"
//     />
//
//     <doc kind="post" id="p_row">
//       Routine returns a pointer to a new WnMatrix__Line structure containing
//       the row.  If the input matrix is invalid, error handling is invoked.
//       It is the caller's responsibility to free the returned WnMatrix__Line
//       structure with WnMatrix__Line__free().
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Get row number 5 from WnMatrix *p_my_matrix:
//       </synopsis>
//
//       <code>
// WnMatrix__Line *p_row;
// p_row = WnMatrix__getRow( p_my_matrix, 5 );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix__Line *WnMatrix__getRow( WnMatrix *, size_t );

/*##############################################################################
// <routine name="WnMatrix__getColumn()">
//
//   <description>
//     <abstract>
//       Gets a column from a WnMatrix structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, column, get
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix__Line *
// WnMatrix__getRow(
//   WnMatrix *self, size_t i_col
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="i_col" 
//       kind="in,positional,required" 
//       doc="i_col"
//     >
//       A size_t integer giving the column to get.
//     </param>
//     <param 
//       name="return"
//       kind="out,positional,required" 
//       doc="p_col"
//     />
//
//     <doc kind="post" id="p_col">
//       Routine returns a pointer to the WnMatrix__Line structure containing
//       the column.  If the input matrix is invalid, error handling is
//       invoked.
//       It is the caller's responsibility to free the returned WnMatrix__Line
//       structure with WnMatrix__Line__free().
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Get column number 3 from WnMatrix *p_my_matrix:
//       </synopsis>
//
//       <code>
// WnMatrix__Line *p_col;
// p_col = WnMatrix__getColumn( p_my_matrix, 3 );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix__Line *WnMatrix__getColumn( WnMatrix *, size_t );

/*##############################################################################
// <routine name="WnMatrix__Line__getNumberOfElements()">
//
//   <description>
//     <abstract>
//       Gets the number of non-zero elements in a row or column
//       in a WnMatrix structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, line, get, number
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t
// WnMatrix__Line__getNumberOfElements(
//   WnMatrix__Line *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_line"
//     >
//       A pointer to a WnMatrix__Line structure.
//     </param>
//     <param 
//       name="return"
//       kind="out,positional,required" 
//       doc="i_num"
//     />
//
//     <doc kind="post" id="i_num">
//       Routine returns the number of non-zero elements in the input row or
//       column.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the number of non-zero elements in row 3 of
//         WnMatrix *p_my_matrix:
//       </synopsis>
//
//       <code>
// p_row = WnMatrix__getRow( p_my_matrix, 3 );
// printf(
//   "%d\n",
//   WnMatrix__getNumberOfElements( p_row )
// );
// WnMatrix__Line__free( p_row );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t WnMatrix__Line__getNumberOfElements( WnMatrix__Line * );

/*##############################################################################
// <routine name="WnMatrix__Line__getNonZeroIndices()">
//
//   <description>
//     <abstract>
//       Gets an array containing the column numbers of the
//       non-zero elements in a row or the row numbers of the non-zero
//       elements in a column in a WnMatrix__Line structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, line, get, indices
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t *
// WnMatrix__Line__getNonZeroIndices(
//   WnMatrix__Line *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_line"
//     >
//       A pointer to a WnMatrix__Line structure.
//     </param>
//     <param 
//       name="return"
//       kind="out,positional,required" 
//       doc="i_index"
//     />
//
//     <doc kind="post" id="i_index">
//       Routine returns a new array of size_ts
//       containing the column numbers of the
//       non-zero elements in a row or the row number of the non-zero elements
//       in a column.  The caller must free the memory for this array.  This
//       is best done by calling WnMatrix__Line__free.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the column numbers of the non-zero elements of row 5 in
//         WnMatrix *p_my_matrix:
//       </synopsis>
//
//       <code>
// size_t i, *a_indices;
// WnMatrix__Line *p_row;
// p_row = WnMatrix__getRow( p_my_matrix, 5 );
// a_indices = WnMatrix__Line__getNonZeroIndices( p_row );
// for(
//      i = 0;
//      i < WnMatrix__Line__getNumberOfElements( p_row );
//      i++
// ) {
//     printf( "%d\n", a_indices[i] );
// }
// WnMatrix__Line__free( p_row );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t *WnMatrix__Line__getNonZeroIndices( WnMatrix__Line * );

/*##############################################################################
// <routine name="WnMatrix__Line__getNonZeroElements()">
//
//   <description>
//     <abstract>
//       Gets an array containing the values of the
//       non-zero elements in a row or of the non-zero
//       elements in a column in a WnMatrix__Line structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, line, get, values
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double *
// WnMatrix__Line__getNonZeroElements(
//   WnMatrix__Line *self
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_line"
//     >
//       A pointer to a WnMatrix__Line structure.
//     </param>
//     <param 
//       name="return"
//       kind="out,positional,required" 
//       doc="i_element"
//     />
//
//     <doc kind="post" id="i_element">
//       Routine returns a new array containing the values of the
//       non-zero elements in a row or the values of the non-zero elements
//       in a column.  The caller must free the memory for this array.  This
//       is best done by calling WnMatrix__Line__free.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the column numbers and values of the non-zero elements of
//         row 5 in WnMatrix *p_my_matrix:
//       </synopsis>
//
//       <code>
// size_t i, *a_indices;
// double *a_values;
// WnMatrix__Line *p_row;
// p_row = WnMatrix__getRow( p_my_matrix, 5 );
// a_indices = WnMatrix__Line__getNonZeroIndices( p_row );
// a_values = WnMatrix__Line__getNonZeroElements( p_row );
// for(
//      i = 0;
//      i < WnMatrix__Line__getNumberOfElements( p_row );
//      i++
// ) {
//     printf( "%d  %e\n", a_indices[i], a_elements[i] );
// }
// WnMatrix__Line__free( p_row );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double *WnMatrix__Line__getNonZeroElements( WnMatrix__Line * );

/*##############################################################################
// <routine name="WnMatrix__extractMatrix()">
//
//   <description>
//     <abstract>
//       Extracts a smaller WnMatrix from a larger one.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, extract
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix *
// WnMatrix__extractMatrix(
//   WnMatrix *self,
//   size_t i_row,
//   size_t i_col,
//   size_t i_row_offset,
//   size_t i_col_offset
// );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="i_row"
//       kind="in,positional,required" 
//       doc="row"
//     >
//       A size_t integer representing the first row where the matrix will be
//       extracted.
//     </param>
//     <param
//       name="i_col"
//       kind="in,positional,required"
//       doc="col"
//     >
//       A size_t inteeger representing the first column where the
//       matrix will be extracted.
//     </param>
//     <param
//       name="i_row_offset"
//       kind="in,positional,required"
//       doc="row_offset"
//     >
//       A size_t integer representing the number of rows to be extracted.
//     </param>
//     <param
//       name="i_col_offset"
//       kind="in,positional,required"
//       doc="col_offset"
//     >
//       A size_t inteeger representing the number of columns to be extracted.
//     </param>
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="pre" id="row">
//       i_row > 0 and i_row <= WnMatrix__getNumberOfRows( self )
//     </doc>
//     <doc kind="pre" id="col">
//       i_col > 0 and i_col <= WnMatrix__getNumberOfColumns( self )
//     </doc>
//     <doc kind="pre" id="row_offset">
//       i_row_offset > 0 and 
//       (i_row + i_row_offset) <= ( WnMatrix__getNumberOfRows( self ) + 1 )
//     </doc>
//     <doc kind="pre" id="col_offset">
//       i_col_offset > 0 and 
//       (i_col + i_col_offset) <= ( WnMatrix__getNumberOfColumns( self ) + 1 )
//     </doc>
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to the new extracted matrix if the routine
//       was successful.  Routine does not free the input
//       matrix, and the caller must free the memory for both the input
//       matrix and the extracted matrix after use (with WnMatrix__free).
//       If the input matrix is invalid, error handling is invoked.
//       If the matrix to be extracted does not fully lie within the parent
//       matrix or if sufficient memory could not be allocated,
//       error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Extract a 2x2 matrix from self beginning at row 3, column 5:
//       </synopsis>
//
//       <code>
// WnMatrix *p_extracted_matrix;
// p_extracted_matrix =
//   WnMatrix__extractMatrix( self, 3, 5, 2, 2 );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix *WnMatrix__extractMatrix(
  WnMatrix *, size_t, size_t, size_t, size_t
);

/*##############################################################################
// <routine name="WnMatrix__getGslMatrix()">
//
//   <description>
//     <abstract>
//       Stores the matrix in gsl_matrix format.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, dense, doubly indexed array
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_matrix *
// WnMatrix__getGslMatrix( WnMatrix *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="mat_out"
//     />
//
//     <doc kind="post" id="mat_out">
//       For valid input, routine returns a new gsl_matrix that contains the
//       elements of WnMatrix *self.  Routine does
//       not free the input matrix.  The user must free the memory for the
//       gsl_matrix matrix.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example1">
//       <synopsis>
//         Store p_matrix in dense format in gsl_matrix p_mat:
//       </synopsis>
//
//       <code>
// gsl_matrix *p_mat;
// p_mat = WnMatrix__getGslMatrix( p_matrix );
// ... (do something with the matrix)
// gsl_matrix_free( p_mat );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_matrix *
WnMatrix__getGslMatrix( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__scaleMatrix()">
//
//   <description>
//     <abstract>
//       Multiplies all the elements in the matrix by a scalar
//       constant.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, scale, multiply, scalar 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void WnMatrix__scaleMatrix( WnMatrix *self, double d_val );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param
//       name="d_val"
//       kind="in,positional,required"
//     >
//       A double containing the factor by which to scale the matrix.
//     </param>
//     <param 
//       name="self"
//       kind="return"
//       doc="matrix_out"
//     />
//
//     <doc kind="post" id="matrix_out">
//       Routine returns matrix with all elements scaled by factor d_val. 
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Multiply all elements of  WnMatrix *p_matrix by a factor of 3.5:
//       </synopsis>
//
//       <code>
// WnMatrix__scaleMatrix( p_matrix, 3.5 );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void WnMatrix__scaleMatrix( WnMatrix *, double );

/*##############################################################################
// <routine name="WnMatrix__Line__free()">
//
//   <description>
//     <abstract>
//       Frees the memory allocated for a WnMatrix__Line structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, free, column, delete, clear, line, row 
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void WnMatrix__Line__free( WnMatrix__Line *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Line structure.
//     </param>
//     <param 
//       name="self"
//       kind="out,positional,required" 
//       doc="line_out"
//     />
//
//     <doc kind="post" id="line_out">
//       No longer has memory allocated for it. 
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Delete WnMatrix__Line *p_row and free the allocated memory:
//       </synopsis>
//
//       <code>
// WnMatrix__Line__free( p_row );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void WnMatrix__Line__free( WnMatrix__Line * );

/*##############################################################################
// <routine name="WnMatrix__getCopy()">
//
//   <description>
//     <abstract>
//       Returns a copy of the matrix.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, copy
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix * WnMatrix__getCopy( WnMatrix *self );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, routine returns a new WnMatrix * that contains
//       a copy of the matrix.  Routine does not free the input matrix.
//       The caller must free the memory for this new matrix with
//       WnMatrix__free.  If the input matrix is invalid, the behavior
//       is undefined, although this error may be caught by the error
//       handler.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve a copy of the WnMatrix *p_matrix and store
//         it in WnMatrix *p_copy:
//       </synopsis>
//       <code>
// p_copy = WnMatrix__getCopy( p_matrix );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix* WnMatrix__getCopy( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__getTranspose()">
//
//   <description>
//     <abstract>
//       Returns the transpose of the matrix.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, copy
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix * WnMatrix__getTranspose( WnMatrix *self );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, routine returns a new WnMatrix * that contains
//       the transpose of the matrix.  Routine does not free the input matrix.
//       The caller must free the memory for this new matrix with
//       WnMatrix__free.  If the input matrix is invalid, the behavior
//       is undefined, although this error may be caught by the error handler.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Retrieve the tranpose of the WnMatrix *p_matrix and store
//         it in WnMatrix *p_transpose_matrix:
//       </synopsis>
//       <code>
// p_tranpose_matrix = WnMatrix__getTranspose( p_matrix );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix* WnMatrix__getTranspose( WnMatrix * );

/*##############################################################################
// <routine name="WnMatrix__Coo__writeToXmlFile()">
//
//   <description>
//     <abstract>
//       Outputs the coordinate matrix to an XML file.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, XML, file
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__Coo_writeToXmlFile(
//   WnMatrix__Coo *self,
//   const char *s_output_xml_file,
//   const char *s_format
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix__Coo structure.
//     </param>
//
//     <param
//       name="s_output_xml_file"
//       kind="in,positional,required"
//       doc="coo_output_file"
//     >
//       A string giving the name of the file to which the XML output is
//       to be written.
//     </param>
//
//     <param
//       name="s_format"
//       kind="in,positional,required"
//       doc="format"
//     >
//       A string giving the C-style format code for output of the matrix
//       element value or NULL for the default (%g).
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, outputs the data for the coordinate matrix in
//       XML format in the designated file for the designated format
//       code.  If the input is not valid,
//       or if the XML output routines fail, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Write the coordinate matrix data in p_my_coo_matrix to the
//         file "coo.xml" using the default format:
//       </synopsis>
//       <code>
// WnMatrix__Coo__writeToXmlFile(
//   p_my_coo_matrix,
//   "coo.xml",
//   NULL
// );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnMatrix__Coo__writeToXmlFile( WnMatrix__Coo *, const char *, const char * );

/*##############################################################################
// <routine name="WnMatrix__Csr__writeToXmlFile()">
//
//   <description>
//     <abstract>
//       Outputs the compressed sparse row matrix to an XML file.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, XML, file
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__Csr_writeToXmlFile(
//   WnMatrix__Csr *self,
//   const char *s_output_xml_file,
//   const char *s_format
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix__Csr structure.
//     </param>
//
//     <param
//       name="s_output_xml_file"
//       kind="in,positional,required"
//       doc="csr_output_file"
//     >
//       A string giving the name of the file to which the XML output is
//       to be written.
//     </param>
//
//     <param
//       name="s_format"
//       kind="in,positional,required"
//       doc="format"
//     >
//       A string giving the C-style format code for output of the matrix
//       element value or NULL for the default (%g).
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, outputs the data for the csr matrix in
//       XML format in the designated file with the designated format
//       code for the matrix element value.  If the input is not valid,
//       or if the XML output routines fail, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Write the compressed sparse row matrix data in p_my_csr_matrix to the
//         file "csr.xml" using exponential format with fifteen decimal
//         place precision:
//       </synopsis>
//       <code>
// WnMatrix__Csr__writeToXmlFile(
//   p_my_csr_matrix,
//   "csr.xml",
//   "%.15e"
// );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnMatrix__Csr__writeToXmlFile( WnMatrix__Csr *, const char *, const char * );

/*##############################################################################
// <routine name="WnMatrix__Yale__writeToXmlFile()">
//
//   <description>
//     <abstract>
//       Outputs the Yale sparse matrix to an XML file.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, XML, file
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__Yale_writeToXmlFile(
//   WnMatrix__Yale *self,
//   const char *s_output_xml_file,
//   const char *s_format
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix__Yale structure.
//     </param>
//     <param
//       name="s_output_xml_file"
//       kind="in,positional,required"
//       doc="yale_output_file"
//     >
//       A string giving the name of the file to which the XML output is
//       to be written.
//     </param>
//     <param
//       name="s_format"
//       kind="in,positional,required"
//       doc="format"
//     >
//       A string giving the C-style format code for output of the matrix
//       element value or NULL for the default (%g).
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, outputs the data for the yale matrix in
//       XML format in the designated file with the designated format
//       code for the matrix element value.  If the input is not valid,
//       or if the XML output routines fail, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Write the yale matrix data in p_my_yale_matrix to the
//         file "yale.xml" using float format with ten decimal place
//         precision for output of the matrix element values:
//       </synopsis>
//       <code>
// WnMatrix__Yale__writeToXmlFile(
//   p_my_yale_matrix,
//   "yale.xml",
//   "%.10f"
// );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnMatrix__Yale__writeToXmlFile( WnMatrix__Yale *, const char *, const char * );

/*##############################################################################
// <routine name="WnMatrix__new_from_xml()">
//
//   <description>
//     <abstract>
//       Creates a new matrix from data in an xml file.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, XML, file
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix *
// WnMatrix__new_from_xml
//   const char *s_input_xml_file, const char *s_xpath_expression
// );
//     </calling_sequence>
//
//     <param
//       name="s_input_xml_file"
//       kind="in,positional,required"
//       doc="input_xml"
//     >
//       A string giving the name of the input xml file.
//     </param>
//     <param
//       name="s_xpath_expression"
//       kind="in,positional,required"
//       doc="xpath_expression"
//     >
//       A string giving an XPath expression to use in selecting out the
//       matrix data.
//     </param>
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, routine returns a pointer to a WnMatrix
//       structure with the input data.  If the input is not valid,
//       or if the XML parser routines fail, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Create a WnMatrix structure from data in the file input.xml:
//       </synopsis>
//       <code>
// p_my_matrix = WnMatrix__new_from_xml( "input.xml", NULL );
//       </code>
//     </doc>
//       
//     <doc kind="example" id="example2">
//       <synopsis>
//         Create a WnMatrix structure from data
//         in the file input.xml such that only non-negative values are
//         included:
//       </synopsis>
//       <code>
// p_my_matrix = WnMatrix__new_from_xml( "input.xml", "[value >= 0]" );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix *WnMatrix__new_from_xml( const char *, const char * );

/*##############################################################################
// <routine name="WnMatrix__is_valid_input_xml()">
//
//   <description>
//     <abstract>
//       Validates an input wn_matrix XML file. 
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, XML, file
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int
// WnMatrix__is_valid_input_xml(
//   const char *s_input_xml_file
// );
//     </calling_sequence>
//
//     <param
//       name="s_input_xml_file"
//       kind="in,positional,required"
//       doc="input_xml"
//     >
//       A string giving the name of the input xml file.
//     </param>
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, routine returns 1 (true) if the XML file
//       is valid or 0 (false) if it is not.  If the file is not valid,
//       the routine also prints out the error message.  If there is
//       a problem with the schema, or the routine cannot find the schema
//       over the web,
//       or if the XML parser routines fail, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Validate the input XML file input.xml:
//       </synopsis>
//       <code>
// if( WnMatrix__is_valid_input_xml( "input.xml" ) ) {
//   printf( "Valid input XML!\n" );
// }
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int
WnMatrix__is_valid_input_xml( const char * );

/*##############################################################################
// <routine name="WnMatrix__is_valid_vector_input_xml()">
//
//   <description>
//     <abstract>
//       Validates an input wn_matrix vector XML file. 
//     </abstract>
//     <keywords>
//       sparse, webnucleo, XML, file, vector
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int
// WnMatrix__is_valid_vector_input_xml(
//   const char *s_input_xml_file
// );
//     </calling_sequence>
//
//     <param
//       name="s_input_xml_file"
//       kind="in,positional,required"
//       doc="input_xml"
//     >
//       A string giving the name of the input xml file.
//     </param>
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, routine returns 1 (true) if the XML file
//       is valid or 0 (false) if it is not.  If the file is not valid,
//       the routine also prints out the error message.  If there is
//       a problem with the schema, or the routine cannot find the schema
//       over the web,
//       or if the XML parser routines fail, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Validate the input XML file input.xml:
//       </synopsis>
//       <code>
// if( WnMatrix__is_valid_vector_input_xml( "input.xml" ) ) {
//   printf( "Valid input XML vector!\n" );
// }
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

int WnMatrix__is_valid_vector_input_xml( const char * );

/*##############################################################################
// <routine name="WnMatrix__new_gsl_vector_from_xml()">
//
//   <description>
//     <abstract>
//       Creates a new vector from data in an xml file.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, vector, XML, file
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnMatrix__new_gsl_vector_from_xml
//   const char *s_input_xml_file, const char *s_xpath_expression
// );
//     </calling_sequence>
//
//     <param
//       name="s_input_xml_file"
//       kind="in,positional,required"
//       doc="input_xml"
//     >
//       A string giving the name of the input xml file.
//     </param>
//     <param
//       name="s_xpath_expression"
//       kind="in,positional,required"
//       doc="xpath_expression"
//     >
//       A string giving an XPath expression to use in selecting out the
//       vector data.
//     </param>
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, routine returns a pointer to a new gsl_vector
//       structure with the input data.  If the input is not valid,
//       or if the XML parser routines fail, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Create a new gsl_vector structure from
//         data in the file input.xml:
//       </synopsis>
//       <code>
// p_my_vector =
//   WnMatrix__new_gsl_vector_from_xml( "input.xml", NULL );
//       </code>
//     </doc>
//       
//     <doc kind="example" id="example2">
//       <synopsis>
//         Create a new gsl_vector structure
//         from data in the file input.xml such that only non-negative values
//         are included:
//       </synopsis>
//       <code>
// p_my_vector =
//   WnMatrix__new_gsl_vector_from_xml( "input.xml", "[. >= 0]" );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
WnMatrix__new_gsl_vector_from_xml(
  const char *,
  const char *
);

/*##############################################################################
// <routine name="WnMatrix__write_gsl_vector_to_xml_file()">
//
//   <description>
//     <abstract>
//       Outputs the vector to an XML file.
//     </abstract>
//     <keywords>
//       webnucleo, vector, XML, file
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__write_gsl_vector_to_xml_file(
//   gsl_vector *self,
//   const char *s_output_xml_file,
//   const char *s_format
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="vector"
//     >
//       A pointer to a gsl_vector.
//     </param>
//     <param
//       name="s_output_xml_file"
//       kind="in,positional,required"
//       doc="vector_output_file"
//     >
//       A string giving the name of the file to which the XML output is
//       to be written.
//     </param>
//
//     <param
//       name="s_format"
//       kind="in,positional,required"
//       doc="format"
//     >
//       A string giving the format code for output of the vector
//       component value or NULL for the default (%g).
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, outputs the data for the vector in
//       XML format in the designated file with the designated format
//       for the vector component value.  If the input is not valid,
//       or if the XML output routines fail, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Write the vector data in p_my_vector to the file "vector.xml"
//         using the default format:
//       </synopsis>
//       <code>
// WnMatrix__write_gsl_vector_to_xml_file(
//   p_my_vector, "vector.xml", NULL
// );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnMatrix__write_gsl_vector_to_xml_file(
  gsl_vector *, const char *, const char *
);

/*##############################################################################
// <routine name="WnMatrix__get_gsl_vector_size()">
//
//   <description>
//     <abstract>
//       Gets the size of a gsl vector.
//     </abstract>
//     <keywords>
//       webnucleo, vector, XML, size
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t
// WnMatrix__get_gsl_vector_size(
//   gsl_vector *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="vector"
//     >
//       A pointer to a gsl_vector.
//     </param>
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, returns the size of the
//       vector.  If the input is not valid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Print the size of the gsl_vector p_my_vector:
//       </synopsis>
//       <code>
// printf(
//   "The size of the vector is %d\n",
//   WnMatrix__get_gsl_vector_size(
//     p_my_vector
//   )
// );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t
WnMatrix__get_gsl_vector_size(
  gsl_vector *
);

/*##############################################################################
// <routine name="WnMatrix__get_gsl_vector_array()">
//
//   <description>
//     <abstract>
//       Gets the data array of a gsl vector.
//     </abstract>
//     <keywords>
//       webnucleo, vector, XML, array
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// double *
// WnMatrix__get_gsl_vector_array(
//   gsl_vector *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="vector"
//     >
//       A pointer to a gsl_vector.
//     </param>
//     <param 
//       kind="return" 
//       doc="result" 
//     />
//
//     <doc kind="post" id="result">
//       For valid input, returns a pointer to a double array containing
//       the data for the vector.  The caller does not need to free this
//       array; it is freed when the user frees the original gsl_vector.
//       If the input is not valid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Get the double array containing the data for the
//         gsl_vector *p_my_vector:
//       </synopsis>
//       <code>
// a_data = WnMatrix__get_gsl_vector_array( p_my_vector );
//       </code>
//       
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

double *
WnMatrix__get_gsl_vector_array(
  gsl_vector *
);

/*##############################################################################
// <routine name="WnMatrix__solve()">
//
//   <description>
//     <abstract>
//       Uses the gsl LU decomposition routines to solve the matrix equation
//       Ax = b for x given A and b.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, equation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnMatrix__solve
//   WnMatrix *self,
//   gsl_vector *p_input_vector
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//     <param 
//       name="p_input_vector" 
//       kind="in,positional,required" 
//       doc="x" 
//     >
//       A pointer to a gsl_vector structure containing the right-hand
//       side vector.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new gsl_vector structure containing
//       the solution vector.  It is the caller's
//       responsibility to free the output vector with gsl_vector_free
//       when done with it.  If the number of columns
//       in the input matrix does not equal the length of the input vector, or
//       if any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Solve the matrix equation for input matrix WnMatrix * p_matrix
//         and right-hand-side vector p_in and return the result as p_out:
//       </synopsis>
//
//       <code>
// p_out =
//   WnMatrix__solve(
//     p_matrix, p_in
//   );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
WnMatrix__solve( WnMatrix *, gsl_vector * );

/*##############################################################################
// <routine name="WnMatrix__Arrow__solve()">
//
//   <description>
//     <abstract>
//       Uses internal wn_matrix routines to solve the matrix equation
//       Ax = b for x given A and b by Gaussian elimination.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, equation, gaussian,
//       elimination, arrow
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnMatrix__Arrow__solve(
//   WnMatrix__Arrow *self,
//   gsl_vector *p_input_vector
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix__Arrow structure.
//     </param>
//     <param 
//       name="p_input_vector" 
//       kind="in,positional,required" 
//       doc="x" 
//     >
//       A pointer to a gsl_vector structure containing the right-hand
//       side vector.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new gsl_vector structure containing
//       solution vector.  The rhs vector and the
//       input matrix return as modified by the row operations.
//       It is the caller's responsibility to free the output vector with
//       gsl_vector_free when done with it.  If any input is invalid,
//       or if the number of rows of the matrix does not equal the rhs vector
//       length, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Solve the matrix equation for the input arrow matrix
//         WnMatrix__Arrow * p_arrow and right-hand-side vector p_in,
//         return the result as p_out, convert the modified arrow matrix
//         to a WnMatrix and output as an ascii file (for diagnostic purposes):
//       </synopsis>
//
//       <code>
// p_out =
//   WnMatrix__Arrow__solve(
//     p_arrow, p_in
//   );
// p_modified_arrow =
//   WnMatrix__Arrow__getWnMatrix( p_arrow );
// WnMatrix__writeMatrixToAsciiFile(
//   p_modified_arrow,
//   "modified_arrow.txt",
//   0.
// );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
WnMatrix__Arrow__solve( WnMatrix__Arrow *, gsl_vector * );

/*##############################################################################
// <routine name="WnMatrix__getArrow()">
//
//   <description>
//     <abstract>
//       Returns an arrow matrix, that is, one in which the matrix is
//       stored in arrow format (a central band with wings along the far
//       right side and bottom).
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, equation, gaussian,
//       elimination, arrow
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix__Arrow *
// WnMatrix__getArrow(
//   WnMatrix *self,
//   size_t i_wing_width
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix.
//     </param>
//     <param 
//       name="i_wing_width"
//       kind="in,positional,required" 
//       doc="wing" 
//     >
//       A size_t giving the width of the arrow wings.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="pre" id="wing">
//       The width of the arrow wings must be less than the number of
//       rows in the matrix.
//     </doc>
//
//     <doc kind="post" id="result">
//       Routine returns a pointer to a new WnMatrix__Arrow structure.
//       It is the caller's responsibility to free the arrow matrix with
//       WnMatrix__Arrow__free() when done with it.  If the number of columns
//       in the input matrix does not equal the length of the input vector, or
//       if any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Get the arrow matrix form of the WnMatrix *p_my_matrix with
//         arrow wing width of 3:
//       </synopsis>
//
//       <code>
// p_arrow = WnMatrix__getArrow( p_my_matrix, 3L );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix__Arrow *
WnMatrix__getArrow( WnMatrix *, size_t i_wing );

/*##############################################################################
// <routine name="WnMatrix__Arrow__free()">
//
//   <description>
//     <abstract>
//       Frees the memory allocated for a WnMatrix__Arrow structure.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, remove, free, clear, arrow
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__Arrow__free( WnMatrix__Arrow *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Arrow structure.
//     </param>
//
//     <doc kind="post" id="matrix_out">
//       On successful return, the matrix and its elements have all
//       been cleared and freed.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Delete WnMatrix__Arrow *p_arrow and free the allocated memory:
//       </synopsis>
//
//       <code>
// WnMatrix__Arrow__free( p_arrow );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnMatrix__Arrow__free( WnMatrix__Arrow * );

/*##############################################################################
// <routine name="WnMatrix__Arrow__getBandWidth()">
//
//   <description>
//     <abstract>
//       Returns the central band width for a WnMatrix__Arrow matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, band, width, arrow
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t
// WnMatrix__Arrow__getBandWidth( WnMatrix__Arrow *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Arrow structure.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the width of the central band is
//       returned.  If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Get the band width of WnMatrix__Arrow *p_arrow:
//       </synopsis>
//
//       <code>
// i_band_width =
//   WnMatrix__Arrow__getBandWidth( p_arrow );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t
WnMatrix__Arrow__getBandWidth( WnMatrix__Arrow * );

/*##############################################################################
// <routine name="WnMatrix__Arrow__getWingWidth()">
//
//   <description>
//     <abstract>
//       Returns the width of the wings for a WnMatrix__Arrow matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, wing, width, arrow
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t
// WnMatrix__Arrow__getWingWidth( WnMatrix__Arrow *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Arrow structure.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the width of the arrow wings is
//       returned.  If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Get the arrow width of WnMatrix__Arrow *p_arrow:
//       </synopsis>
//
//       <code>
// i_wing_width =
//   WnMatrix__Arrow__getWingWidth( p_arrow );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t
WnMatrix__Arrow__getWingWidth( WnMatrix__Arrow * );

/*##############################################################################
// <routine name="WnMatrix__Arrow__getNumberOfRows()">
//
//   <description>
//     <abstract>
//       Returns the number of rows in a WnMatrix__Arrow matrix.
//     </abstract>
//     <keywords>
//        sparse, webnucleo, matrix, nubmer, rows, arrow
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// size_t
// WnMatrix__Arrow__getNumberOfRows( WnMatrix__Arrow *self );
//     </calling_sequence>
//
//     <param 
//       name="self" 
//       kind="in,positional,required" 
//       doc="matrix_in"
//     >
//       A pointer to a WnMatrix__Arrow structure.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the number of rows is
//       returned.  If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Get the number of rows in WnMatrix__Arrow *p_arrow:
//       </synopsis>
//
//       <code>
// i_number_rows =
//   WnMatrix__Arrow__getNumberOfRows( p_arrow );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

size_t
WnMatrix__Arrow__getNumberOfRows( WnMatrix__Arrow * );

/*##############################################################################
// <routine name="WnMatrix__Arrow__getWnMatrix()">
//
//   <description>
//     <abstract>
//       Gets a WnMatrix from an arrow matrix.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, arrow
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnMatrix *
// WnMatrix__Arrow__getWnMatrix(
//   WnMatrix__Arrow *self
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix__Arrow structure.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Routine returns the input arrow matrix as a WnMatrix.
//       If any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Get the WnMatrix form of the WnMatrix__Arrow matrix *p_arrow:
//       </synopsis>
//
//       <code>
// p_matrix = WnMatrix__Arrow__getWnMatrix( p_arrow );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnMatrix *
WnMatrix__Arrow__getWnMatrix( WnMatrix__Arrow * );

/*##############################################################################
// <routine name="WnMatrix__removeRow()">
//
//   <description>
//     <abstract>
//       Removes a row from a matrix.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, remove, row
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__removeRow(
//   WnMatrix *self,
//   size_t i_row
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//
//     <param
//       name="i_row"
//       kind="in,positional,required"
//       doc="row"
//     >
//       A size_t giving the row to remove.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="pre" id="row">
//       The row to be removed must be within the matrix; that is, its
//       index must be less than or equal to the number of rows in the matrix.
//     </doc>
//
//     <doc kind="post" id="result">
//       On successful return, the row has been removed.  The number of
//       rows in the matrix is one less than upon input.
//       If any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Remove row 3 from matrix p_matrix:
//       </synopsis>
//
//       <code>
//  WnMatrix__removeRow( p_matrix, 3L );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnMatrix__removeRow( WnMatrix *, size_t );

/*##############################################################################
// <routine name="WnMatrix__removeColumn()">
//
//   <description>
//     <abstract>
//       Removes a column from a matrix.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, remove, column
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__removeColumn(
//   WnMatrix *self,
//   size_t i_column
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//
//     <param
//       name="i_column"
//       kind="in,positional,required"
//       doc="column"
//     >
//       A size_t giving the column to remove.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="pre" id="column">
//       The column to be removed must be within the matrix; that is, its
//       index must be less than or equal to the number of columns in the
//       matrix.
//     </doc>
//
//     <doc kind="post" id="result">
//       On successful return, the column has been removed.  The number of
//       columns in the matrix is one less than upon input.
//       If any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Remove column 3 from matrix p_matrix:
//       </synopsis>
//
//       <code>
//  WnMatrix__removeColumn( p_matrix, 3L );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnMatrix__removeColumn( WnMatrix *, size_t );

/*##############################################################################
// <routine name="WnMatrix__insertRow()">
//
//   <description>
//     <abstract>
//       Inserts a row into a matrix.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, insert, row
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__insertRow(
//   WnMatrix *self,
//   size_t i_row,
//   gsl_vector *p_row
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//
//     <param
//       name="i_row"
//       kind="in,positional,required"
//       doc="row"
//     >
//       A size_t giving the row index at which to insert the new row.
//     </param>
//
//     <param
//       name="p_row"
//       kind="in,positional,required"
//       doc="vector"
//     >
//       A gsl_vector pointer containing the values for the new row.
//     </param>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="pre" id="row">
//       The row to be inserted must be within the matrix or must be one
//       larger than the number of rows in the matrix; that is, the row
//       will be inserted between two existing rows, inserted before the
//       first row, or tacked on at the end of the matrix.
//     </doc>
//
//     <doc kind="pre" id="vector">
//       The number of elements in the vector must be equal to the number of
//       columns in the matrix.
//     </doc>
//
//     <doc kind="post" id="result">
//       On successful return, the row has been inserted.  The number of
//       rows in the matrix is one more than upon input.
//       If any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Insert vector p_vector into matrix p_matrix at row 3:
//       </synopsis>
//
//       <code>
//  WnMatrix__insertRow( p_matrix, 3L, p_vector );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnMatrix__insertRow( WnMatrix *, size_t, gsl_vector * );

/*##############################################################################
// <routine name="WnMatrix__insertColumn()">
//
//   <description>
//     <abstract>
//       Inserts a column into a matrix.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, insert, column
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnMatrix__insertColumn(
//   WnMatrix *self,
//   size_t i_column,
//   gsl_vector *p_column
// );
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure.
//     </param>
//
//     <param
//       name="i_column"
//       kind="in,positional,required"
//       doc="column"
//     >
//       A size_t giving the column index at which to insert the new column.
//     </param>
//
//     <param
//       name="p_column"
//       kind="in,positional,required"
//       doc="vector"
//     >
//       A gsl_vector pointer containing the values for the new column.
//     </param>
//
//     <doc kind="pre" id="column">
//       The column to be inserted must be within the matrix or must be one
//       larger than the number of columns in the matrix; that is, the column
//       will be inserted between two existing columns, inserted before the
//       first column, or tacked on at the end of the matrix.
//     </doc>
//
//     <doc kind="pre" id="vector">
//       The number of elements in the vector must be equal to the number of
//       rows in the matrix.
//     </doc>
//
//     <param 
//       kind="return" 
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       On successful return, the column has been inserted.  The number of
//       columns in the matrix is one more than upon input.
//       If any input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Insert vector p_vector into matrix p_matrix at column 3:
//       </synopsis>
//
//       <code>
//  WnMatrix__insertColumn( p_matrix, 3L, p_vector );
//       </code>
//         
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnMatrix__insertColumn( WnMatrix *, size_t, gsl_vector * );

/*##############################################################################
// Non-API routines.
//############################################################################*/

void WnMatrix__freeElement( WnMatrix__Element *, xmlChar * );

void
WnMatrix__computeMatrixTimesVectorCallback(
  WnMatrix__Element *, void *, xmlChar *, xmlChar *, xmlChar *
);

void WnMatrix__computeTransposeMatrixTimesVectorCallback(
  WnMatrix__Element *, void *, xmlChar *, xmlChar *, xmlChar *
);

void WnMatrix__writeMatrixToAsciiFileCallback(
  WnMatrix__Element *, WnMatrix *, xmlChar *, xmlChar *, xmlChar *
);

void WnMatrix__getGslMatrixCallback(
  WnMatrix__Element *, gsl_matrix *, xmlChar *, xmlChar *, xmlChar *
);

void WnMatrix__insertMatrixCallback(
  WnMatrix__Element *, void *, xmlChar *, xmlChar *, xmlChar *
);

void WnMatrix__extractMatrixCallback(
  WnMatrix__Element *, void *, xmlChar *, xmlChar *, xmlChar *
);

void WnMatrix__getTransferMatrixCallback(
  WnMatrix__Element *, void *, xmlChar *, xmlChar *, xmlChar *
);

int
WnMatrix__data_compare(
  const void *, const void *
);

WnMatrix__Elements *
WnMatrix__getSortedElements( WnMatrix *, xmlChar *, xmlChar * );

void WnMatrix__countElements(
  WnMatrix__Element *, WnMatrix__Elements *, xmlChar *, xmlChar *, xmlChar *
);

void WnMatrix__getElementsCallback(
  WnMatrix__Element *, WnMatrix__Elements *, xmlChar *, xmlChar *, xmlChar *
);

void WnMatrix__Elements__free( WnMatrix__Elements * );

void
WnMatrix__scaleMatrixCallback(
  WnMatrix__Element *, 
  double *,
  xmlChar *,
  xmlChar *,
  xmlChar *
);

void
WnMatrix__getTransposeCallback(
  WnMatrix__Element *, 
  WnMatrix *,
  xmlChar *,
  xmlChar *,
  xmlChar *
);

WnMatrix__Element *
WnMatrix__Element__copier(
  WnMatrix__Element *,
  xmlChar *
);

xmlDocPtr
WnMatrix__Coo__makeXmlDocument(
  WnMatrix__Coo *,
  xmlChar *
);

void
WnMatrix__get_element_array(
  WnMatrix__Element *,
  void *,
  xmlChar *,
  xmlChar *,
  xmlChar *
);

xmlDocPtr
WnMatrix__Csr__makeXmlDocument(
  WnMatrix__Csr *,
  xmlChar *
);

xmlDocPtr
WnMatrix__Yale__makeXmlDocument(
  WnMatrix__Yale *,
  xmlChar *
);

xmlDocPtr
WnMatrix__make_gsl_vector_xml_document(
  gsl_vector *,
  const xmlChar *
);

void
WnMatrix__get_arrow_width_callback(
  WnMatrix__Element *, WnMatrix__Arrow *, xmlChar *, xmlChar *, xmlChar *
);

void
WnMatrix__arrow_assign_callback(
  WnMatrix__Element *, WnMatrix__Arrow *, xmlChar *, xmlChar *, xmlChar *
);

gsl_vector *
WnMatrix__Arrow__backsub(
  WnMatrix__Arrow *, gsl_vector *
);

gsl_vector *
WnMatrix__backsub(
  WnMatrix *,
  gsl_vector *
);

int
WnMatrix__value_is_zero( double );

size_t
WnMatrix__sx_to_sizet( xmlChar * );

void
WnMatrix__remove_row_callback(
  WnMatrix__Element *,
  void *,
  xmlChar *,
  xmlChar *,
  xmlChar *
);

void
WnMatrix__remove_column_callback(
  WnMatrix__Element *,
  void *,
  xmlChar *,
  xmlChar *,
  xmlChar *
);

void
WnMatrix__insert_row_callback(
  WnMatrix__Element *,
  void *,
  xmlChar *,
  xmlChar *,
  xmlChar *
);

void
WnMatrix__insert_column_callback(
  WnMatrix__Element *,
  void *,
  xmlChar *,
  xmlChar *,
  xmlChar *
);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* WN_MATRIX_H */
