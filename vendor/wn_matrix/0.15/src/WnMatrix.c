/*//////////////////////////////////////////////////////////////////////////////
// <file type = "public">
//   <description>
//     <abstract>
//       Source code for the wn_matrix module.  For documentation,
//       see the header file WnMatrix.h.
//     </abstract>
//   </description>
//   <license>
//      Copyright (c) 2006-2011 Clemson University.
//
//      This file contains the source code for the Clemson Webnucleo group's
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
//   </license>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Includes.
//############################################################################*/

#include <WnMatrix.h>

/*##############################################################################
// WnMatrix__new().
//############################################################################*/

WnMatrix *WnMatrix__new( size_t i_rows, size_t i_columns) {

  WnMatrix *self;
#ifdef WN_USE_MATRIX_LOOKUP
  xmlChar sx_tmp[WNMATRIX_BUF_SIZE];
  size_t i;
#endif

  self = ( WnMatrix * ) malloc( sizeof( WnMatrix ) );

  if( !self )
    return NULL;

  if( i_rows >= INDEX_MAX )
    WNMATRIX__ERROR( "Number of rows larger than can be stored" );

  if( i_columns >= INDEX_MAX )
    WNMATRIX__ERROR( "Number of columns larger than can be stored" );

  self->iRows = i_rows;
  self->iCols = i_columns;
  self->pMatrixHash = xmlHashCreate( 0 );

#ifdef WN_USE_MATRIX_LOOKUP

  if (i_rows < 1 || i_columns < 1) WNMATRIX__ERROR( "Invalid matrix size" );

  self->iMax = GSL_MAX( i_rows, i_columns );
 
  self->pLookup =
    ( xmlChar * ) malloc( self->iMax * WNMATRIX_BUF_SIZE * sizeof( xmlChar ) );

  if( !self->pLookup )
    WNMATRIX__ERROR( "Couldn't allocate look up array" );

  for (i = 0; i< self->iMax; ++i) {
    xmlStrPrintf( sx_tmp, WNMATRIX_BUF_SIZE, (const xmlChar *) "%lu", i+1 );
    memcpy(
      self->pLookup + ( i * WNMATRIX_BUF_SIZE ),
      sx_tmp,
      sizeof( xmlChar ) * WNMATRIX_BUF_SIZE
    );
  }
#endif

  return self;

}

/*##############################################################################
// WnMatrix__assignElement().
//############################################################################*/

int
WnMatrix__assignElement(
  WnMatrix * self,
  size_t i_row,
  size_t i_col,
  double d_val
)
{

#ifndef WN_USE_MATRIX_LOOKUP
  xmlChar sx_row[WNMATRIX_BUF_SIZE], sx_col[WNMATRIX_BUF_SIZE];
#endif
  WnMatrix__Element *p_element;

  if( !self ) {
    WNMATRIX__ERROR( "Invalid matrix pointer" );
  }

  if( i_row < 1 || i_row > self->iRows ) {
    WNMATRIX__ERROR( "Invalid row" );
  }

  if( i_col < 1 || i_col > self->iCols ) {
    WNMATRIX__ERROR( "Invalid column" );
  }

  if ( WnMatrix__value_is_zero( d_val ) )
    return 0;

#ifndef WN_USE_MATRIX_LOOKUP

  xmlStrPrintf( sx_row, WNMATRIX_BUF_SIZE, (const xmlChar *) "%lu", i_row );
  xmlStrPrintf( sx_col, WNMATRIX_BUF_SIZE, (const xmlChar *) "%lu", i_col );

  if( 
      (
        p_element = ( WnMatrix__Element * )
          xmlHashLookup2(
            self->pMatrixHash, sx_row, sx_col
          )
      )
  )
  {
    p_element->dValue += d_val;
    return 0;
  }

#else

  if( 
      (
        p_element = ( WnMatrix__Element * )
          xmlHashLookup2(
            self->pMatrixHash,
            &(self->pLookup[(i_row-1)*WNMATRIX_BUF_SIZE]),
            &(self->pLookup[(i_col-1)*WNMATRIX_BUF_SIZE])
          )
      )
  )
  {
    p_element->dValue += d_val;
    return 0;
  }

#endif

  p_element = (WnMatrix__Element *) malloc( sizeof( WnMatrix__Element ) );

  if( !p_element ) {
      WNMATRIX__ERROR( "Couldn't allocate memory for element" );
  }

  p_element->iRow = i_row;
  p_element->iCol = i_col;
  p_element->dValue = d_val;

#ifndef WN_USE_MATRIX_LOOKUP

  xmlHashAddEntry2(
    self->pMatrixHash,
    sx_row,
    sx_col,
    p_element
  );

#else

  xmlHashAddEntry2(
    self->pMatrixHash,
    &(self->pLookup[(i_row-1)*WNMATRIX_BUF_SIZE]),
    &(self->pLookup[(i_col-1)*WNMATRIX_BUF_SIZE]),
    p_element
  );

#endif

  return 0;

}

/*##############################################################################
/  WnMatrix__updateElement().
/#############################################################################*/

int
WnMatrix__updateElement(
  WnMatrix * self,
  size_t i_row,
  size_t i_col,
  double d_val
)
{

#ifndef WN_USE_MATRIX_LOOKUP
  xmlChar sx_row[WNMATRIX_BUF_SIZE], sx_col[WNMATRIX_BUF_SIZE];
#endif
  WnMatrix__Element *p_element;

  if( !self ) {
    WNMATRIX__ERROR( "Invalid matrix pointer" );
  }

  if( i_row < 1 || i_row > self->iRows ) {
    WNMATRIX__ERROR( "Invalid row" );
  }

  if( i_col < 1 || i_col > self->iCols ) {
    WNMATRIX__ERROR( "Invalid column" );
  }

  if ( WnMatrix__value_is_zero( d_val ) )
    return 0;

  p_element = (WnMatrix__Element *) malloc( sizeof( WnMatrix__Element ) );

  if( !p_element ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  p_element->iRow = i_row;
  p_element->iCol = i_col;
  p_element->dValue = d_val;

#ifndef WN_USE_MATRIX_LOOKUP

  xmlStrPrintf( sx_row, WNMATRIX_BUF_SIZE, (const xmlChar *) "%lu", i_row );
  xmlStrPrintf( sx_col, WNMATRIX_BUF_SIZE, (const xmlChar *) "%lu", i_col );

  return
    xmlHashUpdateEntry2(
      self->pMatrixHash,
      sx_row,
      sx_col,
      p_element,
      (xmlHashDeallocator) WnMatrix__freeElement
    );

#else

  return
    xmlHashUpdateEntry2(
      self->pMatrixHash,
      &(self->pLookup[(i_row-1)*WNMATRIX_BUF_SIZE]),
      &(self->pLookup[(i_col-1)*WNMATRIX_BUF_SIZE]),
      p_element,
      (xmlHashDeallocator) WnMatrix__freeElement
    );

#endif

}

/*##############################################################################
// WnMatrix__removeElement().
//############################################################################*/

int
WnMatrix__removeElement(
  WnMatrix * self, size_t i_row, size_t i_col
) {

  xmlChar sx_row[WNMATRIX_BUF_SIZE], sx_col[WNMATRIX_BUF_SIZE];

  xmlStrPrintf( sx_row, WNMATRIX_BUF_SIZE, (const xmlChar *) "%lu", i_row );
  xmlStrPrintf( sx_col, WNMATRIX_BUF_SIZE, (const xmlChar *) "%lu", i_col );

  return
    xmlHashRemoveEntry2(
      self->pMatrixHash,
      sx_row,
      sx_col,
      (xmlHashDeallocator) WnMatrix__freeElement
    );

}

/*##############################################################################
// WnMatrix__getElement().
//############################################################################*/

double
WnMatrix__getElement(
  WnMatrix *self, size_t i_row, size_t i_col
) {

  WnMatrix__Element *p_element; 
#ifndef WN_USE_MATRIX_LOOKUP
  xmlChar sx_row[WNMATRIX_BUF_SIZE], sx_col[WNMATRIX_BUF_SIZE];
#endif

  if( !self ) {
    WNMATRIX__ERROR( "Invalid matrix pointer" );
  }

  if( i_row < 1 || i_row > self->iRows ) {
    WNMATRIX__ERROR( "Invalid row" );
  }

  if( i_col < 1 || i_col > self->iCols ) {
    WNMATRIX__ERROR( "Invalid column" );
  }

#ifndef WN_USE_MATRIX_LOOKUP
  xmlStrPrintf( sx_row, WNMATRIX_BUF_SIZE, (const xmlChar *) "%lu", i_row );
  xmlStrPrintf( sx_col, WNMATRIX_BUF_SIZE, (const xmlChar *) "%lu", i_col );

  p_element = ( WnMatrix__Element * )
    xmlHashLookup3(
      self->pMatrixHash,
      sx_row,
      sx_col,
      NULL
    );
#else
  p_element =
    ( WnMatrix__Element * )
    xmlHashLookup2(
      self->pMatrixHash,
      &(self->pLookup[(i_row - 1) * WNMATRIX_BUF_SIZE]),
      &(self->pLookup[(i_col - 1) * WNMATRIX_BUF_SIZE])
    );
#endif

  if( p_element ) {
    return p_element->dValue;
  } else {
    return 0;
  }

}

/*##############################################################################
// WnMatrix__freeElement().
//############################################################################*/

void WnMatrix__freeElement(
  WnMatrix__Element *self, xmlChar *sx_row
) {

    if( !sx_row )
      WNMATRIX__ERROR( "Trying to free non-existent element" );

    free( self );

}

/*##############################################################################
// WnMatrix__free().
//############################################################################*/

void WnMatrix__free( WnMatrix *self ) {

  xmlHashFree( self->pMatrixHash, (xmlHashDeallocator) WnMatrix__freeElement );

#ifdef WN_USE_MATRIX_LOOKUP
  free( self->pLookup );
#endif

  free( self );
  self = NULL;

}

/*##############################################################################
// WnMatrix__clear()
//############################################################################*/

void WnMatrix__clear( WnMatrix *self )
{

  xmlHashFree( self->pMatrixHash, (xmlHashDeallocator) WnMatrix__freeElement );
  self->pMatrixHash = xmlHashCreate( 0 );

}

/*##############################################################################
// WnMatrix__getNumberOfRows().
//############################################################################*/

size_t WnMatrix__getNumberOfRows( WnMatrix *self ) {

  return self->iRows;

}

/*##############################################################################
// WnMatrix__getNumberOfColumns().
//############################################################################*/

size_t WnMatrix__getNumberOfColumns( WnMatrix *self ) {

  return self->iCols;

}

/*##############################################################################
// WnMatrix__getNumberOfElements().
//############################################################################*/

size_t WnMatrix__getNumberOfElements( WnMatrix *self ) {

  return (size_t) xmlHashSize( self->pMatrixHash );

}

/*##############################################################################
// WnMatrix__computeMatrixTimesVector()
//############################################################################*/

gsl_vector *
WnMatrix__computeMatrixTimesVector(
  WnMatrix *self,
  gsl_vector *p_input_vector
) {

  typedef struct extra_data {
    gsl_vector *pInputVector;
    gsl_vector *pOutputVector;
  } extra_data;

  gsl_vector *p_output_vector;
  extra_data *p_data;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( !self || !p_input_vector )
    WNMATRIX__ERROR( "Invalid input" );

  /*============================================================================
  // Check number of columns versus length of vector.
  //==========================================================================*/

  if(
    WnMatrix__getNumberOfColumns( self )
    != p_input_vector->size
  )
    WNMATRIX__ERROR(
      "Number of columns of matrix does not equal length of vector"
    );

  /*============================================================================
  // Get extra data structure and assign data.
  //==========================================================================*/

  p_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  if( !p_data )
    WNMATRIX__ERROR( "Couldn't allocate necessary memory" );

  p_data->pOutputVector =
    gsl_vector_calloc( p_input_vector->size );
  
  if( !p_data->pOutputVector )
    WNMATRIX__ERROR( "Couldn't allocate necessary memory" );

  p_data->pInputVector = p_input_vector;

  /*============================================================================
  // Callback to do multiplication.
  //==========================================================================*/

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__computeMatrixTimesVectorCallback,
    p_data
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  p_output_vector = p_data->pOutputVector;
  free( p_data );

  return p_output_vector;

}

/*##############################################################################
// WnMatrix__computeMatrixTimesVectorCallback()
//############################################################################*/

void
WnMatrix__computeMatrixTimesVectorCallback(
  WnMatrix__Element *p_element,
  void *p_extra_data,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  typedef struct extra_data {
    gsl_vector *pInputVector;
    gsl_vector *pOutputVector;
  } extra_data;
    
  extra_data *p_data;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );


  p_data = (extra_data *) p_extra_data;

  p_data->pOutputVector->data[p_element->iRow - 1] +=
    p_element->dValue *
    p_data->pInputVector->data[p_element->iCol - 1];

}

/*##############################################################################
// WnMatrix__computeTransposeMatrixTimesVector()
//############################################################################*/

gsl_vector *
WnMatrix__computeTransposeMatrixTimesVector(
  WnMatrix *self,
  gsl_vector *p_input_vector
) {

  typedef struct extra_data{
    gsl_vector *pInputVector;
    gsl_vector *pOutputVector;
  } extra_data;

  gsl_vector *p_output_vector;
  extra_data *p_data;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( !self || !p_input_vector )
    WNMATRIX__ERROR( "Invalid input" );

  /*============================================================================
  // Check number of rows versus length of vector.
  //==========================================================================*/

  if(
    WnMatrix__getNumberOfRows( self )
    != p_input_vector->size
  )
    WNMATRIX__ERROR(
      "Number of columns of transpose matrix does not equal length of vector"
    );

  /*============================================================================
  // Get extra data structure and assign data.
  //==========================================================================*/

  p_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  if( !p_data )
    WNMATRIX__ERROR( "Couldn't allocate necessary memory" );

  p_data->pOutputVector =
    gsl_vector_calloc( p_input_vector->size );
  
  if( !p_data->pOutputVector )
    WNMATRIX__ERROR( "Couldn't allocate necessary memory" );

  p_data->pInputVector = p_input_vector;

  /*============================================================================
  // Callback to do multiplication.
  //==========================================================================*/

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__computeTransposeMatrixTimesVectorCallback,
    p_data
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  p_output_vector = p_data->pOutputVector;
  free( p_data );

  return p_output_vector;

}

/*##############################################################################
// WnMatrix__computeTransposeMatrixTimesVectorCallback()
//############################################################################*/

void WnMatrix__computeTransposeMatrixTimesVectorCallback(
  WnMatrix__Element *p_element,
  void *p_extra_data,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  typedef struct extra_data {
    gsl_vector *pInputVector;
    gsl_vector *pOutputVector;
  } extra_data;
    
  extra_data *p_data;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );

  p_data = (extra_data *) p_extra_data;

  p_data->pOutputVector->data[p_element->iCol - 1] +=
    p_element->dValue *
    p_data->pInputVector->data[p_element->iRow - 1];

}

/*##############################################################################
// WnMatrix__writeMatrixToAsciiFile()
//############################################################################*/

int WnMatrix__writeMatrixToAsciiFile(
  WnMatrix *self,
  const char *s_ascii_filename,
  double d_cutoff
) {

  FILE *pFile;
  size_t i;
  WnMatrix__Elements *p_elements;

  if( !self ) {
    WNMATRIX__ERROR( "Invalid matrix pointer" );
  }

  if( ( pFile = fopen( s_ascii_filename, "w" ) ) == NULL ) {
      fprintf( stderr, "Can't open file %s\n", s_ascii_filename );
      exit( EXIT_FAILURE );
  }

  p_elements = WnMatrix__getSortedElements( self, NULL, NULL );

  for( i = 0; i < p_elements->iCount; i++ )
    if( fabs( p_elements->pData[i]->dValue ) > d_cutoff ) {
      fprintf(
        pFile,
        "%lu  %lu  %e\n",
        (unsigned long) p_elements->pData[i]->iRow,
        (unsigned long) p_elements->pData[i]->iCol,
        p_elements->pData[i]->dValue
      );
    }

  WnMatrix__Elements__free( p_elements );

  fclose( pFile );

  return 1;

}

/*##############################################################################
// WnMatrix__getGslMatrix()
//############################################################################*/

gsl_matrix *
WnMatrix__getGslMatrix( WnMatrix *self ) {

  gsl_matrix *p_dense;

  p_dense = gsl_matrix_calloc( self->iRows, self->iCols );

  if( !p_dense ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__getGslMatrixCallback,
    p_dense
  );


  return p_dense;

}

/*##############################################################################
// WnMatrix__getGslMatrixCallback()
//############################################################################*/

void WnMatrix__getGslMatrixCallback(
  WnMatrix__Element *p_element,
  gsl_matrix *p_dense,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );


  gsl_matrix_set(
    p_dense,
    p_element->iRow - 1,
    p_element->iCol - 1,
    p_element->dValue
  );

}

/*##############################################################################
// WnMatrix__insertMatrix()
//############################################################################*/

void WnMatrix__insertMatrix(
  WnMatrix *self,
  WnMatrix *p_inserted_matrix,
  size_t i_row,
  size_t i_col
) {

  typedef struct extra_data {
    WnMatrix *pMatrix;
    size_t iRow;
    size_t iCol;
  } extra_data;

  extra_data *p_data;

  if( !self || !p_inserted_matrix ) {
    WNMATRIX__ERROR( "Invalid matrix pointer" );
  }

  if( i_row < 1 || i_row > self->iRows ) {
    WNMATRIX__ERROR( "Invalid row" );
  }

  if( i_col < 1 || i_col > self->iCols ) {
    WNMATRIX__ERROR( "Invalid column" );
  }

  if( i_row + p_inserted_matrix->iRows - 1 > self->iRows ) {
    WNMATRIX__ERROR( "Inserted matrix will overflow final matrix" );
  }

  if( i_col + p_inserted_matrix->iCols - 1 > self->iCols ) {
    WNMATRIX__ERROR( "Inserted matrix will overflow final matrix" );
  }

  p_data = (extra_data *) malloc( sizeof( extra_data ) );

  if( !p_data ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  p_data->pMatrix = self;
  p_data->iRow = i_row;
  p_data->iCol = i_col;

  xmlHashScanFull(
    p_inserted_matrix->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__insertMatrixCallback,
    p_data
  );

  free( p_data );

}

/*##############################################################################
// WnMatrix__insertMatrixCallback()
//############################################################################*/

void WnMatrix__insertMatrixCallback(
  WnMatrix__Element *p_element,
  void *p_extra_data,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  typedef struct extra_data {
    WnMatrix *pMatrix;
    size_t iRow;
    size_t iCol;
  } extra_data;

  extra_data *p_data;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );

  p_data = (extra_data *) p_extra_data;

  WnMatrix__assignElement(
    p_data->pMatrix,
    p_element->iRow + p_data->iRow - 1,
    p_element->iCol + p_data->iCol - 1,
    p_element->dValue
  );

}

/*##############################################################################
// WnMatrix__extractMatrix()
//############################################################################*/

WnMatrix *WnMatrix__extractMatrix(
  WnMatrix *self,
  size_t i_row,
  size_t i_col,
  size_t i_row_offset,
  size_t i_col_offset
) {

  typedef struct extra_data {
    size_t iRow;
    size_t iRowOffset;
    size_t iCol;
    size_t iColOffset;
    WnMatrix *pMatrix;
  } extra_data;
    
  WnMatrix *p_extracted_matrix;
  extra_data *p_data;

  if( !self ) {
    WNMATRIX__ERROR( "Invalid matrix pointer" );
  }

  if( i_row < 1 || i_row > self->iRows ) {
    WNMATRIX__ERROR( "Invalid row" );
  }

  if( i_col < 1 || i_col > self->iCols ) {
    WNMATRIX__ERROR( "Invalid column" );
  }

  if( i_row_offset < 1 || i_row + i_row_offset - 1 > self->iRows ) {
    WNMATRIX__ERROR( "Invalid row offset" );
  }

  if( i_col_offset < 1 || i_col + i_col_offset - 1 > self->iCols ) {
    WNMATRIX__ERROR( "Invalid column offset" );
  }

  p_extracted_matrix = WnMatrix__new( i_row_offset, i_col_offset );

  if( !p_extracted_matrix ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  p_data = (extra_data *) malloc( sizeof( extra_data ) );

  if( !p_data )
      WNMATRIX__ERROR( "Virtual memory exhausted" );

  p_data->iRow = i_row;
  p_data->iRowOffset = i_row_offset;
  p_data->iCol = i_col;
  p_data->iColOffset = i_col_offset;
  p_data->pMatrix = p_extracted_matrix;

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__extractMatrixCallback,
    p_data
  );

  free( p_data );
  return p_extracted_matrix;

}

/*##############################################################################
// WnMatrix__extractMatrixCallback()
//############################################################################*/

void WnMatrix__extractMatrixCallback(
  WnMatrix__Element *p_element,
  void *p_extra_data,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
) {

  typedef struct extra_data {
    size_t iRow;
    size_t iRowOffset;
    size_t iCol;
    size_t iColOffset;
    WnMatrix *pMatrix;
  } extra_data;
    
  double d_value;
  size_t i_row2, i_col2;
  extra_data *p_data;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );


  p_data = (extra_data *) p_extra_data;

  i_row2 = p_element->iRow - p_data->iRow + 1;
  i_col2 = p_element->iCol - p_data->iCol + 1;

  d_value = p_element->dValue;

  if(
      i_row2 > 0 &&
      i_row2 <= p_data->iRowOffset &&
      i_col2 > 0 &&
      i_col2 <= p_data->iColOffset
  ) {

      WnMatrix__assignElement(
        p_data->pMatrix, i_row2, i_col2, d_value
      );
  }

}

/*##############################################################################
// WnMatrix__getDiagonalElements()
//############################################################################*/

gsl_vector *
WnMatrix__getDiagonalElements(
  WnMatrix *self
) {

  size_t i;
  gsl_vector *p_diagonals = gsl_vector_calloc( self->iRows );

  if( !p_diagonals ) {
    WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  for( i = 0; i < self->iRows; i++ ) {

    gsl_vector_set(
      p_diagonals,
      i,
      WnMatrix__getElement( self, i + 1, i + 1 )
    );

  }
  
  return p_diagonals;

}

/*##############################################################################
// WnMatrix__addValueToDiagonals()
//############################################################################*/

void WnMatrix__addValueToDiagonals( WnMatrix *self, double d_val ) {

  size_t i;

  for( i = 1L; i <= self->iRows; i++ ) {

    WnMatrix__assignElement( self, i, i, d_val );

  }

}

/*##############################################################################
// WnMatrix__scaleMatrix()
//############################################################################*/

void
WnMatrix__scaleMatrix( WnMatrix *self, double d_value ) {

  if( !self ) {
      WNMATRIX__ERROR( "Invalid matrix" );
  }

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__scaleMatrixCallback,
    &d_value
  );

}

/*##############################################################################
// WnMatrix__scaleMatrixCallback()
//############################################################################*/

void
WnMatrix__scaleMatrixCallback(
  WnMatrix__Element *p_element,
  double *p_value,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );

  p_element->dValue *= *p_value;

}

/*##############################################################################
// WnMatrix__getTransferMatrix()
//############################################################################*/

WnMatrix *WnMatrix__getTransferMatrix( WnMatrix *self ) {

  typedef struct extra_data {
    WnMatrix *pOriginalMatrix;
    WnMatrix *pTransferMatrix;
  } extra_data;

  WnMatrix *p_transfer_matrix;
  extra_data *p_extra_data;
  
  p_transfer_matrix = WnMatrix__new( self->iRows, self->iCols );

  p_extra_data = (extra_data *) malloc( sizeof( extra_data ) );

  if( !p_transfer_matrix || !p_extra_data ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  p_extra_data->pOriginalMatrix = self;
  p_extra_data->pTransferMatrix = p_transfer_matrix;

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__getTransferMatrixCallback,
    p_extra_data
  );

  free( p_extra_data );
  return p_transfer_matrix;

}

/*##############################################################################
// WnMatrix__getTransferMatrixCallback()
//############################################################################*/

void WnMatrix__getTransferMatrixCallback(
  WnMatrix__Element *p_element,
  void *p_extra_data,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  typedef struct extra_data {
    WnMatrix *pOriginalMatrix;
    WnMatrix *pTransferMatrix;
  } extra_data;

  double d_value;
  extra_data *p_data;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );

  p_data = (extra_data *) p_extra_data;

  d_value =
    p_element->dValue /
    WnMatrix__getElement(
      p_data->pOriginalMatrix,
      p_element->iCol,
      p_element->iCol
    );

  if( p_element->iRow != p_element->iCol ) {
    WnMatrix__assignElement(
      p_data->pTransferMatrix,
      p_element->iRow,
      p_element->iCol,
      d_value
    );
  }

}

/*##############################################################################
// WnMatrix__getRow()
//############################################################################*/

WnMatrix__Line *WnMatrix__getRow(
  WnMatrix *self, size_t i_row_number
){

  size_t i;
  WnMatrix__Elements *p_elements;
  WnMatrix__Line *p_row;
  xmlChar sx_row[WNMATRIX_BUF_SIZE];

  if( !self ) {
      WNMATRIX__ERROR( "Invalid matrix" );
  }

  if( i_row_number < 1 || i_row_number > self->iRows ) {
      WNMATRIX__ERROR( "Invalid row number" );
  }

  xmlStrPrintf(
    sx_row,
    WNMATRIX_BUF_SIZE,
    (const xmlChar *) "%lu",
    i_row_number
  );

  p_elements = WnMatrix__getSortedElements( self, sx_row, NULL );

  p_row = ( WnMatrix__Line * ) malloc( sizeof( WnMatrix__Line ) );

  if( !p_row ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  p_row->iCount = p_elements->iCount;

  p_row->a1 =
    ( size_t * ) malloc( sizeof( size_t ) * p_row->iCount );

  if( !p_row->a1 ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  p_row->a2 =
    ( double * ) malloc( sizeof( double ) * p_row->iCount );

  if( !p_row->a2 ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  for( i = 0; i < p_elements->iCount; i++ )
  {
    p_row->a1[i] = p_elements->pData[i]->iCol;
    p_row->a2[i] = p_elements->pData[i]->dValue;
  }

  WnMatrix__Elements__free( p_elements );

  return p_row;

}

/*##############################################################################
// WnMatrix__getColumn()
//############################################################################*/

WnMatrix__Line *WnMatrix__getColumn(
  WnMatrix *self, size_t i_column_number
){

  size_t i;
  WnMatrix__Elements *p_elements;
  WnMatrix__Line *p_column;
  xmlChar sx_column[WNMATRIX_BUF_SIZE];

  if( !self ) {
      WNMATRIX__ERROR( "Invalid matrix" );
  }

  if( i_column_number < 1 || i_column_number > self->iCols ) {
      WNMATRIX__ERROR( "Invalid row number" );
  }

  xmlStrPrintf(
    sx_column,
    WNMATRIX_BUF_SIZE,
    (const xmlChar *) "%lu",
    i_column_number
  );

  p_elements = WnMatrix__getSortedElements( self, NULL, sx_column );

  p_column = ( WnMatrix__Line * ) malloc( sizeof( WnMatrix__Line ) );

  if( !p_column ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  p_column->iCount = p_elements->iCount;

  p_column->a1 =
    ( size_t * ) malloc( sizeof( size_t ) * p_column->iCount );

  if( !p_column->a1 ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  p_column->a2 =
    ( double * ) malloc( sizeof( double ) * p_column->iCount );

  if( !p_column->a2 ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  for( i = 0; i < p_elements->iCount; i++ )
  {
    p_column->a1[i] = p_elements->pData[i]->iRow;
    p_column->a2[i] = p_elements->pData[i]->dValue;
  }

  WnMatrix__Elements__free( p_elements );

  return p_column;

}

/*##############################################################################
// WnMatrix__getSortedElements().
//############################################################################*/

WnMatrix__Elements *
WnMatrix__getSortedElements(
  WnMatrix *self,
  xmlChar *sx_row,
  xmlChar *sx_column
)
{

  WnMatrix__Elements *p_elements;

  p_elements = ( WnMatrix__Elements * ) malloc( sizeof( WnMatrix__Elements ) );

  if( !p_elements ) WNMATRIX__ERROR( "Couldn't allocate memory" );

  p_elements->iCount = 0;

  xmlHashScanFull3(
    self->pMatrixHash,
    sx_row,
    NULL,
    NULL,
    (xmlHashScannerFull) WnMatrix__countElements,
    p_elements
  );

  p_elements->pData =
    ( WnMatrix__Element ** )
    malloc( p_elements->iCount * sizeof( WnMatrix__Element ) );

  p_elements->iCount = 0;

  xmlHashScanFull3(
    self->pMatrixHash,
    sx_row,
    sx_column,
    NULL,
    (xmlHashScannerFull) WnMatrix__getElementsCallback,
    p_elements
  );

  qsort(
    p_elements->pData,
    p_elements->iCount,
    sizeof( WnMatrix__Element * ),
    WnMatrix__data_compare
  );

  return p_elements;

}

/*##############################################################################
// WnMatrix__countElements().
//############################################################################*/

void WnMatrix__countElements(
  WnMatrix__Element *p_element,
  WnMatrix__Elements *p_elements,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  if( (!sx_row && !sx_col) || sx_tmp || !p_element )
    WNMATRIX__ERROR( "Invalid input" );

  ++(p_elements->iCount);

}

/*##############################################################################
// WnMatrix__getElementsCallback()
//############################################################################*/

void WnMatrix__getElementsCallback(
  WnMatrix__Element *p_element,
  WnMatrix__Elements *p_elements,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  if( !sx_row || !sx_col || sx_tmp ) 
    WNMATRIX__ERROR( "Invalid input" );

  p_elements->pData[p_elements->iCount++] = p_element;

}

/*##############################################################################
// WnMatrix__Elements__free()
//############################################################################*/

void
WnMatrix__Elements__free( WnMatrix__Elements *self )
{

  free( self->pData );
  free( self );

}

/*##############################################################################
// WnMatrix__getTranspose()
//############################################################################*/

WnMatrix *WnMatrix__getTranspose( WnMatrix *self ) {

  WnMatrix *p_transpose_matrix;
  
  if( !self ) {
      WNMATRIX__ERROR( "Invalid matrix" );
  }
  
  p_transpose_matrix = WnMatrix__new( self->iCols, self->iRows );

  if( !p_transpose_matrix ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__getTransposeCallback,
    p_transpose_matrix
  );

  return p_transpose_matrix;

}

/*##############################################################################
// WnMatrix__getTransposeCallback()
//############################################################################*/

void WnMatrix__getTransposeCallback(
  WnMatrix__Element *p_element,
  WnMatrix *p_matrix,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );


  WnMatrix__assignElement(
    p_matrix,
    p_element->iCol,
    p_element->iRow,
    p_element->dValue
  );

}

/*##############################################################################
// WnMatrix__getCopy()
//############################################################################*/

WnMatrix *WnMatrix__getCopy( WnMatrix *self ) {

  WnMatrix *p_copy;

  if( !self ) {
      WNMATRIX__ERROR( "Invalid matrix" );
  }
  
  p_copy =
    ( WnMatrix * ) malloc( sizeof( WnMatrix ) );

  if( !p_copy ) {
      WNMATRIX__ERROR( "Virtual memory exhausted" );
  }

  p_copy->iRows = self->iRows;
  p_copy->iCols = self->iCols;

#ifdef WN_USE_MATRIX_LOOKUP
  p_copy->pLookup = ( xmlChar * ) malloc( self->iMax * WNMATRIX_BUF_SIZE );
  memcpy( p_copy->pLookup, self->pLookup, self->iMax * WNMATRIX_BUF_SIZE );
#endif

  p_copy->pMatrixHash = 
    xmlHashCopy(
      self->pMatrixHash,
      (xmlHashCopier) WnMatrix__Element__copier
    );

  return p_copy;

}

/*##############################################################################
// WnMatrix__Element__copier()
//############################################################################*/

WnMatrix__Element *
WnMatrix__Element__copier(
  WnMatrix__Element *p_element,
  xmlChar *sx_row
)
{

  WnMatrix__Element *p_new_element;

  if( !sx_row ) WNMATRIX__ERROR( "Invalid element" );

  p_new_element =
    ( WnMatrix__Element * ) malloc( sizeof( WnMatrix__Element ) );

  if( !p_new_element ) WNMATRIX__ERROR( "Couldn't allocate memory" );

  p_new_element->iRow = p_element->iRow;
  p_new_element->iCol = p_element->iCol;
  p_new_element->dValue = p_element->dValue;

  return p_new_element;

}

/*##############################################################################
// WnMatrix__getCoo()
//############################################################################*/

WnMatrix__Coo *WnMatrix__getCoo( WnMatrix *self )
{

  size_t i;
  WnMatrix__Coo *p_coo;
  WnMatrix__Elements *p_elements;
 
  p_elements = WnMatrix__getSortedElements( self, NULL, NULL );

  p_coo = ( WnMatrix__Coo * ) malloc( sizeof( WnMatrix__Coo ) );

  p_coo->iCount = WnMatrix__getNumberOfElements( self );

  p_coo->aRow = ( size_t * ) malloc( sizeof( size_t ) * p_coo->iCount );
  p_coo->aCol = ( size_t * ) malloc( sizeof( size_t ) * p_coo->iCount );
  p_coo->aVal = ( double * ) malloc( sizeof( double ) * p_coo->iCount );

  for( i = 0; i < p_coo->iCount; i++ )
  {
    p_coo->aRow[i] = p_elements->pData[i]->iRow;
    p_coo->aCol[i] = p_elements->pData[i]->iCol;
    p_coo->aVal[i] = p_elements->pData[i]->dValue;
  }

  WnMatrix__Elements__free( p_elements );

  return p_coo;

}

/*##############################################################################
// WnMatrix__data_compare().
//############################################################################*/

int
WnMatrix__data_compare(
  const void *p_data1, const void *p_data2
)
{

  const WnMatrix__Element *p_x = *(WnMatrix__Element * const *) p_data1;
  const WnMatrix__Element *p_y = *(WnMatrix__Element * const *) p_data2;

  if( p_x->iRow < p_y->iRow )
    return -1;
  else if( p_x->iRow > p_y->iRow )
    return 1;
  else
    if( p_x->iCol < p_y->iCol )
      return -1;
    else
      return 1;

}

/*##############################################################################
// WnMatrix__getCsr()
//############################################################################*/

WnMatrix__Csr *WnMatrix__getCsr( WnMatrix *self )
{

  WnMatrix__Elements *p_elements;
  WnMatrix__Csr *p_csr;
  size_t i, j, i_current_row;

  p_elements = WnMatrix__getSortedElements( self, NULL, NULL );

  p_csr = ( WnMatrix__Csr * ) malloc( sizeof( WnMatrix__Csr ) );

  p_csr->iCount = WnMatrix__getNumberOfElements( self );
  p_csr->iRowCount = WnMatrix__getNumberOfRows( self );

  p_csr->aRowptr =
    ( size_t * )
    malloc( ( WnMatrix__getNumberOfRows( self ) + 1 ) * sizeof( size_t ) );

  p_csr->aCol =
    ( size_t * )
    malloc( WnMatrix__getNumberOfElements( self ) * sizeof( size_t ) );

  p_csr->aVal =
    ( double * )
    malloc( WnMatrix__getNumberOfElements( self ) * sizeof( double ) );

  for( i = 0; i < p_elements->iCount; i++ )
  {
    p_csr->aCol[i] = p_elements->pData[i]->iCol;
    p_csr->aVal[i] = p_elements->pData[i]->dValue;
  }

  i_current_row = 1L;
  for( j = 0; j < p_elements->pData[0]->iRow; j++ ) p_csr->aRowptr[j] = 0;

  for( i = 0; i < WnMatrix__getNumberOfElements( self ); i++ )
  {
    if( p_elements->pData[i]->iRow != i_current_row )
    {
      for( j = i_current_row; j < p_elements->pData[i]->iRow; j++ )
        p_csr->aRowptr[j] = i;
      i_current_row = p_elements->pData[i]->iRow;
    }
  }

  p_csr->aRowptr[WnMatrix__getNumberOfRows( self )] =
     WnMatrix__getNumberOfElements( self );

  WnMatrix__Elements__free( p_elements );

  return p_csr;

}

/*##############################################################################
// WnMatrix__getYale()
//############################################################################*/

WnMatrix__Yale *WnMatrix__getYale(
  WnMatrix *self
)
{

  size_t i, j, i_terms, i_count, i_current_row;
  WnMatrix__Yale *p_yale;
  WnMatrix__Elements *p_elements;
 
  p_elements = WnMatrix__getSortedElements( self, NULL, NULL );

  i_terms =
    WnMatrix__getNumberOfRows( self ) +
    WnMatrix__getNumberOfElements( self ) +
    1;

  p_yale = ( WnMatrix__Yale * ) malloc( sizeof( WnMatrix__Yale ) ) ;

  if( !p_yale )
    WNMATRIX__ERROR( "Couldn't allocate memory" );

  p_yale->aIja =
    ( size_t * )
    calloc(
      i_terms, sizeof( size_t )
    );
   
  p_yale->aVal =
    ( double * )
    calloc(
      i_terms, sizeof( double )
    );

  if( !p_yale->aIja || !p_yale->aVal )
    WNMATRIX__ERROR( "Couldn't allocate memory" );

  p_yale->aIja[0] = self->iRows + 1;

  i_current_row = 1;
  i_count = self->iRows + 1;
  for( i = 0; i < WnMatrix__getNumberOfElements( self ); i++ )
  {
    if( p_elements->pData[i]->iRow == p_elements->pData[i]->iCol )
      p_yale->aVal[p_elements->pData[i]->iRow - 1] =
        p_elements->pData[i]->dValue;
    else
    {
      if( p_elements->pData[i]->iRow > i_current_row )
      {
        for(
          j = i_current_row; 
          j < p_elements->pData[i]->iRow;
          j++
        )
          p_yale->aIja[j] = i_count;
        i_current_row = p_elements->pData[i]->iRow;
      }
      p_yale->aIja[i_count] = p_elements->pData[i]->iCol;
      p_yale->aVal[i_count++] = p_elements->pData[i]->dValue;
    }
  }

  p_yale->aIja[WnMatrix__getNumberOfRows(self)] = i_count;

  WnMatrix__Elements__free( p_elements );

  return p_yale;

}

/*##############################################################################
// WnMatrix__Line__getNumberOfElements()
//############################################################################*/

size_t WnMatrix__Line__getNumberOfElements( WnMatrix__Line *self ) {

  return self->iCount;

}

/*##############################################################################
// WnMatrix__Line__getNonZeroIndices()
//############################################################################*/

size_t *WnMatrix__Line__getNonZeroIndices( WnMatrix__Line *self ) {

  return self->a1;

}

/*##############################################################################
// WnMatrix__Line__getNonZeroElements()
//############################################################################*/

double *WnMatrix__Line__getNonZeroElements( WnMatrix__Line *self ) {

  return self->a2;

}

/*##############################################################################
// WnMatrix__Line__free()
//############################################################################*/

void WnMatrix__Line__free( WnMatrix__Line *self ) {

  free( self->a1 );
  free( self->a2 );
  free( self );

}

/*##############################################################################
// WnMatrix__Csr__free()
//############################################################################*/

void WnMatrix__Csr__free( WnMatrix__Csr * self ) {
  
  free( self->aRowptr );
  free( self->aCol );
  free( self->aVal );
  free( self );

}

/*##############################################################################
// WnMatrix__Coo__free()
//############################################################################*/

void WnMatrix__Coo__free( WnMatrix__Coo * self ) {
  
  free( self->aRow );
  free( self->aCol );
  free( self->aVal );
  free( self );

}

/*##############################################################################
// WnMatrix__Yale__free()
//############################################################################*/

void WnMatrix__Yale__free( WnMatrix__Yale * self ) {
  
  free( self->aIja );
  free( self->aVal );
  free( self );

}

/*##############################################################################
// WnMatrix__Coo__getRowVector()
//############################################################################*/

size_t *WnMatrix__Coo__getRowVector( WnMatrix__Coo * self ) {

  if( !self ) {
    WNMATRIX__ERROR( "Invalid coordinate matrix" );
  }

  return self->aRow;

}

/*##############################################################################
// WnMatrix__Coo__getColumnVector()
//############################################################################*/

size_t *WnMatrix__Coo__getColumnVector( WnMatrix__Coo * self ) {

  if( !self ) {
    WNMATRIX__ERROR( "Invalid coordinate matrix" );
  }

  return self->aCol;

}

/*##############################################################################
// WnMatrix__Coo__getValueVector()
//############################################################################*/

double *WnMatrix__Coo__getValueVector( WnMatrix__Coo * self ) {

  if( !self ) {
    WNMATRIX__ERROR( "Invalid coordinate matrix" );
  }

  return self->aVal;

}

/*##############################################################################
// WnMatrix__Csr__getRowPointerVector()
//############################################################################*/

size_t *WnMatrix__Csr__getRowPointerVector(
  WnMatrix__Csr * self
) {

  if( !self ) {
    WNMATRIX__ERROR( "Invalid coordinate matrix" );
  }

  return self->aRowptr;

}

/*##############################################################################
// WnMatrix__Csr__getColumnVector()
//############################################################################*/

size_t *WnMatrix__Csr__getColumnVector( WnMatrix__Csr * self ) {

  if( !self ) {
    WNMATRIX__ERROR( "Invalid coordinate matrix" );
  }

  return self->aCol;

}

/*##############################################################################
// WnMatrix__Csr__getValueVector()
//############################################################################*/

double *WnMatrix__Csr__getValueVector( WnMatrix__Csr * self ) {

  if( !self ) {
    WNMATRIX__ERROR( "Invalid coordinate matrix" );
  }

  return self->aVal;

}

/*##############################################################################
// WnMatrix__Yale__getPointerVector()
//############################################################################*/

size_t *WnMatrix__Yale__getPointerVector( WnMatrix__Yale * self ) {

  if( !self ) {
    WNMATRIX__ERROR( "Invalid coordinate matrix" );
  }

  return self->aIja;

}

/*##############################################################################
// WnMatrix__Yale__getValueVector()
//############################################################################*/

double *WnMatrix__Yale__getValueVector( WnMatrix__Yale * self ) {

  if( !self ) {
    WNMATRIX__ERROR( "Invalid coordinate matrix" );
  }

  return self->aVal;

}

/*##############################################################################
// WnMatrix__Coo__writeToXmlFile()
//############################################################################*/

void WnMatrix__Coo__writeToXmlFile(
  WnMatrix__Coo *self,
  const char *s_output_xml_filename,
  const char *s_format
)
{

  xmlDocPtr p_doc = NULL;
  xmlChar *sx_format;

  /*============================================================================
  // Check input structure.
  //==========================================================================*/

  if( !self ) {
     WNMATRIX__ERROR( "Invalid input structure" );
  }

  /*============================================================================
  // Get element value format.
  //==========================================================================*/

  if( s_format )
    sx_format = xmlCharStrdup( s_format );
  else
    sx_format = xmlCharStrdup( VALUE_FORMAT );

  /*============================================================================
  // Get document.
  //==========================================================================*/

  p_doc = WnMatrix__Coo__makeXmlDocument( self, sx_format );

  /*============================================================================
  // Write file.
  //==========================================================================*/

  if(
       xmlSaveFormatFileEnc( s_output_xml_filename, p_doc, "UTF-8", 1 ) == -1
  ) {
      WNMATRIX__ERROR(
        "DOMImplementation.saveDocToFile: failed\n"
      );
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlFree( sx_format );
  xmlFreeDoc( p_doc );
  xmlCleanupParser();

}

/*##############################################################################
// WnMatrix__Coo__makeXmlDocument()
//############################################################################*/

xmlDocPtr WnMatrix__Coo__makeXmlDocument(
  WnMatrix__Coo *self,
  xmlChar *sx_format
)
{

  xmlDocPtr p_doc;
  xmlNodePtr p_root, p_element;
  size_t i;
  xmlChar sx_data[WNMATRIX_BUF_SIZE];
  xmlChar *sx_str, *sx_str2;

  /*============================================================================
  // Check input structure.
  //==========================================================================*/

  if( !self ) {
     WNMATRIX__ERROR( "Invaid input structure" );
  }

  /*============================================================================
  // Get document.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XML_VERSION );
  p_doc = xmlNewDoc( sx_str );
  xmlFree( sx_str );

  if (p_doc == NULL) {
    WNMATRIX__ERROR( "DOMImplementation.createDocument: failed" );
  }

  sx_str = xmlCharStrdup( COO_MATRIX_WITH_PREFIX );
  p_root = xmlNewNode( NULL, sx_str );
  xmlFree( sx_str );

  xmlDocSetRootElement( p_doc, p_root );

  /*============================================================================
  // Create namespace attributes and add.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XSI_DEC );
  sx_str2 = xmlCharStrdup( W3C__NAMESPACE );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( WNMATRIX_COO_DEC );
  sx_str2 = xmlCharStrdup( WNMATRIX__COO__NAMESPACE );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( XSI_LOCATION );
  sx_str2 = xmlCharStrdup( WNMATRIX__COO__SCHEMALOCATION );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  /*============================================================================
  // Loop on elements.
  //==========================================================================*/

  for( i = 0; i < self->iCount; i++ ) {

  /*---------------------------------------------------------------------------
  // Get element.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_ELEMENT );
    p_element = xmlNewChild( p_root, NULL, sx_str, NULL );
    xmlFree( sx_str );

    if( p_element == NULL ) {
      WNMATRIX__ERROR( "Document.createElement: NULL\n\tException" );
    }

  /*---------------------------------------------------------------------------
  // Insert row.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_ROW );

    xmlStrPrintf(
      sx_data,
      WNMATRIX_BUF_SIZE,
      (const xmlChar *) "%lu",
      self->aRow[i]
    );

    xmlNewChild( p_element, NULL, sx_str, sx_data );

    xmlFree( sx_str );

  /*---------------------------------------------------------------------------
  // Insert column.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_COLUMN );

    xmlStrPrintf(
      sx_data,
      WNMATRIX_BUF_SIZE,
      (const xmlChar *) "%lu",
      self->aCol[i]
    );

    xmlNewChild( p_element, NULL, sx_str, sx_data );

    xmlFree( sx_str );

  /*---------------------------------------------------------------------------
  // Insert value
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_VALUE );

    xmlStrPrintf(
      sx_data,
      WNMATRIX_BUF_SIZE,
      sx_format,
      self->aVal[i]
    );

    xmlNewChild( p_element, NULL, sx_str, sx_data );

    xmlFree( sx_str );

  }

  return p_doc;

}

/*##############################################################################
// WnMatrix__Csr__writeToXmlFile()
//############################################################################*/

void
WnMatrix__Csr__writeToXmlFile(
  WnMatrix__Csr *self,
  const char *s_output_xml_filename,
  const char *s_format
)
{

  xmlDocPtr p_doc = NULL;
  xmlChar *sx_format;

  /*============================================================================
  // Check input structure.
  //==========================================================================*/

  if( !self ) {
     WNMATRIX__ERROR( "Invalid input structure" );
  }

  /*============================================================================
  // Get element value format.
  //==========================================================================*/

  if( s_format )
    sx_format = xmlCharStrdup( s_format );
  else
    sx_format = xmlCharStrdup( VALUE_FORMAT );

  /*============================================================================
  // Get document.
  //==========================================================================*/

  p_doc = WnMatrix__Csr__makeXmlDocument( self, sx_format );

  /*============================================================================
  // Write file.
  //==========================================================================*/

  if(
       xmlSaveFormatFileEnc( s_output_xml_filename, p_doc, "UTF-8", 1 ) == -1
  ) {
      WNMATRIX__ERROR(
        "DOMImplementation.saveDocToFile: failed\n"
      );
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlFree( sx_format );
  xmlFreeDoc( p_doc );
  xmlCleanupParser();

}

/*##############################################################################
// WnMatrix__Csr__makeXmlDocument()
//############################################################################*/

xmlDocPtr WnMatrix__Csr__makeXmlDocument(
  WnMatrix__Csr *self,
  xmlChar *sx_format
)
{

  xmlDocPtr p_doc;
  xmlNodePtr p_root, p_row_pointers, p_element, p_elements;
  size_t i;
  xmlChar sx_data[WNMATRIX_BUF_SIZE];
  xmlChar *sx_str, *sx_str2;

  /*============================================================================
  // Check input structure.
  //==========================================================================*/

  if( !self ) {
     WNMATRIX__ERROR( "Invalid input structure" );
  }

  /*============================================================================
  // Get document.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XML_VERSION );
  p_doc = xmlNewDoc( sx_str );
  xmlFree( sx_str );

  if (p_doc == NULL) {
    WNMATRIX__ERROR( "DOMImplementation.createDocument: failed" );
  }

  sx_str = xmlCharStrdup( CSR_MATRIX_WITH_PREFIX );
  p_root = xmlNewNode( NULL, sx_str );
  xmlFree( sx_str );

  if( p_root == NULL ) {
    WNMATRIX__ERROR( "Couldn't create root" );
  }

  xmlDocSetRootElement( p_doc, p_root );

  /*============================================================================
  // Create namespace attributes and add.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XSI_DEC );
  sx_str2 = xmlCharStrdup( W3C__NAMESPACE );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( WNMATRIX_CSR_DEC );
  sx_str2 = xmlCharStrdup( WNMATRIX__CSR__NAMESPACE );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( XSI_LOCATION );
  sx_str2 = xmlCharStrdup( WNMATRIX__CSR__SCHEMALOCATION );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  /*============================================================================
  // Row pointers.
  //==========================================================================*/

  sx_str = xmlCharStrdup( WNMATRIX_ROW_POINTERS );
  p_row_pointers = xmlNewChild( p_root, NULL, sx_str, NULL );
  xmlFree( sx_str );

  if( p_row_pointers == NULL ) {
    WNMATRIX__ERROR( "Document.createElement: NULL\n\tException" );
  }

  for( i = 0; i <= self->iRowCount; i++ ) {

  /*---------------------------------------------------------------------------
  // Insert row pointer.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_ROW_POINTER );

    xmlStrPrintf(
      sx_data,
      WNMATRIX_BUF_SIZE,
      (const xmlChar *) "%lu",
      self->aRowptr[i] + 1L
    );

    xmlNewChild( p_row_pointers, NULL, sx_str, sx_data );

    xmlFree( sx_str );

  }

  /*============================================================================
  // Elements.
  //==========================================================================*/

  sx_str = xmlCharStrdup( WNMATRIX_ELEMENTS );
  p_elements = xmlNewChild( p_root, NULL, sx_str, NULL );
  xmlFree( sx_str );

  if( p_elements == NULL ) {
    WNMATRIX__ERROR( "Document.createElement: NULL\n\tException" );
  }

  for( i = 0; i < self->iCount; i++ ) {

  /*---------------------------------------------------------------------------
  // Get element.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_ELEMENT );
    p_element = xmlNewChild( p_elements, NULL, sx_str, NULL );
    xmlFree( sx_str );

    if( p_element == NULL ) {
      WNMATRIX__ERROR( "Document.createElement: NULL\n\tException" );
    }

  /*---------------------------------------------------------------------------
  // Insert column.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_COLUMN );

    xmlStrPrintf(
      sx_data,
      WNMATRIX_BUF_SIZE,
      (const xmlChar *) "%lu",
      self->aCol[i]
    );

    xmlNewChild( p_element, NULL, sx_str, sx_data );

    xmlFree( sx_str );

  /*---------------------------------------------------------------------------
  // Insert value
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_VALUE );

    xmlStrPrintf(
      sx_data,
      WNMATRIX_BUF_SIZE,
      sx_format,
      self->aVal[i]
    );

    xmlNewChild( p_element, NULL, sx_str, sx_data );

    xmlFree( sx_str );

  }

  return p_doc;

}

/*##############################################################################
// WnMatrix__Yale__writeToXmlFile()
//############################################################################*/

void
WnMatrix__Yale__writeToXmlFile(
  WnMatrix__Yale *self,
  const char *s_output_xml_filename,
  const char *s_format
)
{

  xmlDocPtr p_doc = NULL;
  xmlChar *sx_format;

  /*============================================================================
  // Check input structure.
  //==========================================================================*/

  if( !self ) {
     WNMATRIX__ERROR( "Invalid input structure" );
  }

  /*============================================================================
  // Get element value format.
  //==========================================================================*/

  if( s_format )
    sx_format = xmlCharStrdup( s_format );
  else
    sx_format = xmlCharStrdup( VALUE_FORMAT );

  /*============================================================================
  // Get document.
  //==========================================================================*/

  p_doc = WnMatrix__Yale__makeXmlDocument( self, sx_format );

  /*============================================================================
  // Write file.
  //==========================================================================*/

  if(
       xmlSaveFormatFileEnc( s_output_xml_filename, p_doc, "UTF-8", 1 ) == -1
  ) {
      WNMATRIX__ERROR(
        "DOMImplementation.saveDocToFile: failed\n"
      );
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlFree( sx_format );
  xmlFreeDoc( p_doc );
  xmlCleanupParser();

}

/*##############################################################################
// WnMatrix__Yale__makeXmlDocument()
//############################################################################*/

xmlDocPtr WnMatrix__Yale__makeXmlDocument(
  WnMatrix__Yale *self,
  xmlChar *sx_format
)
{

  xmlDocPtr p_doc;
  xmlNodePtr p_root, p_element;
  size_t i;
  xmlChar *sx_str, *sx_str2;
  xmlChar sx_data[WNMATRIX_BUF_SIZE];

  /*============================================================================
  // Check input structure.
  //==========================================================================*/

  if( !self ) {
     WNMATRIX__ERROR( "Invalid input structure" );
  }

  /*============================================================================
  // Get document.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XML_VERSION );
  p_doc = xmlNewDoc( sx_str );
  xmlFree( sx_str );

  if (p_doc == NULL) {
    WNMATRIX__ERROR( "DOMImplementation.createDocument: failed" );
  }

  sx_str = xmlCharStrdup( YALE_MATRIX_WITH_PREFIX );
  p_root = xmlNewNode( NULL, sx_str );
  xmlFree( sx_str );

  if( p_root == NULL ) {
    WNMATRIX__ERROR( "Couldn't create root" );
  }

  xmlDocSetRootElement( p_doc, p_root );

  /*============================================================================
  // Create namespace attributes and add.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XSI_DEC );
  sx_str2 = xmlCharStrdup( W3C__NAMESPACE );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( WNMATRIX_YALE_DEC );
  sx_str2 = xmlCharStrdup( WNMATRIX__YALE__NAMESPACE );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( XSI_LOCATION );
  sx_str2 = xmlCharStrdup( WNMATRIX__YALE__SCHEMALOCATION );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  /*============================================================================
  // Add entries.
  //==========================================================================*/

  for( i = 0; i < self->aIja[self->aIja[0]-1]; i++ ) {

  /*---------------------------------------------------------------------------
  // Get element.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_ELEMENT );
    p_element =
      xmlNewChild( p_root, NULL, sx_str, NULL );
    xmlFree( sx_str );

    if( p_element == NULL ) {
      WNMATRIX__ERROR( "Document.createElement: NULL\n\tException" );
    }

  /*---------------------------------------------------------------------------
  // Insert ija.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_IJA );

    xmlStrPrintf(
      sx_data,
      WNMATRIX_BUF_SIZE,
      (const xmlChar *) "%lu",
      self->aIja[i]
    );

    xmlNewChild( p_element, NULL, sx_str, sx_data );

    xmlFree( sx_str );

  /*---------------------------------------------------------------------------
  // Insert sa
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_SA );

    xmlStrPrintf(
      sx_data,
      WNMATRIX_BUF_SIZE,
      sx_format,
      self->aVal[i]
    );

    xmlNewChild( p_element, NULL, sx_str, sx_data );

    xmlFree( sx_str );

  }

  return p_doc;

}

/*##############################################################################
// WnMatrix__new_from_xml().
//############################################################################*/

WnMatrix *WnMatrix__new_from_xml(
  const char *s_xml_filename, const char *s_xpath_suffix
 ) {

  xmlDocPtr p_doc;
  xmlXPathContextPtr p_xpathCtx, p_xpathCtx_data;
  xmlXPathObjectPtr p_xpathObj, p_xpathObj_data;
  xmlChar *sx_xpath, *sx_xpath_suffix, *sx_data, *sx_str;

  size_t *a_row, *a_col, i_max, i_row_max, i_col_max, i;
  double *a_val;

  WnMatrix *self;

  /*============================================================================
  // Create document
  //==========================================================================*/

  p_doc = xmlParseFile( s_xml_filename );

  if ( p_doc == NULL ) {

    WNMATRIX__ERROR( "Document could not be opened" );

   }

  /*============================================================================
  // Create main xpath expression.
  //==========================================================================*/

  sx_xpath = xmlCharStrdup( WNMATRIX_XPATH );

  if( !sx_xpath ) {
     WNMATRIX__ERROR( "Couldn't allocate memory for xpath string" );
  }

  sx_xpath_suffix = xmlCharStrdup( s_xpath_suffix );

  if( sx_xpath_suffix ) {

     sx_xpath = xmlStrcat( sx_xpath, sx_xpath_suffix );

     if( !sx_xpath ) {
        WNMATRIX__ERROR( "Couldn't allocate memory for xpath string" );
     }

     xmlFree( sx_xpath_suffix );

  }

  /*============================================================================
  // Get element XPath.
  //==========================================================================*/

  p_xpathCtx = xmlXPathNewContext( p_doc );
  p_xpathObj =
    xmlXPathEvalExpression( sx_xpath, p_xpathCtx );

  xmlFree( sx_xpath );

  if( !p_xpathObj ) {
    WNMATRIX__ERROR( "Invalid xpath expression" );
  }

  if( !p_xpathObj->nodesetval ) {
    WNMATRIX__ERROR( "No matrix data" );
  }

  /*============================================================================
  // Set relevant integers.
  //==========================================================================*/

  i_max = (size_t) p_xpathObj->nodesetval->nodeNr;
  i_row_max = 0;
  i_col_max = 0;

  /*============================================================================
  // Allocate memory for entries.
  //==========================================================================*/

  a_row =
    ( size_t * ) malloc( sizeof( size_t ) * i_max );

  a_col =
    ( size_t * ) malloc( sizeof( size_t ) * i_max );

  a_val =
    ( double * ) malloc( sizeof( double ) * i_max );

  /*============================================================================
  // Loop over elements.
  //==========================================================================*/

  for( i = 0; i < i_max; i++ ) {

    p_xpathCtx_data = xmlXPathNewContext( p_doc );
    p_xpathCtx_data->node = p_xpathObj->nodesetval->nodeTab[i];

  /*---------------------------------------------------------------------------
  // Get rows.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_ROW );
    p_xpathObj_data =
      xmlXPathEvalExpression( sx_str, p_xpathCtx_data );
    xmlFree( sx_str );

    if( !p_xpathObj_data ) {
      WNMATRIX__ERROR( "Invalid xpath expression" );
    }
    if( p_xpathObj_data->nodesetval->nodeNr == 0 ) {
      WNMATRIX__ERROR( "no row found" );
    }

    sx_data =
      xmlNodeGetContent( p_xpathObj_data->nodesetval->nodeTab[0] );

    a_row[i] = WnMatrix__sx_to_sizet( sx_data );
    if( a_row[i] > i_row_max ) i_row_max = a_row[i];
    xmlXPathFreeObject( p_xpathObj_data );
    xmlFree( sx_data );

  /*---------------------------------------------------------------------------
  // Get columns.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_COLUMN );
    p_xpathObj_data =
      xmlXPathEvalExpression( sx_str, p_xpathCtx_data );
    xmlFree( sx_str );

    if( !p_xpathObj_data ) {
      WNMATRIX__ERROR( "Invalid xpath expression" );
    }
    if( p_xpathObj_data->nodesetval->nodeNr == 0 ) {
      WNMATRIX__ERROR( "no column found" );
    }

    sx_data =
      xmlNodeGetContent( p_xpathObj_data->nodesetval->nodeTab[0] );

    a_col[i] = WnMatrix__sx_to_sizet( sx_data );
    if( a_col[i] > i_col_max ) i_col_max = a_col[i];
    xmlXPathFreeObject( p_xpathObj_data );
    xmlFree( sx_data );

  /*---------------------------------------------------------------------------
  // Get values.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_VALUE );
    p_xpathObj_data =
      xmlXPathEvalExpression( sx_str, p_xpathCtx_data );
    xmlFree( sx_str );

    if( !p_xpathObj_data ) {
      WNMATRIX__ERROR( "Invalid xpath expression" );
    }
    if( p_xpathObj_data->nodesetval->nodeNr == 0 ) {
      WNMATRIX__ERROR( "no value found" );
    }

    sx_data =
      xmlNodeGetContent( p_xpathObj_data->nodesetval->nodeTab[0] );

    a_val[i] = atof( (char *) sx_data );
    xmlXPathFreeObject( p_xpathObj_data );
    xmlFree( sx_data );

    xmlXPathFreeContext( p_xpathCtx_data );

  }

  /*============================================================================
  // Now set up matrix.
  //==========================================================================*/

  self = WnMatrix__new( i_row_max, i_col_max );

  for( i = 0; i < i_max; i++ ) {
    WnMatrix__assignElement( self, a_row[i], a_col[i], a_val[i] );
  }
  
  /*============================================================================
  // Clean up temporary arrays.
  //==========================================================================*/

  free( a_row );
  free( a_col );
  free( a_val );

  /*============================================================================
  // Clean up XPath.
  //==========================================================================*/

  xmlXPathFreeObject( p_xpathObj );
  xmlXPathFreeContext( p_xpathCtx );
  
  /*============================================================================
  // Free document.
  //==========================================================================*/
  
  xmlFreeDoc( p_doc );
  xmlCleanupParser();

  /*============================================================================
  // Done.  Return.
  //==========================================================================*/

  return self;

}

/*##############################################################################
// WnMatrix__is_valid_input_xml()
//############################################################################*/

int
WnMatrix__is_valid_input_xml( const char *s_filename ) {

  xmlDocPtr p_doc;
  xmlSchemaParserCtxtPtr p_parser_ctxt;
  xmlSchemaValidCtxtPtr p_valid_ctxt;
  xmlSchemaPtr p_schema;
  int i_valid;

  p_parser_ctxt = xmlSchemaNewParserCtxt( WNMATRIX__COO__SCHEMA );

  if( !(
        p_schema = xmlSchemaParse( p_parser_ctxt )
      )
  ) {
      WNMATRIX__ERROR( "Schema not found or invalid" );
  }

  p_valid_ctxt = xmlSchemaNewValidCtxt( p_schema );

  p_doc = xmlParseFile( s_filename );

  if( xmlSchemaValidateDoc( p_valid_ctxt, p_doc ) == 0 ) {
     i_valid = 1;
  } else {
     i_valid = 0;
  }

  xmlFreeDoc( p_doc );
  xmlSchemaFreeValidCtxt( p_valid_ctxt );
  xmlSchemaFree( p_schema );
  xmlSchemaFreeParserCtxt( p_parser_ctxt );

  xmlCleanupParser();

  return i_valid;

}

/*##############################################################################
// WnMatrix__new_gsl_vector_from_xml().
//############################################################################*/

gsl_vector *
WnMatrix__new_gsl_vector_from_xml(
  const char *s_xml_filename,
  const char *s_xpath_suffix
)
{

  gsl_vector *self;
  xmlDocPtr p_doc;
  size_t i;
  xmlXPathContextPtr p_xpathCtx;
  xmlXPathObjectPtr p_xpathObj;
  xmlChar *sx_xpath, *sx_xpath_suffix, *sx_data;

  /*============================================================================
  // Create document.
  //==========================================================================*/

  p_doc = xmlParseFile( s_xml_filename );

  if ( p_doc == NULL ) {

    WNMATRIX__ERROR( "Document could not be opened" );

   }

  /*============================================================================
  // Create main xpath expression.
  //==========================================================================*/

  sx_xpath = xmlCharStrdup( WNMATRIX_VECTOR_XPATH );

  if( !sx_xpath ) {
     WNMATRIX__ERROR( "Couldn't allocate memory for xpath string" );
  }

  sx_xpath_suffix = xmlCharStrdup( s_xpath_suffix );

  if( sx_xpath_suffix ) {

     sx_xpath = xmlStrcat( sx_xpath, sx_xpath_suffix );

     if( !sx_xpath ) {
        WNMATRIX__ERROR( "Couldn't allocate memory for xpath string" );
     }

     xmlFree( sx_xpath_suffix );

  }

  /*============================================================================
  // Get element XPath.
  //==========================================================================*/

  p_xpathCtx = xmlXPathNewContext( p_doc );
  p_xpathObj =
    xmlXPathEvalExpression( sx_xpath, p_xpathCtx );

  xmlFree( sx_xpath );

  if( !p_xpathObj ) {
    WNMATRIX__ERROR( "Invalid xpath expression" );
  }

  if( !p_xpathObj->nodesetval ) {
    WNMATRIX__ERROR( "No vector data" );
  }

  /*============================================================================
  // Loop over elements and assign to vector.
  //==========================================================================*/

  self =
    gsl_vector_alloc( (size_t) p_xpathObj->nodesetval->nodeNr );

  for( i = 0; i < (size_t) p_xpathObj->nodesetval->nodeNr; i++ ) {
  
    sx_data =
      xmlNodeGetContent( p_xpathObj->nodesetval->nodeTab[i] );

    gsl_vector_set(
      self,
      i,
      atof( (char *) sx_data )
    );

    xmlFree( sx_data );

  }    

  /*============================================================================
  // Clean up XPath.
  //==========================================================================*/

  xmlXPathFreeObject( p_xpathObj );
  xmlXPathFreeContext( p_xpathCtx );
  
  /*============================================================================
  // Free document.
  //==========================================================================*/
  
  xmlFreeDoc( p_doc );
  xmlCleanupParser();

  /*============================================================================
  // Done.  Return.
  //==========================================================================*/

  return self;

}

/*##############################################################################
// WnMatrix__get_gsl_vector_size().
//############################################################################*/

size_t
WnMatrix__get_gsl_vector_size(
  gsl_vector *self
)
{

  if( !self ) WNMATRIX__ERROR( "Invalid input" );

  return self->size;

}

/*##############################################################################
// WnMatrix__get_gsl_vector_array().
//############################################################################*/

double *
WnMatrix__get_gsl_vector_array(
  gsl_vector *self
)
{

  if( !self ) WNMATRIX__ERROR( "Invalid input" );

  return self->data;

}

/*##############################################################################
// WnMatrix__write_gsl_vector_to_xml_file()
//############################################################################*/

void
WnMatrix__write_gsl_vector_to_xml_file(
  gsl_vector *self,
  const char *s_output_xml_filename,
  const char *s_format
)
{

  xmlDocPtr p_doc = NULL;
  xmlChar *sx_format;

  /*============================================================================
  // Check input structure.
  //==========================================================================*/

  if( !self ) {
     WNMATRIX__ERROR( "Invalid input structure" );
  }

  /*============================================================================
  // Get format.
  //==========================================================================*/

  if( s_format )
    sx_format = xmlCharStrdup( s_format );
  else
    sx_format = xmlCharStrdup( VALUE_FORMAT );

  /*============================================================================
  // Get document.
  //==========================================================================*/

  p_doc = WnMatrix__make_gsl_vector_xml_document( self, sx_format );

  /*============================================================================
  // Write file.
  //==========================================================================*/

  if(
       xmlSaveFormatFileEnc( s_output_xml_filename, p_doc, "UTF-8", 1 ) == -1
  ) {
      WNMATRIX__ERROR(
        "DOMImplementation.saveDocToFile: failed\n"
      );
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlFree( sx_format );
  xmlFreeDoc( p_doc );
  xmlCleanupParser();

}

/*##############################################################################
// WnMatrix__make_gsl_vector_xml_document()
//############################################################################*/

xmlDocPtr WnMatrix__make_gsl_vector_xml_document(
  gsl_vector *self, const xmlChar *sx_format
)
{

  xmlDocPtr p_doc;
  xmlNodePtr p_root;
  size_t i;
  xmlChar sx_data[WNMATRIX_BUF_SIZE];
  xmlChar *sx_str, *sx_str2;

  /*============================================================================
  // Check input structure.
  //==========================================================================*/

  if( !self ) {
     WNMATRIX__ERROR( "Invalid input structure" );
  }

  /*============================================================================
  // Get document.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XML_VERSION );
  p_doc = xmlNewDoc( sx_str );
  xmlFree( sx_str );

  if (p_doc == NULL) {
    WNMATRIX__ERROR( "DOMImplementation.createDocument: failed" );
  }

  sx_str = xmlCharStrdup( VECTOR_WITH_PREFIX );
  p_root = xmlNewNode( NULL, sx_str );
  xmlFree( sx_str );

  xmlDocSetRootElement( p_doc, p_root );

  /*============================================================================
  // Create namespace attributes and add.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XSI_DEC );
  sx_str2 = xmlCharStrdup( W3C__NAMESPACE );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( WNMATRIX_VECTOR_DEC );
  sx_str2 = xmlCharStrdup( WNMATRIX__VECTOR__NAMESPACE );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( XSI_LOCATION );
  sx_str2 = xmlCharStrdup( WNMATRIX__VECTOR__SCHEMALOCATION );
  xmlNewProp( p_root, sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  /*============================================================================
  // Loop on elements.
  //==========================================================================*/

  for( i = 0; i < self->size; i++ ) {

  /*---------------------------------------------------------------------------
  // Get element.
  //--------------------------------------------------------------------------*/

    sx_str = xmlCharStrdup( WNMATRIX_VECTOR_ELEMENT );

  /*---------------------------------------------------------------------------
  // Insert value
  //--------------------------------------------------------------------------*/

    xmlStrPrintf(
      sx_data,
      WNMATRIX_BUF_SIZE,
      sx_format,
      gsl_vector_get( self, i )
    );

    xmlNewChild( p_root, NULL, sx_str, sx_data );

    xmlFree( sx_str );

  }

  return p_doc;

}

/*##############################################################################
// WnMatrix__is_valid_vector_input_xml()
//############################################################################*/

int WnMatrix__is_valid_vector_input_xml( const char *s_filename ) {

  xmlDocPtr p_doc;
  xmlSchemaParserCtxtPtr p_parser_ctxt;
  xmlSchemaValidCtxtPtr p_valid_ctxt;
  xmlSchemaPtr p_schema;
  int i_valid;

  p_parser_ctxt = xmlSchemaNewParserCtxt( WNMATRIX__VECTOR__SCHEMA );

  if( !(
        p_schema = xmlSchemaParse( p_parser_ctxt )
      )
  ) {
      WNMATRIX__ERROR( "Schema not found or invalid" );
  }

  p_valid_ctxt = xmlSchemaNewValidCtxt( p_schema );

  p_doc = xmlParseFile( s_filename );

  if( xmlSchemaValidateDoc( p_valid_ctxt, p_doc ) == 0 ) {
     i_valid = 1;
  } else {
     i_valid = 0;
  }

  xmlFreeDoc( p_doc );
  xmlSchemaFreeValidCtxt( p_valid_ctxt );
  xmlSchemaFree( p_schema );
  xmlSchemaFreeParserCtxt( p_parser_ctxt );

  xmlCleanupParser();

  return i_valid;

}

/*##############################################################################
// WnMatrix__solve()
//############################################################################*/

gsl_vector *
WnMatrix__solve( WnMatrix *p_matrix, gsl_vector *p_b ) {

  size_t i_dims;
  int s;
  gsl_matrix *p_m;
  gsl_vector *p_x;
  gsl_permutation *p_perm;
  
  i_dims = WnMatrix__getNumberOfRows( p_matrix );

  /*==========================================================================
  // Allocate memory for gsl vectors.
  //========================================================================*/

  p_x = gsl_vector_alloc( i_dims );
  p_perm = gsl_permutation_alloc( i_dims );

  /*==========================================================================
  // Assign elements.
  //========================================================================*/

  p_m = WnMatrix__getGslMatrix( p_matrix );

  /*==========================================================================
  // Solve matrix.
  //========================================================================*/

  gsl_linalg_LU_decomp( p_m, p_perm, &s );
  gsl_linalg_LU_solve( p_m, p_perm, p_b, p_x );
  
  /*==========================================================================
  // Free allocated memory
  //========================================================================*/

  gsl_permutation_free( p_perm );
  gsl_matrix_free( p_m );

  /*==========================================================================
  // Return.
  //========================================================================*/

  return p_x;

}

/*##############################################################################
// WnMatrix__Arrow__solve().
//############################################################################*/

gsl_vector *
WnMatrix__Arrow__solve( WnMatrix__Arrow *self, gsl_vector *p_u )
{

  size_t i, j, k, i_len;
  double d_gam;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( !self || !p_u )
    WNMATRIX__ERROR( "Invalid input" );

  if( self->iRows != p_u->size )
    WNMATRIX__ERROR( "Arrow matrix and rhs vector must have the same size" );

  /*============================================================================
  // Get length vector.
  //==========================================================================*/

  i_len = self->iRows - self->iWingWidth;

  /*============================================================================
  // Triangularize matrix.
  //==========================================================================*/

  for( i = 0; i < i_len; i++ )
  {
    if( WnMatrix__value_is_zero( self->a[i][self->iBand] ) )
      WNMATRIX__ERROR( "Zero pivot encountered" );
    for( j = 1L; j <= self->iBand; j++ )
    {
      if( i + j < i_len )
      {
        if( !WnMatrix__value_is_zero( self->a[i+j][self->iBand-j] ) )
        {
          d_gam =
            -self->a[i+j][self->iBand - j] /
               self->a[i][self->iBand];
          p_u->data[i+j] += d_gam * p_u->data[i];
          for( k = self->iBand; k < self->iBandWidth; k++ )
            self->a[i+j][k-j] += d_gam * self->a[i][k];
          for( k = 0; k < self->iWingWidth; k++ )
            self->b[k][i+j] += d_gam * self->b[k][i];
        }
      }
    }
    for( j = 0; j < self->iWingWidth; j++ )
    {
      if( !WnMatrix__value_is_zero( self->c[j][i] ) )
      {
        d_gam = -self->c[j][i] / self->a[i][self->iBand];
        p_u->data[i_len+j] += d_gam * p_u->data[i];
        for( k = self->iBand; k < self->iBandWidth; k++ )
        {
          if( i + k - self->iBand < i_len )
            self->c[j][i+k-self->iBand] += d_gam * self->a[i][k];
        }
        for( k = 0; k < self->iWingWidth; k++ )
          self->d[j][k] += d_gam * self->b[k][i];
      }
    }

  }

  if( self->iWingWidth > 1 )
  {
    for( i = 0; i < self->iWingWidth; i++ )
    {
      if( WnMatrix__value_is_zero( self->d[i][i] ) )
        WNMATRIX__ERROR( "Zero pivot encountered" );
      for( j = i + 1L; j < self->iWingWidth; j++ )
      {
        if( !WnMatrix__value_is_zero( self->d[j][i] ) )
        {
          d_gam = -self->d[j][i] / self->d[i][i];
          p_u->data[i_len+j] += d_gam * p_u->data[i_len+i];
          for( k = i; k < self->iWingWidth; k++ )
            self->d[j][k] += d_gam * self->d[i][k];
        }
      }
    }
  }

  /*============================================================================
  // Return back substitution.
  //==========================================================================*/

  return WnMatrix__Arrow__backsub( self, p_u );

}
  
/*##############################################################################
// WnMatrix__get_arrow_width_callback().
//############################################################################*/

void
WnMatrix__get_arrow_width_callback(
  WnMatrix__Element *p_element,
  WnMatrix__Arrow *p_arrow,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  size_t i_row, i_col, i_width;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );

  if( !p_element )
    WNMATRIX__ERROR( "Invalid matrix entry" );

  i_row = p_element->iRow;
  i_col = p_element->iCol;
  i_width = (size_t) labs( (long int) (i_col - i_row ) );

  if( i_row + p_arrow->iWingWidth > p_arrow->iRows )
    return;

  if( i_col + p_arrow->iWingWidth > p_arrow->iRows )
    return;

  if( i_width > p_arrow->iBand )
    p_arrow->iBand = i_width;

}

/*##############################################################################
// WnMatrix__arrow_assign_callback().
//############################################################################*/

void
WnMatrix__arrow_assign_callback(
  WnMatrix__Element *p_element,
  WnMatrix__Arrow *p_arrow,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  size_t i_row, i_col, i_w1, i_w2;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );

  i_row = p_element->iRow;
  i_col = p_element->iCol;

  if( i_row + p_arrow->iWingWidth <= p_arrow->iRows )
  {
    if( ( i_col + p_arrow->iWingWidth ) <= p_arrow->iRows )
    {
      i_w1 = i_col - i_row + p_arrow->iBand;
      p_arrow->a[i_row - 1][i_w1] = p_element->dValue;
    }
    else
    {
      i_w1 = i_col - p_arrow->iRows + p_arrow->iWingWidth - 1;
      p_arrow->b[i_w1][i_row - 1] = p_element->dValue; 
    }
  }
  else
  {  
    if( i_col + p_arrow->iWingWidth <= p_arrow->iRows )
    {
      i_w1 = i_row - p_arrow->iRows + p_arrow->iWingWidth - 1;
      p_arrow->c[i_w1][i_col-1] = p_element->dValue;
    }
    else
    {
      i_w1 = i_row - p_arrow->iRows + p_arrow->iWingWidth - 1;
      i_w2 = i_col - p_arrow->iRows + p_arrow->iWingWidth - 1;
      p_arrow->d[i_w1][i_w2] = p_element->dValue; 
    }
  }

}

/*##############################################################################
// WnMatrix__Arrow__backsub().
//############################################################################*/

gsl_vector *
WnMatrix__Arrow__backsub(
  WnMatrix__Arrow *self,
  gsl_vector *p_b
)
{

  size_t i, j, k, i_len;
  gsl_vector *p_x;

  i_len = self->iRows - self->iWingWidth;

  p_x = gsl_vector_calloc( p_b->size );

  if( self->iWingWidth > 0 )
  {
    for( k = self->iWingWidth; k > 0; k-- )
    {
      i = k - 1;
      p_x->data[i_len + i] = p_b->data[i_len + i];
      for( j = i + 1L; j < self->iWingWidth; j++ )
        p_x->data[i_len+i] -= self->d[i][j] * p_x->data[i_len+j];
      p_x->data[i_len+i] /= self->d[i][i];
    }
  }

  for( k = i_len; k > 0; k-- )
  {
    i = k - 1;
    p_x->data[i] = p_b->data[i];
    for( j = self->iBand + 1L; j < self->iBandWidth; j++ )
      if( i + j < self->iBand + i_len )
        p_x->data[i] -= self->a[i][j] * p_x->data[i+j-self->iBand];
    for( j = 0; j < self->iWingWidth; j++ )
      p_x->data[i] -= self->b[j][i] * p_x->data[i_len+j];
    p_x->data[i] /= self->a[i][self->iBand];
  }

  return p_x;

}

/*##############################################################################
// WnMatrix__getArrow().
//############################################################################*/

WnMatrix__Arrow *
WnMatrix__getArrow( WnMatrix *self, size_t i_wing_width )
{

  WnMatrix__Arrow *p_arrow;
  size_t i, i_len;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( !self )
    WNMATRIX__ERROR( "Invalid input" );

  if( i_wing_width >= WnMatrix__getNumberOfRows( self ) )
    WNMATRIX__ERROR( "Arrow wing width must be < number of rows in matrix" );

  if(
   WnMatrix__getNumberOfRows( self ) != WnMatrix__getNumberOfColumns( self )
  )
    WNMATRIX__ERROR( "Number of rows must equal number of columns" );

  /*============================================================================
  // Get arrow matrix.
  //==========================================================================*/

  p_arrow = (WnMatrix__Arrow *) malloc( sizeof( WnMatrix__Arrow ) );

  if( !p_arrow )
    WNMATRIX__ERROR( "Couldn't allocate arrow matrix" );

  p_arrow->iBand = 0;
  p_arrow->iRows = WnMatrix__getNumberOfRows( self ); 
  p_arrow->iWingWidth = i_wing_width;

  /*============================================================================
  // Get length parameter.
  //==========================================================================*/

  i_len = p_arrow->iRows - p_arrow->iWingWidth;

  /*============================================================================
  // Get band width.
  //==========================================================================*/

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__get_arrow_width_callback,
    p_arrow
  );

  p_arrow->iBandWidth = 2 * p_arrow->iBand + 1;

  /*============================================================================
  // Allocate a matrix.
  //==========================================================================*/

  p_arrow->a = ( double ** ) calloc( i_len, sizeof( double ) );
  
  if( !p_arrow->a )
    WNMATRIX__ERROR( "Couldn't allocate memory" );

  for( i = 0; i < i_len; i++ )
  {
    p_arrow->a[i] =
      ( double * ) calloc( p_arrow->iBandWidth, sizeof( double ) );
    if( !p_arrow->a[i] )
      WNMATRIX__ERROR( "Couldn't allocate memory" );
  }

  /*============================================================================
  // Allocate small matrices.
  //==========================================================================*/

  if( p_arrow->iWingWidth > 0 )
  {

    /*-------------------------------------------------------------------------
    // Allocate b matrix.
    //------------------------------------------------------------------------*/

    p_arrow->b = ( double ** ) calloc( p_arrow->iWingWidth, sizeof( double ) );

    for( i = 0; i < p_arrow->iWingWidth; i++ )
    {
      p_arrow->b[i] = ( double * ) calloc( i_len, sizeof( double ) );
      if( !p_arrow->b[i] )
        WNMATRIX__ERROR( "Couldn't allocate memory" );
    }

    /*-------------------------------------------------------------------------
    // Allocate c matrix.
    //------------------------------------------------------------------------*/

    p_arrow->c = ( double ** ) calloc( p_arrow->iWingWidth, sizeof( double ) );

    for( i = 0; i < p_arrow->iWingWidth; i++ )
    {
      p_arrow->c[i] =
        ( double * ) calloc( i_len, sizeof( double ) );
      if( !p_arrow->c[i] )
        WNMATRIX__ERROR( "Couldn't allocate memory" );
    }

    /*-------------------------------------------------------------------------
    // Allocate d matrix.
    //------------------------------------------------------------------------*/

    p_arrow->d = ( double ** ) calloc( p_arrow->iWingWidth, sizeof( double ) );

    for( i = 0; i < p_arrow->iWingWidth; i++ )
    {
      p_arrow->d[i] =
        ( double * ) calloc( p_arrow->iWingWidth, sizeof( double ) );
      if( !p_arrow->d[i] )
        WNMATRIX__ERROR( "Couldn't allocate memory" );
    }

  }

  /*============================================================================
  // Assign arrow matrix.
  //==========================================================================*/

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__arrow_assign_callback,
    p_arrow
  );

  /*============================================================================
  // Done. Return.
  //==========================================================================*/

  return p_arrow;

}

/*##############################################################################
// WnMatrix__Arrow__free().
//############################################################################*/

void
WnMatrix__Arrow__free( WnMatrix__Arrow *self )
{

  size_t i, i_len;

  i_len = self->iRows - self->iWingWidth;

  for( i = 0; i < i_len; i++ )
    free( self->a[i] );

  free( self->a );

  if( self->iWingWidth > 0 )
  {

    for( i = 0; i < self->iWingWidth; i++ )
      free( self->b[i] );

    free( self->b );

    for( i = 0; i < self->iWingWidth; i++ )
      free( self->c[i] );

    free( self->c );

    for( i = 0; i < self->iWingWidth; i++ )
      free( self->d[i] );

    free( self->d );

  }

  free( self );

}

/*##############################################################################
// WnMatrix__Arrow__getBandWidth().
//############################################################################*/

size_t
WnMatrix__Arrow__getBandWidth( WnMatrix__Arrow *self )
{

  if( !self )
    WNMATRIX__ERROR( "Invalid input" );

  return self->iBandWidth;

}

/*##############################################################################
// WnMatrix__Arrow__getWingWidth().
//############################################################################*/

size_t
WnMatrix__Arrow__getWingWidth( WnMatrix__Arrow *self )
{

  if( !self )
    WNMATRIX__ERROR( "Invalid input" );

  return self->iWingWidth;

}

/*##############################################################################
// WnMatrix__Arrow__getNumberOfRows().
//############################################################################*/

size_t
WnMatrix__Arrow__getNumberOfRows( WnMatrix__Arrow *self )
{

  if( !self )
    WNMATRIX__ERROR( "Invalid input" );

  return self->iRows;

}

/*##############################################################################
// WnMatrix__Arrow__getWnMatrix().
//############################################################################*/

WnMatrix *
WnMatrix__Arrow__getWnMatrix( WnMatrix__Arrow *self )
{

  WnMatrix *p_matrix;
  size_t i, j, i_len;

  if( !self )
    WNMATRIX__ERROR( "Invalid input" );

  i_len = self->iRows - self->iWingWidth;

  p_matrix =
    WnMatrix__new( 
      self->iRows,
      self->iRows
    );

  for( i = 0; i < i_len; i++ )
    for( j = 0; j < self->iBandWidth; j++ )
      if(
        i + j >= self->iBand &&
        i + j < i_len + self->iBand
      )
      {
        WnMatrix__assignElement(
          p_matrix,
          i + 1,
          i + j - self->iBand + 1,
          self->a[i][j]
        );
      }

  for( i = 0; i < self->iWingWidth; i++ )
    for( j = 0; j < i_len; j++ )
      WnMatrix__assignElement(
        p_matrix,
        j + 1,
        i_len + i + 1,
        self->b[i][j]
      );

  for( i = 0; i < self->iWingWidth; i++ )
    for( j = 0; j < i_len; j++ )
      WnMatrix__assignElement(
        p_matrix,
        i_len + i + 1,
        j + 1,
        self->c[i][j]
      );

  for( i = 0; i < self->iWingWidth; i++ )
    for( j = 0; j < self->iWingWidth; j++ )
      WnMatrix__assignElement(
        p_matrix,
        i_len + i + 1,
        i_len + j + 1,
        self->d[i][j]
      );

  return p_matrix;

}

/*##############################################################################
// WnMatrix__value_is_zero().
//############################################################################*/

int
WnMatrix__value_is_zero( double d_x )
{

  if( GSL_SIGN( d_x ) == GSL_SIGN( -d_x ) )
    return 1;
  else
    return 0;

}

/*##############################################################################
// WnMatrix__sx_to_sizet().
//############################################################################*/

size_t
WnMatrix__sx_to_sizet( xmlChar *sx_str )
{

  return (size_t) atol( (char *) sx_str );

}

/*##############################################################################
// WnMatrix__removeRow().
//############################################################################*/

void
WnMatrix__removeRow( WnMatrix *self, size_t i_row_to_remove )
{

  typedef struct {
    size_t iRow;
    WnMatrix *pMatrix;
  } user_data;

  user_data *p_user_data;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( i_row_to_remove > WnMatrix__getNumberOfRows( self ) )
    WNMATRIX__ERROR( "Row to remove is not in matrix" );

  /*============================================================================
  // Copy matrix to new one and shift rows.
  //==========================================================================*/

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data ) WNMATRIX__ERROR( "Couldn't allocate memory" );

  p_user_data->iRow = i_row_to_remove;

  p_user_data->pMatrix =
    WnMatrix__new(
      WnMatrix__getNumberOfRows( self ) - 1,
      WnMatrix__getNumberOfColumns( self )
    );

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__remove_row_callback,
    p_user_data
  );

  /*============================================================================
  // Clear old matrix hash and copy new matrix hash back to old.
  //==========================================================================*/

  xmlHashFree( self->pMatrixHash, (xmlHashDeallocator) WnMatrix__freeElement );

  self->pMatrixHash = p_user_data->pMatrix->pMatrixHash;

  self->iRows--;

#ifdef WN_USE_MATRIX_LOOKUP
  free( self->pLookup );
  self->pLookup = p_user_data->pMatrix->pLookup;
#endif

  /*============================================================================
  // Clean up and done.  Don't do a WnMatrix__free because we don't want to
  // free the matrix hash.
  //==========================================================================*/

  free( p_user_data->pMatrix );
  free( p_user_data );

}

/*##############################################################################
// WnMatrix__remove_row_callback().
//############################################################################*/

void
WnMatrix__remove_row_callback(
  WnMatrix__Element *p_element,
  void *p_data,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  typedef struct {
    size_t iRow;
    WnMatrix *pMatrix;
  } user_data;

  user_data *p_user_data = ( user_data * ) p_data;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );

  if( p_element->iRow == p_user_data->iRow )
    return;

  if( p_element->iRow < p_user_data->iRow )
    WnMatrix__assignElement(
      p_user_data->pMatrix,
      p_element->iRow,
      p_element->iCol,
      p_element->dValue
    );
  else
    WnMatrix__assignElement(
      p_user_data->pMatrix,
      p_element->iRow - 1,
      p_element->iCol,
      p_element->dValue
    );

}
  
/*##############################################################################
// WnMatrix__removeColumn().
//############################################################################*/

void
WnMatrix__removeColumn( WnMatrix *self, size_t i_column_to_remove )
{

  typedef struct {
    size_t iColumn;
    WnMatrix *pMatrix;
  } user_data;

  user_data *p_user_data;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( i_column_to_remove > WnMatrix__getNumberOfColumns( self ) )
    WNMATRIX__ERROR( "Column to remove is not in matrix" );

  /*============================================================================
  // Copy matrix to new one and shift columns.
  //==========================================================================*/

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data ) WNMATRIX__ERROR( "Couldn't allocate memory" );

  p_user_data->iColumn = i_column_to_remove;
  p_user_data->pMatrix =
    WnMatrix__new(
      WnMatrix__getNumberOfRows( self ),
      WnMatrix__getNumberOfColumns( self ) - 1
    );

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__remove_column_callback,
    p_user_data
  );

  /*============================================================================
  // Clear old matrix hash and set old matrix hash to new.
  //==========================================================================*/

  xmlHashFree( self->pMatrixHash, (xmlHashDeallocator) WnMatrix__freeElement );

  self->pMatrixHash = p_user_data->pMatrix->pMatrixHash;

  self->iCols--;

#ifdef WN_USE_MATRIX_LOOKUP
  free( self->pLookup );
  self->pLookup = p_user_data->pMatrix->pLookup;
#endif

  /*============================================================================
  // Clean up and done.  Don't do a WnMatrix__free because we don't want to
  // free the matrix hash.
  //==========================================================================*/

  free( p_user_data->pMatrix );
  free( p_user_data );

}

/*##############################################################################
// WnMatrix__remove_column_callback().
//############################################################################*/

void
WnMatrix__remove_column_callback(
  WnMatrix__Element *p_element,
  void *p_data,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  typedef struct {
    size_t iColumn;
    WnMatrix *pMatrix;
  } user_data;

  user_data *p_user_data = ( user_data * ) p_data;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );

  if( p_element->iCol == p_user_data->iColumn )
    return;

  if( p_element->iCol < p_user_data->iColumn )
    WnMatrix__assignElement(
      p_user_data->pMatrix,
      p_element->iRow,
      p_element->iCol,
      p_element->dValue
    );
  else
    WnMatrix__assignElement(
      p_user_data->pMatrix,
      p_element->iRow,
      p_element->iCol - 1,
      p_element->dValue
    );

}

/*##############################################################################
// WnMatrix__insertRow().
//############################################################################*/

void
WnMatrix__insertRow(
  WnMatrix *self, size_t i_row_to_insert, gsl_vector *p_vector
) {

  typedef struct {
    size_t iRow;
    WnMatrix *pMatrix;
  } user_data;

  user_data *p_user_data;
  size_t i;
  double d_val;

  /*============================================================================
  // Check row to insert.
  //==========================================================================*/

  if( i_row_to_insert > WnMatrix__getNumberOfRows( self ) + 1 )
    WNMATRIX__ERROR( "Row to insert is too large" );

  /*============================================================================
  // Check input vector for appropriate size.
  //==========================================================================*/

  if(
    WnMatrix__get_gsl_vector_size( p_vector ) !=
    WnMatrix__getNumberOfColumns( self ) 
  ) {
    WNMATRIX__ERROR( "Incorrect size for gsl_vector" );
  }

  /*============================================================================
  // Copy matrix to new one and shift rows.
  //==========================================================================*/

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data ) WNMATRIX__ERROR( "Couldn't allocate memory" );

  p_user_data->iRow = i_row_to_insert;

  p_user_data->pMatrix =
    WnMatrix__new(
      WnMatrix__getNumberOfRows( self ) + 1,
      WnMatrix__getNumberOfColumns( self )
    );

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__insert_row_callback,
    p_user_data
  );

  /*============================================================================
  // Clear old matrix hash and set old matrix hash to new.
  //==========================================================================*/

  xmlHashFree( self->pMatrixHash, (xmlHashDeallocator) WnMatrix__freeElement );

  self->pMatrixHash = p_user_data->pMatrix->pMatrixHash;

  self->iRows++;

#ifdef WN_USE_MATRIX_LOOKUP
  free( self->pLookup );
  self->pLookup = p_user_data->pMatrix->pLookup;
#endif

  /*============================================================================
  // Insert new values from vector.
  //==========================================================================*/

  for( i = 1; i <= self->iCols; ++i ) {
    d_val = gsl_vector_get( p_vector, i - 1 );

    if( !WnMatrix__value_is_zero( d_val ) )
    {
      WnMatrix__assignElement(
        self,
        i_row_to_insert,
        i,
        d_val
      );
    }
  }

  /*============================================================================
  // Clean up and done.  Don't do a WnMatrix__free because we don't want to
  // free the matrix hash.
  //==========================================================================*/

  free( p_user_data->pMatrix );
  free( p_user_data );

}

/*##############################################################################
// WnMatrix__insert_row_callback().
//############################################################################*/

void
WnMatrix__insert_row_callback(
  WnMatrix__Element *p_element,
  void *p_data,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  typedef struct {
    size_t iRow;
    WnMatrix *pMatrix;
  } user_data;

  user_data *p_user_data = ( user_data * ) p_data;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );

  if( p_element->iRow < p_user_data->iRow )
    WnMatrix__assignElement(
      p_user_data->pMatrix,
      p_element->iRow,
      p_element->iCol,
      p_element->dValue
    );
  else
    WnMatrix__assignElement(
      p_user_data->pMatrix,
      p_element->iRow + 1,
      p_element->iCol,
      p_element->dValue
    );

}

/*##############################################################################
// WnMatrix__insertColumn().
//############################################################################*/

void
WnMatrix__insertColumn(
  WnMatrix *self, size_t i_column_to_insert, gsl_vector *p_vector
) {

  typedef struct {
    size_t iColumn;
    WnMatrix *pMatrix;
  } user_data;

  user_data *p_user_data;
  size_t i;
  double d_val;

  /*============================================================================
  // Check column to insert.
  //==========================================================================*/

  if( i_column_to_insert > WnMatrix__getNumberOfColumns( self ) + 1 )
    WNMATRIX__ERROR( "Column to insert is too large" );

  /*============================================================================
  // Check input vector for appropriate size.
  //==========================================================================*/

  if(
    WnMatrix__get_gsl_vector_size( p_vector ) !=
    WnMatrix__getNumberOfRows( self ) 
  ) {
    WNMATRIX__ERROR( "Incorrect size for gsl_vector" );
  }

  /*============================================================================
  // Copy matrix to new one and shift columns.
  //==========================================================================*/

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data ) WNMATRIX__ERROR( "Couldn't allocate memory" );

  p_user_data->iColumn = i_column_to_insert;
  p_user_data->pMatrix =
    WnMatrix__new(
      WnMatrix__getNumberOfRows( self ),
      WnMatrix__getNumberOfColumns( self ) + 1
    );

  xmlHashScanFull(
    self->pMatrixHash,
    (xmlHashScannerFull) WnMatrix__insert_column_callback,
    p_user_data
  );

  /*============================================================================
  // Clear old matrix hash and set old matrix hash to new.
  //==========================================================================*/

  xmlHashFree( self->pMatrixHash, (xmlHashDeallocator) WnMatrix__freeElement );

  self->pMatrixHash = p_user_data->pMatrix->pMatrixHash;

  self->iCols++;

#ifdef WN_USE_MATRIX_LOOKUP
  free( self->pLookup );
  self->pLookup = p_user_data->pMatrix->pLookup;
#endif

  /*============================================================================
  // Insert new values from vector.
  //==========================================================================*/

  for( i = 1; i <= self->iRows; ++i ) {
    d_val = gsl_vector_get( p_vector, i - 1 );

    if( !WnMatrix__value_is_zero( d_val ) )
    {
      WnMatrix__assignElement(
        self,
        i,
        i_column_to_insert,
        d_val
      );
    }
  }

  /*============================================================================
  // Clean up and done.  Don't do a WnMatrix__free because we don't want to
  // free the matrix hash.
  //==========================================================================*/

  free( p_user_data->pMatrix );
  free( p_user_data );

}

/*##############################################################################
// WnMatrix__insert_column_callback().
//############################################################################*/

void
WnMatrix__insert_column_callback(
  WnMatrix__Element *p_element,
  void *p_data,
  xmlChar *sx_row,
  xmlChar *sx_col,
  xmlChar *sx_tmp
)
{

  typedef struct {
    size_t iColumn;
    WnMatrix *pMatrix;
  } user_data;

  user_data *p_user_data = ( user_data * ) p_data;

  if( (!sx_row && !sx_col) || sx_tmp || !p_element ) 
    WNMATRIX__ERROR( "Invalid input" );


  if( p_element->iCol < p_user_data->iColumn )
    WnMatrix__assignElement(
      p_user_data->pMatrix,
      p_element->iRow,
      p_element->iCol,
      p_element->dValue
    );
  else
    WnMatrix__assignElement(
      p_user_data->pMatrix,
      p_element->iRow,
      p_element->iCol + 1,
      p_element->dValue
    );

}
