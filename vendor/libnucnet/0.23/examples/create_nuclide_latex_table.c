/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Bradley S. Meyer.
//
//     This is free software; you can redistribute it and/or modify it
//     under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     Please see the src/README.txt file in this distribution for more
//     information.
//   </license>
//   <description>
//     <abstract>
//       Example to demonstrate how to use libnucnet routines to parse
//       in a Libnucnet__Nuc structure of nuclear species from an input
//       xml file, write species data out to a latex file,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet__Nuc.h>

void write_nuclei( Libnucnet__Nuc *, const char *, const char * );

int write_species( Libnucnet__Species *, void * );

int
main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;

  if ( argc != 3 && argc != 4 ) {
      fprintf(
        stderr, "\nUsage: %s filename xpath latex_file\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  xpath = xpath expression (optional).\n\n"
      );
      fprintf(
        stderr, "  latex_file = output latex file.\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Create nuclear species collection.
  //==========================================================================*/

  if( argc == 3 ) {
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );
  } else {
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], argv[2] );
  }
    
  /*============================================================================
  // Print the nuclear data.
  //==========================================================================*/

  if( argc == 3 )
    write_nuclei( p_my_nuclei, argv[0], argv[2] );
  else
    write_nuclei( p_my_nuclei, argv[0], argv[3] );

  /*============================================================================
  // Free the nuclear species collection.
  //==========================================================================*/

  Libnucnet__Nuc__free( p_my_nuclei );
  
  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// write_nuclei()
//############################################################################*/

void
write_nuclei(
  Libnucnet__Nuc * self,
  const char * s_code,
  const char * s_latex_file
)
{

  typedef struct
  {
    FILE * pFile;
    int iCount;
  } work;

  work * p_work;

  /*============================================================================
  // Work structure.
  //==========================================================================*/

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work )
  {
    fprintf( stderr, "Couldn't allocate work structure.\n" );
    exit( EXIT_FAILURE );
  }

  /*============================================================================
  // Open file.
  //==========================================================================*/

  p_work->pFile = fopen( s_latex_file, "w" );

  if( !p_work->pFile )
  {
    fprintf( stderr, "Couldn't open file.\n" );
    exit( EXIT_FAILURE );
  }

  /*============================================================================
  // Write out header.
  //==========================================================================*/

  fprintf(
    p_work->pFile,
    "%% File created by %s.  To create a pdf, run\n%%\n",
    s_code
  );

  fprintf(
    p_work->pFile,
    "%% pdflatex %s\n%%\n",
    s_latex_file
  );

  fprintf(
    p_work->pFile,
    "\\documentclass{article}[12pt]\n"
  );

  fprintf(
    p_work->pFile,
    "\\begin{document}\n"
  );

  fprintf(
    p_work->pFile,
    "This table shows the species in the nuclide collection.\n"
  );

  fprintf(
    p_work->pFile,
    "\\begin{center}\n"
  );

  fprintf(
    p_work->pFile,
    "\\begin{tabular}{ccc}\n"
  );

  fprintf(
    p_work->pFile,
    "Species & Z & A\\\\\n"
  );

  fprintf(
    p_work->pFile,
    "\\hline\\\\\n"
  );

  /*============================================================================
  // Print out species data.
  //==========================================================================*/

  p_work->iCount = 1;

  Libnucnet__Nuc__iterateSpecies(
    self,
    (Libnucnet__Species__iterateFunction) write_species,
    p_work
  );

  /*============================================================================
  // Close latex and file.
  //==========================================================================*/

  fprintf(
    p_work->pFile,
    "\\hline\n"
  );

  fprintf(
    p_work->pFile,
    "\\end{tabular}\n"
  );

  fprintf(
    p_work->pFile,
    "\\end{center}\n"
  );

  fprintf(
    p_work->pFile,
    "\\end{document}\n"
  );

  fclose( p_work->pFile );
  free( p_work );

  return;

}

/*##############################################################################
// write_species()
//############################################################################*/

int
write_species(
  Libnucnet__Species *p_species,
  void * p_data
)
{

  typedef struct
  {
    FILE * pFile;
    int iCount;
  } work;

  char * s_latex_name;

  work * p_work = ( work * ) p_data;

  s_latex_name = Libnucnet__Species__createLatexString( p_species );

  fprintf(
    p_work->pFile,
    "$%s$ & %u & %u\\\\\n",
    s_latex_name,
    Libnucnet__Species__getZ( p_species ),
    Libnucnet__Species__getA( p_species )
  );

  free( s_latex_name );

  if( p_work->iCount++ % 40 == 0 )
  {

    fprintf(
      p_work->pFile,
      "\\hline\n"
    );

    fprintf(
      p_work->pFile,
      "\\end{tabular}\n"
    );

    fprintf(
      p_work->pFile,
      "\\clearpage\n"
    );

    fprintf(
      p_work->pFile,
      "\\begin{tabular}{ccc}\n"
    );

    fprintf(
      p_work->pFile,
      "Species & Z & A\\\\\n"
    );

    fprintf(
      p_work->pFile,
      "\\hline\\\\\n"
    );

  }

  return 1;

}

