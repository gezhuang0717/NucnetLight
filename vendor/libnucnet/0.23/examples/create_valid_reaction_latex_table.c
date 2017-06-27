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
//       in data from an xml file, write out the valid reactions to a latex
//       file, and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet.h>

void write_reactions( Libnucnet__Net *, const char *, const char * );

int write_reaction( Libnucnet__Reaction *, void * );

int
compare_reactions(
  const Libnucnet__Reaction *,
  const Libnucnet__Reaction *
);

int
main( int argc, char **argv ) {

  Libnucnet__Net *p_my_net;

  if ( argc < 3 || argc > 4 ) {
      fprintf(
        stderr, "\nUsage: %s filename latex_file xpath_nuc xpath_reac\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  latex_file = output latex file.\n\n"
      );
      fprintf(
        stderr, "  xpath_nuc = xpath expression for nuclei (optional--required if xpath_reac present).\n\n"
      );
      fprintf(
        stderr, "  xpath_reac = xpath expression (optional).\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Create the network data.
  //==========================================================================*/

  if( argc == 3 )
    p_my_net = Libnucnet__Net__new_from_xml( argv[1], NULL, NULL );
  else if( argc == 4 )
    p_my_net = Libnucnet__Net__new_from_xml( argv[1], argv[3], NULL );
  else
    p_my_net = Libnucnet__Net__new_from_xml( argv[1], argv[3], argv[4] );
    
  /*============================================================================
  // Print the reactions.
  //==========================================================================*/

  write_reactions( p_my_net, argv[0], argv[argc-1] );

  /*============================================================================
  // Free the network.
  //==========================================================================*/

  Libnucnet__Net__free( p_my_net );
  
  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// write_reactions()
//############################################################################*/

void
write_reactions(
  Libnucnet__Net * self,
  const char * s_code,
  const char * s_latex_file
)
{

  typedef struct
  {
    FILE * pFile;
    Libnucnet__Net * pNet;
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

  p_work->pNet = self;

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
    "This table shows the valid reactions in the network.\n"
  );

  fprintf(
    p_work->pFile,
    "\\begin{center}\n"
  );

  fprintf(
    p_work->pFile,
    "\\begin{tabular}{ll}\n"
  );

  fprintf(
    p_work->pFile,
    "Reaction & Source\\\\\n"
  );

  fprintf(
    p_work->pFile,
    "\\hline\\\\\n"
  );

  /*============================================================================
  // Set the reaction compare function for alphabetically sorted
  // reactions.
  //==========================================================================*/

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( self ),
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  /*============================================================================
  // Print out species data.
  //==========================================================================*/

  p_work->iCount = 1;

  Libnucnet__Reac__iterateReactions(
    Libnucnet__Net__getReac( self ),
    (Libnucnet__Reaction__iterateFunction) write_reaction,
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
// write_reaction()
//############################################################################*/

int
write_reaction(
  Libnucnet__Reaction *p_reaction,
  void * p_data
)
{

  typedef struct
  {
    FILE * pFile;
    Libnucnet__Net * pNet;
    int iCount;
  } work;

  char * s_latex_string;

  work * p_work = ( work * ) p_data;

  s_latex_string =
    Libnucnet__Net__createValidReactionLatexString(
      p_work->pNet,
      p_reaction
    );

  if( !s_latex_string ) return 1;

  fprintf(
    p_work->pFile,
    "$%s$ & %s\\\\\n",
    s_latex_string,
    Libnucnet__Reaction__getSource( p_reaction )
  );

  free( s_latex_string );

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
      "\\begin{tabular}{ll}\n"
    );

    fprintf(
      p_work->pFile,
      "Reaction & Source\\\\\n"
    );

    fprintf(
      p_work->pFile,
      "\\hline\\\\\n"
    );

  }

  return 1;

}

/*##############################################################################
// compare_reactions().
//############################################################################*/

int
compare_reactions(
  const Libnucnet__Reaction *p_reaction1,
  const Libnucnet__Reaction *p_reaction2
)
{

  return
    strcmp(
      Libnucnet__Reaction__getString( p_reaction1 ),
      Libnucnet__Reaction__getString( p_reaction2 )
    );

}
