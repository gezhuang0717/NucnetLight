//##############################################################################
// Include's.
//##############################################################################

#include <Libstatmech.h>
#include <Libnucnet.h>
#include <Libnuceq.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

#include "nnt/string_defs.h"
#include "nnt/auxiliary.h"
#include "nnt/wrappers.h"

#include "user/thermo.h"

/*##############################################################################
// Defines.
//############################################################################*/

#define MY_BUF_SIZE  32

/*##############################################################################
// Prototypes.
//############################################################################*/

void
print_equilibrium( nnt::Zone& );

void
print_cluster_data( Libnuceq__Cluster *, void * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;
  nnt::Zone wrapped_zone;
  double d_log10_t9;
  int i_count;
  char s_property[MY_BUF_SIZE];

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( argc < 5 )
  {
    fprintf(
      stderr,
      "Usage: argv[0] input_net_xml rho s/b Ye [xpath constraint]\n\n"
    );
    fprintf(
      stderr,
      "  input_net_xml = input net xml file\n\n"
    );
    fprintf(
      stderr,
      "  rho = input mass density (g/cc)\n\n"
    );
    fprintf(
      stderr,
      "  s/b = input entropy per nucleon\n\n"
    );
    fprintf(
      stderr,
      "  out_xml = output xml file\n\n"
    );
    fprintf(
      stderr,
      "  Ye = input Ye (optional)\n\n"
    );
    fprintf(
      stderr,
      "  [xpath constraint] = xpath for equilibrium cluster and abundance constraint (optional)\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get input.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    argv[1],
    NULL
  );

  p_zone =
    Libnucnet__Zone__new( Libnucnet__getNet( p_my_nucnet ), "0", "0", "0" );

  Libnucnet__addZone( p_my_nucnet, p_zone );

  /*============================================================================
  // Input.
  //==========================================================================*/

  Libnucnet__Zone__updateProperty( p_zone, nnt::s_RHO, NULL, NULL, argv[2] );

  Libnucnet__Zone__updateProperty(
    p_zone, 
    nnt::s_ENTROPY_PER_NUCLEON,
    NULL,
    NULL,
    argv[3]
  );

  Libnucnet__Zone__updateProperty(
    p_zone,
    nnt::s_EQUIL_ONLY,
    NULL,
    NULL,
    "yes"
  );
    
  /*============================================================================
  // Optional input.
  //==========================================================================*/

  if( argc > 5 )
    Libnucnet__Zone__updateProperty( p_zone, nnt::s_YE, NULL, NULL, argv[5] );

  if( argc > 6 )
  {
    i_count = 6;
    while( i_count < argc )
    {

      sprintf( s_property, "%d", ( i_count - 6 ) / 2 );
      Libnucnet__Zone__updateProperty(
        p_zone,
        nnt::s_CLUSTER_XPATH,
        s_property,
        NULL,
        argv[i_count++]
      );
      Libnucnet__Zone__updateProperty(
        p_zone,
        nnt::s_CLUSTER_CONSTRAINT,
        s_property,
        NULL,
        argv[i_count++]
      );
    }

    sprintf( s_property, "%d", ( i_count - 6 ) / 2 );
    Libnucnet__Zone__updateProperty(
      p_zone,
      nnt::s_NUMBER_CLUSTERS,
      NULL,
      NULL,
      s_property
    );

  }

  /*============================================================================
  // Set wrapped zone.
  //==========================================================================*/

  wrapped_zone.setNucnetZone( p_zone );

  /*============================================================================
  // Compute log10(t9) from entropy.
  //==========================================================================*/

  d_log10_t9 = 
    wrapped_zone.computeRootFromQuantity(
      (nnt::quantityFunction)
         user::compute_log10_t9_entropy_root_with_equilibrium,
      &wrapped_zone
    );

  /*============================================================================
  // Print out basic data.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nTemperature (K) = %g\n",
    pow( 10., d_log10_t9 ) * GSL_CONST_NUM_GIGA
  );

  fprintf(
    stdout,
    "Density = %s g/cc, Ye = %g\n\n",
    argv[2],
    Libnucnet__Zone__computeZMoment( p_zone, 1 )
  );

  /*============================================================================
  // Print out equilibrium cluster data.
  //==========================================================================*/

  print_equilibrium( wrapped_zone );
     
  /*============================================================================
  // Print out entropy components.
  //==========================================================================*/

  std::cout <<
    "Baryon entropy per nucleon = " <<
    user::compute_thermo_quantity(
      wrapped_zone, "entropy per nucleon", "baryon"
    ) <<
    std::endl;
     
  std::cout << std::endl <<
    "Electron entropy per nucleon = " <<
    user::compute_thermo_quantity(
      wrapped_zone,
      "entropy per nucleon",
      "electron"
    ) <<
    std::endl;
     
  std::cout << std::endl <<
    "Photon entropy per nucleon = " <<
    user::compute_thermo_quantity(
      wrapped_zone,
      "entropy per nucleon",
      "photon"
    ) <<
    std::endl;

  std::cout << std::endl;
     
  /*============================================================================
  // Write output to xml file.
  //==========================================================================*/

  Libnucnet__writeToXmlFile( p_my_nucnet, argv[4] );

  /*============================================================================
  // Clean up.  Done.
  //==========================================================================*/

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}

/*##############################################################################
// print_equilibrium().
//############################################################################*/

void
print_equilibrium( nnt::Zone& wrapped_zone )
{

  Libnuceq *p_equil;
  
  p_equil = nnt::get_zone_equilibrium( wrapped_zone );

  Libnuceq__computeEquilibrium(
    p_equil,
    boost::lexical_cast<double>( wrapped_zone.getProperty( nnt::s_T9 ) ),
    boost::lexical_cast<double>( wrapped_zone.getProperty( nnt::s_RHO ) )
  );

  Libnuceq__iterateClusters(
    p_equil,
    (Libnuceq__Cluster__iterateFunction) print_cluster_data,
    NULL
  );

  Libnuceq__free( p_equil );

}

/*##############################################################################
// print_cluster_data().
//############################################################################*/

void
print_cluster_data(
  Libnuceq__Cluster *p_cluster,
  void *p_data
)
{

  if( p_data )
  {
    fprintf( stderr, "Routine should receive no extra data.\n" );
    exit( EXIT_FAILURE );
  }

  fprintf(
    stdout,
    "Cluster: %s   mu_cluster/kT = %e\n",
    Libnuceq__Cluster__getXPathString( p_cluster ),
    Libnuceq__Cluster__getMukT( p_cluster )
  );

}

