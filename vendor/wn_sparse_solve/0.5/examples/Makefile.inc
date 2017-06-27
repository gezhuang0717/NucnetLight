#///////////////////////////////////////////////////////////////////////////////
# <file type="public">
#
#   <description>
#     <abstract>
#       Makefile for dealing with Sparskit2.  It is included from Makefile.
#     </abstract>
#   </description>
#   <license>
#     See the src/README.txt file in this distribution for copyright and license
#     information.
#   </license>
#
# </file>
#///////////////////////////////////////////////////////////////////////////////

SPARSEKIT_URL = http://www-users.cs.umn.edu/~saad/software/SPARSKIT/
SPARSEKIT_DATA = SPARSKIT2.tar

#===============================================================================
# Get sparskit
#===============================================================================

sparskit:
	mkdir -p $(SPARSKITDIR)
	wget $(SPARSEKIT_URL)$(SPARSEKIT_DATA).gz
	mv ${SPARSEKIT_DATA}.gz $(SPARSKITDIR)..
	tar xz -C ${SPARSKITDIR}.. -f $(SPARSKITDIR)../${SPARSEKIT_DATA}.gz
	rm ${SPARSKITDIR}../${SPARSEKIT_DATA}.gz

#===============================================================================
# Sparskit objects
#===============================================================================


SPOBJ = $(SPARSKITDIR)BLASSM/blassm.o         \
        $(SPARSKITDIR)BLASSM/matvec.o         \
        $(SPARSKITDIR)FORMATS/formats.o       \
        $(SPARSKITDIR)FORMATS/unary.o         \
        $(SPARSKITDIR)INFO/infofun.o          \
        $(SPARSKITDIR)INOUT/inout.o           \
        $(SPARSKITDIR)ITSOL/ilut.o            \
        $(SPARSKITDIR)ITSOL/iters.o           \
        $(SPARSKITDIR)MATGEN/FDIF/genmat.o    \
        $(SPARSKITDIR)MATGEN/FEM/elmtlib2.o   \
        $(SPARSKITDIR)MATGEN/FEM/femgen.o     \
        $(SPARSKITDIR)MATGEN/FEM/meshes.o     \
        $(SPARSKITDIR)MATGEN/MISC/sobel.o     \
        $(SPARSKITDIR)MATGEN/MISC/zlatev.o    \
        $(SPARSKITDIR)ORDERINGS/ccn.o         \
        $(SPARSKITDIR)ORDERINGS/color.o       \
        $(SPARSKITDIR)ORDERINGS/dsepart.o

lskit libskit.a: $(SPOBJ)
	$(AR) $(SPARSKITDIR)libskit.a $(SPOBJ)

.PHONY: clean_lib

clean_lskit:
	rm -f $(SPOBJ)
	rm -f $(SPARSKITDIR)libskit.a
