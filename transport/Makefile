# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2016-2017
#
#                    Copyright D. Brambilla 2016
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Stefano Brambilla <s.brambilla93@gmail.com>
#   DATE        : Settembre 2016
# ====================================================================

.PHONY: all doc clean distclean library

all: library
	$(MAKE) -C src/1_singlebranch_transp 
	$(MAKE) -C src/2_Ybifurcation_transp 
	$(MAKE) -C src/3_network_transp 
	$(MAKE) -C src/4_curved_singlebranch 

library: 
	$(MAKE) -C include

doc:
	install -d doc
	doxygen Doxyfile
	
pdf: doc
	$(MAKE) -C doc/latex pdf

clean:
	$(RM) -r *~ *.log
	$(MAKE) -C include clean
	$(MAKE) -C src/1_singlebranch_transp clean
	$(MAKE) -C src/2_Ybifurcation_transp clean
	$(MAKE) -C src/3_network_transp clean
	$(MAKE) -C src/4_curved_singlebranch clean
	
distclean: clean
	$(RM) -r doc/*
	$(MAKE) -C include distclean
	$(MAKE) -C src/1_singlebranch_transp distclean
	$(MAKE) -C src/2_Ybifurcation_transp distclean
	$(MAKE) -C src/3_network_transp distclean
	$(MAKE) -C src/4_curved_singlebranch distclean
