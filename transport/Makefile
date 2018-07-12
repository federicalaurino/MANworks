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
	$(MAKE) -C src/1_uncoupled_branch
	$(MAKE) -C src/2_coupled_branch
	$(MAKE) -C src/3_bifurcation 
	$(MAKE) -C src/4_network

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
	$(MAKE) -C src/1_uncoupled_branch clean
	$(MAKE) -C src/2_coupled_branch clean
	$(MAKE) -C src/3_bifurcation  clean
	$(MAKE) -C src/4_network clean

distclean: clean
	$(RM) -r doc/*
	$(MAKE) -C include distclean
	$(MAKE) -C src/1_uncoupled_branch distclean
	$(MAKE) -C src/2_coupled_branch distclean
	$(MAKE) -C src/3_bifurcation  distclean
	$(MAKE) -C src/4_network distclean
