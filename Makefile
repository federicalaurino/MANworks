# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2016-2017
#
#                    Copyright D. Brambilla2017
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Stefano Brambilla <s.brambilla93@gmail.com>
#   DATE        : March 2017
# ====================================================================

.PHONY: all doc clean distclean library

all: 
	$(MAKE) -C fluid
	$(MAKE) -C transport
	
library:
	$(MAKE) -C fluid library
	$(MAKE) -C transport library

fluid: 
	$(MAKE) -C fluid
	
transport: fluid
	$(MAKE) -C transport

doc:
	$(MAKE) -C fluid doc
	$(MAKE) -C transport doc
	
	
pdf: doc
	$(MAKE) pdf -C fluid
	$(MAKE) pdf -C transport	

clean:
	$(RM) -r *~ *.log
	$(MAKE) -C fluid clean
	$(MAKE) -C transport clean

distclean: clean
	$(RM) -r doc/*
	$(MAKE) -C fluid distclean
	$(MAKE) -C transport distclean
