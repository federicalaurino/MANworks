# Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems
#### *Politecnico di Milano* (ITALY)

**Author** : Stefano Brambilla 

**Mailto** : <s.brambilla93@gmail.com>

**Date**   : March 2017


-------------------------------------------------------
## How to install and run the program
-------------------------------------------------------
## THE PACKAGE
- 'fluid'    : The folder containing the code which solves velocity and pressure problem (Notaro code)
- 'transport': The folder containing the code which solves transport problem (Brambilla code)



INSTALLAZIONE:

1) modificare i file config.mk sia in "fluid" che in "transport"

GETFEM_PREFIX=/path/to/Getfem/
PROBLEM_FLUID_PREFIX= /path/to/fluid/code/


Per esempio, in un computer MOX, con i moduli, si puo usare:
GETFEM_PREFIX= $(mkGetfemInc)/../
in tal caso poi da terminale dovrete scrivere:
$ source module.sh


Il codice di notaro può essere sia quello contenuto nella cartella fluid, o un'altra versione più aggiornata del codice, che avete già compilato sul vostro computer.

2) se volete installare tutto, dalla cartella MANworks scrivete su terminale:
$ make
per compilare tutto (anche i casi test)

in alternativa potete andare nelle singole cartelle "fluid" e "transport" e compilare le due librerie separatamente (consigliabile se avete già compilato altrove il codice di notaro)







