# Mixed Finite Element Methods for Coupled 3D/1D Transport Problems
#### *Politecnico di Milano* (ITALY)

**Author** : Stefano Brambilla 

**Mailto** : <s.brambilla93@gmail.com>

**Date**   : March 2018

#### *Previous projects*

**Author** : Domenico Notaro 

**Mailto** : <domenico.not@gmail.com>

**Date**   : March 2016

**Github Page** : https://github.com/domeniconotaro/PACS

-------------------------------------------------------
## How to install and run the program
-------------------------------------------------------
## THE PACKAGE
- 'fluid'    : The folder containing the code which solves velocity and pressure problem (Notaro code)
- 'transport': The folder containing the code which solves transport problem (Brambilla code)

The details of the two libraries are in their respective README.md


##INSTALLATION:
### Prerequisites

You need the open source finite element library "GetFEM++"
See <http://download.gna.org/getfem/html/homepage>
Version >= 5.1 is necessary


BEWARE: 
Recall to add the GetFEM library path to LD_LIBRARY_PATH. Example:
```
$ export LD_LIBRARY_PATH=/home/...path/to.../getfem/lib

```

Finally, modify the files `fluid/config.mk` and `transport/config.mk` as described in, respectively, 
`fluid/README.md` and `transport/README.md`.

======================

###Installation
Build the whole project with:
```
$ make
```
It first build the shared library "libproblem3d1d.so" and his examples by calling the Makefile in `fluid/`:

``` 
$ make -C fluid/
``` 
Then, it build the shared library "libtransport3d1d.so" and his examples by calling the Makefile in `transport/`:

``` 
$ make -C transport/
``` 
It is also possible to build a single example, by building the two libraries with:
``` 
$ make library
``` 
and then calling in the folder example:
``` 
$ make
``` 

Alternatively, you can build the two libraries separately, calling, in the rispective folders,
``` 
$ make
``` 
Remember to modify the `config.mk` files in the folders, following the instructions of the `README.md` files.

This can be helpful if a new release of "libproblem3d1d" is available.
See: <https://github.com/domeniconotaro/PACS> for the latest release.
======================

### Documentation
The documentation is produced by doxygen. The file Doxyfile contains 
the common doxygen configuration for all examples.
Build the code documentation with:
``` 
$ make pdf
``` 
It first fills fluid/doc/ and transport/doc/ with code documentation ($ make doc) and then compile
the .tex files to produce a portable file ($ pdflatex doc/latex/refman.tex).
You can visualize the documentation with
``` 
$ evince fluid/doc/latex/refman.pdf
$ evince transport/doc/latex/refman.pdf
``` 

##  DEV ENVIRONMENT
OS         : CentOS Linux 7 64-bit 

Compiler   : g++-5.2.1

GetFEM lib : 5.2



