/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   transport3d1d.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Declaration of the main class for the 3D/1D coupled transport problem.
 */
 
#ifndef M3D1D_TRANSPORT3D1D_HPP_
#define M3D1D_TRANSPORT3D1D_HPP_
 
// GetFem++ libraries
#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>   
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_superlu.h>
#include <getfem/getfem_partial_mesh_fem.h>
#include <getfem/getfem_interpolated_fem.h>
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_iter_solvers.h>
// Standard libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
// Project headers
#include <defines.hpp>
#include <mesh3d.hpp>       
#include <mesh1d.hpp>
#include <utilities.hpp>
#include <assembling1d.hpp>          
#include <assembling3d.hpp>        
#include <assembling3d1d.hpp>
#include <node.hpp>
#include <dof3d1d.hpp>
#include <descr3d1d.hpp>
#include <param3d1d.hpp>
//#include <defines.hpp>>
#include <problem3d1d.hpp>

 
 namespace getfem {

//!	Main class defining the coupled 3D/1D transport problem.
class transport3d1d: public problem3d1d {

public:
 //! Initialize the problem
	/*!
		1. Read the .param filename from standard input
		2. Import problem descriptors (file paths, GetFEM types, ...)
		3. Import mesh for tissue (3D) and vessel network (1D)
		4. Set finite elements and integration methods
		5. Build problem parameters
		6. Build the list of tissue boundary data
		7. Build the list of vessel boundary (and junction) data
	 */
	void init (int argc, char *argv[]);
	//! Assemble the problem
	/*!
		1. Initialize problem matrices and vectors
		2. Build the monolithic matrix AM
		3. Build the monolithic rhs FM
	 */
	void assembly (void);
	//! Solve the problem
	/*!
		Solve the monolithic system AM*UM=FM (direct or iterative)
	 */
	bool solve (void);

	 protected:
	 
	 vector_type CM;
	
}; //end of class trasport3d1d

}  //end of namespace

#endif
