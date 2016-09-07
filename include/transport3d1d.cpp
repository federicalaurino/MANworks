/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   transport3d1d.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the main class for the 3D/1D coupled transport problem.
 */
 
 #include <transport3d1d.hpp>
 
 namespace getfem {
 void transport3d1d::init (int argc, char *argv[])
 {
 std::cout << "initialize transport problem"<<std::endl;
 
 }; // end of init
 
 void transport3d1d::assembly (void)
 {
 std::cout<<"assemble trasnport problem"<<std::endl;
 
 }; // end of assembly
 
 
 } // end of namespace
