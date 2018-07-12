/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2015-2016
                  
                Copyright (C) 2018 Stefano Brambilla
======================================================================*/
/*! 
  @file   main.cpp  
  @author Stefano Brambilla <s.brambilla93@gmail.com>   
  @date   January 2017. 
  @brief  Main program for test simulations.    
  @details
    We solve the uncoupled 3D/1D problem of concentration trasport in a 1D  
    network \Lambda and in the 3D interstitial tissue \Omega
            
    ***************************************** 
      Benchmark : uncoupled transport problem on single branch  
      Mixed finite elements approximation  
      Monolithic resolution by SuperLU 3.0 
    *****************************************
     
	See Section "Code verification: test-cases" 
 */  
 	 
 	#define M3D1D_VERBOSE_ 
#include <iostream>
#include <problem3d1d.hpp>  
#include <transport3d1d.hpp> 
#include <gmm/gmm.h>

using namespace std;

//! main program
int main(int argc, char *argv[])   
{

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
 
	try {   
		// Declare a new problem 
		getfem::transport3d1d p; 
		double c_start;
		double c_end;  
		/// fluid problem: velocity field and pressure
		
		// Initialize the problem
			c_start =gmm::uclock_sec();
		p.init_fluid(argc, argv);
			c_end =gmm::uclock_sec();
			cout<<"INIT_FLUID TIME: "<<c_end-c_start<<" s;"<<endl;
		// Build the monolithic system 
			c_start =gmm::uclock_sec();
		p.assembly_fluid();
			c_end =gmm::uclock_sec();
			cout<<"ASSEMBLY_FLUID TIME: "<<c_end-c_start<<" s;"<<endl;
		// Solve the problem 
			c_start =gmm::uclock_sec();
		if (!p.solve_fluid()) GMM_ASSERT1(false, "solve procedure has failed");
			c_end =gmm::uclock_sec();
			cout<<"SOLVE_FLUID TIME: "<<c_end-c_start<<" s;"<<endl;
		// Save results in .vtk format
			c_start =gmm::uclock_sec();
		p.export_vtk_fluid();
			c_end =gmm::uclock_sec();
			cout<<"EXPORT_FLUID TIME: "<<c_end-c_start<<" s;"<<endl;
		// Display some global results: mean pressures, total flow rate
		
		std::cout << "--- FINAL RESULTS -------------------------" << std::endl; 
		std::cout << "  Pt average            = " << p.mean_pt()   << std::endl;
		std::cout << "  Pv average            = " << p.mean_pv()   << std::endl;
		std::cout << "  Network-to-Tissue TFR = " << p.flow_rate() << std::endl;
		std::cout << "-------------------------------------------" << std::endl; 	
		
		
		     
		//transport problem: concentration  
		
		//initialize 
			c_start =gmm::uclock_sec();
		p.init_transp(argc, argv);
			c_end =gmm::uclock_sec();
			cout<<"INIT_TRANSP TIME: "<<c_end-c_start<<" s;"<<endl;
		//assemble        
			c_start =gmm::uclock_sec();
		p.assembly_transp(); 
			c_end =gmm::uclock_sec();
			cout<<"ASSEMBLY_TRANSP TIME: "<<c_end-c_start<<" s;"<<endl;   
		//solve     
			c_start =gmm::uclock_sec();
		if (!p.solve_transp()) GMM_ASSERT1(false, "solve procedure has failed");// the export is in the solve at each time step 
			c_end =gmm::uclock_sec();
			cout<<"SOLVE_TRANSP+EXPORT_TRANSP TIME: "<<c_end-c_start<<" s;"<<endl;  
		//check the mass balance in the junctions
			c_start =gmm::uclock_sec();
		p.mass_balance();
			c_end =gmm::uclock_sec();
			cout<<"MASS BALANCE TIME: "<<c_end-c_start<<" s;"<<endl;  	      
		   
		}  
      
	GMM_STANDARD_CATCH_ERROR;     
		 
	   
    
		    
		   
	return 0;    
	   
} /* end of main program */   
    
   
  

