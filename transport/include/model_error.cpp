/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2018 Stefano Brambilla
======================================================================*/
/*! 
  @file   model_error.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   July 2018
  @brief  Methods of the main class for computing model error.
 */
 
 #include <transport3d1d.hpp>
 #include <AMG_Interface.hpp>
 #include <cmath>
 #include "gmm/gmm_inoutput.h"
 #include "getfem/getfem_import.h"


 namespace getfem {

const int A1=1;
const int A2=2;
const int A3=3;

const int PRIMAL=1;
const int DUAL=2;

	void transport3d1d::model_error(int argc, char *argv[]){// A1: U=U(s)
	
		//initialize 
		init_model(argc, argv);

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A1: solve problem for first assumption:"<<endl;
	cout <<"    U(s,r,theta)=U(s)"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A1);    
		//solve     
		if (!solve_model(A1,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A1);    
		//solve     
		if (!solve_model(A1,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D1 AND L1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_error_1();

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A2: solve problem for second assumption:"<<endl;
	cout <<"    Omega+ == Omega"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A2);    
		//solve     
		if (!solve_model(A2,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A2);    
		//solve     
		if (!solve_model(A2,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D2 AND L2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_error_2();

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A3: solve problem for third assumption:"<<endl;
	cout <<"    Neglect fluctuations"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A3);    
		//solve     
		if (!solve_model(A3,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_dual_model(A3);    
		//solve     
		if (!solve_model(A3,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D3 AND L3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_error_3();

};

	//! Initialize problem U1
	void transport3d1d::init_model(int argc, char *argv[]){

	PARAM.read_command_line(argc, argv);
	//1. Import data (algorithm specifications, boundary conditions, ...)	
	import_data_transp();
	//2. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh_transp();
	//3. Set finite elements and integration methods
	set_im_and_fem_model();
 	//4. Build problem parameters
	build_param_transp();
	//5. Build the list of tissue boundary data
	build_tissue_boundary_transp();
	//6. Build the list of tissue boundary (and junction) data
 	build_vessel_boundary_transp();
};
	//3. Set finite elements and integration methods: mf_t should be defined only on Omega_plus
	void transport3d1d::set_im_and_fem_model(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	
	
	pfem pf_Ct = fem_descriptor(descr_transp.FEM_TYPET_C);
	pfem pf_Cv = fem_descriptor(descr_transp.FEM_TYPEV_C); 

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for tissue ..." << endl;
	#endif
		


	mf_Ct.set_finite_element(mesht.convex_index(), pf_Ct);
	mf_Ct_Omega.set_finite_element(mesht.region(descr_transp.OMEGA).index(), pf_Ct); //.region(descr_transp.OMEGA)
	mf_Ct_Sigma.set_finite_element(mesht.region(descr_transp.SIGMA).index(), pf_Ct); //.region(descr_transp.SIGMA)

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif

	mf_Cv.set_finite_element(meshv.convex_index(), pf_Cv);

	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
	#endif

	dof_transp.set(mf_Ct, mf_Cv, mf_Ct_Omega, mf_Ct_Sigma);
	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dof_transp;
	#endif
	
	mimv.clear();
	mf_Uvi.clear();
	mf_Pv.clear();
	mf_coefv.clear();
	mf_coefvi.clear();

	mimt.clear();
	mf_Ut.clear();
	mf_Pt.clear();
	mf_coeft.clear();
	problem3d1d::set_im_and_fem();


};
	//! Assemble problem U1
	void transport3d1d::assembly_model(const size_type ASSUMPTION){


	// define dofs
	size_type dof_t;
	if(ASSUMPTION==A1){
	dof_t  = dof_transp.Ct_Omega();	}
	else{
	dof_t  = dof_transp.Ct();		}
	size_type dof_v  = dof_transp.Cv();
	size_type dof_tot= dof_t+dof_v;
	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM_transp, UM_transp, FM_transp ..." << endl;
	#endif
	gmm::resize(AM_transp, dof_tot, dof_tot);	gmm::clear(AM_transp);
	gmm::resize(UM_transp, dof_tot);		gmm::clear(UM_transp);
	gmm::resize(FM_transp, dof_tot);	 	gmm::clear(FM_transp);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM_transp ..." << endl;
	#endif
	// Mass(time derivative)  matrix for the interstitial problem
	sparse_matrix_type At(dof_t, dof_t);gmm::clear(At);
	// Mass (time derivative)  matrix for the network problem
	sparse_matrix_type Av(dof_v, dof_v);gmm::clear(Av);	

	// Tissue-to-tissue exchange matrix on Gamma
	sparse_matrix_type Mtt(dof_t, dof_t);gmm::clear(Mtt);
	// Tissue-to-vessel exchange matrix on Gamma
	sparse_matrix_type Mtv(dof_t, dof_v);gmm::clear(Mtv);
	// Vessel-to-tissue exchange matrix on Gamma
	sparse_matrix_type Mvt(dof_v, dof_t);gmm::clear(Mtv);
	// Vessel-to-vessel exchange matrix on Gamma
	sparse_matrix_type Mvv(dof_v, dof_v);gmm::clear(Mvv);

	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof_t, dof_t);gmm::clear(Btt);
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof_t, dof_v);gmm::clear(Btv);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof_v, dof_t);gmm::clear(Bvt);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof_v, dof_v);gmm::clear(Bvv);

	// Aux tissue-to-vessel averaging matrix (Circumference)
	sparse_matrix_type Mbar(dof_v, dof_t);gmm::clear(Mbar);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof_v, dof_t);gmm::clear(Mlin);
	// Aux tissue-to-vessel averaging matrix (Section)
	sparse_matrix_type Mbarbar(dof_v, dof_transp.Ct());gmm::clear(Mbarbar); //In any case, this matrix should be defined on the whole 3d mesh
	
	// Aux tissue source vector
	vector_type Ft(dof_t); gmm::clear(Ft);
	// Aux vessels source vector
	vector_type Fv(dof_v); gmm::clear(Fv);


	//Assemble Perimeter, Area, and build Kt and Kv

	//Permeability
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	vector_type Kv(dof.coefv(),k);
	vector_type Kt(dof.coeft(),k);

	//Perimeter = 2*pi*R*k 
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);

	//Area = pi*R^2
	vector_type Area(dof.coefv());
	gmm::copy(param.R(), Area);
	gmm::vscale(param.R(), Area);
	gmm::scale(Area, pi); 

 
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the stiffness matrix At ..." << endl;
	#endif	
	// Assemble At, stiffness matrix for laplacian in tissue
	if(ASSUMPTION==A1){
	getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct_Omega); 	}
	else{
	getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct      );    }
	
	gmm::add(At,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 


	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the stiffness matrix Av ..." << endl;
	#endif	
	// Assemble Av, stiffness matrix for laplacian in vessel
	getfem::asm_stiffness_matrix_for_laplacian(Av,mimv,mf_Cv,mf_coefv,Area); 
	
	gmm::add(Av,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_t, dof_v), 
					gmm::sub_interval(dof_t, dof_v))); 


	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
	if(ASSUMPTION==A1){
	asm_exchange_aux_mat_transp(Mbar, Mlin, 
			mimv, mf_Ct_Omega, mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);}
	else{
	asm_exchange_aux_mat_transp(Mbar, Mlin, 
			mimv, mf_Ct,       mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);}

	bool NEWFORM = 1;//PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");	
	asm_exchange_mat(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_coefv, Mbar, Mlin, Perimeter, NEWFORM);

	// Copying Bvt
	gmm::add(scaled(Bvt,-1),								
			  gmm::sub_matrix(AM_transp, 
			  		gmm::sub_interval(dof_t, dof_v),
					gmm::sub_interval(0, dof_t)));
	// Copying Bvv
	gmm::add(Bvv,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_t, dof_v), 
					gmm::sub_interval(dof_t, dof_v))); 
	

	if(ASSUMPTION==A1 || ASSUMPTION==A2){
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling mass matrices on Gamma..." << endl;
	#endif
	// Assemble Mtt, mass matrix on tissue
	if(ASSUMPTION==A1){
	getfem::asm_mass_matrix_param (Mtt, mimt, mf_Ct_Omega, mf_coeft, Kt, descr_transp.GAMMA);}
	else{
	getfem::asm_mass_matrix_param (Mtt, mimt, mf_Ct,       mf_coeft, Kt, descr_transp.GAMMA);}

	// Assemble Mtv, exchange matrix on gamma
	gmm::vscale(Area, Kv);
	getfem::asm_mass_matrix_param (Mvv, mimv, mf_Cv,mf_coefv, Kv);
	gmm::mult(gmm::transposed(Mbar), Mvv, Mtv); // M1 * M2 ---> M3	
	
	// Copying Mtt
	gmm::add(Mtt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 
	// Copying Mtv

	gmm::add(scaled(Mtv,-1),
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t),
					gmm::sub_interval(dof_t, dof_v))); 

	} else { //ASSUMPTION==A3)
	// Copying Btt
	gmm::add(Btt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 

	gmm::resize(BTT, dof_transp.Ct(), dof_transp.Ct());	gmm::clear(BTT);
	gmm::copy(Btt,BTT);
	// Copying Btv

	gmm::add(scaled(Btv,-1),
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()),
					gmm::sub_interval(dof_transp.Ct(), dof_v))); 

	}

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling source terms ..." << endl;
	#endif

	// Assemble F: source term in tissue
	bool F_CONSTANT = PARAM.int_value("F_CONSTANT", "flag for using constant source term f in tissue");
	if(F_CONSTANT){
	scalar_type f = PARAM.real_value("F", "Value of source for tissue reduced problem");	
	vector_type F(dof.coeft(),f);

	if(ASSUMPTION==A1){
	getfem::asm_source_term(Ft, mimt, mf_Ct_Omega, mf_coeft, F);}
	else{
	getfem::asm_source_term(Ft, mimt, mf_Ct,       mf_coeft, F);}

	}else{ //F is not constant --> use muparser!


	}
	
	
	// Assemble G: source term in vessel
	//! /todo Extend to variable source term: use muParser and build cross section average (operator Mbarbar)
	bool G_CONSTANT = PARAM.int_value("G_CONSTANT", "flag for using constant source term g in vessel");
	if(G_CONSTANT){
	scalar_type g = PARAM.real_value("G", "Value of source for vessel reduced problem");
	vector_type G(dof.coefv(),g);	
	gmm::vscale(Area, G); //G=G*pi*R^2
	getfem::asm_source_term(Fv, mimv, mf_Cv, mf_coefv, G);
	}
	else{ //G is not constant --> use muparser!
	scalar_type g = PARAM.real_value("G", "Value of source for vessel reduced problem");
	vector_type G(dof.coeft(),g);	
	//vector_type Areat(dof.coeft());
	//getfem::interpolation(mf_coefv, mf_coeft, Area, Areat);
	//gmm::vscale(Areat, G); //G=G*pi*R^2

	// Aux tissue source vector
	vector_type Ftt(dof_transp.Ct()); gmm::clear(Ftt);
	getfem::asm_source_term(Ftt, mimt, mf_Ct, mf_coeft, G,descr_transp.SIGMA);

	//Build average matrix Mbarbar
	asm_exchange_aux_mat_bar_bar(Mbarbar,mimt, mf_Ct, mf_Cv,mf_coefv, param.R(), descr.NInt, descr_transp.NIntA, nb_branches);

	gmm::mult(Mbarbar, Ftt, Fv);
	}

	gmm::add(Ft, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(0,dof_t)));

	gmm::add(Fv, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(dof_t, dof_v)));


	// De-allocate memory
	gmm::clear(At);	   gmm::clear(Av); 
	gmm::clear(Mtt);   gmm::clear(Mtv); 
	gmm::clear(Mvt);   gmm::clear(Mvv); 	    
	gmm::clear(Mbar);  gmm::clear(Mlin);   gmm::clear(Mbarbar);  
	gmm::clear(Btt);   gmm::clear(Btv);
	gmm::clear(Bvt);   gmm::clear(Bvv);

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling boundary conditions ..." << endl;
	#endif	
	gmm::resize(AM_temp, dof_tot,
			     dof_tot);	gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_tot); 	gmm::clear(FM_temp);
	gmm::copy(AM_transp,AM_temp);
	gmm::copy(FM_transp,FM_temp);

	//Boundary conditions: 
	//On vessel, we have homogeneous neumann conditions: we simply don't implement the boundary term due to diffusion.

	//Homogeneous Dirichlet on all the faces of tissue	
	vector_type Dirichlet_null(dof_t,0);
	for (size_type bc=0; bc < BCt_transp.size(); ++bc) {

		if(ASSUMPTION==A1){
		getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct_Omega, BCt_transp[bc].rg, Dirichlet_null);	}
		else{
		getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct,       BCt_transp[bc].rg, Dirichlet_null);	}
	} 
	

};

	//! Assemble problem Z1
	void transport3d1d::assembly_dual_model(const size_type ASSUMPTION){
	size_type dof_t;
	if(ASSUMPTION==A1){
	dof_t  = dof_transp.Ct_Omega();	}
	else{
	dof_t  = dof_transp.Ct();	}
	size_type dof_v  = dof_transp.Cv();
	size_type dof_tot= dof_t+dof_v;

	gmm::clear(UM_transp);

	// Aux tissue source vector
	vector_type Jt(dof_t); gmm::clear(Jt);
	// Aux vessels source vector
	vector_type Jv(dof_v); gmm::clear(Jv);


	// prepare AM and Fm temp
	gmm::resize(AM_temp, dof_tot,
			     dof_tot);	gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_tot); 	gmm::clear(FM_temp);
	gmm::copy(AM_transp,AM_temp);
	gmm::copy(FM_transp,FM_temp);

	gmm::scale(gmm::sub_matrix(AM_temp, 
					gmm::sub_interval(0, dof_t),
					gmm::sub_interval(dof_t, dof_v))
		 ,0);
	gmm::scale(gmm::sub_matrix(AM_temp, 
			  		gmm::sub_interval(dof_t, dof_v),
					gmm::sub_interval(0, dof_t))
		 ,0);

		#ifdef M3D1D_VERBOSE_
	cout << "  Assembling source terms ..." << endl;
	#endif
	string functional= PARAM.string_value("FUNCTIONAL", "descriptor for the functional J of the error you wanto to estimate");

	if(functional=="MEAN_VALUE"){
	// Aux tissue source vector
	vector_type Ft(dof.coeft(),1);
	// Aux vessels source vector
	vector_type Fv(dof.coefv(),1);

	// Assemble Jt: source term in tissue
	if(ASSUMPTION==A1){
	getfem::asm_source_term(Jt,mimt, mf_Ct_Omega, mf_coeft, Ft);}
	else{
	getfem::asm_source_term(Jt,mimt, mf_Ct,       mf_coeft, Ft);}

	// Assemble Jv: source term in vessel
	getfem::asm_source_term(Jv,mimv, mf_Cv, mf_coefv, Fv ); 

	//Add to FM_temp

	gmm::add(Jt, 
			gmm::sub_vector(FM_temp,
					gmm::sub_interval(0,dof_t)));

	gmm::add(Jv, 
			gmm::sub_vector(FM_temp,
					gmm::sub_interval(dof_t, dof_v)));
	}
	else if (functional=="MEAN_STRESS"){};

	//Homogeneous Dirichlet on all the faces of tissue
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling boundary conditions ..." << endl;
	#endif	

	vector_type Dirichlet_null(dof_t,0);
	for (size_type bc=0; bc < BCt_transp.size(); ++bc) {
	if(ASSUMPTION==A1){
	getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct_Omega, BCt_transp[bc].rg, Dirichlet_null);}
	else{
	getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct,       BCt_transp[bc].rg, Dirichlet_null);}
			
	} 

/*
	if(functional=="L2NORM")

	// Assemble F: source term in tissue
	sparse_matrix_type JJt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(JJt);
		
	generic_assembly L2normt;
	L2normt.push_mi(mimt);
	L2normt.push_mf(mf_Ct);
	L2normt.push_mat(JJt);
	L2normt.set("V(#1,#1)+=comp(Base(#1).Base(#1))(:,:)");
	L2normt.assembly();
	
	//gmm::add(gmm::mat_trace(JJt),Jt);
	
	for (int i=0; i>dof_transp.Ct(); i++){
		Jt[i]=JJt(i,i);
		Jt[i]=sqrt(Jt[i]);
	};

	// Assemble G: source term in vessel
	//! /todo Extend to variable source term: use muParser and build cross section average (operator Mbarbar)

	sparse_matrix_type JJv(dof_v, dof_v);gmm::clear(JJv);
		
	generic_assembly L2normv;
	L2normv.push_mi(mimv);
	L2normv.push_mf(mf_Cv);
	L2normv.push_mat(JJv);
	L2normv.set("V(#1,#1)+=comp(Base(#1).Base(#1))(:,:)");
	L2normv.assembly();
	
	//gmm::add(gmm::mat_trace(JJv),Jv);
	for (int i=0; i>dof_v; i++){
		Jv[i]=JJv(i,i);
		Jv[i]=sqrt(Jv[i]);
		
	};
*/	


};


	//! Solve problem U1
	bool transport3d1d::solve_model(const size_type ASSUMPTION,const  size_type VERSION){

  	#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	double time = gmm::uclock_sec();
	gmm::clear(UM_transp);
	gmm::clean(AM_temp, 1E-12);
	gmm::clean(FM_temp, 1E-12);
	// Solve the system on AM_temp, UM_transp, FM_temp
	bool solved = solver_transp();
	if (!solved) return false;

	//export solution
	#ifdef M3D1D_VERBOSE_
	std::cout<<"solved! going to export..."<<std::endl;
	#endif		

	size_type dof_t;
	if(ASSUMPTION==A1){
	dof_t = dof_transp.Ct_Omega();	}
	else{
	dof_t  = dof_transp.Ct();	}
	size_type dof_v  = dof_transp.Cv();
	size_type dof_tot= dof_t+dof_v;

	switch(VERSION){

		case PRIMAL:{
		switch(ASSUMPTION){

			case A1: 
			{	gmm::resize(U1, dof_tot);	gmm::clear(U1);
				gmm::copy(UM_transp, U1);
				cout << endl<<"... time to solve U1: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A2: 
			{	gmm::resize(U2, dof_tot);	gmm::clear(U2);
				gmm::copy(UM_transp, U2);	
				cout << endl<<"... time to solve U2: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A3: 
			{	gmm::resize(U3, dof_tot);	gmm::clear(U3);
				gmm::copy(UM_transp, U3);	
				cout << endl<<"... time to solve U3: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
		 }
				break;} 

		case DUAL:{
		switch(ASSUMPTION){

			case A1: 
			{	gmm::resize(Z1, dof_tot);	gmm::clear(Z1);
				gmm::copy(UM_transp, Z1);	
				cout << endl<<"... time to solve Z1: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A2: 
			{	gmm::resize(Z2, dof_tot);	gmm::clear(Z2);
				gmm::copy(UM_transp, Z2);	
				cout << endl<<"... time to solve Z2: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A3: 
			{	gmm::resize(Z3, dof_tot);	gmm::clear(Z3);
				gmm::copy(UM_transp, Z3);	
				cout << endl<<"... time to solve Z3: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
		 }
				break;} 

		}
	return true;

	};
	//! Export problem U1
	void transport3d1d::export_model(const size_type ASSUMPTION,const  size_type VERSION, const string & suff ){

	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif
	size_type dof_t;
	if(ASSUMPTION==A1){
	dof_t = dof_transp.Ct_Omega();	}
	else{
	dof_t  = dof_transp.Ct();	}
	size_type dof_v  = dof_transp.Cv();
	size_type dof_tot= dof_t+dof_v;
	
	// Array of unknown dof of the interstitial velocity
	vector_type Ct(dof_t); 

	// Array of unknown dof of the network velocity
	vector_type Cv(dof_v); 

	//Copy solution


	switch(VERSION){

		case PRIMAL:{
		switch(ASSUMPTION){

			case A1: 
			{	gmm::copy(gmm::sub_vector(U1, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U1, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A2: 
			{	gmm::copy(gmm::sub_vector(U2, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U2, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A3: 
			{	gmm::copy(gmm::sub_vector(U3, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U3, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
		}
				break;} 

		case DUAL:{
		switch(ASSUMPTION){

			case A1: 
			{	gmm::copy(gmm::sub_vector(Z1, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z1, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A2: 
			{	gmm::copy(gmm::sub_vector(Z2, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z2, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A3: 
			{	gmm::copy(gmm::sub_vector(Z3, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z3, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
		}
				break;} 

	}

	string assumption, version;	
	switch(VERSION){
		case PRIMAL:
			version="U";
			break;
		case DUAL:
			version="Z";
			break;
	}
	
	switch(ASSUMPTION){
		case A1:
			assumption="1";
			break;
		case A2:
			assumption="2";
			break;
		case A3:
			assumption="3";
			break;
	}
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting reduced solution in tissue: "<<version<<assumption << endl;
	#endif
	if(ASSUMPTION==A1){
	vtk_export exp_Ct(descr_transp.OUTPUT+version+"t"+assumption+suff+".vtk");
	exp_Ct.exporting(mf_Ct_Omega);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct_Omega, Ct, version+"t"+assumption);
	}else{
	vtk_export exp_Ct(descr_transp.OUTPUT+version+"t"+assumption+suff+".vtk");
	exp_Ct.exporting(mf_Ct);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct, Ct, version+"t"+assumption);
	}


	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting reduced solution in vessel: "<<version<<assumption << endl;
	#endif
	vtk_export exp_Cv(descr_transp.OUTPUT+version+"v"+assumption+suff+".vtk");
	exp_Cv.exporting(mf_Cv);
	exp_Cv.write_mesh();
	exp_Cv.write_point_data(mf_Cv, Cv, version+"v"+assumption);

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif
};


/*	//! Compute model error of A1	
	void transport3d1d::compute_error_1(void){
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A1 ..." << endl;
	#endif
	// Storing error model from bilinear form 
	vector_type d1(dof_transp.Cv()); gmm::clear(d1);
	cout<<"WARNING: d1(u,z)=0!"<<endl;

	// Storing error model from linear form 
	vector_type l1(dof_transp.Cv()); gmm::clear(l1);
	cout<<"WARNING: l1(u,z)=0 if g=const!!"<<endl;

// do nothing: estimators are zero!


	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting d1 and l1 ..." << endl;
	#endif
	vtk_export exp_d(descr_transp.OUTPUT+"d1.vtk");
	exp_d.exporting(mf_Cv);
	exp_d.write_mesh();
	exp_d.write_point_data(mf_Cv, d1, "d1");

	vtk_export exp_l(descr_transp.OUTPUT+"l1.vtk");
	exp_l.exporting(mf_Cv);
	exp_l.write_mesh();
	exp_l.write_point_data(mf_Cv, l1, "l1");
	
};
*/

 
 } // end of namespace
