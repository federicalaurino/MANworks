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
 std::cout << "initialize transport problem..."<<std::endl<<std::endl;
 
 import_data();
 set_im_and_fem();
 build_param();
 
 }; // end of init


 // Aux methods for init
	
	//! Import algorithm specifications
	void transport3d1d::import_data(void)
	{
		std::cout<<"init part 1: import data!......" <<std::endl;
	#ifdef M3D1D_VERBOSE_
	cout << "Importing descriptors for tissue and vessel problems ..." << endl;
	#endif
	descr_transp.import_transp(PARAM);
	#ifdef M3D1D_VERBOSE_
	cout << descr_transp;
	#endif
	

	};
	
	
	
	//! Import mesh for tissue (3D) and vessel (1D)  
	//void build_mesh(void); 
	//! Set finite elements methods and integration methods 
	void transport3d1d::set_im_and_fem(void)
	{
	std::cout<<"init part 2: set fem methods!......" <<std::endl;
	

	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	

		
	pfem pf_Ct = fem_descriptor(descr_transp.FEM_TYPET_C);
	pfem pf_Cv = fem_descriptor(descr_transp.FEM_TYPEV_C);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for tissue ..." << endl;
	#endif
		

	mf_Ct.set_finite_element(mesht.convex_index(), pf_Ct);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif

	mf_Cv.set_finite_element(meshv.convex_index(), pf_Cv);

	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
	#endif

	dof_transp.set(mf_Ct, mf_Cv);
	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dof_transp;
	#endif
	

	};
	
	
	//! Build problem parameters
	void transport3d1d::build_param(void)
	{
	std::cout<<"init part 3: build dimensionless parameters!" <<std::endl;
	
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
	#endif
	param_transp.build(PARAM, mf_coeft, mf_coefv);
	#ifdef M3D1D_VERBOSE_
	cout << param_transp ;
	#endif
	cout<<param_transp;
	};
  
  
  
  void transport3d1d::assembly (void)
 {
 std::cout<<"assemble transport problem"<<std::endl<<std::endl;
 	//1 Build the monolithic matrix AM
	assembly_mat();
	//2 Build the monolithic rhs FM
	assembly_rhs();
 }; // end of assembly
 
void 
transport3d1d::assembly_mat(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM, UM, FM ..." << endl;
	#endif
	gmm::resize(AM_transp, dof_transp.tot(), dof_transp_transp.tot()); gmm::clear(AM);
	gmm::resize(UM_transp, dof_transp.tot()); gmm::clear(UM_transp);
	gmm::resize(FM_transp, dof_transp.tot()); gmm::clear(FM_transp);
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM ..." << endl;
	#endif
	// Mass matrix for the interstitial problem
	sparse_matrix_type Mt(dof_transp.Ct(), dof_transp.Ct());
	// Diffusion matrix for the interstitial problem
	sparse_matrix_type Dt(dof_transp.Ct(), dofv.Ct());
	//Transport matrix for interstitial problem
	sparse_matrix_type Bt(dof_transp.Ct(), dofv.Ct());
	// Mass(time derivative)  matrix for the interstitial problem
	sparse_matrix_type Tt(dof_transp.Ct(), dof_transp.Ct());
		
	// Mass matrix for the network problem
	sparse_matrix_type Mv(dof_transp.Cv(), dof_transp.Cv()); //useless
	// Diffusion matrix for the network problem
	sparse_matrix_type Dv(dof_transp.Cv(), dofv.Cv());
	//Transport matrix for network problem
	sparse_matrix_type Bv(dof_transp.Cv(), dofv.Cv());
	// Mass (time derivative)  matrix for the network problem
	sparse_matrix_type Tv(dof_transp.Cv(), dof_transp.Cv());
		/*
	// Junction compatibility matrix for the network problem
	sparse_matrix_type Jvv(dof.Pv(), dof.Uv());
	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof.Pt(), dof.Pt());
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof.Pt(), dof.Pv());
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof.Pv(), dof.Pt());
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof.Pv(), dof.Pv());
	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar(dof.Pv(), dof.Pt());
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof.Pv(), dof.Pt());
	*/
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mt and Dt ..." << endl;
	#endif
	asm_tissue_darcy_transp(Mt, Dt, Tt, mimt, mf_Ct);
	gmm::scale(Mt, (1.0/param.Dalpha(0)); // Dalpha scalar
	gmm::scale(Tt, (1.0/param.dt()); // dt time step
	
	
	// Copy Mtt
	gmm::add(Mt, 
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 
	// Copy -Dtt^T
	/*gmm::add(gmm::scaled(gmm::transposed(Dtt), -1.0),  
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof.Ut()),
					gmm::sub_interval(dof.Ut(), dof.Pt()))); */
	// Copy Dt
	gmm::add(Dt,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling the tangent versor ..." << endl;
	#endif
	vector_type lambdax; // tangent versor: x component
	vector_type lambday; // tangent versor: y component
	vector_type lambdaz; // tangent versor: z component
	std::ifstream ifs(descr.MESH_FILEV);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr.MESH_FILEV);
	asm_tangent_versor(ifs, lambdax, lambday, lambdaz);
	ifs.close();

/*
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mvv and Dvv ..." << endl;
	#endif
	// Local matrices
	size_type shift = 0;
	for(size_type i=0; i<nb_branches; ++i){

		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		scalar_type Ri = param.R(mimv, i);
		// Coefficient \pi^2*Ri'^4/\kappa_v
		vector_type ci(mf_coefvi[i].nb_dof(), pi*pi*Ri*Ri*Ri*Ri/param.kv(i));
		// Allocate temp local matrices
		sparse_matrix_type Mvvi(mf_Uvi[i].nb_dof(), mf_Uvi[i].nb_dof());
		sparse_matrix_type Dvvi(dof.Pv(), mf_Uvi[i].nb_dof());
		// Allocate temp local tangent versor
		vector_type lambdax_K, lambday_K, lambdaz_K;
		for(size_type j=0; j<mf_coefvi[i].nb_dof(); j++){
			lambdax_K.emplace_back(lambdax[i]);
			lambday_K.emplace_back(lambday[i]);
			lambdaz_K.emplace_back(lambdaz[i]);
		}
		// Build Mvvi and Dvvi
		asm_network_poiseuille(Mvvi, Dvvi, 
			mimv, mf_Uvi[i], mf_Pv, mf_coefvi[i],
			ci, lambdax_K, lambday_K, lambdaz_K, meshv.region(i));
		gmm::scale(Dvvi, pi*Ri*Ri);
		// Copy Mvvi and Dvvi
		gmm::add(Mvvi, 
			gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof()), 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof()))); 
		gmm::add(gmm::scaled(gmm::transposed(Dvvi), -1.0),
			gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift,    mf_Uvi[i].nb_dof()),
				gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 
		gmm::add(Dvvi, 
			gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()),
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift,     mf_Uvi[i].nb_dof()))); 
		gmm::clear(Mvvi); 
		gmm::clear(Dvvi);
	*/	
	} /* end of branches loop */
	
	// Build Mvvi and Dvvi
		asm_network_poiseuille(Tv, Dv, 
			mimv,mf_Cv, lambdax_K, lambday_K, lambdaz_K, meshv.region(i));
		gmm::scale(Tv, (1.0/param.dt());
		// Copy Mvvi and Dvvi
		gmm::add(Tv, 
			gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.tot()), 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.tot()))); 
		gmm::add(Dv, 
			gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.tot()), 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.tot())));
	/*
if (nb_junctions > 0){
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Jvv" << " ..." << endl;
	#endif
	asm_network_junctions(Jvv, mimv, mf_Uvi, mf_Pv, mf_coefv, 
		Jv, param.R());
	#ifdef M3D1D_VERBOSE_
	cout << "  Copying -Jvv^T" << " ..." << endl;
	#endif		
	gmm::add(gmm::scaled(gmm::transposed(Jvv), -1.0),
		gmm::sub_matrix(AM,
			 gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()), 
			 gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())));
	#ifdef M3D1D_VERBOSE_
	cout << "  Copying Jvv" << " ..." << endl;
	#endif		
	gmm::add(Jvv,
		gmm::sub_matrix(AM,
			 gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()),
			 gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
}
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling aux exchange matrices Mbar and Mlin ..." << endl;
	#endif
	asm_exchange_aux_mat(Mbar, Mlin, 
			mimv, mf_Pt, mf_Pv, param.R(), descr.NInt);
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION");
	asm_exchange_mat(Btt, Btv, Bvt, Bvv,
			mimv, mf_Pv, mf_coefv, Mbar, Mlin, param.Q(), NEWFORM);

	// Copying Btt
	gmm::add(Btt, 
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut(), dof.Pt()), 
					gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	// Copying -Btv
	gmm::add(gmm::scaled(Btv, -1),
	 		  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut(), dof.Pt()),
					gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); 
	// Copying -Bvt
	gmm::add(gmm::scaled(Bvt,-1), 
			  gmm::sub_matrix(AM, 
			  		gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv()	, dof.Pv()),
					gmm::sub_interval(dof.Ut(), dof.Pt()))); 
	// Copying Bvv
	gmm::add(Bvv, 
			  gmm::sub_matrix(AM, 
					gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()), 
					gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()))); */

	// De-allocate memory
	gmm::clear(Mt);  gmm::clear(Dt); 
	gmm::clear(Mv); gmm::clear(Dv);
	gmm::clear(Bt);  gmm::clear(Bv);

}

void 
transport3d1d::assembly_rhs(void)
{
 cout<<"assembling rhs vector..."<< endl; 
 
 for(int i = 0; i<dof_transp.tot(); i++)
 FM_transp(i)= 0.005* (i/i+10);
	/*
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic rhs FM ... " << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Initializing RHS for FM ..." << endl;
	#endif
	// Right Hand Side for the interstitial problem 
	vector_type Ft(dof.Ut());
	// Right Hand Side for the network problem 
	vector_type Fv(dof.Uv());

	// Coefficients for tissue BCs
	scalar_type bcoef  = PARAM.real_value("BETA", "Coefficient for mixed BC");
	scalar_type p0coef = PARAM.real_value("P0"); // default: 0

	#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term ..." << endl;
	#endif
	vector_type beta(dof.coeft(), 1.0/bcoef);
	vector_type P0(dof.coeft(), p0coef);
	
	if (PARAM.int_value("TEST_RHS")) {
		#ifdef M3D1D_VERBOSE_
		cout << "  ... as the divergence of exact velocity ... " << endl;
		#endif
		assembly_tissue_test_rhs();
	}
	else {
		sparse_matrix_type Mtt(dof.Ut(), dof.Ut());
		asm_tissue_bc(Mtt, Ft, mimt, mf_Ut, mf_coeft, BCt, P0, beta);
		gmm::add(Mtt, 
			gmm::sub_matrix(AM,
				gmm::sub_interval(0, dof.Ut()),
				gmm::sub_interval(0, dof.Ut())));
		gmm::add(Ft, gmm::sub_vector(FM, gmm::sub_interval(0, dof.Ut())));
		// De-allocate memory
		gmm::clear(Mtt); 
	}
	#ifdef M3D1D_VERBOSE_
	cout << "  Building vessel boundary term ..." << endl;
	#endif
	sparse_matrix_type Mvv(dof.Uv(), dof.Uv());
	asm_network_bc(Mvv, Fv, 
			mimv, mf_Uvi, mf_coefv, BCv, P0, param.R(), bcoef);
	gmm::add(Mvv, 
		gmm::sub_matrix(AM,
			gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()),
			gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
	gmm::add(Fv, gmm::sub_vector(FM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
	// De-allocate memory
	gmm::clear(Ft); gmm::clear(Fv); gmm::clear(Mvv);
	*/
	
}

 bool transport3d1d::solve (void)
 {
  std::cout<<"solve transport problem"<<std::endl<<std::endl;
  #ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	gmm::csc_matrix<scalar_type> A_transp;
	gmm::clean(AM_transp, 1E-12);
	gmm::copy(AM_transp, A_transp);
	gmm::clear(AM_transp); // to be postponed for preconditioner
	double time = gmm::uclock_sec();	
	
	//scalar_type nn_timestep = ceil(param_transp.T() / param_transp.dt());
	
	double time_count = 0;
	string time_suff = "";
	ostringstream convert;
	for(int t=0; t<=param_transp.T() ; t = t + param_transp.dt()){
	time_count++; 
	
	if ( descr.SOLVE_METHOD == "SuperLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(A_transp, UM_transp, FM_transp, cond);
		cout << "  Condition number (transport problem): " << cond << endl;
	}
	else { // Iterative solver //

		// Iterations
		gmm::iteration iter(descr.RES);  // iteration object with the max residu
		iter.set_noisy(1);               // output of iterations (2: sub-iteration)
		iter.set_maxiter(descr.MAXITER); // maximum number of iterations

		// Preconditioners
		//! \todo Add preconditioner choice to param file
		// See \link http://download.gna.org/getfem/html/homepage/gmm/iter.html
		gmm::identity_matrix PM; // no precond
		//gmm::diagonal_precond<sparse_matrix_type> PM(AM); // diagonal preocond
		//gmm::ilu_precond<sparse_matrix_type> PM(AM);
		// ...
		//gmm::clear(AM);
		// See <http://download.gna.org/getfem/doc/gmmuser.pdf>, pag 15
	
		if ( descr.SOLVE_METHOD == "CG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Conjugate Gradient method ... " << endl;
			#endif
			gmm::identity_matrix PS;  // optional scalar product
			gmm::cg(AM_transp, UM_transp, FM_transp, PS, PM, iter);
		}
		else if ( descr.SOLVE_METHOD == "BiCGstab" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
			#endif
			gmm::bicgstab(AM, UM, FM, PM, iter);
		}
		else if ( descr.SOLVE_METHOD == "GMRES" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Generalized Minimum Residual method ... " << endl;
			#endif
			size_type restart = 50;
			gmm::gmres(A, UM, FM, PM, restart, iter);
		}
		else if ( descr.SOLVE_METHOD == "QMR" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
			#endif
			gmm::qmr(AM, UM, FM, PM, iter);
		}
		else if ( descr.SOLVE_METHOD == "LSCG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
			#endif
			gmm::least_squares_cg(AM, UM, FM, iter);
		}
		// Check
		if (iter.converged())
			cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
		else if (iter.get_iteration() == descr.MAXITER)
			cerr << "  ... reached the maximum number of iterations!" << endl;

	}
	
	convert << time_count;
	time_suff = convert.str();
	export_vtk(time_suff);
	gmm::clear(FM_transp);
	gmm::add(UM_transp, FM_transp);
	gmm::clear(UM_transp);
	
	
	} //end of cycle over time 
	
	cout << "... time to solve : " << gmm::uclock_sec() - time << " seconds\n";

	#ifdef M3D1D_VERBOSE_
	cout << "Compute the total flow rate ... " << endl;
	#endif
	// Aux vector
	vector_type Uphi(dof.Pv()); 
	// Extracting matrices Bvt, Bvv
	sparse_matrix_type Bvt(dof.Pv(), dof.Pt());
	sparse_matrix_type Bvv(dof.Pv(), dof.Pv());
	gmm::copy(gmm::sub_matrix(A, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv()	, dof.Pv()),
			gmm::sub_interval(dof.Ut(), dof.Pt())),
				Bvt); 
	gmm::copy(gmm::sub_matrix(A, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()), 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
				Bvv); 
	// Extracting solutions Pt, Pv 
	vector_type Pt(dof.Pt()); 
	vector_type Pv(dof.Pv()); 
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut(), dof.Pt())), Pt);
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())), Pv);
	// Computing Bvv*Pv - Bvt*Pt
	gmm::mult(Bvt, Pt, Uphi);
	gmm::mult_add(Bvv, Pv, Uphi);
	TFR = std::accumulate(Uphi.begin(), Uphi.end(), 0.0);

	// De-allocate memory
	gmm::clear(Bvt); gmm::clear(Bvv);
	gmm::clear(Pt);  gmm::clear(Pv);  
	gmm::clear(Uphi);

	return true;
 }; // end of solve
	
	
 void transport3d1d::export_vtk (const string & time_suff , const string & suff)
 {
  std::cout<<"export transport problem"<<std::endl<<std::endl;
  
  if (PARAM.int_value("VTK_EXPORT"))
  {
	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif

	// Array of unknown dof of the interstitial concentration
	vector_type Ct(dof_transp.Ct());  
	// Array of unknown dof of the network concentration
	vector_type Cv(dof_transp.Cv()); 
	gmm::copy(gmm::sub_vector(UM_transp, 
		gmm::sub_interval(0, dof_transp.Ct())), Ct);
	gmm::copy(gmm::sub_vector(UM_transp, 
		gmm::sub_interval(dof_transp.Ct(), dof_transp.tot())), Cv);

/*
	#ifdef M3D1D_VERBOSE_
	// Save vessel solution for test-cases
	if (nb_branches==1){;
		std::ofstream outPv("Pv.txt");
		outPv << gmm::col_vector(Pv);
		outPv.close();
	}
	
	#endif
	
*/
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ut ..." << endl;
	#endif
//	pfem pf_Ut = fem_descriptor(descr.FEM_TYPET);
//	if(pf_Ut->is_lagrange()==0){ 
		/*
			There is no built-in export for non-lagrangian FEM.
			If this is the case, we need to project before exporting.
		 */
/*		#ifdef M3D1D_VERBOSE_
		cout << "    Projecting Ut on P1 ..." << endl;
		#endif
		mesh_fem mf_P1(mesht);
		mf_P1.set_qdim(bgeot::dim_type(DIMT)); 
		mf_P1.set_classical_finite_element(1);
		sparse_matrix_type M_RT0_P1(mf_P1.nb_dof(), dof.Ut());
		sparse_matrix_type M_P1_P1(mf_P1.nb_dof(), mf_P1.nb_dof());
		vector_type Ut_P1(mf_P1.nb_dof());
		asm_mass_matrix(M_RT0_P1, mimt, mf_P1, mf_Ut);
		asm_mass_matrix(M_P1_P1,  mimt, mf_P1, mf_P1);
		
		vector_type Utt(mf_P1.nb_dof());
		gmm::mult(M_RT0_P1, Ut, Utt);
		double cond1;
		gmm::SuperLU_solve(M_P1_P1, Ut_P1, Utt, cond1);

		vtk_export exp1(descr.OUTPUT+"Ut.vtk");
		exp1.exporting(mf_P1);
		exp1.write_mesh();
		exp1.write_point_data(mf_P1, Ut_P1, "Ut");
	}	
	else {
		vtk_export exp_Ut(descr.OUTPUT+"Ut.vtk");
		exp_Ut.exporting(mf_Ut);
		exp_Ut.write_mesh();
		exp_Ut.write_point_data(mf_Ut, Ut, "Ut");	 
	}
*/	
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ct ..." << endl;
	#endif
	vtk_export exp_Ct(descr_transp.OUTPUT+"Ct"+suff+"_"+time_suff+".vtk");
	exp_Pt.exporting(mf_Ct);
	exp_Pt.write_mesh();
	exp_Pt.write_point_data(mf_Ct, Ct, "Ct");

/*	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Uv ..." << endl;
	#endif
	size_type start = 0;
	size_type length = 0;
	for(size_type i=0; i<nb_branches; ++i){
		if(i>0) start += mf_Uvi[i-1].nb_dof();
		length = mf_Uvi[i].nb_dof();
		vtk_export exp_Uv(descr.OUTPUT+"Uv"+suff+std::to_string(i)+".vtk");
		exp_Uv.exporting(mf_Uvi[i]);
		exp_Uv.write_mesh();
		exp_Uv.write_point_data(mf_Uvi[i], 
			gmm::sub_vector(Uv, gmm::sub_interval(start, length)), "Uv"); 
	}
*/
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv ..." << endl;
	#endif
	vtk_export exp_Cv(descr_transp.OUTPUT+"Cv"+suff+"_"+time_suff+".vtk");
	exp_Pv.exporting(mf_Cv);
	exp_Pv.write_mesh();
	exp_Pv.write_point_data(mf_Cv, Cv, "Cv");

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif
  }
 }; // end of export
  
  
 
 } // end of namespace
