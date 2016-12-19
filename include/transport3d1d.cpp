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
 build_mesh();
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
	descr_transp.import(PARAM);
	#ifdef M3D1D_VERBOSE_
	cout << descr_transp;
	#endif
	 
 
	};
	
	 
	
	//! Import mesh for tissue (3D) and vessel (1D)  
	void transport3d1d::build_mesh(void){
	//no need to build again the  3d mesh
	//but, in order to have the boundary conditions for the nodes
	//we need to build again the 1D mesh from another pts file
	
	#ifdef M3D1D_VERBOSE_
	cout << "Importing the 1D mesh for the vessel (transport problem)... "   << endl;
	#endif
	std::ifstream ifs(descr_transp.MESH_FILEV);
	mesh meshv_transp;
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr_transp.MESH_FILEV);
	import_pts_file(ifs, meshv_transp, BCv_transp, nb_vertices, descr.MESH_TYPEV);
	ifs.close();
	
	

	
	};
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
  
  
  
void
transport3d1d::build_tissue_boundary (void) 
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building tissue boundary ..." << endl;
	#endif
	BCt_transp.clear();
	BCt_transp.reserve(2*DIMT);
	// Parse BC data
	string label_in = PARAM.string_value("BClabel_transp", "Array of tissue boundary labels");
	string value_in = PARAM.string_value("BCvalue_transp", "Array of tissue boundary values");
	vector<string> labels = split(label_in, ' ');
	vector<string> values = split(value_in, ' ');
	GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
	GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
	for (unsigned f=0; f<2*DIMT; ++f) {
		BCt.emplace_back(labels[f], std::stof(values[f]), 0, f);
		#ifdef M3D1D_VERBOSE_
		cout << "  face " << f << " : " << BCt.back() << endl;
		#endif
	}
		
}

void 
transport3d1d::build_vessel_boundary(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building vessel boundary ..." << endl;
	#endif
try {

	dal::bit_vector junctions; // global idx of junctions vertices in meshv
	dal::bit_vector extrema;   // global idx of extreme vertices in meshv

	Jv.clear();
	nb_extrema=0; 
	nb_junctions=0;
	
	size_type fer = nb_branches; // first empty region
	GMM_ASSERT1(meshv.has_region(fer)==0, 
		"Overload in meshv region assembling!");
	
	// List all the convexes
	dal::bit_vector nn = meshv.convex_index();
	bgeot::size_type cv;
	for (cv << nn; cv != bgeot::size_type(-1); cv << nn) {
		
		bgeot::pconvex_structure cvs = meshv.structure_of_convex(cv);
		if (cvs->nb_points()>2) 
			cerr << "Error: convex #" << cv << "has more than 2 vertices!" << endl;
		if (cvs->nb_faces()>2)  
			cerr << "Error: convex #" << cv << "has more than 2 faces!" << endl;

		// Build regions for BCs and junctions
		// Global idx of mesh vertices
		size_type i0 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(1)[0]];
		size_type i1 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(0)[0]];
		// Identify vertex type
		if (meshv.convex_to_point(i0).size()==1){ /* inflow extremum */
			// Update information
			extrema.add(i0);
			nb_extrema++;
			// Build a new region made by a single face
			GMM_ASSERT1(meshv.has_region(fer)==0, 
				"Overload in meshv region assembling!");
			meshv.region(fer).add(cv, 1);
			// Store the current index and then update it
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv.size())) {
				found = (i0 == BCv[bc].idx);
				if (!found) bc++;
			}
			GMM_ASSERT1(found=true, "Miss a boundary node in BCv list!");
			BCv[bc].rg = fer; 
			fer++;
			// Store the containing branch index
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			BCv[bc].branches.emplace_back(branch); 
		}
		else if (meshv.convex_to_point(i0).size()==2){ /* trivial inflow junction */
			// DO NOTHING
		}
		else if (meshv.convex_to_point(i0).size()>=2){ /* non-trivial inflow junction */
			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i0);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i0);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 1); // single-face region
				// Create a new junction node
				Jv.emplace_back("JUN", 0, i0, fer);
				fer++;
			}
			// Search for index of containing branch (\mathcal{P}^{in}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			// Add the inflow branch (to the right junction node)
			size_type jj = 0;
			bool found = false;
			while (!found && jj < nb_junctions){
				found = (i0 == Jv[jj].idx);
				if (!found) jj++;
			}
			//cout << "Branch -" << branch << " added to junction " << jj << endl;
			Jv[jj].value += param.R(mimv, branch);
			Jv[jj].branches.emplace_back(-branch);
			GMM_ASSERT1(branch>0, 
				"Error in network labeling: -0 makes no sense");
		}
		
		if (meshv.convex_to_point(i1).size()==1){ 
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv.size())) {
				found = (i1 == BCv[bc].idx);
				if (!found) bc++;
			}
			if (found){ /* outlow extremum */
				extrema.add(i1); 
				nb_extrema++; 
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				// Store the current index and then update it
				BCv[bc].value *= -1.0;
				BCv[bc].rg = fer; 
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv[bc].branches.emplace_back(branch); 
			}
			else { /* interior -> Mixed point */
				// "MIX" label via post-processing
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				BCv.emplace_back("MIX", 0.0, i1, fer);
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv.back().branches.emplace_back(branch); 
			}
		}
		else if (meshv.convex_to_point(i1).size()==2){ /* trivial outflow junction */

			// Search for index of first containing branch (\mathcal{P}^{out}_j)
			size_type firstbranch = 0; 
			bool contained = false;
			while (!contained && firstbranch<nb_branches ) {
				contained = meshv.region(firstbranch).is_in(cv);
				if (!contained) firstbranch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i1!");

			// Check if i1 is a trivial junction (or a INT point)
			size_type cv1 = meshv.convex_to_point(i1)[0];
			size_type cv2 = meshv.convex_to_point(i1)[1];
			bool is_junc = (meshv.region(firstbranch).is_in(cv1) < 1 ||
							meshv.region(firstbranch).is_in(cv2) < 1 );
							
			if (is_junc){
				cout << "Found a trivial junction at i1 = " << i1 << endl;
				// Check if jucntion has been already stored, 
				// if not add to the junction list (J) and build a new region
				dal::bit_vector tmp; tmp.add(i1);
				if(!junctions.contains(tmp)){
					// Store the junction vertex
					junctions.add(i1);
					nb_junctions++;
					GMM_ASSERT1(meshv.has_region(fer)==0, 
						"Overload in meshv region assembling!");
					// Build a new region with idx "first empty region"
					meshv.region(fer).add(cv, 0);
					// Create a new junction node
					Jv.emplace_back("JUN", 0, i1, fer);
					fer++;
				}
				// Search for index of second containing branch (\mathcal{P}^{out}_j)
				size_type secondbranch = firstbranch+1; 
				size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
				contained = false;
				while (!contained && secondbranch<nb_branches ) {
					contained = meshv.region(secondbranch).is_in(secondcv);
					if (!contained) secondbranch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				// Add the two branches
				Jv.back().branches.emplace_back(+firstbranch);
				Jv.back().branches.emplace_back(-secondbranch);
				Jv.back().value += param.R(mimv, firstbranch);
				Jv.back().value += param.R(mimv, secondbranch);
			}
		}
		else if (meshv.convex_to_point(i1).size()>=2){ /* non-trivial outflow junction */

			// Search for index of containing branch (\mathcal{P}^{out}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");

			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i1);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i1);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 0);
				// Create a new junction node
				Jv.emplace_back("JUN", 0, i1, fer);
				// Add the outflow branch
				Jv.back().branches.emplace_back(+branch);
				Jv.back().value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << i1 << endl;
				fer++;
			}
			else {
				// Add the outflow branch (to the right junction node)
				size_type jj = 0;
				bool found = false;
				while (!found && jj < nb_junctions){
					found = (i1 == Jv[jj].idx);
					if (!found) jj++;
				}
				Jv[jj].branches.emplace_back(+branch);
				Jv[jj].value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << jj << endl;
			}
		}

	} /* end of convexes loop */
	
	// Ckeck network assembly
	#ifdef M3D1D_VERBOSE_
	cout << "--- NETWORK ASSEMBLY ------------------ "   << endl;
	cout << "  Branches:   " << nb_branches << endl
		 << "  Vertices:   " << nn.size()+1 << endl;
	cout << "  Extrema:    " << extrema << endl;	  
	for (size_type i=0; i<BCv.size(); ++i)
		cout << "    -  label=" << BCv[i].label 
			 << ", value=" << BCv[i].value << ", ind=" << BCv[i].idx 
			 << ", rg=" << BCv[i].rg << ", branches=" << BCv[i].branches << endl; 
	cout << "  Junctions: " << junctions << endl;
	for (size_type i=0; i<Jv.size(); ++i)
		cout << "    -  label=" << Jv[i].label 
			 << ", value=" << Jv[i].value << ", ind=" << Jv[i].idx 
			 << ", rg=" << Jv[i].rg << ", branches=" << Jv[i].branches << endl; 
	cout << "---------------------------------------- "   << endl;
	#endif

} 
GMM_STANDARD_CATCH_ERROR; // catches standard errors

} /* end of build_vessel_boundary */


  
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
	gmm::resize(AM_transp, dof_transp.tot(), dof_transp.tot()); gmm::clear(AM_transp);
	gmm::resize(UM_transp, dof_transp.tot()); gmm::clear(UM_transp);
	gmm::resize(FM_transp, dof_transp.tot()); gmm::clear(FM_transp);
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM ..." << endl;
	#endif
	// Reaction matrix for the interstitial problem
	sparse_matrix_type Rt(dof_transp.Ct(), dof_transp.Ct()); gmm::clear(Rt);
	// Diffusion matrix for the interstitial problem
	sparse_matrix_type Dt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Dt);
	//Transport matrix for interstitial problem
	sparse_matrix_type Bt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Bt);
	// Mass(time derivative)  matrix for the interstitial problem
	sparse_matrix_type Mt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Mt);
		
	// Reaction matrix for the network problem
	sparse_matrix_type Rv(dof_transp.Cv(), dof_transp.Cv()); gmm::clear(Rv);
	// Diffusion matrix for the network problem
	sparse_matrix_type Dv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Dv);
	//Transport matrix for network problem
	sparse_matrix_type Bv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bv);
	// Mass (time derivative)  matrix for the network problem
	sparse_matrix_type Mv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mv);
		
	// Junction compatibility matrix for the network problem
	//sparse_matrix_type Jvv(dof.Pv(), dof.Uv());
	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Btt);
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof_transp.Ct(), dof_transp.Cv());gmm::clear(Btv);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Bvt);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bvv);
	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mbar);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mlin);
	
	//matrice identità 
	//sparse_matrix_type Identity(dof_transp.tot(),dof_transp.tot()); gmm::clear(Identity);
	//gmm::copy(gmm::identity_matrix(), Identity);
	//gmm::add(Identity, AM_transp);

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Rt, Mt and Dt ..." << endl;
	#endif
	
	//assemble diffusion and reaction terms
	vector_type mass_coeff(dof.coeft()); gmm::clear(mass_coeff);
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut(), dof.Pt())) ,  mass_coeff);
	gmm::scale (param_transp.Q_pl(), mass_coeff);
	gmm::add(param_transp.Dalpha(), mass_coeff);
	
	asm_tissue_darcy_transp(Rt, Dt, Mt, mimt, mf_Ct, mf_coeft, mass_coeff, param_transp.At() ); //vedi file assembling3d_transp.hpp
	gmm::scale(Mt, (1.0/param_transp.dt())); // dt time step
	
	gmm::MatrixMarket_IO::write("Dt.mm",Dt);
	gmm::MatrixMarket_IO::write("Mt.mm",Mt);
	
	
	// Copy Rt //MATRICE DI REAZIONE PER IL CONSUMO DI OSSIGENO, MOMENTANEAMENTE TOLTA
	/*gmm::add(Rt, 
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 
	*/
	
	// Copy Tt //PASSO TEMPORALE
	
	gmm::add(Mt,  
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 


	// Copy Dt //DIFFUSIONE
	gmm::add(Dt,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 
			 
			
		#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mv and Dv ..." << endl;
	#endif
	
	// Build Mvvi and Dvvi
		asm_network_poiseuille_transp(Dv, Mv, mimv,mf_Cv);
		gmm::scale(Mv, (1.0/param_transp.dt()));
		
		// mass matrix for TIME STEP for vessel
		gmm::add(Mv, 
			gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 
		
		// diffusion matrix for vessels
		gmm::add(Dv, 
			gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
				
	gmm::MatrixMarket_IO::write("Dv.mm",Dv);
	gmm::MatrixMarket_IO::write("Mv.mm",Mv);
		
		
		
	//ADVECTION	
		#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Bt and Bv ..." << endl;
	#endif		
	bool ADVECTION = PARAM.int_value("ADVECTION");
	ADVECTION=0; //momentaneamente non mettiamolo il trasporto
	if(ADVECTION ==1){
	// advection tissue term: 				
	vector_type Ut(dof.Ut()); gmm::clear(Ut);
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())) ,  Ut);
	asm_advection_matrix(Bt, mimt, mf_Ct, mf_Ut,Ut);
	gmm::add(Bt,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 
			
		
					
	/* assembla il termine di trasporto per le rete: 
	    sarà necessario combinare in qualche modo gli mf_Uvi con mf_Cv 
	    (magari costruendo anche per il trasporto un vettore mf_Cvi di fem_methods
	    
	// advection vessel term: 				
	vector_type Uv(dof.Uv()); gmm::clear(Uv);
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(0, dof.Uv())) ,  Uv);
	asm_advection_matrix(Bv, mimv, mf_Cv, mf_Uv,Uv);
	gmm::add(Bv,
			  gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	*/
	
	
	}
	else{
	cout<<"---NO ADVECTION TERM ---"<<endl;
	}
				

	/*  blocco per i termini di junctions... rivedere il modello: serve anche nelle equazioni del trasporto?*/
	


// TERMINI DI ACCOPPIAMENTO 

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling aux exchange matrices Mbar and Mlin ..." << endl;
	#endif
	asm_exchange_aux_mat(Mbar, Mlin, 
			mimv, mf_Ct, mf_Cv, param.R(), descr.NInt);
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION");
	NEWFORM=false;
	vector_type coeff(dof.Pv());
	
	asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_coefv, Mbar, Mlin, param_transp.Y(), NEWFORM);

	// Copying Btt
	gmm::add(Btt, 
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 
	// Copying Btv
	gmm::add(gmm::scaled(Btv, -1.0),
	 		  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()),
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 
	// Copying -Bvt
	gmm::add(gmm::scaled(Bvt, -1.0),  
			  gmm::sub_matrix(AM_transp, 
			  		gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()),
					gmm::sub_interval(0, dof_transp.Ct()))); 
	// Copying Bvv
	gmm::add(Bvv, 
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 

	gmm::MatrixMarket_IO::write("Btt.mm",Btt);
	gmm::MatrixMarket_IO::write("Btv.mm",Btv);
	gmm::MatrixMarket_IO::write("Bvt.mm",Bvt);
	gmm::MatrixMarket_IO::write("Bvv.mm",Bvv);
	
	
	// De-allocate memory
	gmm::clear(Mt);  gmm::clear(Dt); 
	gmm::clear(Mv); gmm::clear(Dv);
	gmm::clear(Bt);  gmm::clear(Bv);
	gmm::clear(Mbar);  gmm::clear(Mlin);
	gmm::clear(Btt);  gmm::clear(Btv);
	gmm::clear(Bvt);  gmm::clear(Bvv);
}

void 
transport3d1d::assembly_rhs(void)
{
 cout<<"assembling rhs vector..."<< endl;  
 
 /*genera una CONDIZIONE INIZIALE NON NULLA nella mesh 3d
 
 for(int i = 0; i<dof_transp.Ct()/4; i++)
 FM_transp[i]= 0.0005* i/(i+3);
 */
 cout<<"assembled rhs vector..."<< endl; 
	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic rhs FM ... " << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Initializing RHS for FM ..." << endl;
	#endif
	// Right Hand Side for the interstitial problem 
	vector_type Ft(dof_transp.Ct());
	// Right Hand Side for the network problem 
	vector_type Fv(dof_transp.Cv());

	// Coefficients for tissue BCs
	scalar_type bcoef  = PARAM.real_value("BETA_transp", "Coefficient for mixed BC for transport problem");
	scalar_type penalty  = PARAM.real_value("PENALTY", "Penalty coefficient for Dirichlet BC");
	if (penalty==0) penalty =1E-30;
	penalty=1./penalty;

	//Al momento le condizioni al contorno di dirichlet in inflow sono fatte con un metodo di penalizzazione:
	//NON CORRETTO, induce errore non accettabile; sistemare con rilavemento.
	//vedi in file d'esempio laplacian.cc la funzione asm_dirichlet_condition 
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term ..." << endl;
	#endif
	vector_type beta(dof.coeft(), 1.0/bcoef);

	
	if (1==0) {
		#ifdef M3D1D_VERBOSE_
		cout << "  ... as the divergence of exact velocity ... " << endl;
		#endif
		assembly_tissue_test_rhs();
	}
	else { cout<<"assemble BC for tissue"<<endl;
		sparse_matrix_type Mtt(dof_transp.Ct(), dof_transp.Ct());
		asm_tissue_bc_transp(Mtt, mimt, mf_Ct, mf_coeft, BCt_transp,beta);
		gmm::add(Mtt, 
			gmm::sub_matrix(AM_transp,
				gmm::sub_interval(0, dof_transp.Ct()),
				gmm::sub_interval(0, dof_transp.Ct()))); 
					
		// De-allocate memory
		gmm::clear(Mtt); 
		
		cout<<"assemble BC for vessel"<<endl;
	#ifdef M3D1D_VERBOSE_
	cout << "  Building vessel boundary term ..." << endl;
	#endif
		sparse_matrix_type Mvv(dof_transp.Cv(), dof_transp.Cv());
		
		vector_type Fv(dof_transp.Cv());

		asm_network_bc_transp(Mvv, Fv, mimv, mf_Cv, mf_coefv, BCv_transp, param_transp.Av(), penalty );
		gmm::add(Mvv, 
			gmm::sub_matrix(AM_transp,
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())));
		gmm::add(Fv, 
			gmm::sub_vector(FM_transp,
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
		// De-allocate memory
		gmm::clear(Mvv);
		gmm::clear(Fv);
	}
	
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
	
	
	vector_type F_transp(gmm::vect_size(FM_transp));
	gmm::clean(FM_transp, 1E-12);
	gmm::copy(FM_transp, F_transp);
	
	//gmm::clear(AM_transp); // to be postponed for preconditioner
	double time = gmm::uclock_sec();	
	
	gmm::MatrixMarket_IO::write("Atransp.mm",A_transp);	
	double time_count = 0;
	for(double t=0;t<=param_transp.T() ; t = t + param_transp.dt()){ 
	time_count++; 
	std::cout<<"iteration number:"<<time_count<<std::endl;
	
	
	if ( descr_transp.SOLVE_METHOD == "SuperLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(A_transp, UM_transp, F_transp, cond);
		cout << "  Condition number (transport problem): " << cond << endl;
	}
	else { // Iterative solver //

		// Iterations
		gmm::iteration iter(descr_transp.RES);  // iteration object with the max residu
		iter.set_noisy(1);               // output of iterations (2: sub-iteration)
		iter.set_maxiter(descr_transp.MAXITER); // maximum number of iterations

		// Preconditioners
		//! \todo Add preconditioner choice to param file
		// See \link http://download.gna.org/getfem/html/homepage/gmm/iter.html
		gmm::identity_matrix PM; // no precond
		//gmm::diagonal_precond<sparse_matrix_type> PM(AM); // diagonal preocond
		//gmm::ilu_precond<sparse_matrix_type> PM(AM);
		// ...
		//gmm::clear(AM);
		// See <http://download.gna.org/getfem/doc/gmmuser.pdf>, pag 15
	
		if ( descr_transp.SOLVE_METHOD == "CG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Conjugate Gradient method ... " << endl;
			#endif
			gmm::identity_matrix PS;  // optional scalar product
			gmm::cg(AM_transp, UM_transp, F_transp, PS, PM, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "BiCGstab" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
			#endif
			gmm::bicgstab(AM, UM, FM, PM, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "GMRES" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Generalized Minimum Residual method ... " << endl;
			#endif
			size_type restart = 50;
			gmm::gmres(A_transp, UM, FM, PM, restart, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "QMR" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
			#endif
			gmm::qmr(AM, UM, FM, PM, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "LSCG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
			#endif
			gmm::least_squares_cg(AM, UM, FM, iter);
		}
		// Check
		if (iter.converged())
			cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
		else if (iter.get_iteration() == descr_transp.MAXITER)
			cerr << "  ... reached the maximum number of iterations!" << endl;

	}
	
	//export solution
	std::cout<<"solved! going to export..."<<std::endl;
	string time_suff = "";
	std::ostringstream convert;
	convert << time_count;
	time_suff = convert.str();
	export_vtk(time_suff); 
	
	std::ofstream outFF(descr_transp.OUTPUT+"FF"+"_t"+time_suff+".txt");
		outFF << gmm::col_vector(F_transp);
		outFF.close();
		
	std::cout<<"exported! now new iteration..."<<std::endl;
	// update rhs (there is the time step mass term)
	gmm::clear(F_transp);
	gmm::copy(FM_transp, F_transp); 
	vector_type TFt(dof_transp.Ct());
	vector_type TFv(dof_transp.Cv());
	//vector_type Ct(dof_transp.Ct());
	//vector_type Cv(dof_transp.Cv());
	//gmm::add(gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())), Ct);
	//gmm::add(gmm::sub_vector(UM_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);
	asm_source_term(TFt,mimt, mf_Ct, mf_Ct,gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct()))); // tentativo con generic_assembly: lo trovi in assembling3d_transp, la funzione è asm_time_rhs_transp
	asm_source_term(TFv,mimv, mf_Cv, mf_Cv,gmm::sub_vector(UM_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	gmm::scale(TFt, (1.0/param_transp.dt())); // dt time step
	gmm::scale(TFv, (1.0/param_transp.dt())); // dt time step
	gmm::add(TFt, gmm::sub_vector(F_transp, gmm::sub_interval(0, dof_transp.Ct())));
	gmm::add(TFv, gmm::sub_vector(F_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	gmm::clear(UM_transp);
	gmm::clear(TFt); gmm::clear(TFv);
	
	} //end of cycle over time 
	
	cout << "... time to solve : " << gmm::uclock_sec() - time << " seconds\n";

	return true;
 }; // end of solve
	
	
 void transport3d1d::export_vtk (const string & time_suff,const string & suff)
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
	
	// Array of unknown dof of the interstitial velocity
	vector_type Ct(dof_transp.Ct()); 

	// Array of unknown dof of the network velocity
	vector_type Cv(dof_transp.Cv()); 

	gmm::copy(gmm::sub_vector(UM_transp, 
		gmm::sub_interval(0, dof_transp.Ct())), Ct);
	gmm::copy(gmm::sub_vector(UM_transp, 
		gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);

	

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ct ..." << endl;
	#endif
	vtk_export exp_Ct(descr_transp.OUTPUT+"Ct"+suff+"_t"+time_suff+".vtk");
	exp_Ct.exporting(mf_Ct);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct, Ct, "Ct");



	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv ..." << endl;
	#endif
	vtk_export exp_Cv(descr_transp.OUTPUT+"Cv"+suff+"_t"+time_suff+".vtk");
	exp_Cv.exporting(mf_Cv);
	exp_Cv.write_mesh();
	exp_Cv.write_point_data(mf_Cv, Cv, "Cv");

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif
  }
 }; // end of export
 
  
  
	void transport3d1d::test(void){

// risolvi il problema: (At + Btt)Ct = Btv Cv
// stazionario
//con Cv considerato costante e uguale a 1 su tutta la mesh
//At è la matrice di stiffness per il laplaciano
//Btt e Btv sono le matrici di accoppiamento tra le mesh 3d - 1d

	//definizioni di matrici e vettori
	vector_type Ct(dof_transp.Ct()); gmm::clear(Ct);
	sparse_matrix_type At(dof_transp.Ct(), dof_transp.Ct());gmm::clear(At);
	sparse_matrix_type Btt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Btt);
	sparse_matrix_type Btv(dof_transp.Ct(), dof_transp.Cv());gmm::clear(Btv);
	sparse_matrix_type Bvt(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Bvt); // non servono, ma devo comunque passarle alla funzione asm_exchange_mat
	sparse_matrix_type Bvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bvv); //	
	vector_type Cv(dof_transp.Cv(), 1.0);
	vector_type F(dof_transp.Ct()); gmm::clear(F);  //F= Btv Cv
	
	//stiffness matrix
	getfem::generic_assembly
	  assem("M$1(#1,#1) += sym(comp(vGrad(#1).vGrad(#1)) (:, i,k, : ,i,k) )");
	  assem.push_mi(mimt);
	  assem.push_mf(mf_Ct);
	  assem.push_mat(At);
	  assem.assembly();
	
	//matrici di accoppiamento
	sparse_matrix_type Mbar(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mbar);
	sparse_matrix_type Mlin(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mlin);
  	asm_exchange_aux_mat(Mbar, Mlin, mimv, mf_Ct, mf_Cv, param.R(), descr.NInt);
	
	bool NEWFORM = true;	
	asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_coefv, Mbar, Mlin, param_transp.Y(), NEWFORM);

	// Copying Btt and add it to At
	gmm::add(Btt, At); 

	// F è il termine noto, F= Btv*Cv
	gmm::mult(Btv, Cv, F);

	//costruisci una matrice del tipo adatto per SuperLU_solver
	gmm::csc_matrix<scalar_type> A;
	gmm::clean(At, 1E-12); 
	gmm::copy(At, A);
	scalar_type cond;
	gmm::SuperLU_solve(A, Ct, F, cond);
	cout << "  Condition number (test diffusion problem): " << cond << endl;

	std::ofstream outUU("Ct.txt");
	outUU << gmm::col_vector(Ct);
	outUU.close();		

		
	vtk_export exp_Ct("test_Ct.vtk");
	exp_Ct.exporting(mf_Ct);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct, Ct, "Ct");
 	

	};//end of test
	
	
	  
	void transport3d1d::test2(void){

	vector_type U(dof_transp.Cv()); gmm::clear(U);
	sparse_matrix_type A(dof_transp.Cv(), dof_transp.Cv());gmm::clear(A);
	vector_type F(dof_transp.Cv());gmm::clear(F);
	vector_type Rg(dof_transp.Cv());	gmm::clear(Rg);
	
	getfem::generic_assembly
	  assem("M$1(#1,#1) += sym(comp(vGrad(#1).vGrad(#1)) (:, i,k, : ,i,k) )");
	  assem.push_mi(mimv);
	  assem.push_mf(mf_Cv);
	  assem.push_mat(A);
	  assem.assembly();
	  
	  gmm::MatrixMarket_IO::write("A.mm",A);
	  
	  vector_type ones(dof_transp.Cv(), 0.01);


	getfem::asm_source_term(F, mimv, mf_Cv, mf_Cv, ones);
	Rg[1]= 1;

	std::ofstream outFF("F.txt");
		outFF << gmm::col_vector(F);
		outFF.close();


		scalar_type cond;
		gmm::SuperLU_solve(A, U, F, cond);
		cout << "  Condition number (test diffusion problem): " << cond << endl;

	std::ofstream outUU("u_tilde.txt");
		outUU << gmm::col_vector(U);
		outUU.close();		

gmm::add(Rg, U);


	std::ofstream outU("U.txt");
		outU << gmm::col_vector(U);
		outU.close();
		
		
	vtk_export exp_U("test2.vtk");
	exp_U.exporting(mf_Cv);
	exp_U.write_mesh();
	exp_U.write_point_data(mf_Cv, U, "test2");
 	
 	
	


	};//end of test
 
 } // end of namespace
