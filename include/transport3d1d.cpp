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
 build_vessel_boundary();
 build_tissue_boundary();
 
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

	//but, in order to have the boundary conditions for the nodes
	//we need to build again the 1D mesh from another pts file
	mesht.clear();
		bool test = 0;
	test = PARAM.int_value("TEST_GEOMETRY");
	if(test==0){
		#ifdef M3D1D_VERBOSE_
		cout << "Importing the 3D mesh for the tissue ...  "   << endl;
		#endif
		 import_msh_file(descr.MESH_FILET, mesht);
	}else{
		#ifdef M3D1D_VERBOSE_
		cout << "Building the regular 3D mesh for the tissue ...  "   << endl;
		#endif
		string st("GT='" + PARAM.string_value("GT_T") + "'; " +
					   "NSUBDIV=" + PARAM.string_value("NSUBDIV_T") + "; " +  
					   "ORG=" + PARAM.string_value("ORG_T") + "; " +  
					   "SIZES=" + PARAM.string_value("SIZES_T") + "; " +  
					   "NOISED=" + PARAM.string_value("NOISED_T")); 
		cout << "mesht description: " << st << endl;
		regular_mesh(mesht, st);
	}
	
	#ifdef M3D1D_VERBOSE_
	cout << "Importing the 1D mesh for the vessel (transport problem)... "   << endl;
	#endif
	std::ifstream ifs(descr_transp.MESH_FILEV);
	meshv.clear();
//	mesh meshv_transp;
//	vector_size_type nb_vertices_transp;
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr_transp.MESH_FILEV);
	import_pts_file(ifs, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV);
	nb_branches = nb_vertices.size();
	ifs.close();
	
	/*
	cout<<"BC per il problema di trasporto"<<endl;
	for (size_type bc=0; bc < BCv_transp.size(); bc++) {
	cout<<"bc:"<<bc<<endl;
	cout <<BCv_transp[bc]<<endl; 
	}
	cout<<"BC per il problema di stokes"<<endl;
	for (size_type bc=0; bc < BCv.size(); bc++) {
	cout<<"bc:"<<bc<<endl;
	cout <<BCv[bc]<<endl;
	}
	*/
	
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
	

//HO DOVUTO CANCELLARE MESHV!!!!!!!!!!!!!!
// quindi metodi fem e di integrazione definiti sulla rete non valgono più! 
// li devo caricare di nuovo!
mimv.clear();
mf_Uvi.clear();
mf_Pv.clear();
mf_coefv.clear();

mimt.clear();
mf_Ut.clear();
mf_Pt.clear();
mf_coeft.clear();


problem3d1d::set_im_and_fem();
/*
#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs for tissue and vessel problems ..." << endl;
	#endif
	pintegration_method pim_v = int_method_descriptor(descr.IM_TYPEV);
	mimv.set_integration_method(meshv.convex_index(), pim_v);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	bgeot::pgeometric_trans pgt_v = bgeot::geometric_trans_descriptor(descr.MESH_TYPEV);
	pfem pf_Uv = fem_descriptor(descr.FEM_TYPEV);
	pfem pf_Pv = fem_descriptor(descr.FEM_TYPEV_P);
	pfem pf_coefv = fem_descriptor(descr.FEM_TYPEV_DATA);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif
	mf_Uvi.reserve(nb_branches);
	mf_coefvi.reserve(nb_branches);
	for(size_type i=0; i<nb_branches; ++i){
		
		mesh_fem mf_tmp(meshv);
		mf_tmp.set_finite_element(meshv.region(i).index(), pf_coefv);
		mf_coefvi.emplace_back(mf_tmp);
		mf_tmp.clear();
		
		mf_tmp.set_finite_element(meshv.region(i).index(), pf_Uv);
		mf_Uvi.emplace_back(mf_tmp);
		mf_tmp.clear();
	}
	mf_Pv.set_finite_element(meshv.convex_index(), pf_Pv);
	mf_coefv.set_finite_element(meshv.convex_index(), pf_coefv);
	
*/


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
		BCt_transp.emplace_back(labels[f], std::stof(values[f]), 0, f);
		#ifdef M3D1D_VERBOSE_
		cout << "  face " << f << " : " << BCt_transp.back() << endl;
		#endif
	} 
	
	for (size_type bc=0; bc < BCt_transp.size(); bc++)
	cout<<BCt_transp[bc]<<endl;
	
	// Build mesht regions
	mesh_region border_faces;
	outer_faces_of_mesh(mesht, border_faces);

	for (mr_visitor i(border_faces); !i.finished(); ++i) {

		assert(i.is_face());

		// Unit outward normal : used to identify faces
		//! \todo Use getfem 5.0's function select_faces_of_normal?
		base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		if (gmm::abs(un[0] + 1.0) < 1.0E-7)      // back
			mesht.region(0).add(i.cv(), i.f());
		else if (gmm::abs(un[0] - 1.0) < 1.0E-7) // front
			mesht.region(1).add(i.cv(), i.f());
		else if (gmm::abs(un[1] + 1.0) < 1.0E-7) // left
			mesht.region(2).add(i.cv(), i.f());
		else if (gmm::abs(un[1] - 1.0) < 1.0E-7) // right
			mesht.region(3).add(i.cv(), i.f());
		else if (gmm::abs(un[2] + 1.0) < 1.0E-7) // bottom
			mesht.region(4).add(i.cv(), i.f());
		else if (gmm::abs(un[2] - 1.0) < 1.0E-7) // top
			mesht.region(5).add(i.cv(), i.f());
				
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
			while (!found && (bc<BCv_transp.size())) {
				found = (i0 == BCv_transp[bc].idx);
				if (!found) bc++;
			}
			GMM_ASSERT1(found=true, "Miss a boundary node in BCv list!");
			BCv_transp[bc].rg = fer; 
			fer++;
			// Store the containing branch index
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			BCv_transp[bc].branches.emplace_back(branch); 
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
			while (!found && (bc<BCv_transp.size())) {
				found = (i1 == BCv_transp[bc].idx);
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
				BCv_transp[bc].value *= +1.0;
				BCv_transp[bc].rg = fer; 
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_transp[bc].branches.emplace_back(branch); 
			}
			else { /* interior -> Mixed point */
				// "MIX" label via post-processing
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				BCv_transp.emplace_back("MIX", 0.0, i1, fer);
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_transp.back().branches.emplace_back(branch); 
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
	for (size_type i=0; i<BCv_transp.size(); ++i)
		cout << "    -  label=" << BCv_transp[i].label 
			 << ", value=" << BCv_transp[i].value << ", ind=" << BCv_transp[i].idx 
			 << ", rg=" << BCv_transp[i].rg << ", branches=" << BCv_transp[i].branches << endl; 
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
	//assembly_rhs();
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
	
	//gmm::MatrixMarket_IO::write("Dt.mm",Dt);
	//gmm::MatrixMarket_IO::write("Mt.mm",Mt);
	
	
	// Copy Rt //MATRICE DI REAZIONE PER IL CONSUMO DI OSSIGENO, MOMENTANEAMENTE TOLTA
		bool REACTION = PARAM.int_value("REACTION", "flag for reaction term");
	if(REACTION ==1){
	gmm::add(Rt, 
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
				 	gmm::sub_interval(0, dof_transp.Ct()))); 
	}
	 
	// Copy Tt //PASSO TEMPORALE
	if(descr_transp.STATIONARY ==0)
	gmm::add(Mt,  
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 


	// Copy Dt //DIFFUSIONE
	bool DIFFUSION_T = PARAM.int_value("DIFFUSION_T", "flag for diffusion term in tissue");
	if(DIFFUSION_T ==1){	
	gmm::add(Dt,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 
			 
		}	
		#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mv and Dv ..." << endl;
	#endif
	
	// Build Mvvi and Dvvi
		asm_network_poiseuille_transp(Dv, Mv, mimv,mf_Cv, mf_coefv, param_transp.Av());
		gmm::scale(Mv, (1.0/param_transp.dt()));
		
		// mass matrix for TIME STEP for vessel
		if(descr_transp.STATIONARY ==0)
		gmm::add(Mv, 
			gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 
		
		// diffusion matrix for vessels
	bool DIFFUSION_V = PARAM.int_value("DIFFUSION_V", "flag for diffusion term in vessel");
	if(DIFFUSION_V ==1){	
			gmm::add(Dv, 
			gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	}			
	//gmm::MatrixMarket_IO::write("Dv.mm",Dv);
	//gmm::MatrixMarket_IO::write("Mv.mm",Mv);
		
		
		
	//ADVECTION	
		#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Bt and Bv ..." << endl;
	#endif		
	bool ADVECTION_T = PARAM.int_value("ADVECTION_T", "flag for advection term in tissue");
	if(ADVECTION_T ==1){
	// advection tissue term: 				
	vector_type Ut(dof.Ut()); gmm::clear(Ut);
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())) ,  Ut);
	asm_advection_tissue(Bt, mimt, mf_Ct, mf_Ut,Ut);
	gmm::add(gmm::scaled(Bt, 1.0),
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 
			
		
	gmm::MatrixMarket_IO::write("Bt.mm",Bt);
	}				
	/* assembla il termine di trasporto per le rete: 
	    sarà necessario combinare in qualche modo gli mf_Uvi con mf_Cv 
	    (magari costruendo anche per il trasporto un vettore mf_Cvi di fem_methods
	        
	*/    
	// advection vessel term: 		
	bool ADVECTION_V = PARAM.int_value("ADVECTION_V", "flag for advection term in vessels");
	if(ADVECTION_V ==1){
	
	//contatore per la posizione iniziale del ramo
	size_type shift = 0;
	
	//costruisci i versori tangenti lambda
	vector_type lambdax; // tangent versor: x component
	vector_type lambday; // tangent versor: y component
	vector_type lambdaz; // tangent versor: z component
	std::ifstream ifs(descr.MESH_FILEV);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr.MESH_FILEV);
	asm_tangent_versor(ifs, lambdax, lambday, lambdaz);
	ifs.close();
	
	//comincia ciclo su tutti i rami
	for (size_type i=0; i<nb_branches; ++i){
	
	//posizionati all'inizio del ramo
	if(i>0) shift += mf_Uvi[i-1].nb_dof();
	//ottieni nel vettore Uvi le velocità del ramo i
 	vector_type Uvi( mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);
	//ottieni i versori tangenti corrispondenti al ramo
	vector_type lambdax_K, lambday_K, lambdaz_K;
	for(size_type j=0; j<mf_coefvi[i].nb_dof(); j++){
			lambdax_K.emplace_back(lambdax[i]);
			lambday_K.emplace_back(lambday[i]);
			lambdaz_K.emplace_back(lambdaz[i]);
	}
	
	//costruisci la matrice del termine di trasporto
	asm_advection_network(Bv, mimv, mf_Cv, mf_coefvi[i], mf_Uvi[i], Uvi, lambdax_K, lambday_K, lambdaz_K, meshv.region(i) );
	}

	//somma la matrice -Bv alla matrice di rigidezza Am_transp
	gmm::add(gmm::scaled(Bv, 1.0),
			  gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	
	
	gmm::MatrixMarket_IO::write("Bv.mm",Bv);
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
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");
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

	//gmm::MatrixMarket_IO::write("Btt.mm",Btt);
	//gmm::MatrixMarket_IO::write("Btv.mm",Btv);
	//gmm::MatrixMarket_IO::write("Bvt.mm",Bvt);
	//gmm::MatrixMarket_IO::write("Bvv.mm",Bvv);
	
	
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



	
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term ..." << endl;
	#endif

		cout<<"assemble BC for tissue"<<endl;

		
		
		
		sparse_matrix_type Att(dof_transp.Ct(), dof_transp.Ct());
		vector_type Ft(dof_transp.Ct());
		gmm::add(	gmm::sub_matrix(AM_temp,
				gmm::sub_interval(0,dof_transp.Ct()),
				gmm::sub_interval(0,dof_transp.Ct()))
				, Att);
		gmm::scale(	gmm::sub_matrix(AM_temp,
				gmm::sub_interval(0,dof_transp.Ct()),
				gmm::sub_interval(0,dof_transp.Ct()))
				, 0.0);	
				
		gmm::add(	 gmm::sub_vector(FM_temp, gmm::sub_interval(0,dof_transp.Ct()))
				,Ft);	 
		gmm::scale(	 gmm::sub_vector(FM_temp, gmm::sub_interval(0,dof_transp.Ct()))
				,0.0);

cout<<"prima di entrare nella funzione asm_tissue_bc: queste sono le facce della mia mesh 3d: "<<endl;
	for (size_type bc=0; bc < BCt_transp.size(); ++bc) {
	cout<< BCt_transp[bc]<<endl;}
	
	scalar_type bcoef  = PARAM.real_value("BETA_transp", "Coefficient for mixed BC for transport problem");
	vector_type beta(dof.coeft(), bcoef);
	asm_tissue_bc_transp(Att, Ft, mimt, mf_Ct, mf_coeft, BCt_transp,beta);
		gmm::add(Att, 
			gmm::sub_matrix(AM_temp,
				gmm::sub_interval(0,dof_transp.Ct()),
				gmm::sub_interval(0,dof_transp.Ct())));
		gmm::add(Ft, 
			gmm::sub_vector(FM_temp,
				gmm::sub_interval(0,dof_transp.Ct())));
		// De-allocate memory
		gmm::clear(Att);
		gmm::clear(Ft);
		cout<<"assemble BC for vessel"<<endl;
	#ifdef M3D1D_VERBOSE_
	cout << "  Building vessel boundary term ..." << endl;
	#endif
		
		//costruisco le condizioni al contorno: influenzeranno AM_temp e Fv

		sparse_matrix_type Avv(dof_transp.Cv(), dof_transp.Cv());
		vector_type Fv(dof_transp.Cv());
		gmm::add(	gmm::sub_matrix(AM_temp,
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()))
				, Avv);
		gmm::scale(	gmm::sub_matrix(AM_temp,
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()))
				, 0.0);	
				
		gmm::add(	 gmm::sub_vector(FM_temp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))
				,Fv);	
		gmm::scale(	 gmm::sub_vector(FM_temp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))
				,0.0);

		asm_network_bc_transp(Avv, Fv, mimv, mf_Cv, mf_coefv, BCv_transp );
		gmm::add(Avv, 
			gmm::sub_matrix(AM_temp,
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())));
		gmm::add(Fv, 
			gmm::sub_vector(FM_temp,
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
		// De-allocate memory
		gmm::clear(Avv);
		gmm::clear(Fv);
		/*
		sparse_matrix_type Avv(dof_transp.Cv(), dof_transp.Cv());
		vector_type Fv(dof_transp.Cv());
		gmm::add(	gmm::sub_matrix(AM_transp,
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()))
				, Avv);
		gmm::scale(	gmm::sub_matrix(AM_transp,
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()))
				, 0.0);		

		asm_network_bc_transp(Avv, Fv, mimv, mf_Cv, mf_coefv, BCv_transp );
		gmm::add(Avv, 
			gmm::sub_matrix(AM_transp,
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
				gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())));
		gmm::add(Fv, 
			gmm::sub_vector(FM_transp,
				gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
		// De-allocate memory
		gmm::clear(Avv);
		gmm::clear(Fv);
		*/
		
		/*sparse_matrix_type Mvv(dof_transp.Cv(), dof_transp.Cv());
		vector_type Fv(dof_transp.Cv());
	

		asm_network_bc_transp(Mvv, Fv, mimv, mf_Cv, mf_coefv, BCv_transp );
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
		*/
		 
		/*
		
		asm_network_bc_transp(  gmm::sub_matrix(AM_transp,
						gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
						gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())),	
					gmm::sub_vector(FM_transp,
						gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())),	
					 mimv, mf_Cv, mf_coefv, BCv_transp);
*/
	
	
}


void transport3d1d::update (void){

// ho la matrice AM assemblata con tutti i termini (da non modificare mai più!)
// ho  il vettore F ancora vuoto (da aggiungere: il termine temporale; le condizioni al contorno 1d)
gmm::copy(AM_transp, AM_temp);
gmm::copy(FM_transp, FM_temp);

// faccio l'update temporale: sommo a F le concentrazioni al passo precedente

// update rhs (there is the time step mass term)
 
	vector_type TFt(dof_transp.Ct());
	vector_type TFv(dof_transp.Cv());
	//vector_type Ct(dof_transp.Ct());
	//vector_type Cv(dof_transp.Cv());
	//gmm::add(gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())), Ct);
	//gmm::add(gmm::sub_vector(UM_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);
	asm_source_term(TFt,mimt, mf_Ct, mf_Ct,gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct()))); // tentativo con generic_assembly: lo trovi in 										assembling3d_transp, la funzione è asm_time_rhs_transp
	asm_source_term(TFv,mimv, mf_Cv, mf_Cv,gmm::sub_vector(UM_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	gmm::scale(TFt, (1.0/param_transp.dt())); // dt time step
	gmm::scale(TFv, (1.0/param_transp.dt())); // dt time step
	gmm::add(TFt, gmm::sub_vector(FM_temp, gmm::sub_interval(0, dof_transp.Ct())));
	gmm::add(TFv, gmm::sub_vector(FM_temp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	gmm::clear(UM_transp);
	gmm::clear(TFt); gmm::clear(TFv);

assembly_rhs();


}
 bool transport3d1d::solve (void)
 {
  std::cout<<"solve transport problem"<<std::endl<<std::endl;
  #ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	gmm::resize(AM_temp, dof_transp.tot(), dof_transp.tot()); gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_transp.tot()); gmm::clear(FM_temp);
	
	double time = gmm::uclock_sec();
	double time_count = 0;	


	for(double t=0;t<=param_transp.T()*(!descr_transp.STATIONARY) ; t = t + param_transp.dt() + (param_transp.dt()==0) ){ 
	time_count++; 
	std::cout<<"iteration number:"<<time_count<<std::endl;
	std::cout<<"time = "<<t<<" s"<<std::endl;	
	
	update();
	
	
	gmm::csc_matrix<scalar_type> A_transp;
	gmm::clean(AM_transp, 1E-12);
	gmm::copy(AM_temp, A_transp);
	
	vector_type F_transp(gmm::vect_size(FM_transp));
	gmm::clean(FM_transp, 1E-12);
	gmm::copy(FM_temp, F_transp);
	
	//gmm::clear(AM_transp); // to be postponed for preconditioner
			
	//gmm::MatrixMarket_IO::write("Atransp.mm",A_transp);	
	
	
	
	
	
	

		
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
	
	
	/*

		 
	*/	
	
	} //end of cycle over time 
	
	cout << endl<<"... time to solve : " << gmm::uclock_sec() - time << " seconds\n";

	return true;
 }; // end of solve
	
	
 void transport3d1d::export_vtk (const string & time_suff,const string & suff)
 {
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

	//Copy solution
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
	
	gmm::MatrixMarket_IO::write("At.mm",At);
	gmm::MatrixMarket_IO::write("Btt.mm",Btt);
	gmm::MatrixMarket_IO::write("Btv.mm",Btv);
	gmm::MatrixMarket_IO::write("Bvt.mm",Bvt);
	gmm::MatrixMarket_IO::write("Bvv.mm",Bvv);	
	std::ofstream outF("F.txt");
	outF << gmm::col_vector(F);
	outF.close();	
	
	int a;
	cout<<"type 0 for superLU; 1 for GMRES"<<endl;
	cin>>a;
	
	if(a==0){
	cout<<"	superLU"<<endl;
	scalar_type cond;
	gmm::SuperLU_solve(At, Ct, F, cond);
	cout << "  Condition number (test diffusion problem): " << cond << endl;}
	else if(a==1){
	cout<<"	GMRES"<<endl;
	gmm::iteration iter (0.0000001, 1,40000);
	//gmm::ilutp_precond<sparse_matrix_type> P(At, 20, 1E-6);
	gmm::identity_matrix P;
	gmm::gmres (At, Ct, F, P, 50, iter);
	
	}



	std::ofstream outUU("Ct.txt");
	outUU << gmm::col_vector(Ct);
	outUU.close();		

	vtk_export exp_Ct("test_Ct.vtk");
	exp_Ct.exporting(mf_Ct);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct, Ct, "Ct");
 	



/*
proviamo a risolvere il problema di diffusione con i gradi di libertà e la mesh 3d definita per il problema dei fluidi.


	//! Monolithic matrix for the coupled problem
	sparse_matrix_type AM;
	//! Monolithic array of unknowns for the coupled problem
	vector_type        UM;
	//! Monolithic right hand side for the coupled problem
	vector_type        FM;

	gmm::resize(AM, dof.Pt(), dof.Pt()); gmm::clear(AM);
	gmm::resize(UM, dof.Pt()); gmm::clear(UM);
	gmm::resize(FM, dof.Pt()); gmm::clear(FM);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM ..." << endl;
	#endif
	// Mass matrix for the interstitial problem
	//getfem::asm_mass_matrix(AM, mimt, mf_Pt);
	  /*getfem::generic_assembly
	  assem("M$1(#1,#1) += sym(comp(vGrad(#1).vGrad(#1)) (:, i,k, : ,i,k) )");
	  assem.push_mi(mimt);
	  assem.push_mf(mf_Pt);
	  assem.push_mat(AM);
	  assem.assembly();
	 getfem::asm_stiffness_matrix_for_homogeneous_laplacian(AM, mimt, mf_Pt);
	  
	vector_type F(dof.coeft(), 1.0);
	getfem::asm_source_term(FM, mimt, mf_Pt, mf_coeft, F);
	gmm::MatrixMarket_IO::write("A_mass.mm",AM);
	std::ofstream outFF("Ct.txt");
	outFF << gmm::col_vector(FM);
	outFF.close();	
	
	gmm::csc_matrix<scalar_type> A;
	gmm::clean(AM, 1E-12);
	gmm::copy(AM, A);
	gmm::clear(AM); // to be postponed for preconditioner
	double time = gmm::uclock_sec();	
	
	// direct solver //
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(A, UM, FM, cond);
		cout << "  Condition number : " << cond << endl;

	std::ofstream outU("U_mass.txt");
	outU << gmm::col_vector(UM);
	outU.close();	
	
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
 	
 */	
	


	};//end of test
 
 } // end of namespace
