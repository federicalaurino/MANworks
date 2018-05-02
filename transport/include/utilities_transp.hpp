/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2015-2016
                  
                Copyright (C) 2016 Domenico Notaro
======================================================================*/
/*! 
  @file   utilities_transp.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   November 2017.
  @brief  Declaration of some miscelleanous auxiliary routines.
 */
#ifndef M3D1D_UTILITIES_TRANSP_HPP_
#define M3D1D_UTILITIES_TRANSP_HPP_

#include <gmm/gmm.h>
#include <defines.hpp>

namespace gmm {

//Calculate the product between two vector, component by component.
//That means: C=A**B   --> C[i]= A[i]* B[i]

//As a matter of fact, it returns: B= A**B
template
<typename VEC>
void scale (const VEC & A, VEC & B){

//calculate the length of the two vectors

int lengthA = gmm::mat_nrows(gmm::col_vector(A));
int lengthB = gmm::mat_nrows(gmm::col_vector(B));
GMM_ASSERT1(lengthA==lengthB, "impossible to scale the vectors: different lengths");

for(int i=0; i<lengthA; i++)
  B[i]= A[i]*B[i];
  
} /* end of scale */

} /* end of gmm namespace */


namespace getfem {

//Calculate the max value of a scalar or a vector field contained in a vector V.
//Aux function for computing Peclet
template
<typename VEC>
scalar_type max_vec (const VEC & V, size_type dim){


GMM_ASSERT1(dim>0, "wrong dimension of field: must be greater or equal to 1");

scalar_type max=0;
size_type size = gmm::vect_size(V);
if(dim>1){  // compute the norm of the vector field, then find the maximum
GMM_ASSERT1((size%dim)==0, "ERROR: the dimension of the field and the size of the vector do not match!");
VEC Vtemp(size/dim, 0);
for(int i=0; i<size/dim; i++){
   for( int j=0; j<dim; j++){
   Vtemp[i]+=(V[i+j])*(V[i+j]);
    }
   Vtemp[i]= sqrt(Vtemp[i]);
   if(Vtemp[i]>max) max= Vtemp[i]; //update maximum
  }
}
else if (dim==1){ //just find the maximum
for(int i=0; i<size; i++){
if(std::abs(V[i])>max) max= std::abs(V[i]);
}

} 
return max;}/* end of max_vec*/


//! Compute the Peclèt Number, defined as: 
//! @f$ Pe =  \frac{U~h}{A} @f$
//! that is the ratio between the advection flux and the diffusion coefficient
template
<typename VEC>
scalar_type peclet(const mesh & mesh,const VEC & U, const scalar_type & A, size_type  dim){


	#ifdef M3D1D_VERBOSE_
	cout <<"computing Peclet number...  "<<endl;
	#endif

// find Umax
scalar_type Umax= max_vec (U, dim);

// find h
scalar_type h=0;
scalar_type temp=0;
for(dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
	if(dim==1)
		temp= estimate_h(mesh, i);
	else if(dim==3)
		temp= 2*mesh.convex_radius_estimate(i);
	if(temp>h) h=temp;
	}

// compute peclet
scalar_type peclet= Umax*h/A;

	#ifdef M3D1D_VERBOSE_
	cout <<"U:   "<<Umax<<endl;
	cout <<"h:   "<<h<<endl;
	cout <<"A:   "<<A<<endl;
	cout <<"Peclet:   "<<peclet<<endl;
	#endif

return peclet;

} /* end of peclet_vessel */

  /**
     Faster (and simpler) assembly of simple Dirichlet conditions (
     u(x) = F(x) on a boundary). 

     @param mf should be Lagrangian.
     @param boundary the boundary number.
     @param DIR the dirichlet condition value.
     @param B,F are modified to enforce the Dirichlet condition. The
     symmetry properties of RM are kept.

     @ingroup asm
  */
 template<typename MATRM, typename VECT1, typename VECT2>
  void assembling_Dirichlet_condition_coupled_tissue
  (MATRM &B, VECT1 &F, const mesh_fem &mf1, const mesh_fem &mf2, size_type boundary,
   const VECT2 &DIR) {
   

    size_type Q1=mf1.get_qdim();
    size_type Q2=mf2.get_qdim();

    size_type nb_dof1=mf1.nb_dof();
    size_type nb_dof2=mf2.nb_dof();    
    
    GMM_ASSERT1(!(mf1.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
    GMM_ASSERT1(!(mf2.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
		
	std::cout<<"inizio uncoupled"<<std::endl;				
    dal::bit_vector nndof = mf1.basic_dof_on_region(boundary);
    pfem pf1;
    
    for (dal::bv_visitor cv(mf1.convex_index()); !cv.finished(); ++cv) {	 	//per tutti i convessi cv della mesh 1
    
      pf1 = mf1.fem_of_element(cv);
      pdof_description ldof = lagrange_dof(pf1->dim());
      size_type nbd = pf1->nb_dof(cv);	      
      for (size_type i = 0; i < nbd; i++) {					//per tutti i dof i del convesso cv
	size_type dof1 = mf1.ind_basic_dof_of_element(cv)[i*Q1];				//trova l'indice delle colonne riferite all
	if (nndof.is_in(dof1) && pf1->dof_types()[i] == ldof) {			//se il dof i del convesso cv è in "boundary"
  
	  for (size_type j = nb_dof1; j < nb_dof1+ nb_dof2; j++) {				//allora per tutti i dof j della mesh 2
		for (size_type l = 0; l < Q1; ++l) {
			F[j] -= B(j, dof1+l) * DIR[dof1+l];
	    		B(j, dof1+l) =  0;
	    		B(dof1+l, j) =  0;
	    	}
	    }
	  } 
	}
     }
   } /* end of assembling_Dirichlet_condition_coupled_tissue*/

  /**
     Faster (and simpler) assembly of simple Dirichlet conditions (
     u(x) = F(x) on a boundary). 

     @param mf should be Lagrangian.
     @param boundary the boundary number.
     @param DIR the dirichlet condition value.
     @param B,F are modified to enforce the Dirichlet condition. The
     symmetry properties of RM are kept.

     @ingroup asm
  */
 template<typename MATRM, typename VECT1, typename VECT2>
  void assembling_Dirichlet_condition_coupled_vessel
  (MATRM &B, VECT1 &F, const mesh_fem &mf1, const mesh_fem &mf2, size_type boundary,
   const VECT2 &DIR) {
   

    size_type Q1=mf1.get_qdim();
    size_type Q2=mf2.get_qdim();

    size_type nb_dof1=mf1.nb_dof();
    size_type nb_dof2=mf2.nb_dof();    
    
    GMM_ASSERT1(!(mf1.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
    GMM_ASSERT1(!(mf2.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
		
	std::cout<<"inizio uncoupled"<<std::endl;				
    dal::bit_vector nndof = mf2.basic_dof_on_region(boundary);
    pfem pf2;
    
    for (dal::bv_visitor cv(mf2.convex_index()); !cv.finished(); ++cv) {	 	//per tutti i convessi cv della mesh 1
    
      pf2 = mf2.fem_of_element(cv);
      pdof_description ldof = lagrange_dof(pf2->dim());
      size_type nbd = pf2->nb_dof(cv);	      
      for (size_type i = 0; i < nbd; i++) {					//per tutti i dof i del convesso cv
	size_type dof2 = mf2.ind_basic_dof_of_element(cv)[i*Q2];				//trova l'indice delle colonne riferite all
	if (nndof.is_in(dof2) && pf2->dof_types()[i] == ldof) {			//se il dof i del convesso cv è in "boundary"
	  
	  for (size_type j = 0; j < nb_dof1; j++) {				//allora per tutti i dof j della mesh 2
		for (size_type l = 0; l < Q2; ++l) {
			F[j] -= B(j, nb_dof1 + dof2+l) * DIR[nb_dof1 + dof2+l];
	    		B(j, nb_dof1 + dof2+l) =  0;
	    		B(nb_dof1 + dof2+l, j) =  0;
	    	}
	    }
	  } 
	}
     }
   } /* end of assembling_Dirichlet_condition_coupled_vessel*/
   
  /*! Build the mixed boundary conditions (weak form) and dirichlet (strong form) for vessels
    @f$ M=\int_{\Gamma_{MIX}} \beta~u~v~d\sigma@f$ and
    @f$ F=\int_{\Gamma_{MIX}} \beta~c0~v~d\sigma@f$
 */
/*!
	@param F        BC contribution to rhs
	@param M        BC contribution to mass matrix
	@param mim      The integration method to be used
	@param mf_c     The finite element method for the concentration @f$u@f$
	@param mf_data  The finite element method for the coefficients
	@param BC       Array of values of network boundary points
	@param beta     The beta value for mix condition @f$p_0@f$
	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void
asm_coupled_bc_transp
	(VEC & F,
	 MAT & M,
	 const mesh_fem & mf_ct,
	 const mesh_fem & mf_cv,
	 const mesh_fem & mf_data_t,
	 const mesh_fem & mf_data_v,
	 const std::vector<getfem::node> & BC_tissue,
	 const std::vector<getfem::node> & BC_vessel
	 )
{

	
	
	GMM_ASSERT1(mf_ct.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");
	GMM_ASSERT1(mf_cv.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");
	GMM_ASSERT1(mf_data_t.get_qdim()==1, "invalid data mesh fem (Qdim=1 required)");
	GMM_ASSERT1(mf_data_v.get_qdim()==1, "invalid data mesh fem (Qdim=1 required)");


//cycle over the tissue boundary nodes
	for (size_type bc=0; bc < BC_tissue.size(); ++bc) {
		GMM_ASSERT1(mf_ct.linked_mesh().has_region(bc), "missed mesh region" << bc);
		if (BC_tissue[bc].label=="DIR") { // Dirichlet BC
			VEC BC_temp(mf_ct.nb_dof(), BC_tissue[bc].value);
			getfem::assembling_Dirichlet_condition_coupled_tissue(M, F, mf_ct, mf_cv, BC_tissue[bc].rg, BC_temp);
			gmm::clear(BC_temp);				
		} 
	}
	
//cycle over the vessels boundary nodes
	for (size_type bc=0; bc < BC_vessel.size(); ++bc) {
		GMM_ASSERT1(mf_cv.linked_mesh().has_region(bc), "missed mesh region" << bc);
		if (BC_vessel[bc].label=="DIR") { // Dirichlet BC
			VEC BC_temp(mf_cv.nb_dof(), BC_vessel[bc].value);
			getfem::assembling_Dirichlet_condition_coupled_vessel(M, F, mf_ct, mf_cv, BC_vessel[bc].rg, BC_temp);
			gmm::clear(BC_temp);				
		} 
	}


} /* end of asm_coupled_bc_transp */

} /* end of getfem namespace */



#endif
