/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2015-2016
                  
                Copyright (C) 2016 Domenico Notaro
======================================================================*/
/*! 
  @file   utilities.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   November 2017.
  @brief  Declaration of some miscelleanous auxiliary routines.
 */
#ifndef M3D1D_UTILITIES_TRANSP_HPP_
#define M3D1D_UTILITIES_TRANSP_HPP_

#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_mesh.h>
#include <gmm/gmm.h>
#include <defines.hpp>
#include  <transport3d1d.hpp>

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
}/* end of max_vec*/




} /* end of getfem namespace */



#endif
