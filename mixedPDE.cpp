#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <assert.h>
#include "mixedPDE.h"
#include "fcns.h"
#include "Mesh/baseObjects.h"
#include "Mesh/polyMesh.h"
#include "directMixed.h"
#include "parameterData.h"
#include "polyQuadrature.h"
#include "Utilities/monitor.h"
#include "Utilities/debug.h"
#include <complex.h>
#include "lapacke.h"
#include <cblas.h>

using namespace directserendipity;
using namespace polymesh;
using namespace polyquadrature;

lapack_int mat_inv(double *A, unsigned n)
{
  // inplace inverse n x n matrix A.
  // matrix A is Column Major (i.e. firts line, second line ... *not* C[][] order)
  // returns:
  //   ret = 0 on success
  //   ret < 0 illegal argument value
  //   ret > 0 singular matrix
  int ipiv[n+1];
  lapack_int ret;

  ret =  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n,n,A,n,ipiv);

  if (ret !=0) return ret;
  ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR,n,A,n,ipiv);
  return ret;
}

int MixedPDE::solve(Monitor& monitor) {
  monitor(0,"Polynomial Degree = ", parameterDataPtr()->dsSpace.degPolyn());

  ParameterData& param = *parameterDataPtr();

  // TEST SINGLE ELEMENT MESH //////////////////////////////////////////////

  if(false) {
    monitor(0,"Test Single Element Mesh and Space");

    polymesh::PolyMesh s_mesh(parameterDataPtr()->mesh.elementPtr(3));
    std::string fileName = parameterDataPtr()->directory_name;
    fileName += "mesh_e4";
    s_mesh.write_matlab(fileName);

    DirectMixed s_dmSpace(8,&s_mesh);
    fileName = parameterDataPtr()->directory_name;
    fileName += "dmSpace_e4";
    s_dmSpace.write_matlab(fileName);
  }
  
  // TEST BASIS FUNCTIONS //////////////////////////////////////////////////

  if(true) {
    monitor(0,"\nTest basis functions for element 0\n");

    DirectMixedArray u(&(parameterDataPtr()->dmSpace), 'f');
    DirectDGArray p(&(parameterDataPtr()->dmSpace), 'f');
    DirectEdgeDGArray l(&(parameterDataPtr()->dmSpace));

    for(int i=0; i<u.size(); i++) {
      u[i]=0;
    }
    u[0]=1;

    for(int i=0; i<p.size(); i++) {
      p[i]=0;
    }
    p[0]=1;

    for(int i=0; i<l.size(); i++) {
      l[i]=0;
    }
    l[0]=1;

    monitor(1,"Write Array");

    std::string fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh_mixed";
    u.write_matlab_mesh(fileName,51,51);

    fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh_DG";
    p.write_matlab_mesh(fileName,51,51);

    fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh_EDG";
    l.write_matlab_mesh(fileName,51,51);
  }

  // TEST QUADRATURE ///////////////////////////////////////////////////////

  if(false) {
    monitor(1,"Test Quadrature Rules");
    testPolyQuadrature(&(parameterDataPtr()->mesh),1e-6);
    return 0;
  }
  
  // SOLVE THE PDE ///////////////////////////////////////////////////////////

  monitor(0,"\nSolve the PDE\n");

  DirectMixedArray solution_u_f(&(param.dmSpace),'f');
  DirectMixedArray solution_u_r(&(param.dmSpace),'r');
  DirectDGArray solution_p_f(&(param.dmSpace),'f');
  DirectDGArray solution_p_r(&(param.dmSpace),'r');
  DirectEdgeDGArray solution_l(&(param.dmSpace));

  // Initialize matrix A for both full and reduced space
  int dimAfull = param.dmSpace.nMixedDoFs('f');
  int dimAreduced = param.dmSpace.nMixedDoFs('r');

  std::vector<double> mat_A_full_vector(dimAfull*dimAfull,0);
  double* mat_A_full = mat_A_full_vector.data();
  std::vector<double> mat_A_reduced_vector(dimAreduced*dimAreduced,0);
  double* mat_A_reduced = mat_A_reduced_vector.data();

  // Initialize matrix B for both full and reduced space
  int rowBfull = dimAfull, rowBreduced = dimAreduced;
  int colBfull = param.dmSpace.nDGDoFs('f');
  int colBreduced = param.dmSpace.nDGDoFs('r');

  std::vector<double> mat_B_full_vector(rowBfull*colBfull,0);
  double* mat_B_full = mat_B_full_vector.data();
  std::vector<double> mat_B_reduced_vector(rowBreduced*colBreduced,0);
  double* mat_B_reduced = mat_B_reduced_vector.data();

  // Initialize matrix L
  int rowLfull = dimAfull, rowLreduced = dimAreduced;
  int colL = param.dmSpace.nIntEdgeDGDoFs();

  std::vector<double> mat_L_full_vector(rowLfull*colL,0);
  double* mat_L_full = mat_L_full_vector.data();
  std::vector<double> mat_L_reduced_vector(rowLreduced*colL,0);
  double* mat_L_reduced = mat_L_reduced_vector.data();

  // Initialize right hand side (W_s coefficients only)
  std::vector<double> rhs_full_vector(colBfull,0);
  double* rhs_full = rhs_full_vector.data();

  std::vector<double> rhs_reduced_vector(colBreduced,0);
  double* rhs_reduced = rhs_reduced_vector.data();


  // quadrature points
  // Note that we update quadRule in each iElement loop
  // But we define quadEdgeRule for all edges at once
  // and get them by interior edge indexing

  polyquadrature::PolyQuadrature quadRule(2*param.dmSpace.degPolyn()+3);
  std::vector<polyquadrature::PolyEdgeQuadrature> quadEdgeRule(param.dmSpace.nInteriorEdges());
  for (int i = 0; i < param.dmSpace.nInteriorEdges(); i++) {
    quadEdgeRule[i].set(param.dmSpace.degPolyn()+1, param.dmSpace.edgeInteriorPtr(i));
  }

  monitor(1,"Matrix and RHS Assembly"); ////////////////////////////////////////

  int starting_Afull = 0, starting_Areduced = 0;
  int starting_colBfull = 0, starting_colBreduced = 0;

  int loc_dimAfull = 0, loc_dimAreduced = 0, loc_colBfull = 0, loc_colBreduced = 0;

  int curr_full_index, curr_reduced_index;
  double evaluation;
  int numEdges = 0, loc_to_int = 0;

  // We construct an array of pointer to all the DirectEdgeDGFE
  // and get them by interior edge indexing

  std::vector<DirectEdgeDGFE*> eePtr(param.dmSpace.nInteriorEdges());
  for (int i = 0; i < param.dmSpace.nInteriorEdges(); i++) {
    eePtr[i] = param.dmSpace.DGEdgeInteriorPtr(i);
    eePtr[i] -> initBasis(quadEdgeRule[i].pts(), quadEdgeRule[i].num());
  }


  for(int iElement=0; iElement<param.mesh.nElements(); iElement++) {
    DirectMixedFE* mePtr = param.dmSpace.MixedElementPtr(iElement);
    DirectDGFE* dgePtr = param.dmSpace.DGElementPtr(iElement);

    quadRule.setElement(mePtr->elementPtr());
    
    mePtr->initBasis(quadRule.pts(), quadRule.num());
    dgePtr->initBasis(quadRule.pts(), quadRule.num());

    loc_dimAfull = mePtr -> dimVFull();
    loc_dimAreduced = mePtr -> dimVReduced();
    loc_colBfull = dgePtr -> dimFull();
    loc_colBreduced = dgePtr -> dimReduced();

    // Matrix A, B and rhs assembly over elements
 
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt)[0];
      double y = quadRule.pt(iPt)[1];

      Tensor2 valD;
      coefD_inv(x,y,valD);
  
      // Assemble matrix A
      for (int j = 0; j < loc_dimAfull; j++) {
        Tensor1 v_j = mePtr -> basis(j,iPt);
        for (int i = 0; i < loc_dimAfull; i++) {
          Tensor1 u_i = mePtr -> basis(i,iPt);
          curr_full_index = dimAfull * (starting_Afull + j) + (starting_Afull + i);
          evaluation = ((valD*u_i)*v_j)*quadRule.wt(iPt);
          mat_A_full[curr_full_index] += evaluation;
//          assert(curr_full_index < dimAfull * dimAfull);
          if (j < loc_dimAreduced && i < loc_dimAreduced) {
            curr_reduced_index = dimAreduced * (starting_Areduced + j) + (starting_Areduced + i);
            mat_A_reduced[curr_reduced_index] += evaluation;
//            assert(curr_reduced_index < dimAreduced * dimAreduced);
          }
        }
      }

      // Assemble matrix B
      // For first j = 0 -> dimCurlPart()-1 rows, div(v_j) = 0, 
      // so we only need to consider divXPo part

      for (int j = mePtr -> dimCurlPart(); j < loc_dimAfull; j++) {
        double divv_j = mePtr -> basisdivXPo(j - mePtr -> dimCurlPart(),iPt);
        for (int i = 0; i < loc_colBfull; i++ ) {
          double p_i = dgePtr -> basis(i,iPt);
          curr_full_index = colBfull * (starting_Afull + j) + (starting_colBfull + i);
          evaluation = divv_j * p_i * quadRule.wt(iPt);
          mat_B_full[curr_full_index] += evaluation;
//          assert(curr_full_index < dimAfull * colBfull);
          if (j < loc_dimAreduced && i < loc_colBreduced) {
            curr_reduced_index = colBreduced * (starting_Areduced + j) + (starting_colBreduced + i);
            mat_B_reduced[curr_reduced_index] += evaluation;
//            assert(curr_reduced_index < dimAreduced * colBreduced);
          }
        }
      }

      // Assemble RHS
      double f_wted = sourceVal(quadRule.pt(iPt)[0],quadRule.pt(iPt)[1]) * quadRule.wt(iPt);
      for (int j = 0; j < loc_colBfull; j++) {
        curr_full_index = starting_colBfull + j;
        evaluation = f_wted * dgePtr->basis(j,iPt);
        rhs_full[curr_full_index] += evaluation;
//        assert(curr_full_index < colBfull);
        if (j < loc_colBreduced) {
          curr_reduced_index = starting_colBreduced + j;
          rhs_reduced[curr_reduced_index] += evaluation;
//          assert(curr_reduced_index < colBreduced);
        }
      }
    }

    // Matrix L assembly over elements

     numEdges = param.mesh.elementPtr(iElement) -> nVertices();

      // Column indexing of L: (global edge indexed by interior)
      // {edge(0),Func(0)}, {edge(0),Func(1)}, ..., {edge(0),Func(eePtr[0]->dim()-1)},
      // {edge(1),Func(0)}, {edge(1),Func(1)}, ..., {edge(1),Func(eePtr[1]->dim()-1)},
      // ..., {edge(nInteriorEdge()-1),Func(eePtr[nInteriorEdge()-1]->dim()-1)} 

      for (int iEdge = 0; iEdge < numEdges; iEdge ++) {
        loc_to_int = param.dmSpace.interiorEdgeIndex(iElement,iEdge);
        if (loc_to_int == -1) continue; // If the edge is on boundary, we skip the loop

        mePtr->initBasis(quadEdgeRule[loc_to_int].pts(), quadEdgeRule[loc_to_int].num());

        for (int iPt = 0; iPt < quadEdgeRule[loc_to_int].num(); iPt++) {
          for (int i = 0; i < eePtr[loc_to_int]->dim(); i++){
            double l_i = eePtr[loc_to_int]->basis(i,iPt);
            for (int j = 0; j < loc_dimAfull; j++) {
              double v_jdotNu = mePtr->basisDotNu(j,iEdge,iPt);
              // Here we use the property that eePtr[loc_to_int]->dim() is the same (degPolyn()+1)for every edge in our mesh
              curr_full_index = colL * (starting_Afull + j) + (loc_to_int*(param.dmSpace.degPolyn()+1)+i);
              evaluation = l_i * v_jdotNu * quadEdgeRule[loc_to_int].wt(iPt);
              mat_L_full[curr_full_index] += evaluation;
//              assert(curr_full_index < dimAfull*colL);
              if (j < loc_dimAreduced) {
                curr_reduced_index = colL * (starting_Areduced + j) + (loc_to_int*(param.dmSpace.degPolyn()+1)+i);
                mat_L_reduced[curr_reduced_index] += evaluation;
//                assert(curr_reduced_index < dimAreduced*colL);
              }
            } 
          }
        }
      
      }

    starting_Afull += loc_dimAfull;
    starting_Areduced += loc_dimAreduced;
    starting_colBreduced += loc_colBreduced;
    starting_colBfull += loc_colBfull;
  }

  // OUTPUT MATRICES AND RHS //////////////////////////////////////////////////

  if (true) {
  
  std::ofstream fout1("test/A_full.txt");
  for(int j = 0; j < dimAfull; j++) {
    for(int i = 0; i < dimAfull; i++) {
      fout1 << mat_A_full[i + dimAfull*j] << "\t";
    }
    if (j < dimAfull - 1) fout1 << "\n";
  }

  std::ofstream fout2("test/A_reduced.txt");
  for(int j = 0; j < dimAreduced; j++) {
    for(int i = 0; i < dimAreduced; i++) {
      fout2 << mat_A_reduced[i + dimAreduced*j] << "\t";
    }
    if (j < dimAreduced - 1) fout2 << "\n";
  }

  std::ofstream fout3("test/B_full.txt");
  for(int j = 0; j < rowBfull; j++) {
    for(int i = 0; i < colBfull; i++) {
      fout3 << mat_B_full[i + colBfull*j] << "\t";
    }
    if (j < rowBfull - 1) fout3 << "\n";
  }

  std::ofstream fout4("test/B_reduced.txt");
  for(int j = 0; j < rowBreduced; j++) {
    for(int i = 0; i < colBreduced; i++) {
      fout4 << mat_B_reduced[i + colBreduced*j] << "\t";
    }
    if (j < rowBreduced - 1) fout4 << "\n";
  }

  std::ofstream fout5("test/L_full.txt");
  for(int j = 0; j < rowLfull; j++) {
    for(int i = 0; i < colL; i++) {
      fout5 << mat_L_full[i + colL*j] << "\t";
    }
    if (j < rowLfull - 1) fout5 << "\n";
  }

  std::ofstream fout6("test/L_reduced.txt");
  for(int j = 0; j < rowLreduced; j++) {
    for(int i = 0; i < colL; i++) {
      fout6 << mat_L_reduced[i + colL*j] << "\t";
    }
    if (j < rowLreduced - 1) fout6 << "\n";
  }

  std::ofstream fout7("test/rhs_full.txt");
  for(int i = 0; i < colBfull; i++) {
    fout7 << rhs_full[i];
    if (i < colBfull - 1) fout7 << "\n";
  }

  std::ofstream fout8("test/rhs_reduced.txt");
  for(int i = 0; i < colBreduced; i++) {
    fout8 << rhs_reduced[i];
    if (i < colBreduced - 1) fout8 << "\n";
  }

  }

  // SOLVE LINEAR SYSTEM //////////////////////////////////////////////////

  monitor(1,"Solution of linear system"); ////////////////////////////////////////

  // Solve the matrix
  // B^{T} A^{-1} [B - L * (L^{T} A^{-1} L)^{-1} * (L^{T} A^{-1} B)] b = rhs


  // Initialize A^{-1} and evaluate as A

  std::vector<double> A_full_inv_vector(dimAfull*dimAfull,0);
  double* A_full_inv = A_full_inv_vector.data();
  std::vector<double> A_reduced_inv_vector(dimAreduced*dimAreduced,0);
  double* A_reduced_inv = A_reduced_inv_vector.data();

  for (int j = 0; j < dimAfull; j++) {
    for (int i = 0; i < dimAfull; i++) {
      A_full_inv[j * dimAfull + i] = mat_A_full[j * dimAfull + i];
      if (i < dimAreduced && j < dimAreduced) {
        A_reduced_inv[j * dimAreduced + i] = mat_A_reduced[j * dimAreduced + i];
      }
    }
  }

  // Update A^{-1} to be the inverse matrix A^{-1} 

  lapack_int ierr = mat_inv(A_full_inv, dimAfull);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }

  ierr = mat_inv(A_reduced_inv, dimAreduced);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }

  std::ofstream fout15("test/A_full_inv.txt");
  for(int j = 0; j < dimAfull; j++) {
    for(int i = 0; i < dimAfull; i++) {
      fout15 << A_full_inv[i + dimAfull*j] << "\t";
    }
    if (j < dimAfull - 1) fout15 << "\n";
  }


  // Calculate L^{T}A^{-1}
  std::vector<double> ltai_full_vector(colL*dimAfull,0);
  double* ltai_full = ltai_full_vector.data();
  std::vector<double> ltai_reduced_vector(colL*dimAreduced,0);
  double* ltai_reduced = ltai_reduced_vector.data();


  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colL, dimAfull,
              dimAfull, 1, mat_L_full, colL, A_full_inv, dimAfull,
              0, ltai_full, dimAfull);

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colL, dimAreduced,
              dimAreduced, 1, mat_L_reduced, colL, A_reduced_inv, dimAreduced,
              0, ltai_reduced, dimAreduced);


  std::ofstream fout11("test/ltai_full.txt");
  for(int j = 0; j < colL; j++) {
    for(int i = 0; i < dimAfull; i++) {
      fout11 << ltai_full[i + dimAfull*j] << "\t";
    }
    if (j < colL - 1) fout11 << "\n";
  }



  // Caluculate (L^{T} A^{-1} L)^{-1} and (L^{T} A^{-1} B)
  // Multiply them together to be mat_b_c
  // Notice that mat_b_c * b = c

  std::vector<double> ltaili_full_vector(colL*colL,0);
  double* ltaili_full = ltaili_full_vector.data();

  std::vector<double> ltaili_reduced_vector(colL*colL,0);
  double* ltaili_reduced = ltaili_reduced_vector.data();

  std::vector<double> mat_b_c_full_vector(colL*colBfull,0);
  double* mat_b_c_full = mat_b_c_full_vector.data();
  std::vector<double> mat_b_c_reduced_vector(colL*colBreduced,0);
  double* mat_b_c_reduced = mat_b_c_reduced_vector.data();

  // Calculate L^{T} A^{-1} * L and store in ltaili

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, colL,
              dimAfull, 1, ltai_full, dimAfull, mat_L_full, colL,
              0, ltaili_full, colL);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, colL,
              dimAreduced, 1, ltai_reduced, dimAfull, mat_L_reduced, colL,
              0, ltaili_reduced, colL);

  // Calculate inverse of L^{T} A^{-1} * L and store in ltaili

  ierr = mat_inv(ltaili_full, colL);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }

  ierr = mat_inv(ltaili_reduced, colL);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }

  std::ofstream fout12("test/ltaili_full.txt");
  for(int j = 0; j < colL; j++) {
    for(int i = 0; i < colL; i++) {
      fout12 << ltaili_full[i + colL*j] << "\t";
    }
    if (j < colL - 1) fout12 << "\n";
  }


  // Calculate L^{T} A^{-1} * B and store in mat_b_c

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, colBfull,
              dimAfull, 1, ltai_full, dimAfull, mat_B_full, colBfull,
              0, mat_b_c_full, colBfull);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, colBreduced,
              dimAreduced, 1, ltai_reduced, dimAfull, mat_B_reduced, colBreduced,
              0, mat_b_c_reduced, colBreduced);





  // Multiply (L^{T} A^{-1} L)^{-1} and (L^{T} A^{-1} B)
  // Store in mat_b_c

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, colBfull,
              colL, 1, ltaili_full, colL, mat_b_c_full, colBfull,
              0, mat_b_c_full, colBfull);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, colBreduced,
              colL, 1, ltaili_reduced, colL, mat_b_c_reduced, colBreduced,
              0, mat_b_c_reduced, colBreduced);


  std::ofstream fout13("test/mat_b_c_full.txt");
  for(int j = 0; j < colL; j++) {
    for(int i = 0; i < colBfull; i++) {
      fout13 << ltaili_full[i + colBfull*j] << "\t";
    }
    if (j < colL - 1) fout13 << "\n";
  }


  // -L * mat_b_c + B -> workmat

  std::vector<double> workmat_full_vector(dimAfull*colBfull,0);
  double* workmat_full = workmat_full_vector.data();
  std::vector<double> workmat_reduced_vector(dimAreduced*colBreduced,0);
  double* workmat_reduced = workmat_reduced_vector.data();

  // We first update workmat to be B

  for (int j = 0; j < dimAfull; j++) {
    for (int i = 0; i < colBfull; i++) {
      workmat_full[j * colBfull + i] = mat_B_full[j * colBfull + i];
      if (j < dimAreduced && i < colBreduced) {
        workmat_reduced[j * colBreduced + i] = mat_B_reduced[j * colBreduced + i];
      }
    }
  }

  // Calculate -L * mat_b_c + B

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, colBfull,
              colL, -1, mat_L_full, colL, mat_b_c_full, colBfull,
              1, workmat_full, colBfull);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, colBreduced,
              colL, -1, mat_L_reduced, colL, mat_b_c_reduced, colBreduced,
              1, workmat_reduced, colBreduced);

  std::ofstream fout14("test/workmat_full.txt");
  for(int j = 0; j < dimAfull; j++) {
    for(int i = 0; i < colBfull; i++) {
      fout14 << workmat_full[i + colBfull*j] << "\t";
    }
    if (j < dimAfull - 1) fout14 << "\n";
  }



  // Do multiplication A^{-1} * workmat and update workmat as this
  // Then workmat = A^{-1} [B - L * (L^{T} A^{-1} L)^{-1} * (L^{T} A^{-1} B)]

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, colBfull,
              dimAfull, 1, A_full_inv, dimAfull, workmat_full, colBfull,
              0, workmat_full, colBfull);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, colBreduced,
              dimAreduced, 1, A_reduced_inv, dimAreduced, workmat_reduced, colBreduced,
              0, workmat_reduced, colBreduced);

  std::ofstream fout16("test/Aiworkmat_full.txt");
  for(int j = 0; j < dimAfull; j++) {
    for(int i = 0; i < colBfull; i++) {
      fout16 << workmat_full[i + colBfull*j] << "\t";
    }
    if (j < dimAfull - 1) fout16 << "\n";
  }


  // Now we initialize mat_b_rhs and evaluate it as B^{T} * workmat
  // Then mat_b_rhs * b = rhs

  std::vector<double> mat_b_rhs_full_vector(colBfull*colBfull,0);
  double* mat_b_rhs_full = mat_b_rhs_full_vector.data();
  std::vector<double> mat_b_rhs_reduced_vector(colBreduced*colBreduced,0);
  double* mat_b_rhs_reduced = mat_b_rhs_reduced_vector.data();

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colBfull, colBfull,
              dimAfull, 1, mat_B_full, colBfull, workmat_full, colBfull,
              0, mat_b_rhs_full, colBfull);  

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colBreduced, colBreduced,
              dimAreduced, 1, mat_B_reduced, colBreduced, workmat_reduced, colBreduced,
              0, mat_b_rhs_reduced, colBreduced);


  std::ofstream fout9("test/mat_b_rhs_full.txt");
  for(int j = 0; j < colBfull; j++) {
    for(int i = 0; i < colBfull; i++) {
      fout9 << mat_b_rhs_full[i + colBfull*j] << "\t";
    }
    if (j < colBfull - 1) fout9 << "\n";
  }

  std::ofstream fout10("test/mat_b_rhs_reduced.txt");
  for(int j = 0; j < colBreduced; j++) {
    for(int i = 0; i < colBreduced; i++) {
      fout10 << mat_b_rhs_reduced[i + colBreduced*j] << "\t";
    }
    if (j < colBreduced - 1) fout10 << "\n";
  }






  // Now we initialize b with rhs and solve the linear system

  std::vector<double> b_full_vector(colBfull,0);
  double* b_full = b_full_vector.data();
  std::vector<double> b_reduced_vector(colBreduced,0);
  double* b_reduced = b_reduced_vector.data();

  for (int i = 0; i < colBfull; i++) {
    b_full[i] = rhs_full[i];
    if (i < colBreduced) b_reduced[i] = rhs_reduced[i];
  }


  lapack_int* ipiv; char norm = 'I'; 
  ipiv = (lapack_int*)malloc(colBfull * sizeof(lapack_int));
  double anorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, colBfull, colBfull, mat_b_rhs_full, colBfull);
  ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, colBfull, 1, mat_b_rhs_full, colBfull, ipiv, b_full, 1); //mat updated to be LU
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  double rcond = 0;
  ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, colBfull, mat_b_rhs_full, colBfull, anorm, &rcond);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  rcond = 1/rcond;

  ipiv = (lapack_int*)malloc(colBreduced * sizeof(lapack_int));
  double anorm_r = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, colBreduced, colBreduced, mat_b_rhs_reduced, colBreduced);
  ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, colBreduced, 1, mat_b_rhs_reduced, colBreduced, ipiv, b_reduced, 1); //mat updated to be LU
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  double rcond_r = 0;
  ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, colBreduced, mat_b_rhs_reduced, colBreduced, anorm_r, &rcond_r);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  rcond_r = 1/rcond_r;

  //Calculate inf condition number
  std::cout << "\tNorm Format:\t" << norm << std::endl;
  std::cout << "\tNorm of mat:\t" << "full:\t" << anorm << "\treduced:\t"<< anorm_r << std::endl;
  std::cout << "\tCond number:\t" << "full:\t" << rcond << "\treduced:\t"<< rcond_r  << std::endl;




/*

  
  
  //Solve the matrix, result would be stored in rhs
  lapack_int* ipiv; char norm = 'I'; 
  ipiv = (lapack_int*)malloc(nn * sizeof(lapack_int));
  double anorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, nn, nn, mat, nn);
  int ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, nn, 1, mat, nn, ipiv, rhs, 1); //mat updated to be LU
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  double rcond = 0;
  ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, nn, mat, nn, anorm, &rcond);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  rcond = 1/rcond;

  //Calculate inf condition number
  std::cout << "\tNorm Format:\t" << norm << std::endl;
  std::cout << "\tNorm of mat:\t" << anorm << std::endl;
  std::cout << "\tCond number:\t" << rcond << std::endl;


  for(int i=0; i<solution.size(); i++) {
    if(index_correction[i] == -1) {
      double x = param.dsSpace.nodePtr(i)->val(0);
      double y = param.dsSpace.nodePtr(i)->val(1);
      solution[i] = bcVal(x,y);
    } else {
      solution[i] = rhs[index_correction[i]];
    }
  }

*/
  return 0;
} 
