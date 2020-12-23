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

using namespace directserendipity;
using namespace polymesh;
using namespace polyquadrature;

lapack_int matInv(double *A, unsigned n)
{
  // inplace inverse n x n matrix A.
  // matrix A is Column Major (i.e. firts line, second line ... *not* C[][] order)
  // returns:
  //   ret = 0 on success
  //   ret < 0 illegal argument value
  //   ret > 0 singular matrix
  int ipiv[n+1];
  lapack_int ret;

  ret =  LAPACKE_dgetrf(LAPACK_COL_MAJOR, n,n,A,n,ipiv);

  if (ret !=0) return ret;
  ret = LAPACKE_dgetri(LAPACK_COL_MAJOR,n,A,n,ipiv);
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
    u[1]=1;

    for(int i=0; i<p.size(); i++) {
      p[i]=0;
    }
    p[1]=1;

    for(int i=0; i<l.size(); i++) {
      l[i]=0;
    }
    l[1]=1;

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
          if (j < loc_dimAreduced && i < loc_dimAreduced) {
            curr_reduced_index = dimAreduced * (starting_Areduced + j) + (starting_Areduced + i);
            mat_A_reduced[curr_reduced_index] += evaluation;
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
          if (j < loc_dimAreduced && i < loc_colBreduced) {
            curr_reduced_index = colBreduced * (starting_Areduced + j) + (starting_colBreduced + i);
            mat_B_reduced[curr_reduced_index] += evaluation;
          }
        }
      }

      // Assemble RHS
      double f_wted = sourceVal(quadRule.pt(iPt)[0],quadRule.pt(iPt)[1]) * quadRule.wt(iPt);
      for (int j = 0; j < loc_colBfull; j++) {
        curr_full_index = starting_colBfull + j;
        evaluation = f_wted * dgePtr->basis(j,iPt);
        rhs_full[curr_full_index] += evaluation;
        if (j < loc_colBreduced) {
          curr_reduced_index = starting_colBreduced + j;
          rhs_reduced[curr_reduced_index] += evaluation;
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
              if (j < loc_dimAreduced) {
                curr_reduced_index = colL * (starting_Areduced + j) + (loc_to_int*(param.dmSpace.degPolyn()+1)+i);
                mat_L_reduced[curr_reduced_index] += evaluation;
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

/*

  monitor(1,"Solution of linear system"); ////////////////////////////////////////
  
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
