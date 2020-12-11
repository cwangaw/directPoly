#include <string>
#include <vector>
#include <cmath>
#include <iostream>

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

    DirectSerendipity s_dsSpace(8,&s_mesh);
    fileName = parameterDataPtr()->directory_name;
    fileName += "dsSpace_e4";
    s_dsSpace.write_matlab(fileName);
  }
  
  // TEST BASIS FUNCTIONS //////////////////////////////////////////////////

  if(false) {
    monitor(0,"\nTest basis functions for element 0\n");

    DirectSerendipityArray u(&(parameterDataPtr()->dsSpace));

    double PI = 3.141592653589793238463;
    
    for(int i=0; i<u.size(); i++) {
      double x = parameterDataPtr()->dsSpace.nodePtr(i)->val(0);
      double y = parameterDataPtr()->dsSpace.nodePtr(i)->val(1);
      u[i] = sin(PI*x)*sin(PI*y); //x*x+y*y;
      //u[i]=0;
      //if (fabs(x-5)<1e-6&&fabs(y-6)<1e-6) {u[i]=1;}
    }
    //u[21]=1;
    monitor(1,"Write Array");

    std::string fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh";
    std::string fileName_grad = parameterDataPtr()->directory_name;
    fileName_grad += "basis_grad_mesh";
    //u.write_matlab_mesh_by_pt(fileName,fileName_grad,3,3);
    u.write_matlab_mesh(fileName,fileName_grad,301,301);
  }

  // TEST QUADRATURE ///////////////////////////////////////////////////////

  if(false) {
    monitor(1,"Test Quadrature Rules");
    testPolyQuadrature(&(parameterDataPtr()->mesh),1e-6);
    return 0;
  }
  
  // SOLVE THE PDE ///////////////////////////////////////////////////////////

  monitor(0,"\nSolve the PDE\n");
  
  DirectMixedArray solution(&(param.dmSpace));

  // Initialize matrix A for both full and reduced space
  int dimAfull = 0, dimAreduced = 0;

  for(int iElement=0; iElement<param.mesh.nElements(); iElement++) { 
    dimAfull += param.dmSpace.MixedElementPtr(iElement)->dimVFull();
    dimAreduced += param.dmSpace.MixedElementPtr(iElement)->dimVReduced();
  }

  std::vector<double> mat_A_full_vector(dimAfull*dimAfull,0);
  double* mat_A_full = mat_A_full_vector.data();
  std::vector<double> mat_A_reduced_vector(dimAreduced*dimAreduced,0);
  double* mat_A_reduced = mat_A_reduced_vector.data();

  // Initialize matrix B for both full and reduced space
  int rowBfull = dimAfull, rowBreduced = dimAreduced;
  int colBfull = 0, colBreduced = 0;

  for (int iElement=0; iElement<param.mesh.nElements(); iElement++) {
    colBfull += param.dmSpace.DGElementPtr(iElement)->dimFull();
    colBreduced += param.dmSpace.DGElementPtr(iElement)->dimReduced();
  }

  std::vector<double> mat_B_full_vector(rowBfull*colBfull,0);
  double* mat_B_full = mat_B_full_vector.data();
  std::vector<double> mat_B_reduced_vector(rowBreduced*colBreduced,0);
  double* mat_B_reduced = mat_B_reduced_vector.data();

  // Initialize matrix L
  int rowLfull = dimAfull, rowLreduced = dimAreduced;
  int colL = param.dmSpace.nInteriorEdges() * (param.dmSpace.degPolyn()+1);

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
  int starting_colL = 0;

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

    // Matrix and rhs assembly over elements
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt)[0];
      double y = quadRule.pt(iPt)[1];

      Tensor2 valD;
      coefD_inv(x,y,valD);

      // Assemble matrix A

      for (int j = 0; j < loc_dimAfull; j++) {
        Tensor1 v_j = mePtr -> basis(j,iPt);
        for (int i = 0; i < loc_dimAfull; i++) {
          Tensor1 u_i = mePtr -> basis(j,iPt);
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
        double divv_j = mePtr -> divXPo(j - mePtr -> dimCurlPart(),iPt);
        for (int i = 0; i < loc_colBfull; i++ ) {
          double p_i = dgePtr -> basis(i,iPt);
          curr_full_index = dimAfull * (starting_Afull + j) + (starting_colBfull + i);
          evaluation = divv_j * p_i * quadRule.wt(iPt);
          mat_B_full[curr_full_index] += evaluation;
          if (j < loc_dimAreduced && i < loc_colBreduced) {
            curr_reduced_index = dimAreduced * (starting_Areduced + j) + (starting_colBreduced + i);
            mat_B_reduced[curr_reduced_index] += evaluation;
          }
        }
      }

      // Assemble matrix L
      numEdges = param.mesh.elementPtr(iElement) -> nVertices();

      // Column indexing of L: (global edge indexed by interior)
      // {edge(0),Func(0)}, {edge(0),Func(1)}, ..., {edge(0),Func(eePtr[0]->dim()-1)},
      // {edge(1),Func(0)}, {edge(1),Func(1)}, ..., {edge(1),Func(eePtr[1]->dim()-1)},
      // ..., {edge(nInteriorEdge()-1),Func(eePtr[nInteriorEdge()-1]->dim()-1)} 

      for (int iEdge = 0; iEdge < numEdges; iEdge ++) {
        loc_to_int = param.dmSpace.interiorEdgeIndex(iElement,iEdge);
        if (loc_to_int == -1) continue; // If the edge is on boundary, we skip the loop
        for (int i = 0; i < eePtr[loc_to_int]->dim(); i++){
          double l_i = eePtr[loc_to_int]->basis(i,iPt);
          for (int j = 0; j < loc_dimAfull; j++) {
            double v_jdotNu = mePtr->basisDotNu(j,iEdge,iPt);
            // Here we use the property that eePtr[loc_to_int]->dim() is the same (degPolyn()+1)for every edge in our mesh
            curr_full_index = dimAfull * (starting_Afull + j) + (loc_to_int*(param.dmSpace.degPolyn()+1)+i);
            evaluation = l_i * v_jdotNu * quadEdgeRule[loc_to_int].wt(iPt);
            mat_L_full[curr_full_index] += evaluation;
            if (j < loc_dimAreduced) {
              curr_reduced_index = dimAreduced * (starting_Areduced + j) + (loc_to_int*(param.dmSpace.degPolyn()+1)+i);
              mat_L_reduced[curr_reduced_index] += evaluation;
            }
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


    starting_Afull += loc_dimAfull;
    starting_Areduced += loc_dimAreduced;
    starting_colBreduced += loc_colBreduced;
    starting_colBfull += loc_colBfull;

  }



  // std::ofstream fout("test/matrix.txt");
  // for(int j=0; j<nn; j++) {
  //   for(int i=0; i<nn; i++) {
  //     fout << mat[i + nn*j] << "\t";
  //   }
  //   if (j < nn - 1) fout << "\n";
  // }
  // std::ofstream rout("test/rhs.txt");
  // for(int i=0; i<nn; i++) {
  //   rout << rhs[i];
  //   if (i < nn - 1) rout << "\n";
  // }

  // monitor(1,"Solution of linear system"); ////////////////////////////////////////
  
  // //Solve the matrix, result would be stored in rhs
  // lapack_int* ipiv; char norm = 'I'; 
  // ipiv = (lapack_int*)malloc(nn * sizeof(lapack_int));
  // double anorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, nn, nn, mat, nn);
  // int ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, nn, 1, mat, nn, ipiv, rhs, 1); //mat updated to be LU
  // if(ierr) { // ?? what should we do ???
  //   std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  // }
  // double rcond = 0;
  // ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, nn, mat, nn, anorm, &rcond);
  // if(ierr) { // ?? what should we do ???
  //   std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  // }
  // rcond = 1/rcond;

  // //Calculate inf condition number
  // std::cout << "\tNorm Format:\t" << norm << std::endl;
  // std::cout << "\tNorm of mat:\t" << anorm << std::endl;
  // std::cout << "\tCond number:\t" << rcond << std::endl;


  // for(int i=0; i<solution.size(); i++) {
  //   if(index_correction[i] == -1) {
  //     double x = param.dsSpace.nodePtr(i)->val(0);
  //     double y = param.dsSpace.nodePtr(i)->val(1);
  //     solution[i] = bcVal(x,y);
  //   } else {
  //     solution[i] = rhs[index_correction[i]];
  //   }
  // }

  // if(param.output_soln_format) {
  

  //   monitor(1,"Write Solution"); //////////////////////////////////////////////////

  //   switch(param.output_soln_format) {
  //   case 1: {
  //     std::string fileName(param.directory_name);
  //     fileName += "solution_raw";
  //     solution.write_raw(fileName);
  //     break;
  //   }
  //   case 2: {
  //     std::string fileName(param.directory_name);
  //     fileName += "solution_mesh";
  //     std::string fileNameGrad(param.directory_name);
  //     fileNameGrad += "solution_grad_mesh";
  //     solution.write_matlab_mesh(fileName,fileNameGrad,
	// 			 param.output_mesh_numPts_x,param.output_mesh_numPts_y);
  //     break;
  //   }
  //   }

  // if(trueSolnKnown()) {
  //   monitor(0,"\nError estimate\n"); ///////////////////////////////////////////////
  
  //   double h = param.dsSpace.mesh()->maxElementDiameter();
    
  //   double l2Error = 0, l2GradError = 0, l2Norm = 0, l2GradNorm = 0;
  //   solution.l2normError(l2Error, l2GradError, l2Norm, l2GradNorm, trueSoln, trueGradSoln);
    
  //   std::cout << "  Max Element Diameter h:  " << h << std::endl;
  //   std::cout << "  L_2 Error:      " << l2Error << std::endl;
  //   std::cout << "  L_2 Grad Error: " << l2GradError << std::endl;
  //   std::cout << std::endl;
  //   std::cout << "  Relative L_2 Error:      " << l2Error/l2Norm << std::endl;
  //   std::cout << "  Relative L_2 Grad Error: " << l2GradError/l2GradNorm << std::endl;
  //   std::cout << std::endl;

  //   if(param.output_soln_format) {
  //     monitor(1,"Write True Solution"); ////////////////////////////////////////////

  //     DirectSerendipityArray u(&(param.dsSpace));

  //     for(int i=0; i<u.size(); i++) {
  //       double x = param.dsSpace.nodePtr(i)->val(0);
  //       double y = param.dsSpace.nodePtr(i)->val(1);
  //       u[i] = trueSoln(x,y);
  //     }

  //     switch(param.output_soln_format) {
  //     case 1: {
  //       std::string fileName(param.directory_name);
  //       fileName += "true_solution_raw";
  //       u.write_raw(fileName);
  //       break;
  //     }
  //     case 2: {
  //       std::string fileName(param.directory_name);
  //       fileName += "true_solution_mesh";
  //       std::string fileNameGrad(param.directory_name);
  //       fileNameGrad += "true_solution_grad_mesh";
  //       u.write_matlab_mesh(fileName,fileNameGrad,
  //       param.output_mesh_numPts_x,param.output_mesh_numPts_y);
  //       break;
  //     }
  //     }
  //   }
  // }
  
  return 0;
} 
