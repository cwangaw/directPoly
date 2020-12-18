#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "ellipticPDE.h"
#include "fcns.h"
#include "Mesh/baseObjects.h"
#include "Mesh/polyMesh.h"
#include "directSerendipity.h"
#include "parameterData.h"
#include "polyQuadrature.h"
#include "Utilities/monitor.h"
#include "Utilities/debug.h"
#include <complex.h>
#include "lapacke.h"

using namespace directserendipity;
using namespace polymesh;
using namespace polyquadrature;


int EllipticPDE::solve(Monitor& monitor) {
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

  if(true) {
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
  
  DirectSerendipityArray solution(&(param.dsSpace));

  // Correct for BCSs: If node i is an interior node, i-th element of index_correction 
  // would give its index in our vector without BC nodes. If node i is a BC node, the  
  // i-th element would give -1.
  std::vector<int> index_correction(param.dsSpace.nNodes());
  int nn = 0;
  for (int i = 0; i < param.dsSpace.nNodes(); i++) {
    if (param.dsSpace.bcType(i) == BCType::interior) {
      index_correction[i] = nn;
      nn++;
    } else {
      index_correction[i] = -1;
    }
  }

  std::vector<double> mat_vector(nn*nn,0);
  double* mat = mat_vector.data();
  std::vector<double> rhs_vector(nn,0);
  double* rhs = rhs_vector.data();

  // quadrature points
  polyquadrature::PolyQuadrature quadRule(2*param.dsSpace.degPolyn()+3);

  monitor(1,"Matrix and RHS Assembly"); ////////////////////////////////////////

  for(int iElement=0; iElement<param.mesh.nElements(); iElement++) {
    DirectSerendipityFE* fePtr = param.dsSpace.finiteElementPtr(iElement);
    
    quadRule.setElement(fePtr->elementPtr());

    fePtr->initBasis(quadRule.pts(), quadRule.num());

    // Local matrix and rhs
    int nn_loc = fePtr->nNodes();
    std::vector<double> mat_loc(nn_loc*nn_loc,0), rhs_loc(nn_loc,0);

    // Determine local to global map
    int node_loc_to_gbl[nn_loc];
    for(int i=0; i<fePtr->nVertexNodes(); i++) {
      node_loc_to_gbl[i] = index_correction[fePtr->vertexNodePtr(i)->nodeIndex()];
    }
    for(int i=0; i<fePtr->nEdgeNodes(); i++) {
      node_loc_to_gbl[i + fePtr->nVertexNodes()] = index_correction[fePtr->edgeNodePtr(i)->nodeIndex()]; 
    }
    for(int i=0; i<fePtr->nCellNodes(); i++) {
      node_loc_to_gbl[i + fePtr->nVertexNodes() + fePtr->nEdgeNodes()]
	      = index_correction[fePtr->cellNodePtr(i)->nodeIndex()];
    }

    // Matrix and rhs assembly over elements
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt)[0];
      double y = quadRule.pt(iPt)[1];

      double valA = coefA(x,y);
      Tensor1 valB, valC;
      Tensor2 valD;
      coefB(x,y,valB); coefC(x,y,valC); coefD(x,y,valD);

      // Local interactions	      
      for(int jNode=0; jNode<nn_loc; jNode++) {
        if (node_loc_to_gbl[jNode] == -1) continue;
	      double valj = fePtr->basis(jNode, iPt); 
        Tensor1 gradValj = fePtr->basisGrad(jNode, iPt); 

        for(int iNode=0; iNode<nn_loc; iNode++) {
          //In case "shape functions" are not delta_{i,j} for BC nodes one day
          //if on BC, we use nodal basis function, otherwise we use shape functions
	        //double vali = (node_loc_to_gbl[jNode] == -1)? fePtr->basis(iNode, iPt) : fePtr->basisSF(iNode, iPt, quadRule.num());
          //Tensor1 gradVali = (node_loc_to_gbl[jNode] == -1)? fePtr->basisGrad(iNode, iPt) : fePtr->basisGradSF(iNode, iPt, quadRule.num());

          double vali = fePtr->basis(iNode, iPt);
          Tensor1 gradVali = fePtr->basisGrad(iNode, iPt);

          if (node_loc_to_gbl[iNode] != -1) {
            // +(a N_i,N_j);
	          mat_loc[iNode + nn_loc*jNode] += valA*vali*valj*quadRule.wt(iPt);
            // +(c dot grad(N_i), N_j)
            mat_loc[iNode + nn_loc*jNode] += (valC*gradVali)*valj*quadRule.wt(iPt);
            // +(b N_i, grad(N_j))
            mat_loc[iNode + nn_loc*jNode] += vali*(valB*gradValj)*quadRule.wt(iPt);
            // +(D grad(N_i), grad(N_j))
            mat_loc[iNode + nn_loc*jNode] += ((valD*gradVali)*gradValj)*quadRule.wt(iPt);
          }
          
          //rhs
          Node node_i = *(fePtr->nodePtr(iNode));
          if (node_loc_to_gbl[iNode] == -1) {
	          double bcVali = bcVal(node_i[0],node_i[1]);
            // -(a g(x_i)*N^{BC}_i, N_j)
            rhs_loc[jNode] -= valA*vali*valj*bcVali*quadRule.wt(iPt); 
            // -(c dot g(x_i)*grad(N^{BC}_i), N_j)
            rhs_loc[jNode] -= (valC*gradVali)*valj*bcVali*quadRule.wt(iPt);
            // -(b g(x_i)*N^{BC}_i, grad(N_j))
            rhs_loc[jNode] -= vali*(valB*gradValj)*bcVali*quadRule.wt(iPt);
            // -(D g(x_i)*grad(N^{BC}_i), grad(N_j))
            rhs_loc[jNode] -= ((valD*gradVali)*gradValj)*bcVali*quadRule.wt(iPt);
          }
	}
        // +(f, N_j)
	rhs_loc[jNode] += sourceVal(quadRule.pt(iPt)[0],quadRule.pt(iPt)[1])*valj*quadRule.wt(iPt);
      }
    }

    // Map local matrix and rhs to global
    
    for(int jNode=0; jNode<nn_loc; jNode++) {
      if (node_loc_to_gbl[jNode] != -1) {
        for(int iNode=0; iNode<nn_loc; iNode++) {
          if (node_loc_to_gbl[iNode] != -1) {
            mat[node_loc_to_gbl[iNode] + nn*node_loc_to_gbl[jNode]]
              += mat_loc[iNode + nn_loc*jNode];
          }
        }
      }
    }

    for(int iNode=0; iNode<nn_loc; iNode++) {
      if (node_loc_to_gbl[iNode] != -1) {
        rhs[node_loc_to_gbl[iNode]] += rhs_loc[iNode];
      }
    }
  }


  std::ofstream fout("test/matrix.txt");
  for(int j=0; j<nn; j++) {
    for(int i=0; i<nn; i++) {
      fout << mat[i + nn*j] << "\t";
    }
    if (j < nn - 1) fout << "\n";
  }
  std::ofstream rout("test/rhs.txt");
  for(int i=0; i<nn; i++) {
    rout << rhs[i];
    if (i < nn - 1) rout << "\n";
  }

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

  if(param.output_soln_DS_format > 0) {
  
    monitor(1,"Write Solution"); //////////////////////////////////////////////////

    switch(param.output_soln_DS_format) {
    case 1: {
      std::string fileName(param.directory_name);
      fileName += "solution_raw";
      solution.write_raw(fileName);
      break;
    }
    case 2: {
      std::string fileName(param.directory_name);
      fileName += "solution_mesh";
      std::string fileNameGrad(param.directory_name);
      fileNameGrad += "solution_grad_mesh";
      solution.write_matlab_mesh(fileName,fileNameGrad,
				 param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y);
      break;
    }
    }

  if(trueSolnKnown()) {
    monitor(0,"\nError estimate\n"); ///////////////////////////////////////////////
  
    double h = param.dsSpace.mesh()->maxElementDiameter();
    
    double l2Error = 0, l2GradError = 0, l2Norm = 0, l2GradNorm = 0;
    solution.l2normError(l2Error, l2GradError, l2Norm, l2GradNorm, trueSoln, trueGradSoln);
    
    std::cout << "  Max Element Diameter h:  " << h << std::endl;
    std::cout << "  L_2 Error:      " << l2Error << std::endl;
    std::cout << "  L_2 Grad Error: " << l2GradError << std::endl;
    std::cout << std::endl;
    std::cout << "  Relative L_2 Error:      " << l2Error/l2Norm << std::endl;
    std::cout << "  Relative L_2 Grad Error: " << l2GradError/l2GradNorm << std::endl;
    std::cout << std::endl;

    if(param.output_soln_DS_format > 0) {
      monitor(1,"Write True Solution"); ////////////////////////////////////////////

      DirectSerendipityArray u(&(param.dsSpace));

      for(int i=0; i<u.size(); i++) {
        double x = param.dsSpace.nodePtr(i)->val(0);
        double y = param.dsSpace.nodePtr(i)->val(1);
        u[i] = trueSoln(x,y);
      }

      switch(param.output_soln_DS_format) {
      case 1: {
        std::string fileName(param.directory_name);
        fileName += "true_solution_raw";
        u.write_raw(fileName);
        break;
      }
      case 2: {
        std::string fileName(param.directory_name);
        fileName += "true_solution_mesh";
        std::string fileNameGrad(param.directory_name);
        fileNameGrad += "true_solution_grad_mesh";
        u.write_matlab_mesh(fileName,fileNameGrad,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y);
        break;
      }
      }
    }
  }
  
  return 0;
} }
