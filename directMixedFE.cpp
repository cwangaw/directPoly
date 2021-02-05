#include <cmath>
#include <iostream>
#include <vector>

#include <complex.h>
#include "lapacke.h"
#include <stdio.h>
#include <assert.h>

using namespace std;

#include "Utilities/debug.h"
#include "Mesh/polyMesh.h"
using namespace polymesh;
#include "directMixed.h"
using namespace directserendipity;

////////////////////////////////////////////////////////////////////////////////
// class DirectMixedFE

void DirectMixedFE::
         set_directmixedfe(DirectMixed* dmSpace, polymesh::PolyElement* element) {
  my_dm_space = dmSpace;
  my_poly_element = element;

  num_vertices = element->nVertices();
  polynomial_degree = my_dm_space->degPolyn();
  ref_origin = Point(my_poly_element->vertexPtr(num_vertices-1)->val(0), my_poly_element->vertexPtr(num_vertices-1)->val(1));

  if (polynomial_degree >= num_vertices - 3) {
    dim_supp = num_vertices * (num_vertices - 3) / 2;
  } else {
    dim_supp = (2*num_vertices-3) * (polynomial_degree+1) - pow(polynomial_degree+1,2) - 2;
    dim_supp /= 2;
  }

  dim_v = (polynomial_degree+3)*(polynomial_degree+1) + dim_supp;
  dim_v_div = (polynomial_degree+2) * (polynomial_degree+1) /2;
  dim_curlpart = dim_v - dim_v_div;
};

DirectMixedFE::~DirectMixedFE() {
  if (high_order_ds_space) delete high_order_ds_space;
  if (one_element_mesh) delete one_element_mesh;
  if (v_value_n) delete[] v_value_n;
  if (v_div_value_n) delete[] v_div_value_n;
  if (v_edge_value_n) delete[] v_edge_value_n;
}

void DirectMixedFE::initBasis(const Point* pt, int num_pts) {
  if(num_pts <= 0) return;
  num_eval_pts = num_pts;

  // Allocate space for the resulting values
  // Note that dimension of v_div is same as dim_w

  if(v_value_n) delete[] v_value_n;
  v_value_n = new Tensor1[num_pts * dim_v];
  if (v_div_value_n) delete[] v_div_value_n;
  v_div_value_n = new double[num_pts * dim_v_div];
  if (v_edge_value_n) delete[] v_edge_value_n;
  v_edge_value_n = new double[num_pts * dim_v * num_vertices];

  int curr_index = 0; // A variable that store current function index
  double x,y; // Store the position of pt[pt_index]
  Tensor1 result; // Store the evaluated result

  //////////////////////////////////////////////////////
  //                                                  //
  //      Construct mixed space V, div(V) and         //
  //          its normal components on each edge      //
  //                                                  //
  //////////////////////////////////////////////////////                        

  /*    curl_\Po_{r+1}   */

  // Order: x^0y^1, x^1y^0, x^0y^2, x^1y^1, x^2y^0, x^0y^3, x^1y^2, ...
  for (int k = 1; k <= polynomial_degree + 1; k++) {
    for (int m = 0; m <= k; m++) {
      for (int pt_index = 0; pt_index < num_pts; pt_index++) {
        x = pt[pt_index].val(0)-ref_origin.val(0);
        y = pt[pt_index].val(1)-ref_origin.val(1);
        if (m == 0) {
          result.set( (k-m)*pow(y,k-m-1) ,0);
        } else if (m == k) {
          result.set( 0, -m*pow(x,m-1) );
        } else {
          result.set( (k-m)*pow(x,m)*pow(y,k-m-1), -m*pow(x,m-1)*pow(y,k-m) );
        }
        v_value_n[pt_index * dim_v + curr_index] = result;
        for ( int nEdge = 0; nEdge < num_vertices; nEdge++ ) {
          v_edge_value_n[pt_index * dim_v * num_vertices + dim_v * nEdge + curr_index] = result * my_poly_element -> edgePtr(nEdge) -> normal();
        }
      }
      curr_index += 1;
    }
  }

  /*    curl_\Supp^{\cDS}_{r+1}(E)    */

  // Set up {\cDS}_{r+1}(E)
  if(one_element_mesh) delete one_element_mesh;
	one_element_mesh = new polymesh::PolyMesh(my_poly_element);
	if(high_order_ds_space) delete high_order_ds_space;
  high_order_ds_space = new DirectSerendipity(polynomial_degree+1,one_element_mesh);
  high_order_ds_space->finiteElementPtr(0)->initBasis(pt, num_pts);
  int higher_order = high_order_ds_space->finiteElementPtr(0)->polynomial_degree;

  if (polynomial_degree < num_vertices - 3) {
    // Get supplemental functions for small r
    // Get gradient of phi_{v,i} in A_\Supp
    for (int i = 0; i < num_vertices; i++) {
      // Exclude vertices in A_\Po
      if ((i >= higher_order - 2) && (i <= higher_order)) continue;
      if ((higher_order == 1) && (i == num_vertices - 1)) continue;

      for (int pt_index = 0; pt_index < num_pts; pt_index++) {
        result.set(high_order_ds_space->finiteElementPtr(0)->basisGrad(i,pt_index).val(1), 
                    -high_order_ds_space->finiteElementPtr(0)->basisGrad(i,pt_index).val(0));
        v_value_n[pt_index * dim_v + curr_index] = result;
        for ( int nEdge = 0; nEdge < num_vertices; nEdge++ ) {
          v_edge_value_n[pt_index * dim_v * num_vertices + dim_v * nEdge + curr_index] = result * my_poly_element -> edgePtr(nEdge) -> normal();
        }
      }
      curr_index += 1;
    } 

    // Get gradient of phi_{e,nEdge,jNode} in A_\Supp
    for (int nEdge = 0; nEdge < num_vertices; nEdge++) {
      for (int jNode = 0; jNode < (higher_order - 1); jNode++) {
        //Exclude phi_{e,nEdge,jNode} in A_\Po
        if ((nEdge <= higher_order - 2) && (jNode <= nEdge)) continue;
        if ((nEdge == higher_order - 1) || (nEdge == higher_order)) continue;

        for (int pt_index = 0; pt_index < num_pts; pt_index++) {
          result.set(high_order_ds_space->finiteElementPtr(0)->gradEdgeBasis(nEdge,jNode,pt_index).val(1), 
                                          -high_order_ds_space->finiteElementPtr(0)->gradEdgeBasis(nEdge,jNode,pt_index).val(0));
          v_value_n[pt_index * dim_v + curr_index] = result;
          for ( int iEdge = 0; iEdge < num_vertices; iEdge++ ) {
            v_edge_value_n[pt_index * dim_v * num_vertices + dim_v * iEdge + curr_index] = result * my_poly_element -> edgePtr(iEdge) -> normal();
          }
        }

        curr_index += 1;
      }
    }

    } 
    else {
      // Get supplemental functions for big r

      // Initialization of variables for calling phi_k_l
      Tensor1 gradresult_k_l,gradresult_l_k;
      double result_of_no_use;

      for (int k=0; k <= num_vertices-3; k++) {
        for (int l=k+2; l <= num_vertices-1; l++) {
          if (k==0 && l==num_vertices-1) { continue; }
          for (int pt_index = 0; pt_index < num_pts; pt_index++) {
            // Update gradresult to evaluate phi_k_l at pt[pt_index]
            high_order_ds_space->finiteElementPtr(0)->phi_k_l(k,l,pt[pt_index],result_of_no_use,gradresult_k_l);
            high_order_ds_space->finiteElementPtr(0)->phi_k_l(l,k,pt[pt_index],result_of_no_use,gradresult_l_k);
            // Store value of curl of phi_k_l
            //result.set(gradresult_k_l.val(1),-gradresult_k_l.val(0));
            result.set(gradresult_k_l.val(1)-gradresult_l_k.val(1),-gradresult_k_l.val(0)+gradresult_l_k.val(0));
            v_value_n[pt_index * dim_v + curr_index] = result;
            for ( int nEdge = 0; nEdge < num_vertices; nEdge++ ) {
              v_edge_value_n[pt_index * dim_v * num_vertices + dim_v * nEdge + curr_index] = result * my_poly_element -> edgePtr(nEdge) -> normal();
            }
          }
          curr_index += 1;
        }
      }
    }

  /*    \x\Po_s(E), where s = r-1, r, and its divergence    */

  // Here we calculate polynomials up to r, 
  // where the first r(r+1)/2 ones are polynomials up to s=r-1
  // and the last r+1 ones are polynomials with EXACTLY degree r
  
  int curr_v_div_index = 0;
  for (int k=0; k<= polynomial_degree; k++) {
    for (int m=0; m <= k; m++) {
      for (int pt_index = 0; pt_index < num_pts; pt_index++) {
        x = pt[pt_index].val(0)-ref_origin.val(0);
        y = pt[pt_index].val(1)-ref_origin.val(1);

        v_div_value_n[pt_index * dim_v_div + curr_v_div_index] = (k+2)*pow(x,m)*pow(y,k-m);
        result.set( pow(x,m+1)*pow(y,k-m), pow(x,m)*pow(y,k-m+1) );
        v_value_n[pt_index * dim_v + curr_index] = result;
        for ( int nEdge = 0; nEdge < num_vertices; nEdge++ ) {
          v_edge_value_n[pt_index * dim_v * num_vertices + dim_v * nEdge + curr_index] = result * my_poly_element -> edgePtr(nEdge) -> normal();
        }
      }
      curr_index += 1;
      curr_v_div_index += 1;
    }
  }

  // Check if our indexing is working as expected
  assert(curr_index == dim_v);
  assert(curr_v_div_index == dim_v_div);
};

// Eval for DirectMixedFE

void DirectMixedFE::eval(const Point* pt, Tensor1* result, int num_pts, char type, double* dofs) {
  initBasis(pt,num_pts);
  int dim = (type == 'f') ? dimVFull() : dimVReduced();
  for(int n=0; n<num_pts; n++) { 
    result[n].set(0,0);
    for(int i=0; i<dim; i++) {
      result[n] += dofs[i]*basis(i,n);
    }
  }
}


void DirectMixedFE::eval(const Point* pt, Tensor1* fullResult, Tensor1* reducedResult, int num_pts, 
              double* full_dofs, double* reduced_dofs) {
  initBasis(pt,num_pts);
  for (int n=0; n<num_pts; n++) {
    fullResult[n].set(0,0);
    reducedResult[n].set(0,0);

    for(int i=0; i<dimVFull(); i++) {
      fullResult[n] += full_dofs[i]*basis(i,n);
      if (i<dimVReduced()) reducedResult[n] += reduced_dofs[i]*basis(i,n);
    }
  }
};

void DirectMixedFE::eval_div(const Point* pt, double* result, int num_pts, char type, double* dofs) {
  initBasis(pt,num_pts);
  int dim = (type == 'f') ? dimXPoFull() : dimXPoReduced();
  for(int n=0; n<num_pts; n++) { 
    result[n] = 0;

    for(int i=0; i<dim; i++) {
      result[n] += dofs[dimCurlPart()+i]*basisdivXPo(i,n);
    }
  }
}


void DirectMixedFE::eval_div(const Point* pt, double* fullResult, double* reducedResult, int num_pts, 
              double* full_dofs, double* reduced_dofs) {
  initBasis(pt,num_pts);
  for (int n=0; n<num_pts; n++) {
    fullResult[n] = 0;
    reducedResult[n] = 0;

    for(int i=0; i<dimXPoFull(); i++) {
      fullResult[n] += full_dofs[dimCurlPart()+i]*basisdivXPo(i,n);
      if (i<dimXPoReduced()) reducedResult[n] += reduced_dofs[dimCurlPart()+i]*basisdivXPo(i,n);
    }
  }
};


// Output for DirectMixedFE

void DirectMixedFE::write_raw(std::ofstream& fout) const {
  fout << "    DIRECT MIXED FE\n";
  fout << "    my_dm_space       = " << my_dm_space << "\n";
  fout << "    my_poly_element   = " << my_poly_element << "\n";
  fout << "    num_vertices      = " << num_vertices << "\n";
  fout << "    polynomial_degree = " << polynomial_degree << "\n";
  fout << "    dimension of V_full = " << dim_v << "\n";
  fout << "    dimension of V_reduced = " << dimVReduced() << "\n";
  fout << "    dimension of Curl part = " << dim_curlpart << "\n";
  fout << "    dimension of xPo_full = " << dim_v_div << "\n";
  fout << "    dimension of xPo_reduced = " << dimXPoReduced() << "\n";

  fout << "    vertex nodes:";
  for(int i=0; i< my_poly_element->nVertices(); i++) {
    fout << " " << my_poly_element->vertexPtr(i)
	 << " (" << i << ")";
  }
  fout << "\n";
}

int DirectMixedFE::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// Class DirectDGFE

void DirectDGFE::
         set_directdgfe(DirectMixed* dmSpace, polymesh::PolyElement* element) {
  my_dm_space = dmSpace;
  my_poly_element = element;
  num_vertices = element->nVertices();
  polynomial_degree = my_dm_space->degPolyn();
  ref_origin = Point(my_poly_element->vertexPtr(num_vertices-1)->val(0), my_poly_element->vertexPtr(num_vertices-1)->val(1));
//my_poly_element -> center();

  dim_w = (polynomial_degree+2) * (polynomial_degree+1) /2;
};

DirectDGFE::~DirectDGFE() {
  if (value_n) delete[] value_n;
}

void DirectDGFE::initBasis(const Point* pt, int num_pts) {
  if(num_pts <= 0) return;
  num_eval_pts = num_pts;

  // Allocate space for the resulting values
  if (value_n) delete[] value_n;
  value_n = new double[num_pts * dim_w];

  /*    W_s = \Po_s(E), where s = r-1, r    */

  // Here we calculate polynomials up to r, 
  // where the first r(r+1)/2 ones are polynomials up to s=r-1
  // and the last r+1 ones are polynomials with EXACTLY degree r
  
  int curr_index = 0;
  double x,y;

  for (int k=0; k<= polynomial_degree; k++) {
    for (int m=0; m <= k; m++) {
      for (int pt_index = 0; pt_index < num_pts; pt_index++) {
        x = pt[pt_index].val(0)-ref_origin.val(0);
        y = pt[pt_index].val(1)-ref_origin.val(1);
        value_n[pt_index * dim_w + curr_index] = pow(x,m)*pow(y,k-m);
      }
      curr_index += 1;
    }
  }

  //Check if our indexing is working as expected
  assert(curr_index == dim_w);
};


// Eval for DirectDGFE

void DirectDGFE::eval(const Point* pt, double* result, int num_pts, char type, double* dofs) {
  initBasis(pt,num_pts);
  int dim = (type == 'f') ? dimFull() : dimReduced();
  for(int n=0; n<num_pts; n++) {
    result[n] = 0;
    Point p(pt[n]);

    for(int i=0; i<dim; i++) {
      result[n] += dofs[i]*basis(i,n);
    }
  }
}

void DirectDGFE::eval(const Point* pt, double* fullResult, double* reducedResult, int num_pts, 
              double* full_dofs, double* reduced_dofs) {
  initBasis(pt,num_pts);
  for (int n=0; n<num_pts; n++) {
    fullResult[n] = 0;
    reducedResult[n] = 0;
    Point p(pt[n]);

    for(int i=0; i<dimFull(); i++) {
      fullResult[n] += full_dofs[i]*basis(i,n);
      if (i<dimReduced()) reducedResult[n] += reduced_dofs[i]*basis(i,n);
    }
  }
};


// Output for DirectDGFE

void DirectDGFE::write_raw(std::ofstream& fout) const {
  fout << "    DIRECT DG FE\n";
  fout << "    my_dm_space       = " << my_dm_space << "\n";
  fout << "    my_poly_element   = " << my_poly_element << "\n";
  fout << "    num_vertices      = " << num_vertices << "\n";
  fout << "    polynomial_degree = " << polynomial_degree << "\n";
  fout << "    dimension of W_full = " << dim_w << "\n";
  fout << "    dimension of W_reduced = " << dimReduced() << "\n";

  fout << "    vertex nodes:";
  for(int i=0; i< my_poly_element->nVertices(); i++) {
    fout << " " << my_poly_element->vertexPtr(i)
	 << " (" << i << ")";
  }
  fout << "\n";
}

int DirectDGFE::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}



////////////////////////////////////////////////////////////////////////////////
// Class DirectEdgeDGFE

void DirectEdgeDGFE::
         set_directedgedgfe(DirectMixed* dmSpace, Edge* edge) {
  my_dm_space = dmSpace;
  polynomial_degree = my_dm_space->degPolyn();
  my_edge = edge;

  dim_l =  polynomial_degree + 1;
};

DirectEdgeDGFE::~DirectEdgeDGFE() {
  if (value_n) delete[] value_n;
}

void DirectEdgeDGFE::initBasis(const Point* pt, int num_pts) {
  if(num_pts <= 0) return;
  num_eval_pts = num_pts;

  // Allocate space for the resulting values
  if (value_n) delete[] value_n;
  value_n = new double[num_pts * dim_l];

  
  /////////////////////////////////////////////////////////////
  //                                                         //
  //      Construct lagrange multiplyer space \Lambda_r      //
  //                                                         //
  /////////////////////////////////////////////////////////////

  // Initialize curr_index to 0
  int curr_index = 0;

  for (int deg = 0; deg <= polynomial_degree; deg++) {
    for (int pt_index = 0; pt_index < num_pts; pt_index++) {
      value_n[pt_index * dim_l + curr_index] = pow(projToEdge(pt[pt_index]), deg);
    }
    curr_index += 1;
  }


  //Check if our indexing is working as expected
  assert(curr_index == dim_l);
};

double DirectEdgeDGFE::projToEdge(const Point& p) const {
  Tensor1 tau(my_edge->tangent());
  return (p[0] - my_edge->vertexPtr(0)->val(0))*tau[0] + (p[1] - (my_edge->vertexPtr(0)->val(1)))*tau[1];
};

// Eval for DirectEdgeDGFE

void DirectEdgeDGFE::eval(const Point* pt, double* result, int num_pts, double* dofs) {
  initBasis(pt,num_pts);
  for(int n=0; n<num_pts; n++) {
    result[n] = 0;
    Point p(pt[n]);

    for(int i=0; i<dim(); i++) {
      result[n] += dofs[i]*basis(i,n);
    }
  }
}

// Output for DirectEdgeDGFE

void DirectEdgeDGFE::write_raw(std::ofstream& fout) const {
  fout << "    DIRECT EDGE DG FE\n";
  fout << "    my_dm_space       = " << my_dm_space << "\n";
  fout << "    my_edge   = " << my_edge << "\n";
  fout << "    polynomial_degree = " << polynomial_degree << "\n";
  fout << "    dimension of lambda space = " << dim_l << "\n";

  fout << "    vertex nodes:";
  for(int i=0; i < 2; i++) {
    fout << " " << my_edge->vertexPtr(i)
	 << " (" << i << ")";
  }
  fout << "\n";
}

int DirectEdgeDGFE::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}