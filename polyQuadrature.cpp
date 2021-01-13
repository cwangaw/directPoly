#include <cmath>

#include "polyQuadrature.h"
#include "debug.h"

using namespace polyquadrature;

////////////////////////////////////////////////////////////////////////////////
// class PolyQuadrature

void PolyQuadrature::set_rule(int desired_dop) {
  my_desired_dop = desired_dop;

  my_rule = ruleForTriangle.size()-1;
  for(unsigned long int i=0; i<ruleForTriangle.size(); i++) {
    if(ruleForTriangle[i].dop >= desired_dop) { my_rule = i; break; }
  }
  my_dop  = ruleForTriangle[my_rule].dop;

  num_pts_ref = ruleForTriangle[my_rule].num;
  my_pts_ref  = ruleForTriangle[my_rule].pts;
  my_wts_ref  = ruleForTriangle[my_rule].wts;
}

void PolyQuadrature::set_element(polymesh::PolyElement* element) {
  my_element = element;
  if(!element) return;
  
  int num_triangles = my_element->nVertices();
  num_pts = num_triangles*num_pts_ref;

  // Quadrature points and weights

  if(my_pts) delete[] my_pts;
  my_pts = new Point[num_pts];

  if(my_wts) delete[] my_wts;
  my_wts = new double[num_pts];

  // Triangles from center of polygon
  Point center = my_element->center();

  for(int i=0; i<num_triangles; i++) {
    // v0 = (0,0)
    Point v1(*(my_element->vertexPtr(i))); v1 -= center;
    Point v2(*(my_element->vertexPtr((i+1) % num_triangles))); v2 -= center;

    Tensor2 mappingMatrix(v1[0], v2[0], v1[1], v2[1]);
    double jacobian = std::abs(mappingMatrix.determinant())/2;

    int kk = i*num_pts_ref;
    for(int j=0; j<num_pts_ref; j++) {
      mappingMatrix.mult(ruleForTriangle[my_rule].pts[j], my_pts[kk+j]);
      my_pts[kk+j] += center;
      my_wts[kk+j] = ruleForTriangle[my_rule].wts[j]*jacobian;
    }
  }
}

PolyQuadrature::~PolyQuadrature() {
  if(my_pts) delete[] my_pts;
  if(my_wts) delete[] my_wts;
}

// Test degree of precision on a mesh over the domain [0,10]^2
void polyquadrature::testPolyQuadrature(polymesh::PolyMesh* mesh, double eps, int toDOP) {
  auto f = [](Point& x, int i, int j) { return std::pow(x[0],i)*std::pow(x[1],j); };
  auto trueIntegF = [](int i, int j) { return std::pow(10,i+1)/(i+1)*std::pow(10,j+1)/(j+1); };

  PolyQuadrature quadrature;

  std::cout << std::endl;
  for(int testDOP=2; testDOP<=toDOP; testDOP++) {
    quadrature.setRule(testDOP);
    if(testDOP != quadrature.dop() ) continue;
    std::cout << "DOP = " << testDOP << "\n";

    for(int i=0; i<=testDOP; i++) {
      for(int j=0; j<=testDOP-i; j++) {

	double full_integ = 0;
	for(int iElem=0; iElem<mesh->nElements(); iElem++) {
	  quadrature.setElement(mesh->elementPtr(iElem));

	  double integ = 0;
	  for(int k=0; k<quadrature.num(); k++) {
	    integ += f(quadrature.pts()[k],i,j)*quadrature.wt(k);
	  }
	  full_integ += integ;
	}
	double true_integ = trueIntegF(i,j);

	double err = std::abs(true_integ - full_integ);
	if(err > eps) std::cout << "  i = " << i << ", j = " << j << ", err = " << err << "\n";
      }
    }
  }
  std::cout << std::endl;
};

////////////////////////////////////////////////////////////////////////////////
// class PolyEdgeQuadrature

void PolyEdgeQuadrature::set_rule(int desired_dop) {
  my_desired_dop = desired_dop;

    my_rule = ruleForEdge.size()-1;

    for(unsigned long int i=0; i<ruleForEdge.size(); i++) {
      if(ruleForEdge[i].num >= desired_dop) { my_rule = i; break; }
    }

    my_dop  = ruleForEdge[my_rule].num;
    num_pts = ruleForEdge[my_rule].num;
    my_pts_ref  = ruleForEdge[my_rule].pts;
    my_wts_ref  = ruleForEdge[my_rule].wts;
}

void PolyEdgeQuadrature::set_edge(polymesh::Edge* edge) {
  my_edge = edge;
  if(!edge) return;

  // Quadrature points and weights

  if (my_pts) delete[] my_pts;
  my_pts = new Point[num_pts];

  if (my_wts) delete[] my_wts;
  my_wts = new double[num_pts];

  Point v0(*(my_edge->vertexPtr(0)));
  Point v1(*(my_edge->vertexPtr(1)));
  Tensor1 vec(v1-v0);

  for (int i = 0; i < num_pts; i++) {
    my_pts[i].set(v0);
    my_pts[i] += (my_pts_ref[i]+1)/2 * (v1 - v0);
    my_wts[i] = my_wts_ref[i] * vec.norm()/ 2;
  }

}

PolyEdgeQuadrature::~PolyEdgeQuadrature() {
  if (my_pts) delete[] my_pts;
  if (my_wts) delete[] my_wts;
}

