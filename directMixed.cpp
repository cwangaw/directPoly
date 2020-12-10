#include <cmath>
#include <iostream>
#include <vector>

#include <stdio.h>
#include <assert.h>


using namespace std;

#include "Utilities/debug.h"
#include "Mesh/polyMesh.h"
using namespace polymesh;
#include "directMixed.h"
using namespace directserendipity;
#include "polyQuadrature.h"
using namespace polyquadrature;


////////////////////////////////////////////////////////////////////////////////
// Class DirectMixed

void DirectMixed::set_directmixed(int polyDeg, PolyMesh* mesh) {
  polynomial_degree = polyDeg;
  my_mesh = mesh;
  num_edges = my_mesh -> nEdges();

  int index = 0;
  for (int i = 0; i < nEdges(); i++) {
    if (!my_mesh->edgePtr(i)->isOnBoundary()) index++;
  }
  num_interior_edges = index;

  // ALLOCATE ELEMENTS

  // Allocate
  
  if(the_dm_edges) delete[] the_dm_edges;
  the_dm_edges = new Edge[num_edges];
  
  if(the_bc_type) delete[] the_bc_type;
  the_bc_type = new EdgeBCType[num_edges];

  if(the_dm_elements) delete[] the_dm_elements;
  the_dm_elements = new DirectMixedFE[my_mesh->nElements()];

  if(the_dg_elements) delete[] the_dg_elements;
  the_dg_elements = new DirectDGFE[my_mesh->nElements()];

  if(the_dg_edge_elements) delete[] the_dg_edge_elements;
  the_dg_edge_elements = new DirectEdgeDGFE[num_edges];

  if(interior_edge_indexing) delete[] interior_edge_indexing;
  interior_edge_indexing = new int[num_edges];

  if(global_edge_indexing) delete[] global_edge_indexing;
  global_edge_indexing = new int[num_interior_edges];

  // DETERMINE EDGES (ORDERING AND TYPES)

  // interior_edge_indexing maps from global edge index to interior edge index
  // global_edge_indexing maps from interior edge index to global edge index
  index = 0;
  for (int i = 0; i < nEdges(); i++) {
    if (my_mesh->edgePtr(i)->isOnBoundary()) {
      the_bc_type[i] = EdgeBCType::boundary;
      interior_edge_indexing[i] = -1;
    } else {
      the_bc_type[i] = EdgeBCType::interior;
      interior_edge_indexing[i] = index;
      global_edge_indexing[index] = i;
      index++;
    }
  }


  // ALLOCATE FINITE ELEMENTS

  for(int iElement=0; iElement<my_mesh->nElements(); iElement++) {
    std::vector<int> edge_index;
    edge_index.clear();
    
    PolyElement* element = my_mesh->elementPtr(iElement);
    int nGon = element->nVertices();

    for (int i = 0; i < nGon; i++) {
      int interior_index = interior_edge_indexing[element->edgePtr(i)->meshIndex()];
      edge_index.push_back(interior_index);
    }

    the_dm_elements[iElement].set(this, element);
    the_dg_elements[iElement].set(this, element);
  }

  for (int iEdge = 0; iEdge < num_edges; iEdge++ ) {
    Edge* edge = my_mesh->edgePtr(iEdge);
    the_dg_edge_elements[iEdge].set(this,edge);
  }

};

DirectMixed::~DirectMixed() {
    if (the_dm_edges) delete[] the_dm_edges;
    if (the_bc_type) delete[] the_bc_type;
    if (interior_edge_indexing) delete[] interior_edge_indexing;
    if (global_edge_indexing) delete[] global_edge_indexing;
    if (the_dm_elements) delete[] the_dm_elements;
    if (the_dg_elements) delete[] the_dg_elements;
    if (the_dg_edge_elements) delete[] the_dg_edge_elements;
}
////////////////////////////////////////////////////////////////////////////////
// class DirectMixedArray

void DirectMixedArray::set_directmixedarray(DirectMixed* dmSpace) {
  my_dm_space = dmSpace;
  num_edges = my_dm_space->nEdges();
  
  if(the_array) delete[] the_array;
  the_array = new double[num_edges];
}

DirectMixedArray::~DirectMixedArray() {
  if(the_array) delete[] the_array;
}