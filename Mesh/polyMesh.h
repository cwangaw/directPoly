#ifndef __polymesh_h_included__
#define __polymesh_h_included__

#include <string>
#include <fstream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// PolyMesh class
//   Includes classes for Edge and PolyElement
//   Vertices are base objects called Point
// Mesh of convex polygons in 2D
//
// Assume base objects: Point and Tensor1 (both two numbers)
////////////////////////////////////////////////////////////////////////////////

#include "baseObjects.h"

namespace polymesh {

  class PolyMesh;

  ////////////////////////////////////////////////////////////////////////////////
  // class Vertex
  ////////////////////////////////////////////////////////////////////////////////

  class Vertex : public base_object::Point
  {
  private:
    PolyMesh* my_mesh;
    int my_mesh_index;

    void set_vertex(double p0, double p1, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      the_point[0] = p0; the_point[1] = p1; my_mesh_index = myIndex; my_mesh = myMesh; };

  public:
    Vertex() { set_vertex(0,0); };
    Vertex(double p0, double p1, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex(p0, p1, myIndex, myMesh); };
    Vertex(const base_object::Point& p, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], myIndex, myMesh); };
    Vertex(const base_object::Point* p, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], myIndex, myMesh); };
    
    void set() { set_vertex(0,0); };
    void set(double p0=0, double p1=0, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex(p0, p1, myIndex, myMesh); };
    void set(const base_object::Point& p, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], myIndex, myMesh); };
    void set(const base_object::Point* p, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], myIndex, myMesh); };

    PolyMesh* mesh() const { return my_mesh; };
    int meshIndex() const { return my_mesh_index; };
     
    void nbrElements(std::vector<int>& theNbrIndices) const;
    void nbrEdges(std::vector<int>& theNbrIndices) const;
    bool isOnBoundary() const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class Edge
  //      X ------------->------- X
  //     v[0]    |    tangent    v[1]
  //             V
  //           normal
  //    Linear function lambda for the edge
  //    Length of the edge
  //    Knows the two adjoining elements (or sets as outside)
  ////////////////////////////////////////////////////////////////////////////////

  class Edge
  {
  private:
    Vertex* the_vertex[2];
    base_object::Tensor1 normal_vector;
    base_object::Tensor1 tangent_vector;
    double edge_length;

    PolyMesh* my_mesh;
    int my_mesh_index;

    void set_edge(Vertex* pt0, Vertex* pt1, int myIndex=-1, PolyMesh* myMesh=nullptr);
    
  public:
    Edge(Vertex* pt0, Vertex* pt1, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_edge(pt0, pt1, myIndex, myMesh); };
    Edge() : normal_vector(), tangent_vector(), edge_length(0), my_mesh(nullptr), my_mesh_index(-1) {
      the_vertex[0] = nullptr; the_vertex[1] = nullptr; };
    
    void set(Vertex* pt0, Vertex* pt1, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_edge(pt0, pt1, myIndex, myMesh); };

    Vertex& vertex(int i) const { return *(the_vertex[i % 2]); }
    Vertex* vertexPtr(int i) const { return the_vertex[i % 2]; }
    double length() const { return edge_length; };

    PolyMesh* mesh() const { return my_mesh; };
    int meshIndex() const { return my_mesh_index; };

    void nbrElements(int theNbrIndices[2]) const;
    bool isOnBoundary() const;
    
    // Decide if a point is in this edge
    bool isInEdge(const base_object::Point& pt) const;

    base_object::Tensor1 normal() const { return normal_vector; };
    double normal(int i) const { return normal_vector[i % 2]; };
    base_object::Tensor1 tangent() const { return tangent_vector; };
    double tangent(int i) const { return tangent_vector[i % 2]; };
    
    double lambda(double x, double y) const;
    double lambda(const base_object::Point& p) const { return lambda(p[0], p[1]); };
    base_object::Tensor1 dLambda() const { return -1*normal_vector; };
    double dLambda(int i) const { return -normal_vector[i % 2]; };
    base_object::Tensor1 curlLambda() const { return tangent_vector; };
    double curlLambda(int i) const { return tangent_vector[i % 2]; };

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    friend class OrientedEdge;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class OrientedEdge
  //   An edge, or its reverse
  ////////////////////////////////////////////////////////////////////////////////

  class OrientedEdge
  {
  private:
    Edge* the_edge;
    double the_orientation; // 1 or -1
    int iOffset; // 0 or 1
    
    void set_oriented_edge(Edge* e, double orient);

  public:
    OrientedEdge(Edge* e, double orient) { set_oriented_edge(e,orient); };
    OrientedEdge() : the_edge(nullptr), the_orientation(0), iOffset(0) {};
    
    void set(Edge* e, double orient) { set_oriented_edge(e,orient); };

    Edge& edge() const { return *the_edge; };
    Edge* edgePtr() const { return the_edge; };
    double orientation() const { return the_orientation; }

    Vertex& vertex(int i) const { return *(the_edge->the_vertex[(i+iOffset) % 2]); }
    Vertex* vertexPtr(int i) const { return the_edge->the_vertex[(i+iOffset) % 2]; }
    double length() const { return the_edge->edge_length; };

    PolyMesh* mesh() const { return the_edge->my_mesh; };
    int meshIndex() const { return the_edge->my_mesh_index; };
    
    void nbrElements(int theNbrIndices[2]);
    bool isOnBoundary() { return the_edge->isOnBoundary(); }

    base_object::Tensor1 normal() const { return the_orientation*the_edge->normal_vector; };
    double normal(int i) const { return the_orientation*the_edge->normal_vector[i % 2]; };
    base_object::Tensor1 tangent() const { return the_orientation*the_edge->tangent_vector; };
    double tangent(int i) const { return the_orientation*the_edge->tangent_vector[i % 2]; };
    
    double lambda(double x, double y) const { return the_orientation*the_edge->lambda(x,y); };
    double lambda(const base_object::Point& p) const {
      return the_orientation*the_edge->lambda(p[0], p[1]); };
    base_object::Tensor1 dLambda() const { return -1*(the_orientation*the_edge->normal_vector); };
    double dLambda(int i) const { return -the_orientation*the_edge->normal_vector[i % 2]; };
    base_object::Tensor1 curLambda() const { return the_orientation*the_edge->tangent_vector; };
    double curlLambda(int i) const { return the_orientation*the_edge->tangent_vector[i % 2]; };
    
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class PolyElement
  //       1      e1       5
  //        X ------------X      nVertices = number of vertices or edges
  //    e2 /              |      Composed of nVertices edges
  //     /                | e5   Convex, and ordered counterclockwise
  //  2 X                 |
  //   e3 \               |      The edges are stored as pointers to en external
  //       X ------------ X         data structure (the PolyMesh).
  //      3       e4       4
  //   Area of element
  ////////////////////////////////////////////////////////////////////////////////

  class PolyElement
  {
  private:
    int num_vertices;
    OrientedEdge* the_oriented_edge = nullptr;
    bool good_element;
    
    double the_area;
    base_object::Point my_center;  // of largest inscribed circle
    base_object::Point my_centroid; // center of mass
    double my_diameter;
    double max_inscribed_radius;

    PolyMesh* my_mesh;
    int my_mesh_index;
    
    bool isConnected();
    bool isConvexCounterclockwise();
    double computeAreaTriangle(const base_object::Point* p0, const base_object::Point* p1,
			       const base_object::Point* p2);
    double computeArea();
    void computeCenterOfTriangle(const OrientedEdge* e0, const OrientedEdge* e1,
				 const OrientedEdge* e2, base_object::Point& center,
				 bool& unique) const;

    void set_polyelement(int nGon, Edge** theEdges, int myIndex=-1, PolyMesh* myMesh=nullptr);

  public:
    PolyElement(int nGon, Edge** theEdges, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_polyelement(nGon, theEdges, myIndex, myMesh); };
    PolyElement() : num_vertices(0), the_oriented_edge(nullptr), good_element(false),
		    the_area(0), my_center(), my_centroid(), my_diameter(0),
		    my_mesh(nullptr), my_mesh_index(-1) {};
    ~PolyElement();

    void set(int nGon, Edge** theEdges, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_polyelement(nGon, theEdges, myIndex, myMesh); };
    
    bool isGood() const { return good_element; };
    int nVertices() const { return num_vertices; };

    int longestEdgeIndex() const; // Return the index of longest edge 
    int nShortEdges(double ratio) const; // ratio is the criteria of defining small edges
    
    OrientedEdge& edge(int i) const { return the_oriented_edge[i % num_vertices]; };
    OrientedEdge* edgePtr(int i) const { return &the_oriented_edge[i % num_vertices]; };
    Vertex& vertex(int i) const { return the_oriented_edge[i % num_vertices].vertex(1); };
    Vertex* vertexPtr(int i) const { return the_oriented_edge[i % num_vertices].vertexPtr(1); };

    PolyMesh* mesh() const { return my_mesh; };
    int meshIndex() const { return my_mesh_index; };
    
    double area() const { return the_area; };
    base_object::Point center()   const { return my_center; };
    base_object::Point centroid() const { return my_centroid; };
    double maxInscribedRadius() const { return max_inscribed_radius; }
    double diameter() const { return my_diameter; };
    double chunkinessParam() const;

    bool isInElement(const base_object::Point& pt) const;
    bool isOnElementBoundary(const base_object::Point& pt) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    friend class PolyMesh;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class PolyMesh
  //   Lists of vertices, edges, and elements
  // Construction
  //    On input there are two arrays:
  //      pts: The points or vertices (x1, y1, x2, y2; ... ).
  //           Gives an indexing of the vertices
  //      elementByVerticesIndex: List of elements, by integer numbers for the vertex
  //             index (p1, p2, p3, p7; p1, p4, p5; ...)
  //           Gives an indexing of the elements
  //    Create the intermediate array
  //      edges_by_vertex_and_elements_index: List of edges
  //             (pt1-1, pt1-2, element1-1, element1-2; ...
  //           Gives an indexing of the edges
  //           Gives the neighboring elements (and whether on the boundary).
  ////////////////////////////////////////////////////////////////////////////////

  class PolyMesh
  {
  public:
    
  private:
    typedef int int2[2];

    int num_vertices;
    int num_edges;
    int num_elements;
    int max_num_edges_of_elements;
    bool mesh_is_good;
    
    double min_x, max_x, min_y, max_y;
    int num_boundary_vertices;
    int num_boundary_edges;
    bool mesh_is_periodic;

    // Mesh quality
    double max_element_diameter;
    double max_chunkiness_parameter;
    double min_chunkiness_parameter;
    double avg_chunkiness_parameter;

    // Polymesh of elements
    Vertex* the_vertices = nullptr; // no order
    Edge* the_edges = nullptr; // no order
    PolyElement* the_elements = nullptr; // no order
    
    // Original mesh (input)
    int* num_vertices_of_element = nullptr;
    int** element_by_vertex_index = nullptr;
    
    // Mesh connectivity
    std::vector<int>* nbr_edges_of_vertex = nullptr; // ordered by edge #
    std::vector<int>* nbr_elements_of_vertex = nullptr; // ordered by element #
    int2* nbr_elements_of_edge = nullptr; // ordered by vertex orientation, LR

    int* vertex_periodic_equiv = nullptr; // -1 if corner, -2 if interior
    int* vertex_periodic_corners = nullptr; // 4 corners
    int* edge_periodic_equiv = nullptr; // -1 if interior

    // Mesh setup
    void set_polymesh(int numVertices, double* pts, int numElements,
	     int* numEdgesOfElement, int** elementByVertexIndex);

  public:
    PolyMesh(int numVertices, double* pts, int numElements,
	     int* numEdgesOfElement, int** elementByVertexIndex) {
      set_polymesh(numVertices, pts, numElements,
		   numEdgesOfElement, elementByVertexIndex); };
    PolyMesh() : num_vertices(0), num_edges(0), num_elements(0),
		 max_num_edges_of_elements(0), mesh_is_good(false),
		 min_x(0), max_x(1), min_y(0), max_y(1),
		 num_boundary_vertices(0), num_boundary_edges(0), mesh_is_periodic(false),
		 max_element_diameter(0), max_chunkiness_parameter(0),
		 min_chunkiness_parameter(0), avg_chunkiness_parameter(0),
		 the_vertices(nullptr), the_edges(nullptr), the_elements(nullptr),
		 num_vertices_of_element(nullptr), element_by_vertex_index(nullptr),
		 nbr_edges_of_vertex(nullptr), nbr_elements_of_vertex(nullptr),
		 nbr_elements_of_edge(nullptr) {};
    PolyMesh(PolyElement* single_element); // Create a one element mesh from an element
    ~PolyMesh();

    void set(int numVertices, double* pts, int numElements,
	     int* numEdgesOfElement, int** elementByVertexIndex) {
      set_polymesh(numVertices, pts, numElements,
		   numEdgesOfElement, elementByVertexIndex); };

    int createMesh(char meshType='q', int nx=1, int ny=1, double xMin=0, double xMax=1,
		   double yMin=0, double yMax=1, double distortionFactor=0,
		   bool allowPeriodic=false); // Create simple mesh

    int removeShortEdges(double ratio); // Remove short mesh edges, return # edges removed
    int removeShortEdges(double ratio, std::vector<int> verticesRemoved);

    int makePeriodic(); // return error if not possible
    void removePeriodicity();

    // Access functions
    int nVertices() const { return num_vertices; };
    int nEdges() const { return num_edges; };
    int nElements() const { return num_elements; };
    int maxNEdgesOfElements() const { return max_num_edges_of_elements; };
    int maxNVerticesOfElements() const { return max_num_edges_of_elements; };
    double minX() const { return min_x; }
    double maxX() const { return max_x; }
    double minY() const { return min_y; }
    double maxY() const { return max_y; }
    bool isPeriodic() const { return mesh_is_periodic; }
    
    double maxElementDiameter() const { return max_element_diameter; };
    double maxChunkinessParam() const { return max_chunkiness_parameter; }
    double minChunkinessParam() const { return min_chunkiness_parameter; }
    double avgChunkinessParam() const { return avg_chunkiness_parameter; }
    
    bool isGood() const { return mesh_is_good; };
    
    PolyElement& element(int i) const { return the_elements[i % num_elements]; };
    PolyElement* elementPtr(int i) const { return &the_elements[i % num_elements]; };
    Edge& edge(int i) const { return the_edges[i % num_edges]; };
    Edge* edgePtr(int i) const { return &the_edges[i % num_edges]; };
    Vertex& vertex(int i) const { return the_vertices[i % num_vertices]; };
    Vertex* vertexPtr(int i) const { return &the_vertices[i % num_vertices]; };
    
    // Mesh connectivity
    int nVerticesOfElement(int i) const { return num_vertices_of_element[i]; }
    
    void nbrEdgesOfVertex(int i, std::vector<int>& theNbrIndices) const {
      theNbrIndices = nbr_edges_of_vertex[i]; }
    void nbrElementsOfVertex(int i, std::vector<int>& theNbrIndices) const {
      theNbrIndices = nbr_elements_of_vertex[i]; }
    
    void nbrElementsOfEdge(int i, int theNbrIndices[2]) const {
      theNbrIndices[0] = nbr_elements_of_edge[i][0];
      theNbrIndices[1] = nbr_elements_of_edge[i][1]; }

    int nShortEdges(double ratio); // ratio is the criteria of defining small edges
    
    bool isVertexOnBoundary(int i) const;
    bool isEdgeOnBoundary(int i) const;

    // Functions related to periodicity
    int nVerticesPeriodic() const { return (num_boundary_vertices - 4)/2 + 1; }
    int nEdgesPeriodic() const { return num_boundary_edges/2; }

    void pointOnBoundingBox(const base_object::Point& p,
			    bool& onL, bool& onR, bool& onB, bool& onT);
    void edgeOnBoundingBox(const Edge& e, bool& onL, bool& onR, bool& onB, bool& onT);
      
    int vertexPeriodicEquiv(int i) const {
      return mesh_is_periodic ? vertex_periodic_equiv[i] : -2; } // -1 corner, -2 interior
    int vertexPeriodicCorners(int c) const {
      return mesh_is_periodic ? vertex_periodic_corners[c % 4] : -1; } // order LB,RB,LT,RT
    int edgePeriodicEquiv(int i) const {
      return mesh_is_periodic ? edge_periodic_equiv[i] : -1; }

    // Points in the Mesh

    int elementIndex(const base_object::Point& pt) const;
    
    // Output
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    int write_matlab(std::string& filename) const;

    friend class PolyElement;
    friend class OrientedEdge;
    friend class Edge;
    friend class Vertex;
  };

};

#endif
