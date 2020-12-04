#ifndef __directserendipity_h_included__
#define __directserendipity_h_included__

////////////////////////////////////////////////////////////////////////////////
// DirectSerendipity class
//   Includes classes for DirectSerendipityFE
// Direct Serendipity Elements on convex polygons
//
// Uses the PolyMesh classes for the mesh and elements
// Assume base objects: Point and Tensor1 (both two numbers)
////////////////////////////////////////////////////////////////////////////////

#include "Mesh/baseObjects.h"
#include "Mesh/polyMesh.h"
#include "fcns.h"

namespace directserendipity {

  enum class NodeType { vertex, edge, cell };
  enum class BCType { interior, dirichlet, neumann, robin };
  enum class EdgeBCType { interior, boundary };
  class DirectSerendipity;
  class DirectMixed;
 
  ////////////////////////////////////////////////////////////////////////////////
  // class Node
  //   Stores node type and BC type
  ////////////////////////////////////////////////////////////////////////////////

  class Node : public Point
  {
  private:
    DirectSerendipity* my_ds_space;
    int my_node_index;
    int my_mesh_index; // to vertex, edge, or element index (depending on type)
    
    void set_node(double p0, double p1, int myMeshIndex=-1, int myNodeIndex=-1,
		  DirectSerendipity* myDSSpace=nullptr) {
      the_point[0] = p0; the_point[1] = p1; my_mesh_index = myMeshIndex;
      my_node_index = myNodeIndex; my_ds_space = myDSSpace; };

  public:
    Node() { set_node(0,0); };
    Node(double p0, double p1, int myMeshIndex=-1, int myNodeIndex=-1,
	 DirectSerendipity* myDSSpace=nullptr) {
      set_node(p0, p1, myMeshIndex, myNodeIndex, myDSSpace); };
    Node(const Point& p, int myMeshIndex=-1, int myNodeIndex=-1,
	 DirectSerendipity* myDSSpace=nullptr) {
      set_node(p[0], p[1], myMeshIndex, myNodeIndex, myDSSpace); };
    
    void set() { set_node(0,0); };
    void set(double p0=0, double p1=0, int myMeshIndex=-1, int myNodeIndex=-1,
	     DirectSerendipity* myDSSpace=nullptr) {
      set_node(p0, p1, myMeshIndex, myNodeIndex, myDSSpace); };
    void set(const Point& p, int myMeshIndex=-1, int myNodeIndex=-1,
	     DirectSerendipity* myDSSpace=nullptr) {
      set_node(p[0], p[1], myMeshIndex, myNodeIndex, myDSSpace); };
    void set(const Point* p, int myMeshIndex=-1, int myNodeIndex=-1,
	     DirectSerendipity* myDSSpace=nullptr) {
      set_node((*p)[0], (*p)[1], myMeshIndex, myNodeIndex, myDSSpace); };

    DirectSerendipity* dsSpace() const { return my_ds_space; };
    int nodeIndex() const { return my_node_index; };

    int meshIndex() const { return my_mesh_index; };
    int meshVertexIndex() const { return (nodeType() == NodeType::vertex) ? my_node_index : -1; };
    int meshEdgeIndex() const   { return (nodeType() == NodeType::edge) ? my_node_index : -1; };
    int meshElementIndex() const {return (nodeType() == NodeType::cell) ? my_node_index : -1; };

    NodeType nodeType() const;
    BCType bcType() const;

    bool isVertex() const { return nodeType() == NodeType::vertex; };
    bool isEdge()   const { return nodeType() == NodeType::edge; };
    bool isCell()   const { return nodeType() == NodeType::cell; };
      
    bool isInterior() const { return bcType() == BCType::interior; };
    bool isDirchlet() const { return bcType() == BCType::dirichlet; };
    bool isNeumann()  const { return bcType() == BCType::neumann; };
    bool isRobin()    const { return bcType() == BCType::robin; };
      
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  // class DirectSerendipityArray
  //    Gives values for each nodal point
  //    Also gives the gradient
  ////////////////////////////////////////////////////////////////////////////////

  class DirectSerendipityArray
  {
  private:
    int num_nodes;
    double* the_array = nullptr;

    DirectSerendipity* my_ds_space;

    void set_directserendipityarray(DirectSerendipity* dsSpace);

  public:
    DirectSerendipityArray() : num_nodes(0), the_array(nullptr), my_ds_space(nullptr) {};
    DirectSerendipityArray(DirectSerendipity* dsSpace) {
      set_directserendipityarray(dsSpace); };
    DirectSerendipityArray(const DirectSerendipityArray& a) : the_array(nullptr) {
      set_directserendipityarray(a.dsSpace()); };
    ~DirectSerendipityArray();
    
    void set(DirectSerendipity* dsSpace) { set_directserendipityarray(dsSpace); };

    DirectSerendipity* dsSpace() const { return my_ds_space; };
    int size() const { return num_nodes; };

    double& operator() (int i)       { return the_array[i]; }
    double  operator() (int i) const { return the_array[i]; }
    double& operator[] (int i)       { return the_array[i]; }
    double  operator[] (int i) const { return the_array[i]; }

    void eval(const Point* pts, double* result, Tensor1* gradResult, int num_pts) const;
    void eval(const Point& pt, double& result, Tensor1& gradResult) const;
    double eval(const Point& pt) const;

    void l2normError(double& l2Error, double& l2GradError, double& l2Norm, double& l2GradNorm,
		     double (*referenceFcn)(double,double) = nullptr, 
		     Tensor1 (*referenceGradFcn)(double,double) = nullptr);
    void l2norm(double& l2Norm, double& l2GradNorm) {
      double null1,null2; l2normError(l2Norm,l2GradNorm,null1,null2); };

    void write_matlab_mesh(std::ofstream* fout, std::ofstream* fout_grad,
			   int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh(std::ofstream& fout, std::ofstream& fout_grad,
			   int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, &fout_grad, num_pts_x, num_pts_y); };
    void write_matlab_mesh(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, nullptr, num_pts_x, num_pts_y); };
    int write_matlab_mesh(std::string& filename, std::string& filename_grad,
			  int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_matlab_mesh_by_pt(std::ofstream& fout, std::ofstream& fout_grad,
				 int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, std::string& filename_grad,
				int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  // class DirectSerendipityFE
  //    Defined on a poly-element (class PolyElement)
  //    
  //    Gives nodal basis functions
  //      First call initBasis to evaluate all basis functions at a given set of points
  //      Then access the basis functions by node number, or type and number, and pt number
  //
  // Store the result in such a way that the evaluations of all nodal basis functions at
  // one point are put together. For example, if r = 3, N = 3, and we have two points,
  // then value_n = (Phi_1_1 at Point1, Phi_1_2 at Point 1, Phi_2_1 at Point 1, ...,
  //                  Phi_3_1 at Point 2,  Phi_3_2 at Point 2).
  ////////////////////////////////////////////////////////////////////////////////

  class DirectSerendipityFE 
  {
  private:
    DirectSerendipity* my_ds_space;
    polymesh::PolyElement* my_poly_element;

    int num_vertices; // Redundant with my_poly_element
    int polynomial_degree; // Redundant with my_ds_space
    int deg_cell_poly; // Redundant with above
    int num_cell_nodes; // Redundant with above
    int num_nodes; // Redundant with above

    std::vector<Node*> the_vertex_nodes; // CCW ordered
    std::vector<Node*> the_edge_nodes; // grouped by edge (i) to corresponding vertex (i)
    std::vector<Node*> the_cell_nodes; // ordered lower left to right, bottom to top

    // Evaluation storage
    int num_eval_pts;
    double* value_n = nullptr;
    Tensor1* gradvalue_n = nullptr;

    // If necessary (degPolyn small), the bigger space within which we construct the basis
    polymesh::PolyMesh* one_element_mesh = nullptr;
    DirectSerendipity* high_order_ds_space = nullptr;

    void set_directserendipityfe(DirectSerendipity* dsSpace, polymesh::PolyElement* element,
				 const std::vector<Node*>& vertexNodes,
				 const std::vector<Node*>& edgeNodes,
				 const std::vector<Node*>& cellNodes);

    // Private helper functions
    
    void projToEdge(int iEdge, double x, double y, double& result, Tensor1& gradresult) const;
    void projToEdge(int iEdge, const Point& p, double& result, Tensor1& gradresult) const { projToEdge(iEdge,p[0],p[1],result,gradresult); }
    double projToEdge(int iEdge, double x, double y) const;
    double projToEdge(int iEdge, const Point& p) const { return projToEdge(iEdge,p[0],p[1]); }

    //r_supp is 1 on e_i, 0 on e_j
    void r_supp(int i, int j, const Point& p, double& result, Tensor1& gradresult) const;
    double r_supp(int i, int j, const Point& p) const;
    double r_supp(int i, int j, double x, double y) const { return r_supp(i, j, Point(x,y)); }  


    void lambda_supp (int i, int j, const Point& p, double& result, Tensor1& gradresult) const;
    double lambda_supp (int i, int j, const Point& p) const;
    double lambda_supp(int i, int j, double x, double y) const {
      return lambda_supp(i,j,Point(x,y)); }

    void phi_k_l (int k, int l, const Point& p, double& result, Tensor1& gradresult) const;

    double phi_k_l (int k, int l, const Point& p) const;
    double phi_k_l (int i, int j, double x, double y) const { return phi_k_l(i, j, Point(x,y)); }

    // Reconcile triangular indexing
    int mapI(int k, int deg_cell_poly) const;
    int mapJ(int k, int deg_cell_poly) const;
    int mapK(int i, int j, int deg_cell_poly) const;


    //Matrix Computation Result for edge nodes
    std::vector<double> x_vector_for_all;
    const double* get_result_of_coefficients(int i, int j) const;

    //Lagrange polynomial that is 1 on x_nEdge_jNode,
    //    0 on other nodes including two vertices on edge_nEdge
    void lagrange(int nEdge, int jNode, const Point& pt, double& result, Tensor1& gradresult) const;
    double lagrange(int nEdge, int jNode, const Point& pt) const;

    //Lagrange polynomial that is 1 on x_{v,i}, 0 on other nodes including vertex on edge_nEdge
    void lagrange_v(int i, int nEdge, const Point& pt, double& result, Tensor1& gradresult) const;
    double lagrange_v(int i, int nEdge, const Point& pt) const;

    //varphi before normalization for small r and edge nodes in A_po
    void varphi_edge_nodes_A_Po(int nEdge, int jNode, const Point& pt, double& result, Tensor1& gradresult) const;
    double varphi_edge_nodes_A_Po(int nEdge, int jNode, const Point& pt) const;

    //varphi before normalization for small r and vertex nodes in A_po
    void varphi_vertex_nodes_A_Po(int n, const Point& pt, double& result, Tensor1& gradresult) const;
    double varphi_vertex_nodes_A_Po(int n, const Point& pt) const;

    // Helper data
    //std::vector<double> lambda_at_pt;
    //Point* usual_order_basis_eval_pts = nullptr;
    //Point* high_order_basis_eval_pts = nullptr;
    //Point* basis_eval_pts = nullptr;
    //int num_basis_eval_pts;
    
  public:
    DirectSerendipityFE() : my_ds_space(nullptr), my_poly_element(nullptr), num_vertices(0),
			    polynomial_degree(0), deg_cell_poly(-1), num_cell_nodes(0), num_nodes(0) {};
    DirectSerendipityFE(DirectSerendipity* dsSpace, polymesh::PolyElement* element,
			std::vector<Node*> vertexNodes, std::vector<Node*> edgeNodes,
			std::vector<Node*> cellNodes) {
      set_directserendipityfe(dsSpace, element, vertexNodes, edgeNodes, cellNodes); };
    ~DirectSerendipityFE();
      
    void set(DirectSerendipity* dsSpace, polymesh::PolyElement* element,
	     std::vector<Node*> vertexNodes, std::vector<Node*> edgeNodes,
	     std::vector<Node*> cellNodes) {
      set_directserendipityfe(dsSpace, element, vertexNodes, edgeNodes, cellNodes); };

    // Access functions
    DirectSerendipity* dsSpace() const { return my_ds_space; };
    polymesh::PolyElement* elementPtr() const { return my_poly_element; };
    int nVertices() const { return num_vertices; }
    int degPolyn() const { return polynomial_degree; }
    int degCellPolyn() const { return deg_cell_poly; }
    
    int nVertexNodes() const { return num_vertices; }
    int nEdgeNodes() const { return num_vertices*(polynomial_degree-1); }
    int nCellNodes() const { return num_cell_nodes; }
    int nNodes() const { return num_nodes; }

    Node* vertexNodePtr(int i) const { return the_vertex_nodes[i]; };
    Node* edgeNodePtr(int k)   const { return the_edge_nodes[k];   };
    Node* edgeNodePtr(int i, int j) const { return the_edge_nodes[i*(polynomial_degree - 1) + j]; };
    Node* cellNodePtr(int k)   const { return the_cell_nodes[k];   };
    Node* cellNodePtr(int i, int j) const { return the_cell_nodes[mapK(i, j, deg_cell_poly)]; };
    Node* nodePtr(int i) const;

    // FE basis functions
    double lambda(int iEdge, const Point& pt) const {
      return my_poly_element->edgePtr(iEdge)->lambda(pt); }
    Tensor1 dLambda(int iEdge) const { return my_poly_element->edgePtr(iEdge)->dLambda(); }


    /*
    void vertexBasis(double* result, int num_pts) const;
    void vertexBasisGrad(Tensor1* result, int num_pts) const;
    */
    void initBasis(const Point* pt, int num_pts); // (re)evaluate all basis fcns at the points

    double basis(int iNode, int iPt) const { return value_n[iNode + iPt*num_nodes]; };
    Tensor1 basisGrad(int iNode, int iPt) const { return gradvalue_n[iNode + iPt*num_nodes]; };
    
    double vertexBasis(int iVNode, int iPt) const { return value_n[iVNode + iPt*num_nodes]; };
    Tensor1 gradVertexBasis(int iVNode, int iPt) const { return gradvalue_n[iVNode + iPt*num_nodes]; };


    double edgeBasis(int iENode, int iPt) const {
      return value_n[iENode + num_vertices + iPt*num_nodes]; };
    Tensor1 gradEdgeBasis(int iENode, int iPt) const { return gradvalue_n[iENode + num_vertices + iPt*num_nodes]; }
    double edgeBasis(int iEdge, int j, int iPt) const  {
      return edgeBasis( iEdge*(polynomial_degree-1) + j, iPt); }
    Tensor1 gradEdgeBasis(int iEdge, int j, int iPt) const  {
      return gradEdgeBasis( iEdge*(polynomial_degree-1) + j, iPt ); }
    

    double cellBasis(int iCNode, int iPt) const {
      return value_n[iCNode + num_vertices*polynomial_degree + iPt*num_nodes]; };
    Tensor1 gradCellBasis(int iCNode, int iPt) const { return gradvalue_n[iCNode + num_vertices*polynomial_degree + iPt*num_nodes]; };  
    
    //mode:0,1. If mode = 0, we defaultly use nodal basis functions; if mode = 1, we use shape functions 
    void eval(const Point* pt, double* result, Tensor1* gradResult, int num_pts,
	      double* vertex_dofs, double* edge_dofs=nullptr, double* cell_dofs=nullptr);
    void eval(const Point& pt, double& result, Tensor1& gradResult, double* vertex_dofs,
              double* edge_dofs=nullptr, double* cell_dofs=nullptr) {
      eval(&pt, &result, &gradResult, 1, vertex_dofs, edge_dofs, cell_dofs); };

    void eval(const Point* pt, double* result, int num_pts, double* vertex_dofs,
              double* edge_dofs=nullptr, double* cell_dofs=nullptr);
    double eval(const Point& pt, double* vertex_dofs,
		double* edge_dofs=nullptr, double* cell_dofs=nullptr) { double result;
      eval(&pt, &result, 1, vertex_dofs, edge_dofs, cell_dofs); return result; }

    // Output functions
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    void write_matlab_vertices(std::ofstream& fout, double* vertex_dofs) const;
    int write_matlab_vertices(std::string& filename, double* vertex_dofs) const;

    friend class DirectMixedFE;
  };

  void lambda_for_two_points(const Point& pt0, const Point& pt1, double x, double y, double& result, Tensor1& gradresult);

  inline void lambda_for_two_points(const Point& pt0, const Point& pt1, const Point& p, double& result, Tensor1& gradresult) {
    lambda_for_two_points(pt0, pt1, p[0], p[1],result,gradresult);
  };
  double lambda_for_two_points(const Point& pt0, const Point& pt1, double x, double y);

  inline double lambda_for_two_points(const Point& pt0, const Point& pt1, const Point& p) {
    return lambda_for_two_points(pt0, pt1, p[0], p[1]);


  };

 
  ////////////////////////////////////////////////////////////////////////////////
  // class DirectMixedFE
  //    Defined on a poly-element (class PolyElement)
  //    
  //    Gives basis functions separated into curl(\cDS_{r+1}) part and \x\Po_s(E) part
  //      First call initBasis to evaluate all basis functions at a given set of points
  //      Then access the basis functions by function index and pt number
  //
  // Store the result in such a way that the evaluations of all basis functions at
  // one point are put together. The basis functions are stored in the following order:
  //  curl(\cDS_{r+1}) -> \x\Po_{r-1}(E) -> \x\tilde\Po_r(E) ( polynomials of order r only )
  ////////////////////////////////////////////////////////////////////////////////

  class DirectMixedFE
  {
  private:
    DirectMixed* my_dm_space;
    polymesh::PolyElement* my_poly_element;

    int dim_v;
    int dim_curlpart;
    int dim_v_div;

    int num_vertices; // Redundant with my_poly_element
    int polynomial_degree; // Redundant with my_ds_space

    // If necessary (degPolyn small), the bigger space within which we construct the basis
    polymesh::PolyMesh* one_element_mesh = nullptr;
    DirectSerendipity* high_order_ds_space = nullptr;
    
    // Pointer to Evaluation storage
    int num_eval_pts;
    Tensor1* v_value_n = nullptr;

    double* v_div_value_n = nullptr;
    double* v_edge_value_n = nullptr;


    double* w_value_n = nullptr;
    double* lag_value_n = nullptr;

    void set_directmixedfe(DirectMixed* dmSpace, polymesh::PolyElement* element);

  public:
    DirectMixedFE() : my_dm_space(nullptr), my_poly_element(nullptr) {};
    DirectMixedFE( DirectMixed* dmSpace, polymesh::PolyElement* element) { 
      set_directmixedfe(dmSpace, element); };
    
    ~DirectMixedFE();

    void set(DirectMixed* dmSpace, polymesh::PolyElement* element) {
      set_directmixedfe(dmSpace, element); };

    void initBasis(const Point* pt, int num_pts); // (re)evaluate all basis fcns at the points

    polymesh::PolyElement* elementPtr() const { return my_poly_element; };

    // Access basis functions evaluated at pt[iPt]
    Tensor1 basis(int iFunc, int iPt) const {
      return v_value_n[iPt * dim_v + iFunc];
    };

    double basisDotNu(int iFunc, int nEdge, int iPt) const {
      return v_edge_value_n[iPt * dim_v * num_vertices + dim_v * nEdge + iFunc];
    }

    Tensor1 xPo(int iFunc, int iPt) const {
      return v_value_n[iPt * dim_v + dim_curlpart + iFunc];
    };

    double divXPo(int iFunc, int iPt) const {
      return v_div_value_n[iPt * dim_v_div + iFunc];
    }

    // Get dimensions of spaces
    int dimVFull() { return dim_v; };
    int dimVReduced() { return dim_v - polynomial_degree - 1; };

    int dimCurlPart() { return dim_curlpart; }

    int dimXPoFull() { return dim_v_div; };
    int dimXPoReduced() { return polynomial_degree * (polynomial_degree + 1)/2; };
    
  };

  
  ////////////////////////////////////////////////////////////////////////////////
  // class DirectDGFE
  //    Defined on a poly-element (class PolyElement)
  //    
  //    Gives polynomials of order up to s on the element E, s = r-1, r
  //    Serve as paired space of mixed space that approximate scalar functions
  //      First call initBasis to evaluate all basis functions at a given set of points
  //      Then access the basis functions by function index and pt number
  //
  // Store the result in such a way that the evaluations of all basis functions at
  // one point are put together. The basis functions are stored in the following order:
  //  \Po_{r-1}(E) -> \tilde\Po_r(E) ( polynomials of order r only )
  ////////////////////////////////////////////////////////////////////////////////

    class DirectDGFE
  {
  private:
    DirectMixed* my_dm_space;
    polymesh::PolyElement* my_poly_element;

    int dim_w;

    int num_vertices; // Redundant with my_poly_element
    int polynomial_degree; // Redundant with my_ds_space
    
    // Pointer to Evaluation storage
    int num_eval_pts;
    double* value_n = nullptr;

    void set_directdgfe(DirectMixed* dmSpace, polymesh::PolyElement* element);

  public:
    DirectDGFE() : my_dm_space(nullptr), my_poly_element(nullptr) {};
    DirectDGFE( DirectMixed* dmSpace, polymesh::PolyElement* element ) { 
      set_directdgfe(dmSpace, element); };
    
    ~DirectDGFE();

    void set(DirectMixed* dmSpace, polymesh::PolyElement* element) {
      set_directdgfe(dmSpace, element); };

    void initBasis(const Point* pt, int num_pts); // (re)evaluate all basis fcns at the points

    // Access basis functions evaluated at pt[iPt]
    double basis(int iFunc, int iPt) const {
      return value_n[iPt * dim_w + iFunc];
    };

    // Get dimensions of spaces
    int dimFull() { return dim_w; };
    int dimReduced() { return (polynomial_degree+1) * polynomial_degree /2;}
  };


  ////////////////////////////////////////////////////////////////////////////////
  // class DirectEdgeDGFE
  //
  //    Defined on an interior edge
  //    
  //    Gives polynomials of order up to r on each interior edge of the global mesh
  //
  //    Serve as Lagrange multipliers of mixed space
  //      First call initBasis to evaluate all basis functions at a given set of points
  //      Then access the basis functions by function index and pt number
  //
  // Store the result in such a way that the evaluations of all basis functions at
  // one point are put together. 
  ////////////////////////////////////////////////////////////////////////////////

  class DirectEdgeDGFE
  {
  private:
    DirectMixed* my_dm_space;
    polymesh::Edge* my_edge;

    int dim_l;

    int polynomial_degree; // Redundant with my_dm_space
    
    // Private helper function
    double projToEdge(const Point& p) const;

    // Pointer to Evaluation storage
    int num_eval_pts;
    double* value_n = nullptr;

    void set_directedgedgfe(DirectMixed* dmSpace, polymesh::Edge* edge);

  public:
    DirectEdgeDGFE() : my_dm_space(nullptr), my_edge(nullptr) {};
    DirectEdgeDGFE( DirectMixed* dmSpace, polymesh::Edge* edge ) { 
      set_directedgedgfe(dmSpace, edge); };
    
    ~DirectEdgeDGFE();

    void set(DirectMixed* dmSpace, polymesh::Edge* edge) { set_directedgedgfe(dmSpace, edge); };

    void initBasis(const Point* pt, int num_pts); // (re)evaluate all basis fcns at the points

    // Access basis functions evaluated at pt[iPt]
    double basis(int iFunc, int iPt) const {
      return value_n[iPt * dim_l + iFunc];
    };

    // Get dimensions of spaces
    int dim() { return dim_l; };
  };



  ////////////////////////////////////////////////////////////////////////////////
  // class DirectSerendipity
  //   Defined on a poly-mesh (class PolyMesh), consisting of direct serendipity
  //     finite elements (class DirectSerendipityFE)
  //    
  //   Gives function evaluation from global DoFs for the basis functions 
  //
  // NODES
  // Local numbering:
  //   Local node number is counted from vertex 0 (i.e., mesh edge 0, v[1])
  //   CCW counting all edge and vertex nodes on the boundary, followed
  //     by the interior nodes (if any).
  //   Local numbering is also given by node type
  //
  //             edge0
  //     V0_0--E01_8--E00_7--V2_6
  //         \              /
  //       E10_1   C0_9   E21_5
  //            \        /
  //          E11_2    E20_4
  //               \  /
  //               V1_3
  //   Local indexing for cell nodes: (i,j) and k
  //      3  9   degCellPoly=3
  //         | \                .
  //      2  7--8    k
  //    j    |  | \             .
  //      1  4--5--6
  //         |  |  | \          .
  //      0  0--1--2--3
  //             i
  //    Number of nodes = N = (degCellPoly+2)*(degCellPoly+1)/2
  //    k = mapK(i,j,degCellPoly) = N - (degCellPoly+2-j)*(degCellPoly+1-j)/2 + i
  //    j = mapJ(k,degCellPoly) = extract from an array, my_ds_space->map_j_array[degCellPoly][k],
  //                                                      e.g., [0 0 0 0 1 1 1 2 2 3]
  //    i = mapI(k,degCellPoly) = k - mapK(0,j,degCellPoly)
  //
  // Global numbering:
  //   By element, counting shared nodes only on lowest numbered element.
  //
  ////////////////////////////////////////////////////////////////////////////////

  class DirectSerendipity
  {
  private:
    int polynomial_degree;
    polymesh::PolyMesh* my_mesh;
    //DirectSerendipity my_high_order_ds_space;

    int num_nodes;
    Node* the_ds_nodes = nullptr;
    NodeType* the_node_type = nullptr;
    BCType* the_bc_type = nullptr;

    DirectSerendipityFE* the_ds_elements = nullptr;
    std::vector<int>* map_j_array = nullptr;

    // Connectivity
    int* mesh_vertex_to_node_index = nullptr;
    int* mesh_edge_to_first_node_index = nullptr; // first edge node index only
    int* mesh_element_to_first_node_index = nullptr; // first cell node index only

    void set_directserendipity(int polyDeg, polymesh::PolyMesh* mesh);
    
  public:
    DirectSerendipity() : polynomial_degree(0), my_mesh(nullptr), num_nodes(0), the_ds_nodes(nullptr),
			  the_node_type(nullptr), the_bc_type(nullptr),
			  the_ds_elements(nullptr), map_j_array(nullptr),
			  mesh_vertex_to_node_index(nullptr), mesh_edge_to_first_node_index(nullptr),
			  mesh_element_to_first_node_index(nullptr) {};
    DirectSerendipity(int polyDeg, polymesh::PolyMesh* mesh) {
      set_directserendipity(polyDeg, mesh); };
    ~DirectSerendipity();

    void set(int polyDeg, polymesh::PolyMesh* mesh) { set_directserendipity(polyDeg, mesh); };
    
    int degPolyn() const { return polynomial_degree; };
    polymesh::PolyMesh* mesh() const { return my_mesh; };
    
    int nNodes() const { return num_nodes; };
    int nVertexNodes() const { return my_mesh->nVertices(); };
    int nInteriorEdgeNodes() const { return my_mesh->nEdges() * (polynomial_degree - 1); };
    int nEdgeNodes() const { return nVertexNodes() + nInteriorEdgeNodes(); };
    int nCellNodes() const { return num_nodes - nInteriorEdgeNodes() - nVertexNodes(); };
    int nInteriorElementNodes() const { return nCellNodes(); };

    Node* nodePtr(int i) const { return &the_ds_nodes[i]; }
    DirectSerendipityFE* finiteElementPtr(int i) const { return &the_ds_elements[i]; }


    NodeType nodeType(int i) const { return the_node_type[i]; };
    BCType bcType(int i) const { return the_bc_type[i]; };

    int meshVertexToNodeIndex(int i) const { return mesh_vertex_to_node_index[i]; };
    int meshEdgeToFirstNodeIndex(int i) const { return mesh_edge_to_first_node_index[i]; };
    int meshElementToFirstNodeIndex(int i) const { return mesh_element_to_first_node_index[i]; };

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    int write_matlab(std::string& filename) const;

    friend class DirectSerendipityFE;
    friend class DirectMixedFE;
    friend class DirectSerendipityArray;
    friend class DirectSerendipityTensor1;
    friend class Node;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class DirectMixed
  //   Defined on a poly-mesh (class PolyMesh), consisting of direct mixed
  //     finite elements (class DirectMixedFE)
  //    
  //   Gives function evaluation from global coefficients for the basis functions
  //
  ////////////////////////////////////////////////////////////////////////////////

  class DirectMixed
  {
  private:
    int polynomial_degree;
    polymesh::PolyMesh* my_mesh;
    //DirectSerendipity my_high_order_ds_space;

    int num_edges; //redundant with my_mesh
    int num_interior_edges;
    polymesh::Edge* the_dm_edges = nullptr;
    EdgeBCType* the_bc_type = nullptr;
    int* interior_edge_indexing = nullptr;
    int* global_edge_indexing = nullptr;

    DirectMixedFE* the_dm_elements = nullptr;
    DirectDGFE* the_dg_elements = nullptr;
    DirectEdgeDGFE* the_dg_edge_elements = nullptr;

    void set_directmixed(int polyDeg, polymesh::PolyMesh* mesh);
    
  public:
    DirectMixed() : polynomial_degree(0), my_mesh(nullptr), num_edges(0),
    the_dm_edges(nullptr), the_bc_type(nullptr),
			  the_dm_elements(nullptr) {};
    DirectMixed(int polyDeg, polymesh::PolyMesh* mesh) {
      set_directmixed(polyDeg, mesh); };
    ~DirectMixed();

    void set(int polyDeg, polymesh::PolyMesh* mesh) { set_directmixed(polyDeg, mesh); };
    
    int nEdges() const { return num_edges; };
    int nInteriorEdges() const { return num_interior_edges; };
    int degPolyn() const { return polynomial_degree; };
    polymesh::PolyMesh* mesh() const { return my_mesh; };
    polymesh::PolyElement* elementPtr(int i) const {return my_mesh->elementPtr(i); };

    DirectMixedFE* MixedElementPtr(int i) const { return &the_dm_elements[i]; };
    DirectDGFE* DGElementPtr(int i) const { return &the_dg_elements[i]; };
    DirectEdgeDGFE* DGEdgePtr(int i) const { return &the_dg_edge_elements[i]; };
    DirectEdgeDGFE* DGEdgeInteriorPtr(int i) const { return &the_dg_edge_elements[int_to_glob(i)]; };
    EdgeBCType bcType(int i) const { return the_bc_type[i]; }; //Boundary type for each edge
    
    // Return the interior edge indexing from global edge indexing
    int glob_to_int(int i) const { return interior_edge_indexing[i]; };

    // Return the global edge indexing from interior edge indexing
    int int_to_glob(int i) const {return global_edge_indexing[i]; }; 

    // Return the interior edge indexing of the iEdge-th edge of the iElement-th element
    int interiorEdgeIndex(int iElement, int iEdge) const { 
      return interior_edge_indexing[my_mesh->elementPtr(iElement)->edgePtr(iEdge)->meshIndex()]; 
    };

    // Return the edge from interior edge index
    polymesh::Edge* edgeInteriorPtr(int i) const {
      return my_mesh -> edgePtr(int_to_glob(i));
    };

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    int write_matlab(std::string& filename) const;


    friend class DirectMixedFE;
    friend class DirectDGFE;
    friend class DirectEdgeDGFE;
    friend class DirectMixedArray;
  };


  ////////////////////////////////////////////////////////////////////////////////
  // class DirectMixedArray
  //    Gives values for each coefficient of basis function
  ////////////////////////////////////////////////////////////////////////////////

  class DirectMixedArray
  {
  private:
    int num_edges;
    double* the_array = nullptr;

    DirectMixed* my_dm_space;

    void set_directmixedarray(DirectMixed* dmSpace);

  public:
    DirectMixedArray() : num_edges(0), the_array(nullptr), my_dm_space(nullptr) {};
    DirectMixedArray(DirectMixed* dmSpace) {
      set_directmixedarray(dmSpace); };
    DirectMixedArray(const DirectMixedArray& a) : the_array(nullptr) {
      set_directmixedarray(a.dmSpace()); };
    ~DirectMixedArray();
    
    void set(DirectMixed* dsSpace) { set_directmixedarray(dsSpace); };

    DirectMixed* dmSpace() const { return my_dm_space; };
    int size() const { return num_edges; };

    double& operator() (int i)       { return the_array[i]; }
    double  operator() (int i) const { return the_array[i]; }
    double& operator[] (int i)       { return the_array[i]; }
    double  operator[] (int i) const { return the_array[i]; }

    void eval(const Point* pts, double* result, Tensor1* gradResult, int num_pts) const;
    void eval(const Point& pt, double& result, Tensor1& gradResult) const;
    double eval(const Point& pt) const;

    void l2normError(double& l2Error, double& l2GradError, double& l2Norm, double& l2GradNorm,
		     double (*referenceFcn)(double,double) = nullptr, 
		     Tensor1 (*referenceGradFcn)(double,double) = nullptr);
    void l2norm(double& l2Norm, double& l2GradNorm) {
      double null1,null2; l2normError(l2Norm,l2GradNorm,null1,null2); };

    void write_matlab_mesh(std::ofstream* fout, std::ofstream* fout_grad,
			   int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh(std::ofstream& fout, std::ofstream& fout_grad,
			   int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, &fout_grad, num_pts_x, num_pts_y); };
    void write_matlab_mesh(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, nullptr, num_pts_x, num_pts_y); };
    int write_matlab_mesh(std::string& filename, std::string& filename_grad,
			  int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_matlab_mesh_by_pt(std::ofstream& fout, std::ofstream& fout_grad,
				 int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, std::string& filename_grad,
				int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };
};
#endif
