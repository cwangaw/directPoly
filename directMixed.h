#ifndef __directmixed_h_included__
#define __directmixed_h_included__

#include "directSerendipity.h"

namespace directserendipity {

  enum class EdgeBCType { interior, boundary };
  class DirectSerendipity;
  class DirectMixed;
  
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

    // If the functions are indexed including curl part, 
    // we need to pass iFunc = index of the function - dimCurlPart
    double basisdivXPo(int iFunc, int iPt) const {
      return v_div_value_n[iPt * dim_v_div + iFunc];
    }

    // Get dimensions of spaces
    int dimVFull() const { return dim_v; };
    int dimVReduced() const { return dim_v - polynomial_degree - 1; };

    int dimCurlPart() const { return dim_curlpart; }

    int dimXPoFull() const { return dim_v_div; };
    int dimXPoReduced() const { return polynomial_degree * (polynomial_degree + 1)/2; };

    // Evaluation \u in full or reduced space or both at a point
    void eval(const Point* pt, Tensor1* result, int num_pts, char type, double* dofs=nullptr);
    void eval(const Point& pt, Tensor1& result, char type, double* dofs=nullptr) {
      eval(&pt, &result, 1, type, dofs); };

    Tensor1 eval(const Point& pt, char type, double* dofs=nullptr) { 
      Tensor1 result;
      eval(&pt, &result, 1, type, dofs); 
      return result; 
    }

    void eval(const Point* pt, Tensor1* fullResult, Tensor1* reducedResult, int num_pts, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr);
    void eval(const Point& pt, Tensor1& fullResult, Tensor1& reducedResult, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr) {
      eval(&pt, &fullResult, &reducedResult, 1, full_dofs, reduced_dofs); };
    
    // Output functions
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
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

    polymesh::PolyElement* elementPtr() const { return my_poly_element; };

    void initBasis(const Point* pt, int num_pts); // (re)evaluate all basis fcns at the points

    // Access basis functions evaluated at pt[iPt]
    double basis(int iFunc, int iPt) const {
      return value_n[iPt * dim_w + iFunc];
    };

    // Get dimensions of spaces
    int dimFull() const { return dim_w; };
    int dimReduced() const { return (polynomial_degree+1) * polynomial_degree /2;}

    // Evaluation \p of full or reduced space at a point
    void eval(const Point* pt, double* result, int num_pts, char type, double* dofs=nullptr);
    void eval(const Point& pt, double& result,char type, double* dofs=nullptr) {
      eval(&pt, &result, 1, type, dofs); };

    double eval(const Point& pt, char type, double* dofs=nullptr) { 
      double result;
      eval(&pt, &result, 1, type, dofs); 
      return result; 
    }

    void eval(const Point* pt, double* fullResult, double* reducedResult, int num_pts, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr);
    void eval(const Point& pt, double& fullResult, double& reducedResult, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr) {
      eval(&pt, &fullResult, &reducedResult, 1, full_dofs, reduced_dofs); };
    
    // Output functions
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
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

    polymesh::Edge* edgePtr() const { return my_edge; };

    void initBasis(const Point* pt, int num_pts); // (re)evaluate all basis fcns at the points

    // Access basis functions evaluated at pt[iPt]
    double basis(int iFunc, int iPt) const {
      return value_n[iPt * dim_l + iFunc];
    };

    // Get dimensions of spaces
    int dim() { return dim_l; };

    // Evaluation \lambda at a point
    void eval(const Point* pt, double* result, int num_pts, double* dofs=nullptr);
    void eval(const Point& pt, double& result, double* dofs=nullptr) {
      eval(&pt, &result, 1, dofs); };

    double eval(const Point& pt, double* dofs=nullptr) { 
      double result;
      eval(&pt, &result, 1, dofs); 
      return result; 
    }
    
    // Output functions
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
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

    int mixed_dofs_full;
    int mixed_dofs_reduced;
    int dg_dofs_full;
    int dg_dofs_reduced;
    int dg_edge_dofs; // Here we also include DoFs on the boundary

    void set_directmixed(int polyDeg, polymesh::PolyMesh* mesh);
    
  public:
    DirectMixed() : polynomial_degree(0), my_mesh(nullptr), num_edges(0), num_interior_edges(0),
                    the_dm_edges(nullptr), the_bc_type(nullptr), interior_edge_indexing(nullptr),
                    global_edge_indexing(nullptr), the_dm_elements(nullptr), the_dg_elements(nullptr),
                    the_dg_edge_elements(nullptr), mixed_dofs_full(0), mixed_dofs_reduced(0),
                    dg_dofs_full(0), dg_dofs_reduced(0), dg_edge_dofs(0) {};
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

    // Return the number of global dofs of spaces
    int nMixedDoFs(char type) const;
    int nDGDoFs(char type) const;
    int nEdgeDGDoFs() const;

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
    friend class DirectDGArray;
    friend class DirectEdgeDGArray;
  };


  ////////////////////////////////////////////////////////////////////////////////
  // class DirectMixedArray
  //    Gives values for each coefficient of basis function
  ////////////////////////////////////////////////////////////////////////////////

  class DirectMixedArray
  {
  private:
    int num_dofs;
    int num_elements; // Redundant with my_dm_space
    char space_type; // either 'f' or 'r'

    double* the_array = nullptr;

    DirectMixed* my_dm_space;

    void set_directmixedarray(DirectMixed* dmSpace, char spacetype);

  public:
    DirectMixedArray() : num_dofs(0), num_elements(0), space_type('f'), the_array(nullptr), my_dm_space(nullptr) {};
    DirectMixedArray(DirectMixed* dmSpace, char spacetype) {
      set_directmixedarray(dmSpace, spacetype); };
    DirectMixedArray(const DirectMixedArray& a) : the_array(nullptr) {
      set_directmixedarray(a.dmSpace(), a.spaceType()); };
    ~DirectMixedArray();
    
    void set(DirectMixed* dmSpace, char spacetype) { set_directmixedarray(dmSpace, spacetype); };

    DirectMixed* dmSpace() const { return my_dm_space; };
    char spaceType() const { return space_type; }
    int size() const { return num_dofs; };

    double& operator() (int i)       { return the_array[i]; }
    double  operator() (int i) const { return the_array[i]; }
    double& operator[] (int i)       { return the_array[i]; }
    double  operator[] (int i) const { return the_array[i]; }

    void eval(const Point* pts, Tensor1* result, int num_pts) const;
    void eval(const Point& pt, Tensor1& result) const;
    Tensor1 eval(const Point& pt) const {
      Tensor1 result; eval(pt, result); return result;
    };

    void l2normError(double& l2Error, double& l2Norm, Tensor1 (*referenceFcn)(double,double) = nullptr);
    void l2norm(double& l2Norm) {
      double null; l2normError(l2Norm,null); };

    void write_matlab_mesh(std::ofstream* fout, int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, num_pts_x, num_pts_y); };

    int write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };


  ////////////////////////////////////////////////////////////////////////////////
  // class DirectDGArray
  //    Gives values for each coefficient of basis function
  ////////////////////////////////////////////////////////////////////////////////

  class DirectDGArray
  {
  private:
    int num_dofs;
    int num_elements; // Redundant with my_dm_space
    char space_type; // either 'f' or 'r'

    double* the_array = nullptr;

    DirectMixed* my_dm_space;

    void set_directdgarray(DirectMixed* dmSpace, char spacetype);

  public:
    DirectDGArray() : num_dofs(0), num_elements(0), space_type('f'), the_array(nullptr), my_dm_space(nullptr) {};
    DirectDGArray(DirectMixed* dmSpace, char spacetype) {
      set_directdgarray(dmSpace, spacetype); };
    DirectDGArray(const DirectDGArray& a) : the_array(nullptr) {
      set_directdgarray(a.dmSpace(), a.spaceType()); };
    ~DirectDGArray();
    
    void set(DirectMixed* dmSpace, char spacetype) { set_directdgarray(dmSpace, spacetype); };

    DirectMixed* dmSpace() const { return my_dm_space; };
    char spaceType() const { return space_type; }
    int size() const { return num_dofs; };

    double& operator() (int i)       { return the_array[i]; }
    double  operator() (int i) const { return the_array[i]; }
    double& operator[] (int i)       { return the_array[i]; }
    double  operator[] (int i) const { return the_array[i]; }

    void eval(const Point* pts, double* result, int num_pts) const;
    void eval(const Point& pt, double& result) const;
    double eval(const Point& pt) const {
      double result; eval(pt, result); return result;
    };

    void l2normError(double& l2Error, double& l2Norm, double (*referenceFcn)(double,double) = nullptr);
    void l2norm(double& l2Norm) {
      double null; l2normError(l2Norm,null); };

    void write_matlab_mesh(std::ofstream* fout, int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, num_pts_x, num_pts_y); };

    int write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class DirectEdgeDGArray
  //    Gives values for each coefficient of basis function
  //    Here we do not distinguish boundary and
  ////////////////////////////////////////////////////////////////////////////////

  class DirectEdgeDGArray
  {
  private:
    int num_dofs;
    int num_edges; // Redundant with my_dm_space

    double* the_array = nullptr;

    DirectMixed* my_dm_space;

    void set_directedgedgarray(DirectMixed* dmSpace);

  public:
    DirectEdgeDGArray() : num_dofs(0), num_edges(0), the_array(nullptr), my_dm_space(nullptr) {};
    DirectEdgeDGArray(DirectMixed* dmSpace) {
      set_directedgedgarray(dmSpace); };
    DirectEdgeDGArray(const DirectDGArray& a) : the_array(nullptr) {
      set_directedgedgarray(a.dmSpace()); };
    ~DirectEdgeDGArray();
    
    void set(DirectMixed* dmSpace) { set_directedgedgarray(dmSpace); };

    DirectMixed* dmSpace() const { return my_dm_space; };
    int size() const { return num_dofs; };

    double& operator() (int i)       { return the_array[i]; }
    double  operator() (int i) const { return the_array[i]; }
    double& operator[] (int i)       { return the_array[i]; }
    double  operator[] (int i) const { return the_array[i]; }

    void eval(const Point* pts, double* result, int num_pts) const;
    void eval(const Point& pt, double& result) const;
    double eval(const Point& pt) const {
      double result; eval(pt, result); return result;
    };

    void l2normError(double& l2Error, double& l2Norm, double (*referenceFcn)(double,double) = nullptr);
    void l2norm(double& l2Norm) {
      double null; l2normError(l2Norm,null); };

    void write_matlab_mesh(std::ofstream* fout, int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, num_pts_x, num_pts_y); };

    int write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };
};
#endif