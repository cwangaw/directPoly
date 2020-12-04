#include <fstream>

#include "debug.h"
#include "parameterData.h"
#include "Reader/reader.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// PARAMETERDATA

ParameterData::ParameterData()
  : infile_name("infile"),
    echofile_name("echo"),
    echo(0) {}

ParameterData::ParameterData(const std::string& infileName)
  : infile_name(infileName),
    echofile_name("echo"),
    echo(0) {}

ParameterData::ParameterData(const std::string& infileName,
				  const std::string& echofileName, bool echo0)
  : infile_name(infileName),
    echofile_name(echofileName),
    echo(echo0) {}

// CASES =========================================================================

void ParameterData::print_cases(ostream& fout) {
  int n;

  // Soln output
  fout<< "Soln output formats (" << n_case_soln_output << " total):\n";
  n=0;
  fout<< "  " << case_soln_output_none
      << " = no soln output\n";
  n++;
  fout<< "  " << case_soln_output_raw
      << " = raw text file\n";
  n++;
  fout<< "  " << case_soln_output_matlab
      << " = matlab file\n";
  n++;
  if(n != n_case_soln_output) {
    fout<< "  " << "other cases exist" << endl;
  }

  // Mesh output
  fout<< "Mesh output formats (" << n_case_mesh_output << " total):\n";
  n=0;
  fout<< "  " << case_mesh_output_none
      << " = no mesh output\n";
  n++;
  fout<< "  " << case_mesh_output_raw
      << " = raw text file\n";
  n++;
  fout<< "  " << case_mesh_output_matlab
      << " = matlab file\n";
  n++;
  if(n != n_case_mesh_output) {
    fout<< "  " << "other cases exist" << endl;
  }

  // DSSpace output
  fout<< "DS space output formats (" << n_case_dsSpace_output << " total):\n";
  n=0;
  fout<< "  " << case_dsSpace_output_none
      << " = no DS space output\n";
  n++;
  fout<< "  " << case_dsSpace_output_raw
      << " = raw text file\n";
  n++;
  fout<< "  " << case_dsSpace_output_matlab
      << " = matlab file\n";
  n++;
  if(n != n_case_dsSpace_output) {
    fout<< "  " << "other cases exist" << endl;
  }
}

// ===================================================================================

#define ERRCHK(s) {error=s; if(error){processReaderError(error);return error;}}
#define ERRRET(s) {error=s; if(error) return error;}

// MESH AND FINITE ELEMENTS ==========================================================

// Mesh files. Vertices (x,y) in a list: x0,y0,x1,y1, ...
//             Elements by vertices in CCW direction in a bracked list: [v00,v01,...][v10, v11,...] ...
int ParameterData::readMesh() {
  int error = 0;
  int dimen[NUMBER_DIMEN] = {0,0,0};
  double unitVal;
  char c;

  const int offset = 1;

  double vertices[2*nVertices];
  ERRCHK(readArray(vertices, 2*nVertices, dimen, unitVal));

  int nVerticesPerElement[nElements];
  int* elements[nElements];
  vector<int> elementList;

  ERRCHK(beginList());
  for(int i=0; i<nElements; i++) {
    int ii;
    
    ERRCHK(skipWS(0,0));
    c = getc();
    if(c!='[') { error = ERR_IN_DATA_BLOCK; processReaderError(error); return error; }
    reader_echo(c);
    
    ERRCHK(skipWS(0,0));
    c = getc();
    ungetc(c);
    
    int nElementVertices = 0;
    while(c!=']') {
      ERRCHK(readScalar(ii));
      elementList.push_back(ii);
      nElementVertices++;
      ERRCHK(skipWS(0,0));
      c = getc();
      ungetc(c);
    }
    c = getc();
    reader_echo(c);

    nVerticesPerElement[i] = nElementVertices;
    elements[i] = new int[nElementVertices];
    for(int j=0; j<nElementVertices; j++) {
      elements[i][j] = elementList[j] - offset;
    }
    elementList.clear();
  }
  ERRCHK(endList());

  mesh.set(nVertices,vertices,nElements,nVerticesPerElement,elements);
  if(!mesh.isGood()) { error = ERR_BAD_DATA; processReaderError(error); return ERR_IN_GRID_ARRAY; }

  for(int i=0; i<nElements; i++) delete[] elements[i];
  return error;
}

// Output mesh
int ParameterData::writeMesh() {
  if(output_mesh_format == case_mesh_output_none) return 0;

  string fileName = directory_name;
  fileName += "meshOut";

  switch(output_mesh_format) {
  case case_mesh_output_raw: {
    fileName += ".txt";
    if(mesh.write_raw(fileName)) return ERR_FILESYSTEM;
    break;
  }
  case case_mesh_output_matlab: {
    if(mesh.write_matlab(fileName)) return ERR_FILESYSTEM;
    break;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
  return 0;
}

// Output DS space
int ParameterData::writeDSSpace() {
  if(output_dsSpace_format == case_dsSpace_output_none) return 0;

  string fileName = directory_name;
  fileName += "dsSpaceOut";

  switch(output_dsSpace_format) {
  case case_dsSpace_output_raw: {
    fileName += ".txt";
    if(dsSpace.write_raw(fileName)) return ERR_FILESYSTEM;
    break;
  }
  case case_dsSpace_output_matlab: {
    if(dsSpace.write_matlab(fileName)) return ERR_FILESYSTEM;
    break;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
  return 0;
}

// READ ==========================================================================

#define SIZE_STRING 15

int ParameterData::read() {
  int error = 0;

  ERRCHK(openReader(infile_name, echofile_name, 121, 1, echo));

  // OUTPUT DIRECTORY

  ERRCHK(readString(directory_name));

  // MESH PARAMETERS

  ERRCHK(readString(mesh_type));
  char meshTypeC;
  switch(mesh_type[0]) {
  case 'q': case 'Q': meshTypeC = 'q'; break;
  case 'd': case 'D': meshTypeC = 'd'; break;
  case 't': case 'T': meshTypeC = 't'; break;
  case 'g': case 'G': meshTypeC = 'g'; break;
  default: return processReaderError(ERR_UNSUPPORTED_OPTION);
  }	       

  switch(meshTypeC) {
  case 'q':
  case 'd':
  case 't': {
    double xMin, xMax;
    ERRCHK(readScalar(xMin));
    ERRCHK(readScalar(xMax));
    if(xMin >= xMax) return processReaderError(ERR_BAD_DATA);
    
    double yMin, yMax;
    ERRCHK(readScalar(yMin));
    ERRCHK(readScalar(yMax));
    if(yMin >= yMax) return processReaderError(ERR_BAD_DATA);

    int nx, ny;
    ERRCHK(readScalar(nx));
    ERRCHK(readScalar(ny));
    if(nx < 1) nx = 1;
    if(ny < 1) ny = 1;

    ERRCHK(readScalar(distortion_factor));

    int err = mesh.createMesh(meshTypeC, nx,ny, xMin,xMax, yMin,yMax, distortion_factor);
    switch(err) {
    case 1: return processReaderError(ERR_UNSUPPORTED_OPTION);
    case 2: return processReaderError(ERR_BAD_DATA);
    }
    break;
  }

  case 'g': {
    ERRCHK(readScalar(nVertices));
    ERRCHK(readScalar(nElements));
    ERRCHK(readMesh());
    break;
  }

  default:
    return processReaderError(ERR_UNSUPPORTED_OPTION);
  }

  // FINITE ELEMENTS

  ERRCHK(readScalar(polynomial_degree));
  if(polynomial_degree < 1) polynomial_degree = 1;

  dsSpace.set(polynomial_degree,&mesh);
  dmSpace.set(polynomial_degree,&mesh);
  
  // ALGORITHM PARAMETERS

  /*
  ERRCHK(readScalar(maximumIterations));
  if(maximumIterations < 0) maximumIterations = 0;

  ERRCHK(readScalar(absoluteTolerance));
  if(absoluteTolerance < 0) absoluteTolerance = 0;

  ERRCHK(readScalar(relativeTolerance));
  if(relativeTolerance < 0) relativeTolerance = 0;
  if(absoluteTolerance == 0 && relativeTolerance == 0)
    return processReaderError(ERR_OUT_OF_RANGE);
  */
  
  // OUTPUT PARAMETERS

  ERRCHK(readScalar(output_soln_format));
  if(output_soln_format < 0 || output_soln_format >= n_case_soln_output)
    return processReaderError(ERR_OUT_OF_RANGE);

  ERRCHK(readScalar(output_mesh_numPts_x));
  if(output_mesh_numPts_x < 2) output_mesh_numPts_x = 2;
  ERRCHK(readScalar(output_mesh_numPts_y));
  if(output_mesh_numPts_y < 2) output_mesh_numPts_y = 2;

  ERRCHK(readScalar(output_mesh_format));
  if(output_mesh_format < 0 || output_mesh_format >= n_case_mesh_output)
    return processReaderError(ERR_OUT_OF_RANGE);
  ERRCHK(writeMesh());

  ERRCHK(readScalar(output_dsSpace_format));
  if(output_dsSpace_format < 0 || output_dsSpace_format >= n_case_dsSpace_output)
    return processReaderError(ERR_OUT_OF_RANGE);
  ERRCHK(writeDSSpace());

  ERRCHK(readScalar(monitor_to_level));

  return error;
}
