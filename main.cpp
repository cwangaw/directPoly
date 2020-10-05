///////////////////////////////////////////////////////////////////////////////
//
//           DIRECT SERENDIPITY ELEMENTS ON POLYGONS
//
//  Todd Arbogast <arbogast@ices.utexas.edu>
//  Chuning Wang <cwangaw@utexas.edu>
//
//  Center for Subsurface Modeling
//  Oden Institute for Computational Engineering and Sciences
//  The University of Texas at Austin
//  Austin, Texas  78712
//
//  Begun: July 2020
//
//=============================================================================
//  EQUATIONS
//
//    - div D grad p + div(b p) + c.grad p + a p = f
//                       in Omega (bounded, connected, Lipschitz domain in R^2)
//    p = g              on the boundary of Omega
//    where D is a tensor (matrix), b and c are vectors.
//
//    Find p in DS_{r,g} such that
//      (D grad p , grad q) + (b p , grad q) + (c.grad p , q) + (a p , q)
//          = (f , q) for all q in DS_{r,0}
//    where DS_r are the direct serendipity spaces and DS_{r,g} has the
//      boundary nodes set to g. We use a nodal basis for the space.
//      
//  ELEMENTS E
//    E is a convex, nondegenerate polygon (of at least 3 sides)
//    Polynomials of degree r >= 1
//    The mesh is conforrming
//
//=============================================================================
//
//  ACKNOWLEDGMENTS
//    The development of this code has been supported by
//    the U.S. National Science Foundation.
//
//  COPYRIGHT 2020
//    No portion of this document or program may be reproduced,
//    transmitted or otherwise copied without prior written permission of
//    the author, except for the internal use of the Center for Subsurface
//    Modeling (CSM), The University of Texas at Austin. The authors make
//    no representations or warranties about the correctness of any portion
//    of the program code, supplementary codes, or associated documentation,
//    nor about the suitability of this program for any purpose, and is in
//    no way liable for any damages resulting from its use or misuse.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <sstream>

#include "debug.h"
#include "fcns.h"
#include "Utilities/monitor.h"
#include "parameterData.h"
#include "Utilities/version.h"
#include "Reader/expr.h"
#include "ellipticPDE.h"
//#include "baseObjects.h"
//#include "polyMesh.h"
//#include "directSerendipity.h"

////////////////////////////////////////////////////////////////////////////////
// MAIN ROUTINES

// OUTPUT HEADER INFORMATION
static void write_header (const Version& code, const ParameterData& param)
{
  std::cout << "\n";
  std::cout << "       DIRECT SERENDIPITY ELEMENTS ON POLYGONS\n\n";

  std::cout << "           Name:       " << code.name << "\n";
  std::cout << "           Version:    " << code.version[0];
  for(int i=1; i<code.version_depth; i++) std::cout << "." << code.version[i];
  std::cout << "\n";
  std::cout << "           Date:       " << code.date << "\n";
  std::cout << code.description << "\n\n";

  std::cout << "           Input file: " << param.infile_name << "\n";
  if(param.echo)
  std::cout << "           Echo file:  " << param.echofile_name << "\n";
  std::cout << std::endl;
  return;
}

// PARSE COMMAND LINE ARGUMENTS
static int parse_flags (int argc, char** argv, ParameterData& param, Version& code)
{
  for(int i = 1; i < argc; i++) {
    if(argv[i][0] == '-')
      {
	switch(argv[i][1])
	  {
	  case 'h':
	    std::cout << "Usage: " << argv[0] << " [command-options]\n";
	    std::cout << "  Where the optional [command-options] are:\n";
	    std::cout << "    -i [<file name>]   The initial general input <file"
		      << " name>, or \"" << param.infile_name << "\" if omitted\n";
	    std::cout << "    -e [<file name>]   Echo input to <file name>, or \""
		      << param.echofile_name << "\" if omitted\n";
	    std::cout << "    -E                 Echo input to standard out\n";
	    std::cout << "    -c                 Just display input case information\n";
	    std::cout << "    -u                 Just display input physical units information\n";
	    std::cout << "    -v                 Just display version information\n";
	    std::cout << "    -h                 Print this usage information only";
	    std::cout << std::endl;
	    return 1;
	  case 'i':
	    i++;
	    if(i < argc && argv[i][0] != '-') param.infile_name = argv[i]; else i--;
	    break;
	  case 'e':
	    param.echo = 1;
	    i++;
	    if(i < argc && argv[i][0] != '-') param.echofile_name = argv[i]; else i--;
	    break;
	  case 'E':
	    param.echo = 1;
	    param.echofile_name = '\0';
	    break;
	  case 'c':
	    param.print_cases();
	    return 1;
	  case 'u':
	    expr_listUnits(std::cout);
	    return 1;
	  case 'v':
	    write_header (code, param);
	    return 1;
	  default:
	    std::cerr << "Unrecognized command line option: " << argv[i] << std::endl;
	    return 1;
	  }
      }
    else
      {
	std::cerr << "Unrecognized item on command line: " << argv[i] << "\n";
	std::cerr << "For command line arguments, try: " << argv[0] << " -h"<< std::endl;
	return 1;
      }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// MAIN

int main(int argc, char* argv[]) {
  // Initialize
  Version code("directpoly","September 2020",3);
  code.set_version(1,0,1);
  code.set_description("Solves a second order PDE using DS spaces on polygons\nAllows simple mesh construction");

  ParameterData param;
  Monitor monitor;
  
  try {
    if ( parse_flags(argc, argv, param, code) ) return 0;

    write_header (code, param);

    if ( param.read () ) return 1;
    monitor.reset (0, param.monitor_to_level);
    monitor (0, "THE DATA WAS READ");

    monitor (0, "\nBEGIN COMPUTATION\n");

    EllipticPDE ellipticPDE(param);
    ellipticPDE.solve(monitor);
    
  }
  catch (std::exception &exc) {
    std::cerr << "\n\n----------------------------------------------------\n";
    std::cerr << "ERROR: " << exc.what() << std::endl;
    std::cerr << "----------------------------------------------------" << std::endl;
    return 1;
  }
  catch (int error) {
    std::cerr << "\n\n----------------------------------------------------\n";
    std::cerr << "ERROR, CODE: " << error << std::endl;
    std::cerr << "----------------------------------------------------" << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << "\n\n----------------------------------------------------\n";
    std::cerr << "ERROR, UNSPECIFIED" << std::endl;
    std::cerr << "----------------------------------------------------" << std::endl;
    return 1;
  }

  monitor(0,"\nEND COMPUTATION");
  return 0;
}
