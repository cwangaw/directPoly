#include <cmath>
using namespace std;

#include "debug.h"
#include "fcns.h"

////////////////////////////////////////////////////////////////////////////////
// SPATIALLY DEPENDENT PDE FUNCTIONS

// Coefficients
static const double PI = 3.141592653589793238463;

double coefA(double x, double y) {
  return 1; //0;
}

void coefB(double x, double y, Tensor1& b) {
  b.set(0,0);
}

void coefC(double x, double y, Tensor1& c) {
  c.set(0,0);
}

void coefD(double x, double y, Tensor2& d) {
  d.set(0,0,0,0);
  //d.set(1,0,0,1);
}

// Source f
double sourceVal(double x, double y) {
  //return 2*pow(PI,2)*sin(PI*x)*sin(PI*y); 
  //return -4;
  return 1;
}

// BC values g
double bcVal(double x, double y) {
  //return 0; 
  //return pow(x,2)+pow(y,2);
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
// TRUE SOLUTION (IF KNOWN)

bool trueSolnKnown() { return true; }

// Real solution
double trueSoln(double x, double y) {
  //return sin(PI*x)*sin(PI*y); 
  //return pow(x,2)+pow(y,2);
  return 1;
}

Tensor1 trueGradSoln(double x, double y) {
  //return Tensor1(PI*cos(PI*x)*sin(PI*y),PI*sin(PI*x)*cos(PI*y));
  //return Tensor1(2*x,2*y);
  return Tensor1(0,0);
}
