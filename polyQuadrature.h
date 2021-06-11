#ifndef __polyquadrature_h_included__
#define __polyquadrature_h_included__

////////////////////////////////////////////////////////////////////////////////
// PolyQuadrature class
//  Computes quasrature points and weights for a PolyElement
//
// Uses the PolyElement class
// Assume base objects: Point, Tensor1, and Tensor2
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include "Mesh/baseObjects.h"
#include "Mesh/polyMesh.h"

namespace polyquadrature {

  ////////////////////////////////////////////////////////////////////////////////
  // Quadrature rules on reference right triangle (0,0) (1,0) (0,1)
  //    num = number of points in the rule
  //    dop = degree of precision of the rule
  //    pts = quadrature points (using Point)
  //    wts = quadrature weights
  ////////////////////////////////////////////////////////////////////////////////

  static struct {
    const int num = 3;
    const int dop = 2;
    Point pts[3] = { Point(0.66666666666666666667, 0.16666666666666666667),
		     Point(0.16666666666666666667, 0.66666666666666666667),
		     Point(0.16666666666666666667, 0.16666666666666666667) };
    double wts[3] = { 0.33333333333333333333, 0.33333333333333333333,
		      0.33333333333333333333 };
  } strang1;

  static struct {
    const int num = 3;
    const int dop = 2;
    Point pts[3] = { Point(0.50000000000000000000, 0.00000000000000000000),
		     Point(0.50000000000000000000, 0.50000000000000000000),
		     Point(0.00000000000000000000, 0.50000000000000000000) };
    double wts[3] = { 0.33333333333333333333, 0.33333333333333333333,
		      0.33333333333333333333 };
  } strang2;

  static struct {
    const int num = 4;
    const int dop = 3;
    Point pts[4] = { Point(0.33333333333333333333, 0.33333333333333333333),
		     Point(0.60000000000000000000, 0.20000000000000000000),
		     Point(0.20000000000000000000, 0.60000000000000000000),
		     Point(0.20000000000000000000, 0.20000000000000000000) };
    double wts[4] = { -0.56250000000000000000, 0.52083333333333333333,
		      0.52083333333333333333, 0.52083333333333333333 };
  } strang3;

  static struct {
    const int num = 6;
    const int dop = 4;
    Point pts[6] = { Point(0.816847572980459, 0.091576213509771),
		     Point(0.091576213509771, 0.816847572980459),
		     Point(0.091576213509771, 0.091576213509771),
		     Point(0.108103018168070, 0.445948490915965),
		     Point(0.445948490915965, 0.108103018168070),
		     Point(0.445948490915965, 0.445948490915965) };
    double wts[6] = { 0.109951743655322, 0.109951743655322,
		      0.109951743655322, 0.223381589678011,
		      0.223381589678011, 0.223381589678011 };
  } strang5;

  static struct {
    const int num = 7;
    const int dop = 5;
    Point pts[7] = { Point(0.33333333333333333, 0.33333333333333333),
		     Point(0.79742698535308720, 0.10128650732345633),
		     Point(0.10128650732345633, 0.79742698535308720),
		     Point(0.10128650732345633, 0.10128650732345633),
		     Point(0.05971587178976981, 0.47014206410511505),
		     Point(0.47014206410511505, 0.05971587178976981),
		     Point(0.47014206410511505, 0.47014206410511505) };
    double wts[7] = { 0.22500000000000000, 0.12593918054482717,
		      0.12593918054482717, 0.12593918054482717,
		      0.13239415278850616, 0.13239415278850616,
		      0.13239415278850616 };
  } strang7;

  // Strang8 is incorrect!

  static struct {
    const int num = 12;
    const int dop = 6;
    Point pts[12] = { Point(0.873821971016996, 0.063089014491502),
		      Point(0.063089014491502, 0.873821971016996),
		      Point(0.063089014491502, 0.063089014491502),
		      Point(0.501426509658179, 0.249286745170910),
		      Point(0.249286745170910, 0.501426509658179),
		      Point(0.249286745170910, 0.249286745170910),
		      Point(0.636502499121399, 0.310352451033785),
		      Point(0.636502499121399, 0.053145049844816),
		      Point(0.310352451033785, 0.636502499121399),
		      Point(0.310352451033785, 0.053145049844816),
		      Point(0.053145049844816, 0.636502499121399),
		      Point(0.053145049844816, 0.310352451033785) };
    double wts[12] = { 0.050844906370207, 0.050844906370207,
		       0.050844906370207, 0.116786275726379,
		       0.116786275726379, 0.116786275726379,
		       0.082851075618374, 0.082851075618374,
		       0.082851075618374, 0.082851075618374,
		       0.082851075618374, 0.082851075618374 };
  } strang9;

  static struct {
    const int num = 13;
    const int dop = 7;
    Point pts[13] = { Point(0.333333333333333, 0.333333333333333),
		      Point(0.479308067841923, 0.260345966079038),
		      Point(0.260345966079038, 0.479308067841923),
		      Point(0.260345966079038, 0.260345966079038),
		      Point(0.869739794195568, 0.065130102902216),
		      Point(0.065130102902216, 0.869739794195568),
		      Point(0.065130102902216, 0.065130102902216),
		      Point(0.638444188569809, 0.312865496004875),
		      Point(0.638444188569809, 0.048690315425316),
		      Point(0.312865496004875, 0.638444188569809),
		      Point(0.312865496004875, 0.048690315425316),
		      Point(0.048690315425316, 0.638444188569809),
		      Point(0.048690315425316, 0.312865496004875) };
    double wts[13] = { -0.149570044467670, 0.175615257433204,
		       0.175615257433204, 0.175615257433204,
		       0.053347235608839, 0.053347235608839,
		       0.053347235608839, 0.077113760890257,
		       0.077113760890257, 0.077113760890257,
		       0.077113760890257, 0.077113760890257,
		       0.077113760890257 };
  } strang10;

  static struct {
    const int num = 19;
    const int dop = 8;
    Point pts[19] = { Point(0.3333333333333333, 0.3333333333333333),
		      Point(0.7974269853530872, 0.1012865073234563),
		      Point(0.1012865073234563, 0.7974269853530872),
		      Point(0.1012865073234563, 0.1012865073234563),
		      Point(0.0597158717897698, 0.4701420641051151),
		      Point(0.4701420641051151, 0.0597158717897698),
		      Point(0.4701420641051151, 0.4701420641051151),
		      Point(0.5357953464498992, 0.2321023267750504),
		      Point(0.2321023267750504, 0.5357953464498992),
		      Point(0.2321023267750504, 0.2321023267750504),
		      Point(0.9410382782311209, 0.0294808608844396),
		      Point(0.0294808608844396, 0.9410382782311209),
		      Point(0.0294808608844396, 0.0294808608844396),
		      Point(0.7384168123405100, 0.2321023267750504),
		      Point(0.7384168123405100, 0.0294808608844396),
		      Point(0.2321023267750504, 0.7384168123405100),
		      Point(0.2321023267750504, 0.0294808608844396),
		      Point(0.0294808608844396, 0.7384168123405100),
		      Point(0.0294808608844396, 0.2321023267750504) };
    double wts[19] = { 0.0378610912003147, 0.0376204254131829,
		       0.0376204254131829, 0.0376204254131829,
		       0.0783573522441174, 0.0783573522441174,
		       0.0783573522441174, 0.1162714796569659,
		       0.1162714796569659, 0.1162714796569659,
		       0.0134442673751655, 0.0134442673751655,
		       0.0134442673751655, 0.0375097224552317,
		       0.0375097224552317, 0.0375097224552317,
		       0.0375097224552317, 0.0375097224552317,
		       0.0375097224552317 };
    } toms584_19;

  static struct {
    const int num = 19;
    const int dop = 9;
    Point pts[19] = { Point(0.33333333333333331     , 0.33333333333333331     ),  
		      Point(2.06349616025259287E-002, 0.48968251919873701     ),
		      Point(0.48968251919873701     , 2.06349616025259287E-002),
		      Point(0.48968251919873701     , 0.48968251919873701     ),
		      Point(0.12582081701412900     , 0.43708959149293553     ),
		      Point(0.43708959149293553     , 0.12582081701412900     ),
		      Point(0.43708959149293553     , 0.43708959149293553     ),
		      Point(0.62359292876193562     , 0.18820353561903219     ),
		      Point(0.18820353561903219     , 0.62359292876193562     ),
		      Point(0.18820353561903219     , 0.18820353561903219     ),
		      Point(0.91054097321109406     , 4.47295133944529688E-002),
		      Point(4.47295133944529688E-002, 0.91054097321109406     ),
		      Point(4.47295133944529688E-002, 4.47295133944529688E-002),
		      Point(0.74119859878449801     , 3.68384120547362581E-002),
		      Point(0.74119859878449801     , 0.22196298916076573     ),
		      Point(3.68384120547362581E-002, 0.74119859878449801     ),
		      Point(3.68384120547362581E-002, 0.22196298916076573     ),
		      Point(0.22196298916076573     , 0.74119859878449801     ),
		      Point(0.22196298916076573     , 3.68384120547362581E-002) };
    double wts[19] = { 9.71357962827961025E-002, 3.13347002271398278E-002,
		       3.13347002271398278E-002, 3.13347002271398278E-002,
		       7.78275410047754301E-002, 7.78275410047754301E-002,
		       7.78275410047754301E-002, 7.96477389272090969E-002,
		       7.96477389272090969E-002, 7.96477389272090969E-002,
		       2.55776756586981006E-002, 2.55776756586981006E-002,
		       2.55776756586981006E-002, 4.32835393772893970E-002,
		       4.32835393772893970E-002, 4.32835393772893970E-002,
		       4.32835393772893970E-002, 4.32835393772893970E-002,
		       4.32835393772893970E-002 };
  } toms612_19;

  static struct {
    const int num = 28;
    const int dop = 11;
    Point pts[28] = { Point(0.33333333333333333, 0.333333333333333333),
		      Point(0.9480217181434233 , 0.02598914092828833),
		      Point(0.02598914092828833, 0.9480217181434233),
		      Point(0.02598914092828833, 0.02598914092828833),
		      Point(0.8114249947041546 , 0.09428750264792270),
		      Point(0.09428750264792270, 0.8114249947041546),
		      Point(0.09428750264792270, 0.09428750264792270),
		      Point(0.01072644996557060, 0.4946367750172147),
		      Point(0.4946367750172147 , 0.01072644996557060),
		      Point(0.4946367750172147 , 0.4946367750172147),
		      Point(0.5853132347709715 , 0.2073433826145142),
		      Point(0.2073433826145142 , 0.5853132347709715),
		      Point(0.2073433826145142 , 0.2073433826145142),
		      Point(0.1221843885990187 , 0.4389078057004907),
		      Point(0.4389078057004907 , 0.1221843885990187),
		      Point(0.4389078057004907 , 0.4389078057004907),
		      Point(0.6779376548825902 , 0.04484167758913055),
		      Point(0.6779376548825902 , 0.27722066752827925),
		      Point(0.04484167758913055, 0.6779376548825902),
		      Point(0.04484167758913055, 0.27722066752827925),
		      Point(0.27722066752827925, 0.6779376548825902),
		      Point(0.27722066752827925, 0.04484167758913055),
		      Point(0.8588702812826364 , 0.00000000000000000),
		      Point(0.8588702812826364 , 0.1411297187173636),
		      Point(0.0000000000000000 , 0.8588702812826364),
		      Point(0.0000000000000000 , 0.1411297187173636),
		      Point(0.1411297187173636 , 0.8588702812826364),
		      Point(0.1411297187173636 , 0.0000000000000000) };
    double wts[28] = { 0.08797730116222190, 0.008744311553736190,
		       0.008744311553736190, 0.008744311553736190,
		       0.03808157199393533, 0.03808157199393533,
		       0.03808157199393533, 0.01885544805613125,
		       0.01885544805613125, 0.01885544805613125,
		       0.07215969754474100, 0.07215969754474100,
		       0.07215969754474100, 0.06932913870553720,
		       0.06932913870553720, 0.06932913870553720,
		       0.04105631542928860, 0.04105631542928860,
		       0.04105631542928860, 0.04105631542928860,
		       0.04105631542928860, 0.04105631542928860,
		       0.007362383783300573, 0.007362383783300573,
		       0.007362383783300573, 0.007362383783300573,
		       0.007362383783300573, 0.007362383783300573 };
      
  } toms612_28;

  static struct {
    const int num = 37;
    const int dop = 13;
    Point pts[37] = {
      Point(0.333333333333333333333333333333, 0.333333333333333333333333333333),
      Point(0.950275662924105565450352089520, 0.024862168537947217274823955239),
      Point(0.024862168537947217274823955239, 0.950275662924105565450352089520),
      Point(0.024862168537947217274823955239, 0.024862168537947217274823955239),
      Point(0.171614914923835347556304795551, 0.414192542538082326221847602214),
      Point(0.414192542538082326221847602214, 0.171614914923835347556304795551),
      Point(0.414192542538082326221847602214, 0.414192542538082326221847602214),
      Point(0.539412243677190440263092985511, 0.230293878161404779868453507244),
      Point(0.230293878161404779868453507244, 0.539412243677190440263092985511),
      Point(0.230293878161404779868453507244, 0.230293878161404779868453507244),
      Point(0.772160036676532561750285570113, 0.113919981661733719124857214943),
      Point(0.113919981661733719124857214943, 0.772160036676532561750285570113),
      Point(0.113919981661733719124857214943, 0.113919981661733719124857214943),
      Point(0.009085399949835353883572964740, 0.495457300025082323058213517632),
      Point(0.495457300025082323058213517632, 0.009085399949835353883572964740),
      Point(0.495457300025082323058213517632, 0.495457300025082323058213517632),
      Point(0.062277290305886993497083640527, 0.468861354847056503251458179727),
      Point(0.468861354847056503251458179727, 0.062277290305886993497083640527),
      Point(0.468861354847056503251458179727, 0.468861354847056503251458179727),
      Point(0.022076289653624405142446876931, 0.851306504174348550389457672223),
      Point(0.022076289653624405142446876931, 0.126617206172027096933163647918),
      Point(0.851306504174348550389457672223, 0.022076289653624405142446876931),
      Point(0.851306504174348550389457672223, 0.126617206172027096933163647918),
      Point(0.126617206172027096933163647918, 0.022076289653624405142446876931),
      Point(0.126617206172027096933163647918, 0.851306504174348550389457672223),
      Point(0.018620522802520968955913511549, 0.689441970728591295496647976487),
      Point(0.018620522802520968955913511549, 0.291937506468887771754472382212),
      Point(0.689441970728591295496647976487, 0.018620522802520968955913511549),
      Point(0.689441970728591295496647976487, 0.291937506468887771754472382212),
      Point(0.291937506468887771754472382212, 0.018620522802520968955913511549),
      Point(0.291937506468887771754472382212, 0.689441970728591295496647976487),
      Point(0.096506481292159228736516560903, 0.635867859433872768286976979827),
      Point(0.096506481292159228736516560903, 0.267625659273967961282458816185),
      Point(0.635867859433872768286976979827, 0.096506481292159228736516560903),
      Point(0.635867859433872768286976979827, 0.267625659273967961282458816185),
      Point(0.267625659273967961282458816185, 0.096506481292159228736516560903),
      Point(0.267625659273967961282458816185, 0.635867859433872768286976979827) };
    double wts[37] = {
      0.051739766065744133555179145422, 0.008007799555564801597804123460,
      0.008007799555564801597804123460, 0.008007799555564801597804123460,
      0.046868898981821644823226732071, 0.046868898981821644823226732071,
      0.046868898981821644823226732071, 0.046590940183976487960361770070,
      0.046590940183976487960361770070, 0.046590940183976487960361770070,
      0.031016943313796381407646220131, 0.031016943313796381407646220131,
      0.031016943313796381407646220131, 0.010791612736631273623178240136,
      0.010791612736631273623178240136, 0.010791612736631273623178240136,
      0.032195534242431618819414482205, 0.032195534242431618819414482205,
      0.032195534242431618819414482205, 0.015445834210701583817692900053,
      0.015445834210701583817692900053, 0.015445834210701583817692900053,
      0.015445834210701583817692900053, 0.015445834210701583817692900053,
      0.015445834210701583817692900053, 0.017822989923178661888748319485,
      0.017822989923178661888748319485, 0.017822989923178661888748319485,
      0.017822989923178661888748319485, 0.017822989923178661888748319485,
      0.017822989923178661888748319485, 0.037038683681384627918546472190,
      0.037038683681384627918546472190, 0.037038683681384627918546472190,
      0.037038683681384627918546472190, 0.037038683681384627918546472190,
      0.037038683681384627918546472190 };
  } toms706_37;

  class RuleForTriangle {
  public:
    int num;
    int dop;
    Point* pts;
    double* wts;

    RuleForTriangle(int num_in, int dop_in, Point* pts_in, double* wts_in) :
      num(num_in), dop(dop_in), pts(pts_in), wts(wts_in) {};
  };

  static const std::vector<RuleForTriangle> ruleForTriangle
    = { RuleForTriangle(strang1.num,strang1.dop,strang1.pts,strang1.wts),
	RuleForTriangle(strang2.num,strang2.dop,strang2.pts,strang2.wts),
	RuleForTriangle(strang3.num,strang3.dop,strang3.pts,strang3.wts),
	RuleForTriangle(strang5.num,strang5.dop,strang5.pts,strang5.wts),
	RuleForTriangle(strang7.num,strang7.dop,strang7.pts,strang7.wts),
	RuleForTriangle(strang9.num,strang9.dop,strang9.pts,strang9.wts),
	RuleForTriangle(strang10.num,strang10.dop,strang10.pts,strang10.wts),
	RuleForTriangle(toms584_19.num,toms584_19.dop,toms584_19.pts,toms584_19.wts),
	RuleForTriangle(toms612_19.num,toms612_19.dop,toms612_19.pts,toms612_19.wts),
	RuleForTriangle(toms612_28.num,toms612_28.dop,toms612_28.pts,toms612_28.wts),
	RuleForTriangle(toms706_37.num,toms706_37.dop,toms706_37.pts,toms706_37.wts) };

  ////////////////////////////////////////////////////////////////////////////////
  // PolyQuadrature class
  //  Decomposes an nGon into nGon-2 triangles by criss-crossing the element
  //  One sets the rule type on creation
  //  
  //  One first sets the element, and then the points are computed for use.
  ////////////////////////////////////////////////////////////////////////////////

  class PolyQuadrature
  {
  private:
    polymesh::PolyElement* my_element = nullptr;

    int my_desired_dop; // desired degree of precision
    int my_dop; // actual degree of precision
    int my_refinement_level;
    int my_rule;

    // Reference triangle
    int num_pts_ref;
    Point* my_pts_ref;
    double* my_wts_ref;

    // Polygon
    int num_pts;
    Point* my_pts = nullptr;
    double* my_wts = nullptr;

    void set_rule(int desired_dop);
    void set_element(int refinement_level, polymesh::PolyElement* element);
    
  public:
    PolyQuadrature(int desired_dop=2, int refinement_level=0, polymesh::PolyElement* element=nullptr) {
      set_rule(desired_dop); set_element(refinement_level, element); };
    ~PolyQuadrature();

    void set(int desired_dop=2, int refinement_level=0, polymesh::PolyElement* element=nullptr) {
      set_rule(desired_dop); set_element(refinement_level, element); };
    void setRule(int desired_dop) { set_rule(desired_dop); set_element(my_refinement_level, my_element); };
    void setElement(int refinement_level=0, polymesh::PolyElement* element = nullptr) { set_element(refinement_level, element); };

    polymesh::PolyElement* elementPtr() const { return my_element; };
    int quadratureRule() const { return my_rule; };

    int num() const { return num_pts; };
    Point& pt(int i) const { return my_pts[i]; };
    double wt(int i) const { return my_wts[i]; };
    Point* pts() const { return my_pts; };
    double* wts() const { return my_wts; };

    int desiredDOP() const { return my_desired_dop; };
    int dop() const { return my_dop; };
    bool isDesiredDOP() { return my_desired_dop == my_dop; };
    bool isAtLeastDesiredDOP() { return my_desired_dop <= my_dop; };
  }; 
 
 
  ////////////////////////////////////////////////////////////////////////////////
  // Quadrature rules on reference edge -1 1
  //    num = number of points in the rule = degree of precision of the rule
  //    pts = quadrature points (using double)
  //    wts = quadrature weights
  //    Note that pts are on (-1,1), but wts are for pts scaling to (0,1)
  ////////////////////////////////////////////////////////////////////////////////

  static struct {
    const int num = 1;
    double pts[1] = { 0 };
    double wts[1] = { 2 };
  } gauleg1;

  static struct {
    const int num = 2;
    double pts[2] = { -1/sqrt(3), 1/sqrt(3) };
    double wts[2] = { 1, 1 };
  } gauleg2;

  static struct {
    const int num = 3;
    double pts[3] = { -sqrt(0.6), 0, sqrt(0.6) };
    double wts[3] = { double(5)/double(9), double(8)/double(9), 
                      double(5)/double(9) };
  } gauleg3;

  static struct {
    const int num = 4;
    double pts[4] = { -0.861136311594052575224, -0.3399810435848562648027,
                      0.3399810435848562648027, 0.861136311594052575224 };
    double wts[4] = { 0.3478548451374538573731, 0.6521451548625461426269, 
                      0.6521451548625461426269, 0.3478548451374538573731 };
  } gauleg4;

  static struct {
    const int num = 5;
    double pts[5] = { -0.9061798459386639927976, -0.5384693101056830910363,
                      0, 
                      0.5384693101056830910363, 0.9061798459386639927976 };
    double wts[5] = { 0.2369268850561890875143, 0.4786286704993664680413, 
                      0.5688888888888888888889, 
                      0.4786286704993664680413, 0.2369268850561890875143 };
  } gauleg5;

  static struct {
    const int num = 6;
    double pts[6] = { -0.9324695142031520278123, -0.661209386466264513661,
                      -0.2386191860831969086305, 0.238619186083196908631, 
                      0.661209386466264513661, 0.9324695142031520278123 };
    double wts[6] = { 0.1713244923791703450403, 0.3607615730481386075698, 
                      0.4679139345726910473899, 0.4679139345726910473899,
                      0.3607615730481386075698, 0.1713244923791703450403 };
  } gauleg6;

  static struct {
    const int num = 7;
    double pts[7] = { -0.9491079123427585245262, -0.7415311855993944398639,
                      -0.4058451513773971669066, 0, 
                      0.4058451513773971669066, 0.7415311855993944398639,
                      0.9491079123427585245262 };
    double wts[7] = { 0.1294849661688696932706, 0.2797053914892766679015, 
                      0.38183005050511894495, 0.417959183673469387755,
                      0.38183005050511894495, 0.279705391489276667901,
                      0.129484966168869693271 };
  } gauleg7;


  static struct {
    const int num = 8;
    double pts[8] = { -0.9602898564975362316836, -0.7966664774136267395916,
                      -0.5255324099163289858177, -0.1834346424956498049395, 
                      0.1834346424956498049395, 0.5255324099163289858177,
                      0.7966664774136267395916, 0.9602898564975362316836 };
    double wts[8] = { 0.1012285362903762591525, 0.2223810344533744705444, 
                      0.313706645877887287338, 0.3626837833783619829652,
                      0.3626837833783619829652, 0.313706645877887287338,
                      0.222381034453374470544, 0.1012285362903762591525 };
  } gauleg8;

  static struct {
    const int num = 9;
    double pts[9] = { -0.9681602395076260898356, -0.8360311073266357942994,
                      -0.6133714327005903973087, -0.3242534234038089290385,
                      0,
                      0.3242534234038089290385, 0.6133714327005903973087,
                      0.8360311073266357942994, 0.9681602395076260898356 };
    double wts[9] = { 0.0812743883615744119719, 0.1806481606948574040585, 
                      0.2606106964029354623187, 0.312347077040002840069,
                      0.330239355001259763165,
                      0.312347077040002840069, 0.260610696402935462319,
                      0.1806481606948574040585, 0.081274388361574411972 };
  } gauleg9;

  static struct {
    const int num = 10;
    double pts[10] = { -0.973906528517171720078, -0.8650633666889845107321,
                       -0.6794095682990244062343, -0.4333953941292471907993,
                       -0.1488743389816312108848, 0.1488743389816312108848,
                        0.4333953941292471907993, 0.6794095682990244062343,
                        0.8650633666889845107321, 0.973906528517171720078 };
    double wts[10] = { 0.0666713443086881375936, 0.149451349150580593146,
                       0.219086362515982043996, 0.2692667193099963550912,
                       0.2955242247147528701739, 0.295524224714752870174,
                       0.269266719309996355091, 0.2190863625159820439955,
                       0.1494513491505805931458, 0.0666713443086881375936 };
  } gauleg10;

  static struct {
    const int num = 11;
    double pts[11] = { -0.9782286581460569928039, -0.8870625997680952990752,
                       -0.7301520055740493240934, -0.5190961292068118159257,
                       -0.2695431559523449723315, 0,
                       0.269543155952344972332, 0.5190961292068118159257,
                       0.7301520055740493240934, 0.887062599768095299075,
                       0.9782286581460569928039 };
    double wts[11] = { 0.0556685671161736664828, 0.1255803694649046246347,
                       0.1862902109277342514261, 0.2331937645919904799185,
                       0.2628045445102466621807, 0.2729250867779006307145,
                       0.262804544510246662181, 0.2331937645919904799185,
                       0.1862902109277342514261, 0.1255803694649046246347,
                       0.055668567116173666483 };
  } gauleg11;

  static struct {
    const int num = 12;
    double pts[12] = { -0.9815606342467192506906, -0.9041172563704748566785,
                       -0.769902674194304687037, -0.5873179542866174472967,
                       -0.3678314989981801937527, -0.1252334085114689154724,
                       0.1252334085114689154724, 0.3678314989981801937527,
                       0.5873179542866174472967, 0.7699026741943046870369,
                       0.9041172563704748566785, 0.9815606342467192506906 };
    double wts[12] = { 0.0471753363865118271946, 0.1069393259953184309603,
                       0.1600783285433462263347, 0.2031674267230659217491,
                       0.233492536538354808761, 0.2491470458134027850006,
                       0.2491470458134027850006, 0.233492536538354808761,
                       0.203167426723065921749, 0.160078328543346226335,
                       0.1069393259953184309603, 0.0471753363865118271946 };
  } gauleg12;

  static struct {
    const int num = 13;
    double pts[13] = { -0.9841830547185881494728, -0.9175983992229779652066,
                       -0.8015780907333099127942, -0.642349339440340220644,
                       -0.4484927510364468528779, -0.2304583159551347940655,
                       0,
                       0.2304583159551347940655, 0.448492751036446852878,
                       0.642349339440340220644, 0.8015780907333099127942,
                       0.9175983992229779652066, 0.9841830547185881494728 };
    double wts[13] = { 0.04048400476531587952, 0.0921214998377284479144,
                       0.1388735102197872384636, 0.1781459807619457382801,
                       0.2078160475368885023125, 0.2262831802628972384121,
                       0.2325515532308739101946, 
                       0.2262831802628972384121, 0.2078160475368885023125, 
                       0.17814598076194573828, 0.138873510219787238464, 
                       0.0921214998377284479144, 0.04048400476531587952 };
  } gauleg13;

  class RuleForEdge {
  public:
    int num;
    double* pts;
    double* wts;

    RuleForEdge(int num_in,  double* pts_in, double* wts_in) :
      num(num_in), pts(pts_in), wts(wts_in) {};
  };

  static const std::vector<RuleForEdge> ruleForEdge
    = { RuleForEdge(gauleg1.num, gauleg1.pts, gauleg1.wts),
	      RuleForEdge(gauleg2.num, gauleg2.pts, gauleg2.wts),
        RuleForEdge(gauleg3.num, gauleg3.pts, gauleg3.wts),
        RuleForEdge(gauleg4.num, gauleg4.pts, gauleg4.wts),
        RuleForEdge(gauleg5.num, gauleg5.pts, gauleg5.wts),
        RuleForEdge(gauleg6.num, gauleg6.pts, gauleg6.wts),
        RuleForEdge(gauleg7.num, gauleg7.pts, gauleg7.wts),
        RuleForEdge(gauleg8.num, gauleg8.pts, gauleg8.wts),
        RuleForEdge(gauleg9.num, gauleg9.pts, gauleg9.wts),
        RuleForEdge(gauleg10.num, gauleg10.pts, gauleg10.wts),
        RuleForEdge(gauleg11.num, gauleg11.pts, gauleg11.wts),
        RuleForEdge(gauleg12.num, gauleg12.pts, gauleg12.wts),
        RuleForEdge(gauleg13.num, gauleg13.pts, gauleg13.wts),
 };

  ////////////////////////////////////////////////////////////////////////////////
  // PolyEdgeQuadrature class
  //
  //  One sets the rule type on creation
  //  
  //  One first sets the edge, and then the points are computed for use.
  ////////////////////////////////////////////////////////////////////////////////

  class PolyEdgeQuadrature
  {
  private:
    polymesh::Edge* my_edge = nullptr;

    int my_desired_dop; // desired degree of precision
    int my_dop; // actual degree of precision
    int my_rule;

    // Reference edge
    int num_pts;
    double* my_pts_ref;
    double* my_wts_ref;

    // Edge
    Point* my_pts = nullptr;
    double* my_wts = nullptr;
    void set_rule(int desired_dop);
    void set_edge(polymesh::Edge* edge);
    
  public:
    PolyEdgeQuadrature(int desired_dop=2, polymesh::Edge* edge=nullptr) {
      set_rule(desired_dop); set_edge(edge); };
    ~PolyEdgeQuadrature();

    void set(int desired_dop=2, polymesh::Edge* edge=nullptr) {
      set_rule(desired_dop); set_edge(edge); };
    void setRule(int desired_dop) { set_rule(desired_dop); set_edge(my_edge); };
    void setEdge(polymesh::Edge* edge) { set_edge(edge); };

    polymesh::Edge* edgePtr() const { return my_edge; };
    int edgeQuadratureRule() const { return my_rule; };

    int num() const { return num_pts; };
    Point& pt(int i) const { return my_pts[i]; };
    double wt(int i) const { return my_wts[i]; };
    Point* pts() const { return my_pts; };
    double* wts() const { return my_wts; };

    int desiredDOP() const { return my_desired_dop; };
    int dop() const { return my_dop; };
    bool isDesiredDOP() { return my_desired_dop == my_dop; };
    bool isAtLeastDesiredDOP() { return my_desired_dop <= my_dop; };
  }; 


  void testPolyQuadrature(polymesh::PolyMesh* mesh, int refinement_level, double eps=1e-6,
			  int toDOP=ruleForTriangle[ruleForTriangle.size()-1].dop);

};

#endif
