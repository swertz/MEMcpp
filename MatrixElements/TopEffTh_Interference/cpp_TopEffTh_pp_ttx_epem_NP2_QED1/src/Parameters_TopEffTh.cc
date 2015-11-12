//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "Parameters_TopEffTh.h"

using namespace std; 

// Constructor
Parameters_TopEffTh::Parameters_TopEffTh(const string param_card, const bool
    verbose):
slha(SLHAReader(param_card))
{
  setIndependentParameters(); 
  setIndependentCouplings(); 

  if(verbose)
  {
    printIndependentParameters(); 
    printIndependentCouplings(); 
  }

}

void Parameters_TopEffTh::setIndependentParameters()
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  mdl_WH = slha.get_block_entry("decay", 25, 4.070000e-03); 
  mdl_WW = slha.get_block_entry("decay", 24, 2.085000e+00); 
  mdl_WZ = slha.get_block_entry("decay", 23, 2.495200e+00); 
  mdl_WT = slha.get_block_entry("decay", 6, 1.508336e+00); 
  mdl_ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  mdl_ymt = slha.get_block_entry("yukawa", 6, 1.720000e+02); 
  mdl_ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00); 
  aS = slha.get_block_entry("sminputs", 3, 1.184000e-01); 
  mdl_Gf = slha.get_block_entry("sminputs", 2, 1.166370e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.279000e+02); 
  mdl_MH = slha.get_block_entry("mass", 25, 1.250000e+02); 
  mdl_MZ = slha.get_block_entry("mass", 23, 9.118760e+01); 
  mdl_MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  mdl_MT = slha.get_block_entry("mass", 6, 1.720000e+02); 
  mdl_MB = slha.get_block_entry("mass", 5, 4.700000e+00); 
  mdl_C1qt = slha.get_block_entry("fourfermion", 8, 1.000000e+00); 
  mdl_C1qd = slha.get_block_entry("fourfermion", 7, 1.000000e+00); 
  mdl_C1qu = slha.get_block_entry("fourfermion", 6, 1.000000e+00); 
  mdl_C8dt = slha.get_block_entry("fourfermion", 5, 1.000000e+00); 
  mdl_C8ut = slha.get_block_entry("fourfermion", 4, 1.000000e+00); 
  mdl_C83qq = slha.get_block_entry("fourfermion", 3, 1.000000e+00); 
  mdl_C81qq = slha.get_block_entry("fourfermion", 2, 1.000000e+00); 
  mdl_C13qq = slha.get_block_entry("fourfermion", 1, 1.000000e+00); 
  mdl_CphiG = slha.get_block_entry("dim6", 9, 1.000000e+00); 
  mdl_CG = slha.get_block_entry("dim6", 8, 1.000000e+00); 
  mdl_ICtG = slha.get_block_entry("dim6", 7, 1.000000e+00); 
  mdl_RCtG = slha.get_block_entry("dim6", 6, 1.000000e+00); 
  mdl_ICtW = slha.get_block_entry("dim6", 5, 1.000000e+00); 
  mdl_RCtW = slha.get_block_entry("dim6", 4, 1.000000e+00); 
  mdl_IC3phiq = slha.get_block_entry("dim6", 3, 1.000000e+00); 
  mdl_RC3phiq = slha.get_block_entry("dim6", 2, 1.000000e+00); 
  mdl_Lambda = slha.get_block_entry("dim6", 1, 1.000000e+03); 
  mdl_complexi = std::complex<double> (0., 1.); 
  mdl_C3phiq = mdl_complexi * mdl_IC3phiq + mdl_RC3phiq; 
  mdl_CtG = mdl_complexi * mdl_ICtG + mdl_RCtG; 
  mdl_CtW = mdl_complexi * mdl_ICtW + mdl_RCtW; 
  mdl_MZ__exp__2 = ((mdl_MZ) * (mdl_MZ)); 
  mdl_MZ__exp__4 = ((mdl_MZ) * (mdl_MZ) * (mdl_MZ) * (mdl_MZ)); 
  mdl_sqrt__2 = sqrt(2.); 
  mdl_MH__exp__2 = ((mdl_MH) * (mdl_MH)); 
  mdl_Lambda__exp__2 = ((mdl_Lambda) * (mdl_Lambda)); 
  mdl_conjg__C3phiq = conj(mdl_C3phiq); 
  mdl_conjg__CtG = conj(mdl_CtG); 
  mdl_conjg__CtW = conj(mdl_CtW); 
  mdl_aEW = 1./aEWM1; 
  mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
      mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
  mdl_sqrt__aEW = sqrt(mdl_aEW); 
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt(M_PI); 
  mdl_MW__exp__2 = ((mdl_MW) * (mdl_MW)); 
  mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2; 
  mdl_cw = sqrt(1. - mdl_sw2); 
  mdl_sqrt__sw2 = sqrt(mdl_sw2); 
  mdl_sw = mdl_sqrt__sw2; 
  mdl_g1 = mdl_ee/mdl_cw; 
  mdl_gw = mdl_ee/mdl_sw; 
  mdl_vev = (2. * mdl_MW * mdl_sw)/mdl_ee; 
  mdl_vev__exp__2 = ((mdl_vev) * (mdl_vev)); 
  mdl_lam = mdl_MH__exp__2/(2. * mdl_vev__exp__2); 
  mdl_yb = (mdl_ymb * mdl_sqrt__2)/mdl_vev; 
  mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_vev; 
  mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_vev; 
  mdl_muH = sqrt(mdl_lam * mdl_vev__exp__2); 
  mdl_ee__exp__2 = ((mdl_ee) * (mdl_ee)); 
  mdl_sw__exp__2 = ((mdl_sw) * (mdl_sw)); 
  mdl_cw__exp__2 = ((mdl_cw) * (mdl_cw)); 
}
void Parameters_TopEffTh::setIndependentCouplings()
{
  GC_1 = -(mdl_ee * mdl_complexi)/3.; 
  GC_2 = (2. * mdl_ee * mdl_complexi)/3.; 
  GC_57 = (mdl_ee * mdl_complexi)/(mdl_sw * mdl_sqrt__2); 
  GC_58 = -(mdl_cw * mdl_ee * mdl_complexi)/(2. * mdl_sw); 
  GC_59 = (mdl_cw * mdl_ee * mdl_complexi)/(2. * mdl_sw); 
  GC_69 = -(mdl_ee * mdl_complexi * mdl_sw)/(6. * mdl_cw); 
  GC_78 = (2. * mdl_CphiG * mdl_complexi * mdl_vev)/mdl_Lambda__exp__2; 
  GC_79 = (mdl_CtG * mdl_complexi * mdl_vev * mdl_sqrt__2)/mdl_Lambda__exp__2; 
  GC_81 = (mdl_CtW * mdl_cw * mdl_complexi * mdl_vev)/(mdl_Lambda__exp__2 *
      mdl_sqrt__2);
  GC_91 = (mdl_CtW * mdl_complexi * mdl_sw * mdl_vev)/(mdl_Lambda__exp__2 *
      mdl_sqrt__2);
  GC_95 = -((mdl_complexi * mdl_yt)/mdl_sqrt__2); 
  GC_138 = (mdl_C3phiq * mdl_cw * mdl_ee * mdl_complexi * mdl_vev__exp__2)/(4.
      * mdl_Lambda__exp__2 * mdl_sw) + (mdl_C3phiq * mdl_ee * mdl_complexi *
      mdl_sw * mdl_vev__exp__2)/(4. * mdl_cw * mdl_Lambda__exp__2) + (mdl_cw *
      mdl_ee * mdl_complexi * mdl_vev__exp__2 * mdl_conjg__C3phiq)/(4. *
      mdl_Lambda__exp__2 * mdl_sw) + (mdl_ee * mdl_complexi * mdl_sw *
      mdl_vev__exp__2 * mdl_conjg__C3phiq)/(4. * mdl_cw * mdl_Lambda__exp__2);
  GC_145 = (mdl_complexi * mdl_vev * mdl_conjg__CtG *
      mdl_sqrt__2)/mdl_Lambda__exp__2;
  GC_166 = (mdl_cw * mdl_complexi * mdl_vev *
      mdl_conjg__CtW)/(mdl_Lambda__exp__2 * mdl_sqrt__2);
  GC_170 = (mdl_complexi * mdl_sw * mdl_vev *
      mdl_conjg__CtW)/(mdl_Lambda__exp__2 * mdl_sqrt__2);
  GC_10_81 = (mdl_C81qq * mdl_complexi)/mdl_Lambda__exp__2;
  GC_11_81 = (mdl_C81qq * mdl_complexi)/mdl_Lambda__exp__2;
  GC_10_83 = - (mdl_C83qq * mdl_complexi)/mdl_Lambda__exp__2;
  GC_11_83 = (mdl_C83qq * mdl_complexi)/mdl_Lambda__exp__2;
  GC_12 = -((mdl_C13qq * mdl_complexi)/mdl_Lambda__exp__2); 
  GC_13 = (mdl_C13qq * mdl_complexi)/mdl_Lambda__exp__2; 
  GC_15 = (mdl_C1qd * mdl_complexi)/mdl_Lambda__exp__2; 
  GC_16 = (mdl_C1qt * mdl_complexi)/mdl_Lambda__exp__2; 
  GC_17 = (mdl_C1qu * mdl_complexi)/mdl_Lambda__exp__2; 
  GC_23 = (mdl_C81qq * mdl_complexi)/mdl_Lambda__exp__2; 
  GC_24 = -((mdl_C83qq * mdl_complexi)/mdl_Lambda__exp__2); 
  GC_26 = (mdl_C8dt * mdl_complexi)/mdl_Lambda__exp__2; 
  GC_27 = (mdl_C8ut * mdl_complexi)/mdl_Lambda__exp__2; 
  GC_28 = (-6. * mdl_CG)/mdl_Lambda__exp__2; 
}
void Parameters_TopEffTh::setDependentParameters()
{
  mdl_sqrt__aS = sqrt(aS); 
  G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
  mdl_G__exp__2 = ((G) * (G)); 
  mdl_G__exp__3 = ((G) * (G) * (G)); 
}
void Parameters_TopEffTh::setDependentCouplings()
{
  GC_146 = -((G * mdl_vev * mdl_conjg__CtG * mdl_sqrt__2)/mdl_Lambda__exp__2); 
  GC_7 = mdl_complexi * G; 
  GC_6 = -G; 
  GC_85 = -((mdl_CtG * G * mdl_vev * mdl_sqrt__2)/mdl_Lambda__exp__2); 
}

// Routines for printing out parameters
void Parameters_TopEffTh::printIndependentParameters()
{
  cout <<  "TopEffTh model parameters independent of event kinematics:" <<
      endl;
  cout << setw(20) <<  "mdl_WH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WH << endl;
  cout << setw(20) <<  "mdl_WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WW << endl;
  cout << setw(20) <<  "mdl_WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WZ << endl;
  cout << setw(20) <<  "mdl_WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WT << endl;
  cout << setw(20) <<  "mdl_ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymtau << endl;
  cout << setw(20) <<  "mdl_ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymt << endl;
  cout << setw(20) <<  "mdl_ymb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymb << endl;
  cout << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS << endl;
  cout << setw(20) <<  "mdl_Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Gf << endl;
  cout << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1 << endl;
  cout << setw(20) <<  "mdl_MH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MH << endl;
  cout << setw(20) <<  "mdl_MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MZ << endl;
  cout << setw(20) <<  "mdl_MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MTA << endl;
  cout << setw(20) <<  "mdl_MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MT << endl;
  cout << setw(20) <<  "mdl_MB " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MB << endl;
  cout << setw(20) <<  "mdl_C1qt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_C1qt << endl;
  cout << setw(20) <<  "mdl_C1qd " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_C1qd << endl;
  cout << setw(20) <<  "mdl_C1qu " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_C1qu << endl;
  cout << setw(20) <<  "mdl_C8dt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_C8dt << endl;
  cout << setw(20) <<  "mdl_C8ut " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_C8ut << endl;
  cout << setw(20) <<  "mdl_C83qq " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_C83qq << endl;
  cout << setw(20) <<  "mdl_C81qq " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_C81qq << endl;
  cout << setw(20) <<  "mdl_C13qq " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_C13qq << endl;
  cout << setw(20) <<  "mdl_CphiG " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_CphiG << endl;
  cout << setw(20) <<  "mdl_CG " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_CG << endl;
  cout << setw(20) <<  "mdl_ICtG " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ICtG << endl;
  cout << setw(20) <<  "mdl_RCtG " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_RCtG << endl;
  cout << setw(20) <<  "mdl_ICtW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ICtW << endl;
  cout << setw(20) <<  "mdl_RCtW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_RCtW << endl;
  cout << setw(20) <<  "mdl_IC3phiq " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_IC3phiq << endl;
  cout << setw(20) <<  "mdl_RC3phiq " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_RC3phiq << endl;
  cout << setw(20) <<  "mdl_Lambda " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_Lambda << endl;
  cout << setw(20) <<  "mdl_complexi " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_complexi << endl;
  cout << setw(20) <<  "mdl_C3phiq " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_C3phiq << endl;
  cout << setw(20) <<  "mdl_CtG " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_CtG << endl;
  cout << setw(20) <<  "mdl_CtW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_CtW << endl;
  cout << setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__2 << endl;
  cout << setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__4 << endl;
  cout << setw(20) <<  "mdl_sqrt__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__2 << endl;
  cout << setw(20) <<  "mdl_MH__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__2 << endl;
  cout << setw(20) <<  "mdl_Lambda__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_Lambda__exp__2 << endl;
  cout << setw(20) <<  "mdl_conjg__C3phiq " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__C3phiq << endl;
  cout << setw(20) <<  "mdl_conjg__CtG " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CtG << endl;
  cout << setw(20) <<  "mdl_conjg__CtW " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CtW << endl;
  cout << setw(20) <<  "mdl_aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_aEW << endl;
  cout << setw(20) <<  "mdl_MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MW << endl;
  cout << setw(20) <<  "mdl_sqrt__aEW " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__aEW << endl;
  cout << setw(20) <<  "mdl_ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ee << endl;
  cout << setw(20) <<  "mdl_MW__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw2 << endl;
  cout << setw(20) <<  "mdl_cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cw << endl;
  cout << setw(20) <<  "mdl_sqrt__sw2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__sw2 << endl;
  cout << setw(20) <<  "mdl_sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw << endl;
  cout << setw(20) <<  "mdl_g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_g1 << endl;
  cout << setw(20) <<  "mdl_gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gw << endl;
  cout << setw(20) <<  "mdl_vev " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_vev << endl;
  cout << setw(20) <<  "mdl_vev__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_vev__exp__2 << endl;
  cout << setw(20) <<  "mdl_lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_lam << endl;
  cout << setw(20) <<  "mdl_yb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yb << endl;
  cout << setw(20) <<  "mdl_yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yt << endl;
  cout << setw(20) <<  "mdl_ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ytau << endl;
  cout << setw(20) <<  "mdl_muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_muH << endl;
  cout << setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ee__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sw__exp__2 << endl;
  cout << setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cw__exp__2 << endl;
}
void Parameters_TopEffTh::printIndependentCouplings()
{
  cout <<  "TopEffTh model couplings independent of event kinematics:" << endl; 
  cout << setw(20) <<  "GC_1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_1 << endl;
  cout << setw(20) <<  "GC_2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_2 << endl;
  cout << setw(20) <<  "GC_57 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_57 << endl;
  cout << setw(20) <<  "GC_58 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_58 << endl;
  cout << setw(20) <<  "GC_59 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_59 << endl;
  cout << setw(20) <<  "GC_69 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_69 << endl;
  cout << setw(20) <<  "GC_78 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_78 << endl;
  cout << setw(20) <<  "GC_79 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_79 << endl;
  cout << setw(20) <<  "GC_81 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_81 << endl;
  cout << setw(20) <<  "GC_91 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_91 << endl;
  cout << setw(20) <<  "GC_95 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_95 << endl;
  cout << setw(20) <<  "GC_138 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_138 << endl;
  cout << setw(20) <<  "GC_145 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_145 << endl;
  cout << setw(20) <<  "GC_166 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_166 << endl;
  cout << setw(20) <<  "GC_170 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_170 << endl;
  cout << setw(20) <<  "GC_10_81 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_10_81 << endl;
  cout << setw(20) <<  "GC_11_81 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_11_81 << endl;
  cout << setw(20) <<  "GC_10_83 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_10_83 << endl;
  cout << setw(20) <<  "GC_11_83 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_11_83 << endl;
  cout << setw(20) <<  "GC_12 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_12 << endl;
  cout << setw(20) <<  "GC_13 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_13 << endl;
  cout << setw(20) <<  "GC_15 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_15 << endl;
  cout << setw(20) <<  "GC_16 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_16 << endl;
  cout << setw(20) <<  "GC_17 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_17 << endl;
  cout << setw(20) <<  "GC_23 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_23 << endl;
  cout << setw(20) <<  "GC_24 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_24 << endl;
  cout << setw(20) <<  "GC_26 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_26 << endl;
  cout << setw(20) <<  "GC_27 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_27 << endl;
  cout << setw(20) <<  "GC_28 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_28 << endl;
}
void Parameters_TopEffTh::printDependentParameters()
{
  cout <<  "TopEffTh model parameters dependent on event kinematics:" << endl; 
  cout << setw(20) <<  "mdl_sqrt__aS " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__aS << endl;
  cout << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  cout << setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__2 << endl;
  cout << setw(20) <<  "mdl_G__exp__3 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__3 << endl;
}
void Parameters_TopEffTh::printDependentCouplings()
{
  cout <<  "TopEffTh model couplings dependent on event kinematics:" << endl; 
  cout << setw(20) <<  "GC_146 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_146 << endl;
  cout << setw(20) <<  "GC_7 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_7 << endl;
  cout << setw(20) <<  "GC_6 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_6 << endl;
  cout << setw(20) <<  "GC_85 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_85 << endl;
}


