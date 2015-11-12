//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <string> 
#include <utility> 
#include <vector> 
#include <map> 

#include "cpp_TopEffTh_pp_ttx_epem_NP2_QED1.h"
#include "HelAmps_TopEffTh.h"
#include "process_base_classes.h"

using namespace MG5_TopEffTh; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t t~ NP=2 QED=1 @1
// *   Decay: t > w+ b NP=0
// *     Decay: w+ > e+ ve NP=0
// *   Decay: t~ > w- b~ NP=0
// *     Decay: w- > e- ve~ NP=0
// Process: u u~ > t t~ NP=2 QED=1 @1
// *   Decay: t > w+ b NP=0
// *     Decay: w+ > e+ ve NP=0
// *   Decay: t~ > w- b~ NP=0
// *     Decay: w- > e- ve~ NP=0
// Process: c c~ > t t~ NP=2 QED=1 @1
// *   Decay: t > w+ b NP=0
// *     Decay: w+ > e+ ve NP=0
// *   Decay: t~ > w- b~ NP=0
// *     Decay: w- > e- ve~ NP=0
// Process: d d~ > t t~ NP=2 QED=1 @1
// *   Decay: t > w+ b NP=0
// *     Decay: w+ > e+ ve NP=0
// *   Decay: t~ > w- b~ NP=0
// *     Decay: w- > e- ve~ NP=0
// Process: s s~ > t t~ NP=2 QED=1 @1
// *   Decay: t > w+ b NP=0
// *     Decay: w+ > e+ ve NP=0
// *   Decay: t~ > w- b~ NP=0
// *     Decay: w- > e- ve~ NP=0

//--------------------------------------------------------------------------
// Initialize process.

cpp_TopEffTh_pp_ttx_epem_NP2_QED1::cpp_TopEffTh_pp_ttx_epem_NP2_QED1(Parameters_TopEffTh &params):
params(params)
{
  // Set external particle masses for this matrix element
  mME.push_back(params.ZERO); 
  mME.push_back(params.ZERO); 
  mME.push_back(params.ZERO); 
  mME.push_back(params.ZERO); 
  mME.push_back(params.mdl_MB); 
  mME.push_back(params.ZERO); 
  mME.push_back(params.ZERO); 
  mME.push_back(params.mdl_MB); 

  mapFinalStates[{-11, 12, 5, 11, -12, -5}] = 
  {
    {
      &cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_gg_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex, 
      false, 
      {
        std::make_pair(21, 21)
      }
      , 
      256, 
      256
    }
    , 
    {
      &cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_uux_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex, 
      true, 
      {
        std::make_pair(2, -2)
      }
      , 
      256, 
      36
    }
    , 
    {
      &cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_ccx_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex, 
      true, 
      {
        std::make_pair(4, -4)
      }
      , 
      256, 
      36
    }
    , 
    {
      &cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_ddx_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex, 
      true, 
      {
        std::make_pair(1, -1)
      }
      , 
      256, 
      36
    }
    , 
    {
      &cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_ssx_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex, 
      true, 
      {
        std::make_pair(3, -3)
      }
      , 
      256, 
      36
    }
  }; 

}

//--------------------------------------------------------------------------
// Evaluate |M|^2, return a map of final states

std::map < std::pair < int, int > , double >
    cpp_TopEffTh_pp_ttx_epem_NP2_QED1::sigmaKin(const std::vector <
    std::vector<double> > &initialMomenta, const std::vector < std::pair < int,
    std::vector<double> > > &finalState)
{

  // Set initial particle momenta
  momenta[0] = (double * ) (&initialMomenta[0][0]); 
  momenta[1] = (double * ) (&initialMomenta[1][0]); 

  // Suppose final particles are passed in the "correct" order
  std::vector<int> selectedFinalState(8 - 2); 
  size_t index = 2; 
  for (auto const &finalPart: finalState)
  {
    selectedFinalState[index - 2] = finalPart.first; 
    momenta[index] = (double * ) (&finalPart.second[0]); 
    index++; 
  }

  // Set the event specific parameters
  params.setDependentParameters(); 
  params.setDependentCouplings(); 

  // Initialise result object
  std::map < std::pair < int, int > , double > result; 

  // Define permutation
  int perm[8]; 
  for(int i = 0; i < 8; i++ )
  {
    perm[i] = i; 
  }

  for(auto &me: mapFinalStates[selectedFinalState])
  {
    double me_sum = 0; 
    double me_mirror_sum = 0; 
    for(int ihel = 0; ihel < 256; ihel++ )
    {
      if(me.goodHel[ihel])
      {
        double sum = 0.; 
        calculate_wavefunctions(perm, helicities[ihel]); 
        double meTemp = (this->* (me.meCall))(); 
        sum += meTemp; 
        me_sum += meTemp/me.denominator; 

        if(me.hasMirrorProcess)
        {
          perm[0] = 1; 
          perm[1] = 0; 
          // Calculate wavefunctions
          calculate_wavefunctions(perm, helicities[ihel]); 
          // Mirror back
          perm[0] = 0; 
          perm[1] = 1; 
          meTemp = (this->* (me.meCall))(); 
          sum += meTemp; 
          me_mirror_sum += meTemp/me.denominator; 
        }

        if( !sum)
          me.goodHel[ihel] = false;
      }
    }

    for (auto const &initialState: me.initialStates)
    {
      result[initialState] = me_sum; 
      if (me.hasMirrorProcess)
        result[std::make_pair(initialState.second, initialState.first)] =
            me_mirror_sum;
    }
  }


  return result; 
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void cpp_TopEffTh_pp_ttx_epem_NP2_QED1::calculate_wavefunctions(const int
    perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  // Calculate all wavefunctions
  static std::complex<double> w[26][18]; 

  vxxxxx(&momenta[perm[0]][0], mME[0], hel[0], -1, w[0]); 
  vxxxxx(&momenta[perm[1]][0], mME[1], hel[1], -1, w[1]); 
  ixxxxx(&momenta[perm[2]][0], mME[2], hel[2], -1, w[2]); 
  oxxxxx(&momenta[perm[3]][0], mME[3], hel[3], +1, w[3]); 
  FFV2_3(w[2], w[3], params.GC_57, params.mdl_MW, params.mdl_WW, w[4]); 
  oxxxxx(&momenta[perm[4]][0], mME[4], hel[4], +1, w[5]); 
  FFV2_1(w[5], w[4], params.GC_57, params.mdl_MT, params.mdl_WT, w[6]); 
  oxxxxx(&momenta[perm[5]][0], mME[5], hel[5], +1, w[7]); 
  ixxxxx(&momenta[perm[6]][0], mME[6], hel[6], -1, w[8]); 
  FFV2_3(w[8], w[7], params.GC_57, params.mdl_MW, params.mdl_WW, w[9]); 
  ixxxxx(&momenta[perm[7]][0], mME[7], hel[7], -1, w[10]); 
  FFV2_2(w[10], w[9], params.GC_57, params.mdl_MT, params.mdl_WT, w[11]); 
  VVS2_3(w[0], w[1], params.GC_78, params.mdl_MH, params.mdl_WH, w[12]); // phiG
  VVV2P0_1(w[0], w[1], params.GC_28, params.ZERO, params.ZERO, w[13]); // OG 
  VVV1P0_1(w[0], w[1], params.GC_6, params.ZERO, params.ZERO, w[14]); 
  FFV5_1(w[6], w[0], params.GC_7, params.mdl_MT, params.mdl_WT, w[15]); 
  FFV3_8_1(w[6], w[0], params.GC_145, params.GC_79, params.mdl_MT, // OtG
      params.mdl_WT, w[16]);
  FFV5_2(w[11], w[0], params.GC_7, params.mdl_MT, params.mdl_WT, w[17]); 
  FFV3_8_2(w[11], w[0], params.GC_145, params.GC_79, params.mdl_MT, // OtG
      params.mdl_WT, w[18]);
  ixxxxx(&momenta[perm[0]][0], mME[0], hel[0], +1, w[19]); 
  oxxxxx(&momenta[perm[1]][0], mME[1], hel[1], -1, w[20]); 
  FFV1P0_3(w[19], w[20], params.GC_2, params.ZERO, params.ZERO, w[21]); 
  FFV1P0_3(w[19], w[20], params.GC_7, params.ZERO, params.ZERO, w[22]); 
  FFV2_7_3(w[19], w[20], params.GC_59, params.GC_69, params.mdl_MZ,
      params.mdl_WZ, w[23]);
  FFV1P0_3(w[19], w[20], params.GC_1, params.ZERO, params.ZERO, w[24]); 
  FFV2_4_3(w[19], w[20], params.GC_58, params.GC_69, params.mdl_MZ,
      params.mdl_WZ, w[25]);

  // Calculate all amplitudes
  
  // gg
  FFVV1_2_0(w[11], w[6], w[0], w[1], params.GC_146, params.GC_85, amp[0]); 
  FFS2_0(w[11], w[6], w[12], params.GC_95, amp[1]); 
  FFV5_0(w[11], w[6], w[13], params.GC_7, amp[2]);
  FFV5_0(w[11], w[6], w[14], params.GC_7, amp[3]); 
  FFV3_8_0(w[11], w[6], w[14], params.GC_145, params.GC_79, amp[4]); 
  FFV5_0(w[11], w[15], w[1], params.GC_7, amp[5]); 
  FFV3_8_0(w[11], w[15], w[1], params.GC_145, params.GC_79, amp[6]); 
  FFV5_0(w[11], w[16], w[1], params.GC_7, amp[7]); 
  FFV5_0(w[17], w[6], w[1], params.GC_7, amp[8]); 
  FFV3_8_0(w[17], w[6], w[1], params.GC_145, params.GC_79, amp[9]); 
  FFV5_0(w[18], w[6], w[1], params.GC_7, amp[10]); 
  
  // uux
  FFFF4_5_0(w[11], w[20], w[19], w[6], params.GC_16, params.GC_17, amp[11]); 
  FFFF1_0(w[11], w[20], w[19], w[6], -params.GC_13, amp[12]); 
  FFFF1_6_0(w[11], w[20], w[19], w[6], -params.GC_11, -params.GC_27, amp[13]); 
  FFV5_0(w[11], w[6], w[22], params.GC_7, amp[15]); 
  FFV3_8_0(w[11], w[6], w[22], params.GC_145, params.GC_79, amp[16]); 
  
  // ccx
  FFFF4_5_0(w[19], w[6], w[11], w[20], params.GC_17, params.GC_16, amp[19]); 
  FFFF3_0(w[19], w[6], w[11], w[20], params.GC_13, amp[20]); 
  FFFF3_7_0(w[19], w[6], w[11], w[20], params.GC_11, params.GC_27, amp[21]); 
  amp[23] = amp[15]; //FFV5_0(w[11], w[6], w[22], params.GC_7, amp[23]); 
  amp[24] = amp[16]; //FFV3_8_0(w[11], w[6], w[22], params.GC_145, params.GC_79, amp[24]); 
  
  // ddx
  amp[27] = amp[19]; //FFFF4_5_0(w[19], w[6], w[11], w[20], params.GC_15, params.GC_16, amp[27]); 
  FFFF1_0(w[19], w[6], w[11], w[20], -params.GC_12, amp[28]); 
  FFFF1_2_6_0(w[19], w[6], w[11], w[20], -params.GC_24, -params.GC_23,
      -params.GC_26, amp[29]);
  amp[31] = amp[23]; //FFV5_0(w[11], w[6], w[22], params.GC_7, amp[31]); 
  amp[32] = amp[16]; //FFV3_8_0(w[11], w[6], w[22], params.GC_145, params.GC_79, amp[32]); 
  
  // ssx
  amp[35] = amp[27]; //FFFF4_5_0(w[19], w[6], w[11], w[20], params.GC_15, params.GC_16, amp[35]); 
  amp[36] = amp[28]; //FFFF1_0(w[19], w[6], w[11], w[20], -params.GC_12, amp[36]); 
  FFFF1_6_0(w[19], w[6], w[11], w[20], -params.GC_10, -params.GC_26, amp[37]); 
  amp[39] = amp[31]; //FFV5_0(w[11], w[6], w[22], params.GC_7, amp[39]); 
  amp[40] = amp[32]; //FFV3_8_0(w[11], w[6], w[22], params.GC_145, params.GC_79, amp[40]); 

}
double cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_gg_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[3][2]; 
  // The color matrix
  static const double denom[3] = {3, 3, 1}; 
  static const double cf[3][3] = {{16, -2, 6}, {-2, 16, 6}, {2, 2, 6}}; 

  // Calculate color flows
  static const std::complex<double> cI(0., 1.); 
  // NP=0
  jamp[0][0] = - cI * amp[3] + amp[5];
  jamp[1][0] = + cI * amp[3] + amp[8];
  jamp[2][0] = 0.; 
  // NP=2
  jamp[0][1] = - cI * amp[2] - cI * amp[0] - cI * amp[4] + amp[6] + amp[7];
  jamp[1][1] = + cI * amp[2] + cI * amp[0] + cI * amp[4] + amp[9] + amp[10];
  jamp[2][1] = +2. * (+amp[1]); 

  // Sum and square the color flows to get the matrix element
  double matrix = 0;
  for(int i = 0; i < 3; i++ )
  {
    for(int sq1 = 0; sq1 < 2; sq1++ )
    {
      ztemp = 0.; 
      for(int j = 0; j < 3; j++ )
      {
        for(int sq2 = 0; sq2 < 2; sq2++ )
        {
          if(sq1 + sq2 == 1)
            ztemp = ztemp + cf[i][j] * jamp[j][sq2];
        }
      }
      matrix = matrix + real(ztemp * conj(jamp[i][sq1]))/denom[i]; 
    }
  }

  return matrix; 
}

double cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_uux_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[2][2]; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{9, 3}, {3, 9}}; 

  // Calculate color flows
  // NP=0
  jamp[0][0] = -1./6. * amp[15];
  jamp[1][0] = 1./2. * amp[15]; 
  // NP=2
  jamp[0][1] = +1./6. * amp[13] - amp[12] - 1./6. * amp[16];
  jamp[1][1] = -1./2. * amp[13] - amp[11] + 1./2. * amp[16]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(int i = 0; i < 2; i++ )
  {
    for(int sq1 = 0; sq1 < 2; sq1++ )
    {
      ztemp = 0.; 
      for(int j = 0; j < 2; j++ )
      {
        for(int sq2 = 0; sq2 < 2; sq2++ )
        {
          if(sq1 + sq2 == 1)
            ztemp = ztemp + cf[i][j] * jamp[j][sq2];
        }
      }
      matrix = matrix + real(ztemp * conj(jamp[i][sq1]))/denom[i]; 
    }
  }

  return matrix; 
}

double cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_ccx_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[2][2]; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{9, 3}, {3, 9}}; 

  // Calculate color flows
  // NP=0
  jamp[0][0] = - 1./6. * amp[23];
  jamp[1][0] = 1./2. * amp[23]; 
  // NP=2
  jamp[0][1] = +1./6. * amp[21] - amp[20] - 1./6. * amp[24];
  jamp[1][1] = -1./2. * amp[21] - amp[19] + 1./2. * amp[24]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(int i = 0; i < 2; i++ )
  {
    for(int sq1 = 0; sq1 < 2; sq1++ )
    {
      ztemp = 0.; 
      for(int j = 0; j < 2; j++ )
      {
        for(int sq2 = 0; sq2 < 2; sq2++ )
        {
          if(sq1 + sq2 == 1)
            ztemp = ztemp + cf[i][j] * jamp[j][sq2];
        }
      }
      matrix = matrix + real(ztemp * conj(jamp[i][sq1]))/denom[i]; 
    }
  }

  return matrix; 
}

double cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_ddx_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[2][2]; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{9, 3}, {3, 9}}; 

  // Calculate color flows
  // NP=0
  jamp[0][0] = - 1./6. * amp[31];
  jamp[1][0] = 1./2. * amp[31]; 
  // NP=2
  jamp[0][1] = +1./6. * amp[29] - amp[28] - 1./6. * amp[32];
  jamp[1][1] = -1./2. * amp[29] - amp[27] + 1./2. * amp[32]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(int i = 0; i < 2; i++ )
  {
    for(int sq1 = 0; sq1 < 2; sq1++ )
    {
      ztemp = 0.; 
      for(int j = 0; j < 2; j++ )
      {
        for(int sq2 = 0; sq2 < 2; sq2++ )
        {
          if(sq1 + sq2 == 1)
            ztemp = ztemp + cf[i][j] * jamp[j][sq2];
        }
      }
      matrix = matrix + real(ztemp * conj(jamp[i][sq1]))/denom[i]; 
    }
  }

  return matrix; 
}

double cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_ssx_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[2][2]; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{9, 3}, {3, 9}}; 

  // Calculate color flows
  // NP=0
  jamp[0][0] = - 1./6. * amp[39];
  jamp[1][0] = 1./2. * amp[39]; 
  // NP=2
  jamp[0][1] = +1./6. * amp[37] - amp[36] - 1./6. * amp[40];
  jamp[1][1] = -1./2. * amp[37] - amp[35] + 1./2. * amp[40]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(int i = 0; i < 2; i++ )
  {
    for(int sq1 = 0; sq1 < 2; sq1++ )
    {
      ztemp = 0.; 
      for(int j = 0; j < 2; j++ )
      {
        for(int sq2 = 0; sq2 < 2; sq2++ )
        {
          if(sq1 + sq2 == 1)
            ztemp = ztemp + cf[i][j] * jamp[j][sq2];
        }
      }
      matrix = matrix + real(ztemp * conj(jamp[i][sq1]))/denom[i]; 
    }
  }

  return matrix; 
}



