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
#include <iostream>
#include <iomanip>

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

  // Set the event specific parameters => here, do not change => set here once and for all
  params.setDependentParameters(); 
  params.setDependentCouplings();
  params.printIndependentParameters();
  params.printDependentParameters();

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

std::map < std::pair < int, int > , std::vector<double> >
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

  // Initialise result object
  std::map < std::pair < int, int > , std::vector<double> > result; 

  // Define permutation
  static int perm[8]; 
  for(int i = 0; i < 8; i++ )
  {
    perm[i] = i; 
  }

  for(auto &me: mapFinalStates[selectedFinalState])
  {
    std::vector<double> me_sum(MAXOPS, 0); 
    std::vector<double> me_mirror_sum(MAXOPS, 0); 
    for(int ihel = 0; ihel < 256; ihel++ )
    {
      if(me.goodHel[ihel])
      {
        double sum = 0.; 
        calculate_wavefunctions(perm, helicities[ihel]); 
        std::vector<double> meTemp = (this->* (me.meCall))();
        sum += firstNonZero(meTemp);
        addAlphaV2toV1(me_sum, meTemp, 1./me.denominator);

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
          sum += firstNonZero(meTemp);
          addAlphaV2toV1(me_mirror_sum, meTemp, 1./me.denominator);
        }

        if( !sum){
          me.goodHel[ihel] = false;
        }
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
  FFVV1_2_0(w[11], w[6], w[0], w[1], params.GC_146, params.GC_85, amp[0]); //OtG 
  FFS2_0(w[11], w[6], w[12], params.GC_95, amp[1]); 
  FFV5_0(w[11], w[6], w[13], params.GC_7, amp[2]);
  FFV5_0(w[11], w[6], w[14], params.GC_7, amp[3]); 
  FFV3_8_0(w[11], w[6], w[14], params.GC_145, params.GC_79, amp[4]); //OtG 
  FFV5_0(w[11], w[15], w[1], params.GC_7, amp[5]); 
  FFV3_8_0(w[11], w[15], w[1], params.GC_145, params.GC_79, amp[6]); //OtG 
  FFV5_0(w[11], w[16], w[1], params.GC_7, amp[7]); 
  FFV5_0(w[17], w[6], w[1], params.GC_7, amp[8]); 
  FFV3_8_0(w[17], w[6], w[1], params.GC_145, params.GC_79, amp[9]); //OtG 
  FFV5_0(w[18], w[6], w[1], params.GC_7, amp[10]); 
  
  // uux
  // O1qX FFFF4_5_0(w[11], w[20], w[19], w[6], params.GC_16, params.GC_17, amp[11]); 
  FFFF1_6_0(w[11], w[20], w[19], w[6], -params.GC_11_81, 0., amp[13]); // O81qq
  FFFF1_6_0(w[11], w[20], w[19], w[6], -params.GC_11_83, 0., amp[14]); // O83qq
  FFFF1_6_0(w[11], w[20], w[19], w[6], 0., -params.GC_27, amp[17]); // O8ut
  FFV5_0(w[11], w[6], w[22], params.GC_7, amp[15]); 
  FFV3_8_0(w[11], w[6], w[22], params.GC_145, params.GC_79, amp[16]); //OtG
  
  // ccx
  // O1qX FFFF4_5_0(w[19], w[6], w[11], w[20], params.GC_17, params.GC_16, amp[19]); 
  FFFF3_7_0(w[19], w[6], w[11], w[20], params.GC_11_81, 0., amp[21]); // O81qq
  FFFF3_7_0(w[19], w[6], w[11], w[20], params.GC_11_83, 0., amp[22]); // O83qq
  FFFF3_7_0(w[19], w[6], w[11], w[20], 0., params.GC_27, amp[25]); // O8ut
  amp[23] = amp[15]; //FFV5_0(w[11], w[6], w[22], params.GC_7, amp[23]); 
  amp[24] = amp[16]; //OtG //FFV3_8_0(w[11], w[6], w[22], params.GC_145, params.GC_79, amp[24]); 
  
  // ddx
  // O1qX amp[27] = amp[19]; //FFFF4_5_0(w[19], w[6], w[11], w[20], params.GC_15, params.GC_16, amp[27]); 
  FFFF1_2_6_0(w[19], w[6], w[11], w[20], 0., -params.GC_23, 0., amp[29]); // O81qq
  FFFF1_2_6_0(w[19], w[6], w[11], w[20], -params.GC_24, 0., 0., amp[30]); // O83qq
  FFFF1_2_6_0(w[19], w[6], w[11], w[20], 0., 0., -params.GC_26, amp[33]); // O8dt
  amp[31] = amp[23]; //FFV5_0(w[11], w[6], w[22], params.GC_7, amp[31]); 
  amp[32] = amp[16]; //OtG //FFV3_8_0(w[11], w[6], w[22], params.GC_145, params.GC_79, amp[32]); 
  
  // ssx
  // O1qX amp[35] = amp[27]; //FFFF4_5_0(w[19], w[6], w[11], w[20], params.GC_15, params.GC_16, amp[35]); 
  FFFF1_6_0(w[19], w[6], w[11], w[20], -params.GC_10_81, 0., amp[37]); // O81qq
  FFFF1_6_0(w[19], w[6], w[11], w[20], -params.GC_10_83, 0., amp[38]); // O83qq
  FFFF1_6_0(w[19], w[6], w[11], w[20], 0., -params.GC_26, amp[41]); // O8dt
  amp[39] = amp[31]; //FFV5_0(w[11], w[6], w[22], params.GC_7, amp[39]); 
  amp[40] = amp[32]; //OtG //FFV3_8_0(w[11], w[6], w[22], params.GC_145, params.GC_79, amp[40]);

}

std::vector<double> cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_gg_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
{

  static const unsigned int nOps = OPHIG;
  static std::complex<double> ztemp; 
  static std::complex<double> jamp[3][nOps+1] {}; 
  // The color matrix
  static const double denom[3] = {3, 3, 1}; 
  static const double cf[3][3] = {{16, -2, 6}, {-2, 16, 6}, {2, 2, 6}}; 

  // Calculate color flows
  static const std::complex<double> cI(0., 1.); 
  // SM
  jamp[0][SM] = - cI * amp[3] + amp[5];
  jamp[1][SM] = + cI * amp[3] + amp[8];
  jamp[2][SM] = 0.;
  // OtG
  jamp[0][OTG] = - cI * (amp[0] + amp[4]) + amp[6] + amp[7];
  jamp[1][OTG] = + cI * (amp[0] + amp[4]) + amp[9] + amp[10];
  jamp[2][OTG] = 0.;
  // OG
  jamp[0][OG] = - cI * amp[2];
  jamp[1][OG] = + cI * amp[2];
  jamp[2][OG] = 0.;
  // OphiG
  jamp[0][OPHIG] = 0.;
  jamp[1][OPHIG] = 0.;
  jamp[2][OPHIG] = +2. * (+amp[1]); 

  // Sum and square the color flows to get the matrix element
  std::vector<double> matrix(MAXOPS, 0.);
  for(int sq1 = 0; sq1 < nOps+1; sq1++ )
  {
    for(int i = 0; i < 3; i++ )
    {
      ztemp = 0;
      for(int j = 0; j < 3; j++ )
        ztemp += cf[i][j] * jamp[j][sq1];
      for(int sq2 = 0; sq2 < nOps+1; sq2++)
      {
        if(sq1 == SM && sq2 > 0)
          matrix[sq2] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq2 == SM && sq1 > 0)
          matrix[sq1] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq1 == SM && sq2 == SM)
          matrix[SM] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
      }
    }
  }

  return matrix; 
}

std::vector<double> cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_uux_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
{

  static const unsigned int nOps = O8UT;
  static std::complex<double> ztemp; 
  static std::complex<double> jamp[2][nOps+1] {}; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{9, 3}, {3, 9}}; 

  // Calculate color flows
  // SM
  jamp[0][SM] = -1./6. * amp[15];
  jamp[1][SM] = 1./2. * amp[15]; 
  // OtG
  jamp[0][OTG] = - 1./6. * amp[16];
  jamp[1][OTG] = 1./2. * amp[16]; 
  // O81qq
  jamp[0][O81QQ] = +1./6. * amp[13];
  jamp[1][O81QQ] = -1./2. * amp[13]; 
  // O83qq
  jamp[0][O83QQ] = +1./6. * amp[14];
  jamp[1][O83QQ] = -1./2. * amp[14]; 
  // O8ut
  jamp[0][O8UT] = +1./6. * amp[17];
  jamp[1][O8UT] = -1./2. * amp[17]; 
  // O1qt
  /*jamp[0][1] = 0.;
  jamp[1][1] = - amp[11]; 
  // O1qu
  jamp[0][1] = 0.;
  jamp[1][1] = - amp[11]; */

  // Sum and square the color flows to get the matrix element
  std::vector<double> matrix(MAXOPS, 0.);
  for(int sq1 = 0; sq1 < nOps+1; sq1++ )
  {
    for(int i = 0; i < 2; i++ )
    {
      ztemp = 0;
      for(int j = 0; j < 2; j++ )
        ztemp += cf[i][j] * jamp[j][sq1];
      for(int sq2 = 0; sq2 < nOps+1; sq2++)
      {
        if(sq1 == SM && sq2 > 0)
          matrix[sq2] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq2 == SM && sq1 > 0)
          matrix[sq1] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq1 == SM && sq2 == SM)
          matrix[SM] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
      }
    }
  }

  return matrix; 
}

std::vector<double> cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_ccx_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
{

  static const unsigned int nOps = O8UT;
  static std::complex<double> ztemp; 
  static std::complex<double> jamp[2][nOps+1] {}; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{9, 3}, {3, 9}}; 

  // Calculate color flows
  // SM
  jamp[0][SM] = - 1./6. * amp[23];
  jamp[1][SM] = 1./2. * amp[23]; 
  // OtG
  jamp[0][OTG] = - 1./6. * amp[24];
  jamp[1][OTG] = 1./2. * amp[24]; 
  // O81qq
  jamp[0][O81QQ] = +1./6. * amp[21];
  jamp[1][O81QQ] = -1./2. * amp[21]; 
  // O83qq
  jamp[0][O83QQ] = +1./6. * amp[22];
  jamp[1][O83QQ] = -1./2. * amp[22]; 
  // O8ut
  jamp[0][O8UT] = +1./6. * amp[25];
  jamp[1][O8UT] = -1./2. * amp[25]; 
  /*// O1qt
  jamp[0][1] = 0.;
  jamp[1][1] = - amp[19]; 
  // O1qu
  jamp[0][1] = 0.;
  jamp[1][1] = - amp[19];*/

  // Sum and square the color flows to get the matrix element
  std::vector<double> matrix(MAXOPS, 0.);
  for(int sq1 = 0; sq1 < nOps+1; sq1++ )
  {
    for(int i = 0; i < 2; i++ )
    {
      ztemp = 0;
      for(int j = 0; j < 2; j++ )
        ztemp += cf[i][j] * jamp[j][sq1];
      for(int sq2 = 0; sq2 < nOps+1; sq2++)
      {
        if(sq1 == SM && sq2 > 0)
          matrix[sq2] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq2 == SM && sq1 > 0)
          matrix[sq1] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq1 == SM && sq2 == SM)
          matrix[SM] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
      }
    }
  }

  return matrix; 
}

std::vector<double> cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_ddx_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
{

  static const unsigned int nOps = O8DT;
  static std::complex<double> ztemp; 
  static std::complex<double> jamp[2][nOps+1] {}; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{9, 3}, {3, 9}}; 

  // Calculate color flows
  // SM
  jamp[0][SM] = - 1./6. * amp[31];
  jamp[1][SM] = 1./2. * amp[31]; 
  // OTG
  jamp[0][OTG] = - 1./6. * amp[32];
  jamp[1][OTG] = 1./2. * amp[32]; 
  // O81QQ
  jamp[0][O81QQ] = 1./6. * amp[29];
  jamp[1][O81QQ] = -1./2. * amp[29]; 
  // O83QQ
  jamp[0][O83QQ] = 1./6. * amp[30];
  jamp[1][O83QQ] = -1./2. * amp[30]; 
  // O8DT
  jamp[0][O8DT] = 1./6. * amp[33];
  jamp[1][O8DT] = -1./2. * amp[33]; 
  /*// O1QD
  jamp[0][1] = 0.; 
  jamp[1][1] = - amp[27]; 
  // O1QT
  jamp[0][1] = 0.; 
  jamp[1][1] = - amp[27];*/

  // Sum and square the color flows to get the matrix element
  std::vector<double> matrix(MAXOPS, 0.);
  for(int sq1 = 0; sq1 < nOps+1; sq1++ )
  {
    for(int i = 0; i < 2; i++ )
    {
      ztemp = 0;
      for(int j = 0; j < 2; j++ )
        ztemp += cf[i][j] * jamp[j][sq1];
      for(int sq2 = 0; sq2 < nOps+1; sq2++)
      {
        if(sq1 == SM && sq2 > 0)
          matrix[sq2] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq2 == SM && sq1 > 0)
          matrix[sq1] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq1 == SM && sq2 == SM)
          matrix[SM] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
      }
    }
  }

  return matrix; 
}

std::vector<double> cpp_TopEffTh_pp_ttx_epem_NP2_QED1::matrix_1_ssx_ttx_t_wpb_wp_epve_tx_wmbx_wm_emvex() 
 {
 
   static const unsigned int nOps = O8DT;
   static std::complex<double> ztemp; 
  static std::complex<double> jamp[2][nOps+1] {}; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{9, 3}, {3, 9}}; 

  // Calculate color flows
  // SM
  jamp[0][SM] = - 1./6. * amp[39];
  jamp[1][SM] = 1./2. * amp[39]; 
  // OTG
  jamp[0][OTG] = - 1./6. * amp[40];
  jamp[1][OTG] = 1./2. * amp[40]; 
  // O81QQ
  jamp[0][O81QQ] = 1./6. * amp[37];
  jamp[1][O81QQ] = -1./2. * amp[37]; 
  // O83QQ
  jamp[0][O83QQ] = 1./6. * amp[38];
  jamp[1][O83QQ] = -1./2. * amp[38]; 
  // O8DT
  jamp[0][O8DT] = 1./6. * amp[41];
  jamp[1][O8DT] = -1./2. * amp[41]; 
  /*// O1QD
  jamp[0][1] = 0.; 
  jamp[1][1] = - amp[35]; 
  // O1QT
  jamp[0][1] = 0.; 
  jamp[1][1] = - amp[35];*/

  // Sum and square the color flows to get the matrix element
  std::vector<double> matrix(MAXOPS, 0.);
  for(int sq1 = 0; sq1 < nOps+1; sq1++ )
  {
    for(int i = 0; i < 2; i++ )
    {
      ztemp = 0;
      for(int j = 0; j < 2; j++ )
        ztemp += cf[i][j] * jamp[j][sq1];
      for(int sq2 = 0; sq2 < nOps+1; sq2++)
      {
        if(sq1 == SM && sq2 > 0)
          matrix[sq2] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq2 == SM && sq1 > 0)
          matrix[sq1] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
        if(sq1 == SM && sq2 == SM)
          matrix[SM] += real(ztemp * conj(jamp[i][sq2]))/denom[i];
      }
    }
  }

  return matrix; 
}



