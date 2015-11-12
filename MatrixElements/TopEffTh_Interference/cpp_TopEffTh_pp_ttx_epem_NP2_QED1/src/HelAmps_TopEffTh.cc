//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "HelAmps_TopEffTh.h"
#include <complex> 
#include <cmath> 
#include <iostream> 
#include <cstdlib> 
using namespace std; 

namespace MG5_TopEffTh 
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}

void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = complex<double> (-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3])); 
    if (pp == 0.0)
    {
      sqm[0] = sqrt(abs(fmass)); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = sqrt(p[0] + pp); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (sqrt(pp3 * 0.5/pp), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/sqrt(2.0 * pp * pp3); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = sqrt(max(p[0] + p[3], 0.0)) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * sqrt(2.0 * p[0]), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = complex<double> (0.0, 0.0); 
      fi[5] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}


void txxxxx(double p[4], double tmass, int nhel, int nst, complex<double>
    tc[18])
{
  complex<double> ft[6][4], ep[4], em[4], e0[4]; 
  double pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  sqh = sqrt(0.5); 
  sqs = sqrt(0.5/3); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], sqrt(pt2 + p[3] * p[3])); 
  pt = min(pp, sqrt(pt2)); 

  ft[4][0] = complex<double> (p[0] * nst, p[3] * nst); 
  ft[5][0] = complex<double> (p[1] * nst, p[2] * nst); 

  // construct eps+
  if(nhel >= 0)
  {
    if(pp == 0)
    {
      ep[0] = complex<double> (0, 0); 
      ep[1] = complex<double> (-sqh, 0); 
      ep[2] = complex<double> (0, nst * sqh); 
      ep[3] = complex<double> (0, 0); 
    }
    else
    {
      ep[0] = complex<double> (0, 0); 
      ep[3] = complex<double> (pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = complex<double> (-sqh, 0); 
        ep[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }

  }

  // construct eps-
  if(nhel <= 0)
  {
    if(pp == 0)
    {
      em[0] = complex<double> (0, 0); 
      em[1] = complex<double> (sqh, 0); 
      em[2] = complex<double> (0, nst * sqh); 
      em[3] = complex<double> (0, 0); 
    }
    else
    {
      em[0] = complex<double> (0, 0); 
      em[3] = complex<double> (-pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = complex<double> (sqh, 0); 
        em[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }
  }

  // construct eps0
  if(fabs(nhel) <= 1)
  {
    if(pp == 0)
    {
      e0[0] = complex<double> (0, 0); 
      e0[1] = complex<double> (0, 0); 
      e0[2] = complex<double> (0, 0); 
      e0[3] = complex<double> (1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = complex<double> (pp/tmass, 0); 
      e0[3] = complex<double> (p[3] * emp, 0); 

      if(pt != 0)
      {
        e0[1] = complex<double> (p[1] * emp, 0); 
        e0[2] = complex<double> (p[2] * emp, 0); 
      }
      else
      {
        e0[1] = complex<double> (0, 0); 
        e0[2] = complex<double> (0, 0); 
      }
    }
  }

  if(nhel == 2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = ep[i] * ep[j]; 
    }
  }
  else if(nhel == -2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = em[i] * em[j]; 
    }
  }
  else if(tmass == 0)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = 0; 
    }
  }
  else if(tmass != 0)
  {
    if(nhel == 1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]); 
      }
    }
    else if(nhel == 0)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqs * (ep[i] * em[j] + em[i] * ep[j]
         + 2.0 * e0[i] * e0[j]); 
      }
    }
    else if(nhel == -1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]); 
      }
    }
    else
    {
      std::cerr <<  "Invalid helicity in txxxxx.\n"; 
      std::exit(1); 
    }
  }

  tc[0] = ft[4][0]; 
  tc[1] = ft[5][0]; 

  for(j = 0; j < 4; j++ )
  {
    for(i = 0; i < 4; i++ )
      tc[j * 4 + i + 2] = ft[j][i]; 
  }
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = sqrt(0.5); 
  hel = double(nhel); 
  nsvahl = nsv * abs(hel); 
  pt2 = (p[1] * p[1]) + (p[2] * p[2]); 
  pp = min(p[0], sqrt(pt2 + (p[3] * p[3]))); 
  pt = min(pp, sqrt(pt2)); 
  vc[0] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - abs(hel); 
    if(pp == 0.0)
    {
      vc[2] = complex<double> (0.0, 0.0); 
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsvahl * sqh); 
      vc[5] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = complex<double> (-hel * sqh, 0.0); 
        vc[4] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = sqrt((p[1] * p[1]) + (p[2] * p[2])); 
    vc[2] = complex<double> (0.0, 0.0); 
    vc[5] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3]))); 
    if (pp == 0.000)
    {
      sqm[0] = sqrt(abs(fmass)); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[abs(ip)]; 
      fo[3] = ip * nsf * sqm[abs(ip)]; 
      fo[4] = im * nsf * sqm[abs(im)]; 
      fo[5] = ip * sqm[abs(im)]; 
    }
    else
    {
      pp = min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3]))); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = sqrt(p[0] + pp); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (sqrt(pp3 * 0.5/pp), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/sqrt(2.0 * pp * pp3); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = sqrt(max(p[0] + p[3], 0.00)) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * sqrt(2.0 * p[0]); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = complex<double> (0.00, 0.00); 
      fo[5] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[2] = complex<double> (1.00, 0.00); 
  sc[0] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}

void FFS2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP7; 
  std::complex<double> TMP13; 
  TMP7 = (F1[4] * F2[4] + F1[5] * F2[5]); 
  TMP13 = (F1[2] * F2[2] + F1[3] * F2[3]); 
  vertex = COUP * - S3[2] * (+cI * (TMP13 + TMP7)); 
}


void FFV1P0_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5]
      * F2[3]);
  V3[3] = denom * - cI * (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3]
      * F2[4]);
  V3[4] = denom * - cI * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] *
      F2[4] + F1[4] * F2[3]));
  V3[5] = denom * - cI * (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5]
      * F2[3]);
}


void FFV3_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P3[4]; 
  std::complex<double> denom; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * - 2. * cI * M1 * (F2[3] * (P3[0] * (V3[3] + cI * (V3[4])) +
      (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5]))
      + P3[3] * (V3[3] + cI * (V3[4]))))) + F2[2] * (V3[5] * P3[0] - cI *
      (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]));
  F1[3] = denom * - 2. * cI * M1 * (F2[2] * (P3[0] * (V3[3] - cI * (V3[4])) +
      (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + F2[3] * (V3[2] * P3[3] - cI * (V3[3]
      * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]));
  F1[4] = denom * 2. * cI * (F2[2] * (P1[1] * (P3[0] * (+cI * (V3[4]) - V3[3])
      + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[2]) + cI * (V3[5])) +
      P3[3] * (V3[3] - cI * (V3[4]))))) + (P1[2] * (P3[0] * - 1. * (V3[4] + cI
      * (V3[3])) + (P3[1] * (-cI * (V3[5]) + cI * (V3[2])) + (P3[2] * (V3[2] -
      V3[5]) + P3[3] * (V3[4] + cI * (V3[3]))))) + (P1[0] * (V3[5] * P3[0] - cI
      * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]) + P1[3] *
      (V3[2] * P3[3] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] *
      P3[0])))) + F2[3] * (P1[0] * (P3[0] * (V3[3] + cI * (V3[4])) + (P3[1] * -
      1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) + P3[3] *
      (V3[3] + cI * (V3[4]))))) + (P1[3] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) -
      P3[3] * (V3[3] + cI * (V3[4]))))) + (P1[1] * (V3[5] * P3[0] - cI * (V3[4]
      * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]) + P1[2] * (V3[4] * P3[1]
      - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[3] * P3[2])))));
  F1[5] = denom * 2. * cI * (F2[2] * (P1[0] * (P3[0] * (V3[3] - cI * (V3[4])) +
      (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + (P1[3] * (P3[0] * (V3[3] - cI *
      (V3[4])) + (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (+cI * (V3[4]) - V3[3])))) + (P1[1] * (V3[2] * P3[3] -
      cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]) + P1[2] *
      (V3[4] * P3[1] - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[3] *
      P3[2])))) + F2[3] * (P1[1] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P3[1] * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) - P3[3] *
      (V3[3] + cI * (V3[4]))))) + (P1[2] * (P3[0] * (+cI * (V3[3]) - V3[4]) +
      (P3[1] * - 1. * (+cI * (V3[2] + V3[5])) + (P3[2] * (V3[2] + V3[5]) +
      P3[3] * (+cI * (V3[3]) - V3[4])))) + (P1[0] * (V3[2] * P3[3] - cI *
      (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]) + P1[3] * (V3[2]
      * P3[3] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] *
      P3[0])))));
}

void FFV3_8_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P3[4]; 
  std::complex<double> denom; 
  int i; 
  std::complex<double> Ftmp[6]; 
  FFV3_1(F2, V3, COUP1, M1, W1, F1); 
  FFV8_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFV8_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP1; 
  std::complex<double> TMP0; 
  double P3[4]; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP1 = (F1[4] * (F2[4] * (P3[0] * (V3[2] + V3[5]) + (P3[1] * (+cI * (V3[4]) -
      V3[3]) + (P3[2] * - 1. * (V3[4] + cI * (V3[3])) - P3[3] * (V3[2] +
      V3[5])))) + F2[5] * (P3[0] * (V3[3] + cI * (V3[4])) + (P3[1] * (V3[5] -
      V3[2]) + (P3[2] * (-cI * (V3[2]) + cI * (V3[5])) - P3[3] * (V3[3] + cI *
      (V3[4])))))) + F1[5] * (F2[4] * (P3[0] * (V3[3] - cI * (V3[4])) + (P3[1]
      * - 1. * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) + P3[3] *
      (V3[3] - cI * (V3[4]))))) + F2[5] * (P3[0] * (V3[2] - V3[5]) + (P3[1] * -
      1. * (V3[3] + cI * (V3[4])) + (P3[2] * (+cI * (V3[3]) - V3[4]) + P3[3] *
      (V3[2] - V3[5]))))));
  TMP0 = (F1[4] * (F2[4] * (P3[0] * (V3[2] - V3[5]) + (P3[1] * - 1. * (V3[3] +
      cI * (V3[4])) + (P3[2] * (+cI * (V3[3]) - V3[4]) + P3[3] * (V3[2] -
      V3[5])))) + F2[5] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) + (P3[1] *
      (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) + P3[3] *
      (V3[3] + cI * (V3[4])))))) + F1[5] * (F2[4] * (P3[0] * (+cI * (V3[4]) -
      V3[3]) + (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] +
      V3[5])) + P3[3] * (+cI * (V3[4]) - V3[3])))) + F2[5] * (P3[0] * (V3[2] +
      V3[5]) + (P3[1] * (+cI * (V3[4]) - V3[3]) + (P3[2] * - 1. * (V3[4] + cI *
      (V3[3])) - P3[3] * (V3[2] + V3[5]))))));
  vertex = COUP * (-cI * (TMP1) + cI * (TMP0)); 
}


void FFVV2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP19; 
  std::complex<double> TMP18; 
  TMP19 = (F1[4] * (F2[4] * (V3[2] * (V4[2] + V4[5]) + (V3[3] * (+cI * (V4[4])
      - V4[3]) + (V3[4] * - 1. * (V4[4] + cI * (V4[3])) - V3[5] * (V4[2] +
      V4[5])))) + F2[5] * (V3[2] * (V4[3] + cI * (V4[4])) + (V3[3] * (V4[5] -
      V4[2]) + (V3[4] * (-cI * (V4[2]) + cI * (V4[5])) - V3[5] * (V4[3] + cI *
      (V4[4])))))) + F1[5] * (F2[4] * (V3[2] * (V4[3] - cI * (V4[4])) + (V3[3]
      * - 1. * (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2] + V4[5])) + V3[5] *
      (V4[3] - cI * (V4[4]))))) + F2[5] * (V3[2] * (V4[2] - V4[5]) + (V3[3] * -
      1. * (V4[3] + cI * (V4[4])) + (V3[4] * (+cI * (V4[3]) - V4[4]) + V3[5] *
      (V4[2] - V4[5]))))));
  TMP18 = (F1[4] * (F2[4] * (V3[2] * (V4[2] - V4[5]) + (V3[3] * - 1. * (V4[3] +
      cI * (V4[4])) + (V3[4] * (+cI * (V4[3]) - V4[4]) + V3[5] * (V4[2] -
      V4[5])))) + F2[5] * (V3[2] * - 1. * (V4[3] + cI * (V4[4])) + (V3[3] *
      (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[5]) + cI * (V4[2])) + V3[5] *
      (V4[3] + cI * (V4[4])))))) + F1[5] * (F2[4] * (V3[2] * (+cI * (V4[4]) -
      V4[3]) + (V3[3] * (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] +
      V4[5])) + V3[5] * (+cI * (V4[4]) - V4[3])))) + F2[5] * (V3[2] * (V4[2] +
      V4[5]) + (V3[3] * (+cI * (V4[4]) - V4[3]) + (V3[4] * - 1. * (V4[4] + cI *
      (V4[3])) - V3[5] * (V4[2] + V4[5]))))));
  vertex = COUP * (-cI * (TMP18) + cI * (TMP19)); 
}


void FFV2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP5; 
  TMP5 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - cI * TMP5; 
}


void FFFF4_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP7; 
  std::complex<double> TMP6; 
  TMP7 = (F1[4] * F2[4] + F1[5] * F2[5]); 
  TMP6 = (F4[2] * F3[2] + F4[3] * F3[3]); 
  vertex = COUP * - cI * TMP6 * TMP7; 
}

void FFFF4_5_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> tmp; 
  FFFF4_0(F1, F2, F3, F4, COUP1, vertex); 
  FFFF5_0(F1, F2, F3, F4, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFV8_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P3[4]; 
  std::complex<double> denom; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * 2. * cI * (F2[4] * (P1[1] * (P3[0] * (+cI * (V3[4]) - V3[3])
      + (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + (P1[2] * (P3[0] * - 1. * (V3[4] + cI
      * (V3[3])) + (P3[1] * (+cI * (V3[2] + V3[5])) + (P3[2] * (V3[2] + V3[5])
      - P3[3] * (V3[4] + cI * (V3[3]))))) + (P1[0] * (V3[2] * P3[3] - cI *
      (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]) + P1[3] * (V3[2]
      * P3[3] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]))))
      + F2[5] * (P1[0] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) + (P3[1] *
      (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) + P3[3] *
      (V3[3] + cI * (V3[4]))))) + (P1[3] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P1[1] * (V3[5] * P3[0] -
      cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]) + P1[2] *
      (V3[3] * P3[2] - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[4] *
      P3[1])))));
  F1[3] = denom * 2. * cI * (F2[4] * (P1[0] * (P3[0] * (+cI * (V3[4]) - V3[3])
      + (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + (P1[3] * (P3[0] * (V3[3] - cI *
      (V3[4])) + (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] +
      V3[5])) + P3[3] * (V3[3] - cI * (V3[4]))))) + (P1[1] * (V3[2] * P3[3] -
      cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]) + P1[2] *
      (V3[3] * P3[2] - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[4] *
      P3[1])))) + F2[5] * (P1[1] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (V3[3] + cI * (V3[4]))))) + (P1[2] * (P3[0] * (+cI * (V3[3]) -
      V3[4]) + (P3[1] * (-cI * (V3[2]) + cI * (V3[5])) + (P3[2] * (V3[2] -
      V3[5]) + P3[3] * (V3[4] - cI * (V3[3]))))) + (P1[0] * (V3[5] * P3[0] - cI
      * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]) + P1[3] *
      (V3[2] * P3[3] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] *
      P3[0])))));
  F1[4] = denom * - 2. * cI * M1 * (F2[5] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + F2[4] * (V3[2] * P3[3] -
      cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]));
  F1[5] = denom * - 2. * cI * M1 * (F2[4] * (P3[0] * (+cI * (V3[4]) - V3[3]) +
      (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + F2[5] * (V3[5] * P3[0] - cI * (V3[3]
      * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]));
}


void FFV5_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3] - cI
      * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] - V3[2]))))
      + (F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * - 1. * (V3[2] +
      V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] + cI *
      (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])))));
  F1[3] = denom * - cI * (F2[2] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[2] + V3[5]) + (P1[1] *
      - 1. * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3]
      * (V3[2] + V3[5])))) + M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] *
      (V3[5] - V3[2]))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (V3[5] - V3[2]) + F2[3] *
      (V3[3] + cI * (V3[4])))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] +
      V3[5]))));
}


void FFFF6_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP2; 
  TMP2 = (F1[4] * F3[5] * (F2[3] * F4[2] - F2[2] * F4[3]) + F1[5] * F3[4] *
      (F2[2] * F4[3] - F2[3] * F4[2]));
  vertex = COUP * - 2. * cI * TMP2; 
}


void FFV7_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> denom; 
  double P3[4]; 
  double OM3; 
  std::complex<double> TMP4; 
  std::complex<double> TMP3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP4 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP3 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * - 4. * cI * (OM3 * - 1./4. * P3[0] * (TMP3 + 4. * (TMP4)) +
      (+1./4. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] *
      F2[3]));
  V3[3] = denom * - 4. * cI * (OM3 * - 1./4. * P3[1] * (TMP3 + 4. * (TMP4)) +
      (-1./4. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] *
      F2[2]));
  V3[4] = denom * 4. * cI * (OM3 * 1./4. * P3[2] * (TMP3 + 4. * (TMP4)) +
      (+1./4. * cI * (F1[2] * F2[5]) - 1./4. * cI * (F1[3] * F2[4]) - cI *
      (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
  V3[5] = denom * 4. * cI * (OM3 * 1./4. * P3[3] * (TMP3 + 4. * (TMP4)) +
      (+1./4. * (F1[2] * F2[4]) - 1./4. * (F1[3] * F2[5]) - F1[4] * F2[2] +
      F1[5] * F2[3]));
}


void FFFF3_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP20; 
  TMP20 = (F1[2] * F3[3] * (F2[4] * F4[5] - F2[5] * F4[4]) + F1[3] * F3[2] *
      (F2[5] * F4[4] - F2[4] * F4[5]));
  vertex = COUP * - 2. * cI * TMP20; 
}

void FFFF3_7_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> tmp; 
  FFFF3_0(F1, F2, F3, F4, COUP1, vertex); 
  FFFF7_0(F1, F2, F3, F4, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void VVS2_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> denom; 
  double P1[4]; 
  double P2[4]; 
  std::complex<double> TMP17; 
  double P3[4]; 
  std::complex<double> TMP16; 
  std::complex<double> TMP15; 
  std::complex<double> TMP14; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  S3[0] = +V1[0] + V2[0]; 
  S3[1] = +V1[1] + V2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP15 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP14 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP17 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP16 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * (-cI * (TMP16 * TMP17) + cI * (TMP14 * TMP15)); 
}


void FFV2_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1.
      * (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + F1[3] * (P2[0] * (V3[2] - V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] - V3[5])))));
  F2[4] = denom * - cI * M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]));
  F2[5] = denom * cI * M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] * (V3[2] -
      V3[5]));
}


void FFV3_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP26; 
  std::complex<double> TMP25; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP25 = (F1[2] * (F2[2] * (P3[0] * (V3[2] + V3[5]) + (P3[1] * - 1. * (V3[3] +
      cI * (V3[4])) + (P3[2] * (+cI * (V3[3]) - V3[4]) - P3[3] * (V3[2] +
      V3[5])))) + F2[3] * (P3[0] * (V3[3] + cI * (V3[4])) + (P3[1] * - 1. *
      (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) + P3[3] *
      (V3[3] + cI * (V3[4])))))) + F1[3] * (F2[2] * (P3[0] * (V3[3] - cI *
      (V3[4])) + (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (+cI * (V3[4]) - V3[3])))) + F2[3] * (P3[0] * (V3[2] -
      V3[5]) + (P3[1] * (+cI * (V3[4]) - V3[3]) + (P3[2] * - 1. * (V3[4] + cI *
      (V3[3])) + P3[3] * (V3[2] - V3[5]))))));
  TMP26 = (F1[2] * (F2[2] * (P3[0] * (V3[2] - V3[5]) + (P3[1] * (+cI * (V3[4])
      - V3[3]) + (P3[2] * - 1. * (V3[4] + cI * (V3[3])) + P3[3] * (V3[2] -
      V3[5])))) + F2[3] * (P3[0] * - 1. * (V3[3] + cI * (V3[4])) + (P3[1] *
      (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) - P3[3] * (V3[3] + cI
      * (V3[4])))))) + F1[3] * (F2[2] * (P3[0] * (+cI * (V3[4]) - V3[3]) +
      (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[2]) + cI * (V3[5])) +
      P3[3] * (V3[3] - cI * (V3[4]))))) + F2[3] * (P3[0] * (V3[2] + V3[5]) +
      (P3[1] * - 1. * (V3[3] + cI * (V3[4])) + (P3[2] * (+cI * (V3[3]) - V3[4])
      - P3[3] * (V3[2] + V3[5]))))));
  vertex = COUP * (-cI * (TMP26) + cI * (TMP25)); 
}

void FFV3_8_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> tmp; 
  FFV3_0(F1, F2, V3, COUP1, vertex); 
  FFV8_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFV2_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> denom; 
  double P3[4]; 
  double OM3; 
  std::complex<double> TMP3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP3 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP3); 
  V3[3] = denom * - cI * (-F1[2] * F2[5] - F1[3] * F2[4] - P3[1] * OM3 * TMP3); 
  V3[4] = denom * - cI * (-cI * (F1[2] * F2[5]) + cI * (F1[3] * F2[4]) - P3[2]
      * OM3 * TMP3);
  V3[5] = denom * - cI * (F1[3] * F2[5] - F1[2] * F2[4] - P3[3] * OM3 * TMP3); 
}

void FFV2_7_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> denom; 
  double P3[4]; 
  double OM3; 
  int i; 
  std::complex<double> Vtmp[6]; 
  FFV2_3(F1, F2, COUP1, M3, W3, V3); 
  FFV7_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}
void FFV2_4_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> denom; 
  double P3[4]; 
  double OM3; 
  int i; 
  std::complex<double> Vtmp[6]; 
  FFV2_3(F1, F2, COUP1, M3, W3, V3); 
  FFV4_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void FFV5_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] * (+cI *
      (V3[4]) - V3[3]))));
  F2[3] = denom * - cI * (F1[2] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] * (V3[2] +
      V3[5]))));
  F2[4] = denom * - cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] +
      cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * - 1. *
      (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) + P2[3] * (V3[3] - cI
      * (V3[4]))))) + M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) - P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] *
      (V3[2] - V3[5]))));
}


void FFFF5_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP12; 
  std::complex<double> TMP13; 
  TMP13 = (F1[2] * F2[2] + F1[3] * F2[3]); 
  TMP12 = (F4[4] * F3[4] + F4[5] * F3[5]); 
  vertex = COUP * - cI * TMP12 * TMP13; 
}


void FFFF1_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP9; 
  TMP9 = (F1[2] * F3[3] * (F2[5] * F4[4] - F2[4] * F4[5]) + F1[3] * F3[2] *
      (F2[4] * F4[5] - F2[5] * F4[4]));
  vertex = COUP * - 2. * cI * TMP9; 
}

void FFFF1_2_6_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> COUP3,
    std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> tmp; 
  FFFF1_0(F1, F2, F3, F4, COUP1, vertex); 
  FFFF2_0(F1, F2, F3, F4, COUP2, tmp); 
  vertex = vertex + tmp; 
  FFFF6_0(F1, F2, F3, F4, COUP3, tmp); 
  vertex = vertex + tmp; 
}
void FFFF1_6_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> tmp; 
  std::complex<double> COUP3; 
  FFFF1_0(F1, F2, F3, F4, COUP1, vertex); 
  FFFF6_0(F1, F2, F3, F4, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFV2_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])));
  F1[3] = denom * - cI * M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] * (V3[5]
      - V3[2]));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * - cI * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * -
      1. * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3]
      - cI * (V3[4]))))) + F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))));
}


void FFVV1_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP28; 
  std::complex<double> TMP27; 
  TMP27 = (F1[2] * (F2[2] * (V3[2] * (V4[2] + V4[5]) + (V3[3] * - 1. * (V4[3] +
      cI * (V4[4])) + (V3[4] * (+cI * (V4[3]) - V4[4]) - V3[5] * (V4[2] +
      V4[5])))) + F2[3] * (V3[2] * (V4[3] + cI * (V4[4])) + (V3[3] * - 1. *
      (V4[2] + V4[5]) + (V3[4] * - 1. * (+cI * (V4[2] + V4[5])) + V3[5] *
      (V4[3] + cI * (V4[4])))))) + F1[3] * (F2[2] * (V3[2] * (V4[3] - cI *
      (V4[4])) + (V3[3] * (V4[5] - V4[2]) + (V3[4] * (-cI * (V4[5]) + cI *
      (V4[2])) + V3[5] * (+cI * (V4[4]) - V4[3])))) + F2[3] * (V3[2] * (V4[2] -
      V4[5]) + (V3[3] * (+cI * (V4[4]) - V4[3]) + (V3[4] * - 1. * (V4[4] + cI *
      (V4[3])) + V3[5] * (V4[2] - V4[5]))))));
  TMP28 = (F1[2] * (F2[2] * (V3[2] * (V4[2] - V4[5]) + (V3[3] * (+cI * (V4[4])
      - V4[3]) + (V3[4] * - 1. * (V4[4] + cI * (V4[3])) + V3[5] * (V4[2] -
      V4[5])))) + F2[3] * (V3[2] * - 1. * (V4[3] + cI * (V4[4])) + (V3[3] *
      (V4[2] + V4[5]) + (V3[4] * (+cI * (V4[2] + V4[5])) - V3[5] * (V4[3] + cI
      * (V4[4])))))) + F1[3] * (F2[2] * (V3[2] * (+cI * (V4[4]) - V4[3]) +
      (V3[3] * (V4[2] - V4[5]) + (V3[4] * (-cI * (V4[2]) + cI * (V4[5])) +
      V3[5] * (V4[3] - cI * (V4[4]))))) + F2[3] * (V3[2] * (V4[2] + V4[5]) +
      (V3[3] * - 1. * (V4[3] + cI * (V4[4])) + (V3[4] * (+cI * (V4[3]) - V4[4])
      - V3[5] * (V4[2] + V4[5]))))));
  vertex = COUP * (-cI * (TMP27) + cI * (TMP28)); 
}

void FFVV1_2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> tmp; 
  FFVV1_0(F1, F2, V3, V4, COUP1, vertex); 
  FFVV2_0(F1, F2, V3, V4, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFV5_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP11; 
  std::complex<double> TMP10; 
  TMP11 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP10 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - 1. * (+cI * (TMP10 + TMP11)); 
}


void FFV3_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> denom; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * - 2. * cI * M2 * (F1[3] * (P3[0] * (V3[3] - cI * (V3[4])) +
      (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + F1[2] * (V3[5] * P3[0] - cI * (V3[4]
      * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]));
  F2[3] = denom * - 2. * cI * M2 * (F1[2] * (P3[0] * (V3[3] + cI * (V3[4])) +
      (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5]))
      + P3[3] * (V3[3] + cI * (V3[4]))))) + F1[3] * (V3[2] * P3[3] - cI *
      (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]));
  F2[4] = denom * - 2. * cI * (F1[2] * (P2[1] * (P3[0] * (V3[3] + cI * (V3[4]))
      + (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] +
      V3[5])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[2] * (P3[0] * (V3[4] -
      cI * (V3[3])) + (P3[1] * (+cI * (V3[2] + V3[5])) + (P3[2] * - 1. * (V3[2]
      + V3[5]) + P3[3] * (V3[4] - cI * (V3[3]))))) + (P2[0] * (V3[5] * P3[0] -
      cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]) + P2[3] *
      (V3[5] * P3[0] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] *
      P3[3])))) + F1[3] * (P2[0] * (P3[0] * (V3[3] - cI * (V3[4])) + (P3[1] *
      (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) + P3[3] * (+cI
      * (V3[4]) - V3[3])))) + (P2[3] * (P3[0] * (V3[3] - cI * (V3[4])) + (P3[1]
      * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) + P3[3] *
      (+cI * (V3[4]) - V3[3])))) + (P2[1] * (V3[2] * P3[3] - cI * (V3[3] *
      P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]) + P2[2] * (V3[4] * P3[1] -
      cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[3] * P3[2])))));
  F2[5] = denom * - 2. * cI * (F1[2] * (P2[0] * (P3[0] * (V3[3] + cI * (V3[4]))
      + (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] +
      V3[5])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[3] * (P3[0] * - 1. *
      (V3[3] + cI * (V3[4])) + (P3[1] * (V3[2] + V3[5]) + (P3[2] * (+cI *
      (V3[2] + V3[5])) - P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[1] * (V3[5] *
      P3[0] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] * P3[3]) +
      P2[2] * (V3[4] * P3[1] - cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) -
      V3[3] * P3[2])))) + F1[3] * (P2[1] * (P3[0] * (V3[3] - cI * (V3[4])) +
      (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + (P2[2] * (P3[0] * (V3[4] + cI *
      (V3[3])) + (P3[1] * (-cI * (V3[2]) + cI * (V3[5])) + (P3[2] * (V3[5] -
      V3[2]) - P3[3] * (V3[4] + cI * (V3[3]))))) + (P2[0] * (V3[2] * P3[3] - cI
      * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[5] * P3[0]) + P2[3] *
      (V3[5] * P3[0] - cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[2] *
      P3[3])))));
}

void FFV3_8_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> Ftmp[6]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> denom; 
  int i; 
  FFV3_2(F1, V3, COUP1, M2, W2, F2); 
  FFV8_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFFF2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP9; 
  TMP9 = (F1[2] * F3[3] * (F2[5] * F4[4] - F2[4] * F4[5]) + F1[3] * F3[2] *
      (F2[4] * F4[5] - F2[5] * F4[4]));
  vertex = COUP * - 2. * cI * TMP9; 
}


void VVV2P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP22; 
  double P2[4]; 
  std::complex<double> TMP23; 
  std::complex<double> TMP17; 
  double P3[4]; 
  std::complex<double> TMP30; 
  std::complex<double> TMP21; 
  std::complex<double> denom; 
  std::complex<double> TMP14; 
  std::complex<double> TMP29; 
  std::complex<double> TMP24; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP24 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP21 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP22 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP29 = (P3[0] * P2[0] - P3[1] * P2[1] - P3[2] * P2[2] - P3[3] * P2[3]); 
  TMP14 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP17 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP30 = (P3[0] * P1[0] - P3[1] * P1[1] - P3[2] * P1[2] - P3[3] * P1[3]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (P2[0] * (-cI * (TMP24 * TMP30) + cI * (TMP21 * TMP23)) +
      (P3[0] * (-cI * (TMP14 * TMP22) + cI * (TMP17 * TMP24)) + (TMP29 * (-cI *
      (V2[2] * TMP21) + cI * (V3[2] * TMP14)) + (-cI * (V3[2] * TMP17 * TMP23)
      + cI * (V2[2] * TMP22 * TMP30)))));
  V1[3] = denom * (P2[1] * (-cI * (TMP24 * TMP30) + cI * (TMP21 * TMP23)) +
      (P3[1] * (-cI * (TMP14 * TMP22) + cI * (TMP17 * TMP24)) + (TMP29 * (-cI *
      (V2[3] * TMP21) + cI * (V3[3] * TMP14)) + (-cI * (V3[3] * TMP17 * TMP23)
      + cI * (V2[3] * TMP22 * TMP30)))));
  V1[4] = denom * (P2[2] * (-cI * (TMP24 * TMP30) + cI * (TMP21 * TMP23)) +
      (P3[2] * (-cI * (TMP14 * TMP22) + cI * (TMP17 * TMP24)) + (TMP29 * (-cI *
      (V2[4] * TMP21) + cI * (V3[4] * TMP14)) + (-cI * (V3[4] * TMP17 * TMP23)
      + cI * (V2[4] * TMP22 * TMP30)))));
  V1[5] = denom * (P2[3] * (-cI * (TMP24 * TMP30) + cI * (TMP21 * TMP23)) +
      (P3[3] * (-cI * (TMP14 * TMP22) + cI * (TMP17 * TMP24)) + (TMP29 * (-cI *
      (V2[5] * TMP21) + cI * (V3[5] * TMP14)) + (-cI * (V3[5] * TMP17 * TMP23)
      + cI * (V2[5] * TMP22 * TMP30)))));
}


void FFV8_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> denom; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * - 2. * cI * (F1[4] * (P2[1] * (P3[0] * (V3[3] + cI * (V3[4]))
      + (P3[1] * (V3[5] - V3[2]) + (P3[2] * (-cI * (V3[2]) + cI * (V3[5])) -
      P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[2] * (P3[0] * (V3[4] - cI *
      (V3[3])) + (P3[1] * (-cI * (V3[5]) + cI * (V3[2])) + (P3[2] * (V3[5] -
      V3[2]) + P3[3] * (+cI * (V3[3]) - V3[4])))) + (P2[0] * (V3[2] * P3[3] -
      cI * (V3[4] * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]) + P2[3] *
      (V3[5] * P3[0] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] *
      P3[3])))) + F1[5] * (P2[0] * (P3[0] * (+cI * (V3[4]) - V3[3]) + (P3[1] *
      (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) + P3[3] * (+cI
      * (V3[4]) - V3[3])))) + (P2[3] * (P3[0] * (V3[3] - cI * (V3[4])) + (P3[1]
      * - 1. * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] + V3[5])) + P3[3] *
      (V3[3] - cI * (V3[4]))))) + (P2[1] * (V3[2] * P3[3] - cI * (V3[4] *
      P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]) + P2[2] * (V3[3] * P3[2] -
      cI * (V3[2] * P3[3]) + cI * (V3[5] * P3[0]) - V3[4] * P3[1])))));
  F2[3] = denom * - 2. * cI * (F1[4] * (P2[0] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[3] * (P3[0] * - 1. *
      (V3[3] + cI * (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI *
      (V3[5]) + cI * (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + (P2[1] *
      (V3[5] * P3[0] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] *
      P3[3]) + P2[2] * (V3[3] * P3[2] - cI * (V3[2] * P3[3]) + cI * (V3[5] *
      P3[0]) - V3[4] * P3[1])))) + F1[5] * (P2[1] * (P3[0] * (V3[3] - cI *
      (V3[4])) + (P3[1] * - 1. * (V3[2] + V3[5]) + (P3[2] * (+cI * (V3[2] +
      V3[5])) + P3[3] * (V3[3] - cI * (V3[4]))))) + (P2[2] * (P3[0] * (V3[4] +
      cI * (V3[3])) + (P3[1] * - 1. * (+cI * (V3[2] + V3[5])) + (P3[2] * - 1. *
      (V3[2] + V3[5]) + P3[3] * (V3[4] + cI * (V3[3]))))) + (P2[0] * (V3[5] *
      P3[0] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]) +
      P2[3] * (V3[5] * P3[0] - cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) -
      V3[2] * P3[3])))));
  F2[4] = denom * - 2. * cI * M2 * (F1[5] * (P3[0] * (+cI * (V3[4]) - V3[3]) +
      (P3[1] * (V3[2] + V3[5]) + (P3[2] * - 1. * (+cI * (V3[2] + V3[5])) +
      P3[3] * (+cI * (V3[4]) - V3[3])))) + F1[4] * (V3[2] * P3[3] - cI * (V3[4]
      * P3[1]) + cI * (V3[3] * P3[2]) - V3[5] * P3[0]));
  F2[5] = denom * - 2. * cI * M2 * (F1[4] * (P3[0] * - 1. * (V3[3] + cI *
      (V3[4])) + (P3[1] * (V3[2] - V3[5]) + (P3[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P3[3] * (V3[3] + cI * (V3[4]))))) + F1[5] * (V3[5] * P3[0] -
      cI * (V3[3] * P3[2]) + cI * (V3[4] * P3[1]) - V3[2] * P3[3]));
}


void FFV4_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> denom; 
  double P3[4]; 
  double OM3; 
  std::complex<double> TMP4; 
  std::complex<double> TMP3; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP4 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP3 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * 2. * cI * (OM3 * 1./2. * P3[0] * (TMP3 - 2. * (TMP4)) +
      (-1./2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] *
      F2[3]));
  V3[3] = denom * 2. * cI * (OM3 * 1./2. * P3[1] * (TMP3 - 2. * (TMP4)) +
      (+1./2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] *
      F2[2]));
  V3[4] = denom * - 2. * cI * (OM3 * 1./2. * P3[2] * (+2. * (TMP4) - TMP3) +
      (-1./2. * cI * (F1[2] * F2[5]) + 1./2. * cI * (F1[3] * F2[4]) - cI *
      (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
  V3[5] = denom * - 2. * cI * (OM3 * 1./2. * P3[3] * (+2. * (TMP4) - TMP3) +
      (-1./2. * (F1[2] * F2[4]) + 1./2. * (F1[3] * F2[5]) - F1[4] * F2[2] +
      F1[5] * F2[3]));
}


void FFFF7_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP8; 
  TMP8 = (F1[4] * F3[5] * (F2[2] * F4[3] - F2[3] * F4[2]) + F1[5] * F3[4] *
      (F2[3] * F4[2] - F2[2] * F4[3]));
  vertex = COUP * - 2. * cI * TMP8; 
}


void VVV1P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP22; 
  double P2[4]; 
  std::complex<double> TMP23; 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> denom; 
  std::complex<double> TMP14; 
  std::complex<double> TMP24; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP24 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP14 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP21 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP22 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP23 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP24 * (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] * (-cI *
      (TMP21) + cI * (TMP22)) + V3[2] * (-cI * (TMP23) + cI * (TMP14))));
  V1[3] = denom * (TMP24 * (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] * (-cI *
      (TMP21) + cI * (TMP22)) + V3[3] * (-cI * (TMP23) + cI * (TMP14))));
  V1[4] = denom * (TMP24 * (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] * (-cI *
      (TMP21) + cI * (TMP22)) + V3[4] * (-cI * (TMP23) + cI * (TMP14))));
  V1[5] = denom * (TMP24 * (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] * (-cI *
      (TMP21) + cI * (TMP22)) + V3[5] * (-cI * (TMP23) + cI * (TMP14))));
}


}  // end namespace $(namespace)s_TopEffTh

