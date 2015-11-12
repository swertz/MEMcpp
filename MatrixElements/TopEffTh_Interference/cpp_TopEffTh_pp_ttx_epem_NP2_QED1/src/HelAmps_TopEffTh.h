//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_TopEffTh_H
#define HelAmps_TopEffTh_H

#include <cmath> 
#include <complex> 

namespace MG5_TopEffTh 
{
double Sgn(double e, double f); 

void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void txxxxx(double p[4], double tmass, int nhel, int nst, std::complex<double>
    fi[18]);

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void FFS2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex);

void FFV1P0_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[]);

void FFV3_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[]);
void FFV3_8_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[]);

void FFV8_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex);

void FFVV2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex);

void FFV2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex);

void FFFF4_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex);
void FFFF4_5_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> & vertex);

void FFV8_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[]);

void FFV5_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[]);

void FFFF6_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex);

void FFV7_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[]);

void FFFF3_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex);
void FFFF3_7_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> & vertex);

void VVS2_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[]);

void FFV2_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[]);

void FFV3_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex);
void FFV3_8_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex);

void FFV2_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[]);
void FFV2_7_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[]);
void FFV2_4_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[]);

void FFV5_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[]);

void FFFF5_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex);

void FFFF1_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex);
void FFFF1_2_6_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> COUP3,
    std::complex<double> & vertex);
void FFFF1_6_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> & vertex);

void FFV2_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[]);

void FFVV1_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex);
void FFVV1_2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> & vertex);

void FFV5_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex);

void FFV3_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[]);
void FFV3_8_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[]);

void FFFF2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex);

void VVV2P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[]);

void FFV8_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[]);

void FFV4_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[]);

void FFFF7_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> F3[], std::complex<double> F4[], std::complex<double>
    COUP, std::complex<double> & vertex);

void VVV1P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[]);

}  // end namespace MG5_TopEffTh

#endif  // HelAmps_TopEffTh_H

