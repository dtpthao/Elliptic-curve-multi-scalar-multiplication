#ifndef _OPTION_BIN_H
#define _OPTION_BIN_H

#include "Duration.h"

typedef epoint* pepoint;

/* option_Bin.cpp */
void ScalarMul_Bin_L2R(big k, pepoint P, pepoint R);

/* option_Bin2.cpp */
inline void PreMul_Bin(pepoint P, pepoint Q, pepoint *plist);
// use this function
void ShamirMul_Bin(big a, pepoint P, big b, pepoint Q, pepoint R);
void TestShrMul_Bin(csprng &Rng, pepoint P, big n, std::string &msg);

/* option_Bin3.cpp */
void ShamirMul_Bin3(big k1, big k2, big k3, pepoint P1, pepoint P2, pepoint P3, pepoint R);
void test_bin3(csprng &Rng, pepoint P, big n, std::string msg);

/* option_Bin_n.cpp */
void ShamirMul_Bin_n(int n, PL *shrBin, big *k, pepoint *P, pepoint R);
void test_bin_n(int d, csprng &Rng, pepoint P, big n, std::string msg);


/*******************************************************************************/
/**_________________________The Very First Option_____________________________**/
/**                                                                           **/
void ScalarMul_Bin_R2L(big k, pepoint P, pepoint R);
void PreMul_Bin_old(pepoint P, pepoint Q, pepoint *plist);
void ShamirMul_Bin_old(big a, pepoint P, big b, pepoint Q, pepoint R);
/**_____________________End "The very First Option"___________________________**/
/*******************************************************************************/


#endif
