#ifndef _OPTION_BIN_H
#define _OPTION_BIN_H

extern "C" {
#include "miracl.h"
}
#include <Windows.h>
#include "ShamirMul.h"

typedef epoint* pepoint;


//ScalarMul
void ScalarMul_Bin_L2R(big k, pepoint P, pepoint R);

void ScalarMul_Bin_R2L(big k, pepoint P, pepoint R);

//ShamirMul
inline void PreMul_Bin(pepoint P, pepoint Q, pepoint *plist);


void SingleShrMulBin(PL *shrBin, big k, pepoint P, pepoint R);

void ShamirMul_Bin_ptr(PL *shrBin, big a, pepoint P, big b, pepoint Q, pepoint R);

void ShamirMul_Bin_add(PL &shrBin, big a, pepoint P, big b, pepoint Q, pepoint R);

void ShamirMul_Bin_val(PL shrBin, big a, pepoint P, big b, pepoint Q, pepoint R);

void ShamirMul_Bin(big a, pepoint P, big b, pepoint Q, pepoint R);

void ShamirMul_Bin1(big a, pepoint P, big b, pepoint Q, pepoint R);

//void TestShrMul_Bin(csprng &Rng, pepoint P, big n);

#endif
