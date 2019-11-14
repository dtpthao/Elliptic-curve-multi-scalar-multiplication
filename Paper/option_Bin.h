#ifndef _OPTION_BIN_H
#define _OPTION_BIN_H

#include "Duration.h"

typedef epoint* pepoint;

void ScalarMul_Bin_L2R(big k, pepoint P, pepoint R);
inline void PreMul_Bin(pepoint P, pepoint Q, pepoint *plist);
void ShamirMul_Bin_ptr(PL *shrBin, big a,
	pepoint P, big b, pepoint Q, pepoint R);	//in use
void ShamirMul_Bin(big a, pepoint P, big b, pepoint Q, pepoint R);
void TestShrMul_Bin(csprng &Rng, pepoint P, big n, std::string &msg);


/*******************************************************************************/
/**_________________________The Very First Option_____________________________**/
/**                                                                           **/
void ScalarMul_Bin_R2L(big k, pepoint P, pepoint R);
void PreMul_Bin_old(pepoint P, pepoint Q, pepoint *plist);
void ShamirMul_Bin_old(big a, pepoint P, big b, pepoint Q, pepoint R);
/**_____________________End "The very First Option"___________________________**/
/*******************************************************************************/


#endif
