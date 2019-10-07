#ifndef _SHR_DECOMPOSIT_H
#define _SHR_DECOMPOSIT_H

#include "Ellipse.h"

/*
Compute a, b, Q and R = aP + bQ
*/
//void ShamirMul(big k, pepoint P, pepoint R,
//	void(*func) (big, pepoint, big, pepoint, pepoint));
void ShamirMul(big k, pepoint P, pepoint R, PL *opt,
	void(*func) (PL *, big, pepoint, big, pepoint, pepoint));

/*
Compute a, b and Q
(Copied from function above but no computing R = aP + bQ)
*/
void ShamirDecomposit(big k, pepoint P, big a, pepoint Q, big b);
void ShamirDecomposit1(big k, pepoint P, big a, pepoint Q, big b);

/*
Compute a, b
(Copied from function above but no computing Q and R)
(Used in deleted test but I keep it as it doesn't matter anyway)
*/
void ShamirDecomposit(big k, big a, big b);

#endif

