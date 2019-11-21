#ifndef _SHR_DECOMP_H
#define _SHR_DECOMP_H

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
void ShamirDecompose(big k, pepoint P, big a, pepoint Q, big b);

/*
Compute a, b
(Copied from function above but no computing Q and R)
(Used in deleted test but I keep it as it doesn't matter anyway)
*/
void ShamirDecompose(big k, big a, big b);


void ShamirDecompose3(big k, big k1, big k2, big k3,
	pepoint P, pepoint P1, pepoint P2, pepoint P3);

void ShamirDecompose_n(int n, big k, big *kx, pepoint P, pepoint *Px);

#endif

