#include "option_JSF.h"

DWORD GenJSF(big R, big S, char *JSFr, char *JSFs)
{
	big R1 = mirvar(0), S1 = mirvar(0);
	DWORD lenJSF = 0;

	copy(R, R1); copy(S, S1);
	while (S1->len > 0 || R1->len > 0) {
		lenJSF++;
		JSFr[lenJSF] = R1->w[0] & 1;
		JSFs[lenJSF] = S1->w[0] & 1;
		if (JSFr[lenJSF] & JSFs[lenJSF]) {
			if (R1->w[0] & 2) JSFr[lenJSF] = -JSFr[lenJSF];
			if (S1->w[0] & 2) JSFs[lenJSF] = -JSFs[lenJSF];
		}
		else if (JSFr[lenJSF] ^ JSFs[lenJSF]) {
			if ((R1->w[0] & 2) ^ (S1->w[0] & 2)) {
				JSFr[lenJSF] = -JSFr[lenJSF];
				JSFs[lenJSF] = -JSFs[lenJSF];
			}
		}
		sftbit(R1, -1, R1);
		sftbit(S1, -1, S1);
		if (JSFr[lenJSF] == -1) incr(R1, 1, R1);
		if (JSFs[lenJSF] == -1) incr(S1, 1, S1);
	}
	mirkill(R1); mirkill(S1);
	return lenJSF;
}

// {P+Q, P, P-Q, Q, 0, -Q, Q-P, -P, -P-Q}
inline void PreMul_JSF(pepoint P, pepoint Q, pepoint *plist)
{
	pepoint tmp = epoint_init();
	epoint2_copy(P, plist[0]);
	ecurve2_padd(Q, plist[0]);
	epoint2_norm(plist[0]);

	epoint2_copy(P, plist[1]);
	epoint2_copy(Q, plist[3]);

	epoint2_copy(Q, plist[5]);
	epoint2_negate(plist[5]);		//plist[5] is normalized

	epoint2_copy(P, plist[7]);
	epoint2_negate(plist[7]);

	epoint2_copy(P, plist[2]);
	ecurve2_padd(plist[5], plist[2]);
	epoint2_norm(plist[2]);

	epoint2_copy(Q, plist[6]);
	ecurve2_padd(plist[7], plist[6]);
	epoint2_norm(plist[6]);

	epoint2_copy(plist[0], plist[8]);
	epoint2_negate(plist[8]);
}

void ShamirMul_JSF(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char JSFa[1000] = { 0 };
	char JSFb[1000] = { 0 };
	DWORD lenJSF;
	int index;

	PreMul_JSF(P, Q, glob_epoints);
	lenJSF = GenJSF(a, b, JSFa, JSFb);
	index = 4 - 3 * JSFa[lenJSF] - JSFb[lenJSF];
	epoint_copy(glob_epoints[index], R);
	for (int i = lenJSF - 1; i > 0; i--) {
		index = 4 - 3 * JSFa[i] - JSFb[i];
		ecurve2_double(R);
		if (index != 4) ecurve2_padd(glob_epoints[index], R);
	}
}

/*************************************************************************************************/
/**__________________________________Another Option_____________________________________________**/
/**                                                                                             **/

inline void subGenJSF(big Var1, big Var2, char &JSFi, bool &d, big RS)
{
	DWORD v1 = Var1->w[0], v2 = Var2->w[0];
	if (!(v1 & 1)) JSFi = 0;
	else {
		JSFi = (v1 & 2) ? -1 : 1;
		if ((!(v1 & 7 ^ 3) || !(v1 & 7 ^ 5)) && ((v2 & 3) == 2))
			JSFi = -JSFi;
	}
	if (((int)d << 1) == (JSFi + 1)) d = 1 - d;
	sftbit(RS, -1, RS);
}

DWORD GenJSF2(big R, big S, char *JSFr, char *JSFs)
{
	big L1 = mirvar(1), L2 = mirvar(1),
		R1 = mirvar(0), S1 = mirvar(0);
	bool d1 = 0, d2 = 0;
	DWORD lenJSF = 0;

	copy(R, L1); copy(S, L2); copy(R, R1); copy(S, S1);
	while (L1->len > 0 || L2->len > 0) {
		lenJSF++;
		subGenJSF(L1, L2, JSFr[lenJSF], d1, R1);
		subGenJSF(L2, L1, JSFs[lenJSF], d2, S1);
		incr(R1, d1, L1);
		incr(S1, d2, L2);
	}
	mirkill(L1); mirkill(L2);
	mirkill(R1); mirkill(S1);
	return lenJSF;
}

/**_______________________________End "Another Option"__________________________________________**/
/*************************************************************************************************/


