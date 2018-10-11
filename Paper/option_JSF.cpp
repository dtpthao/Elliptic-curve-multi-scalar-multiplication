#include "option_JSF.h"

inline void subGenJSF(big Var1, big Var2, char *JSF, DWORD &lenJSF, bool &d, big RS) {
	if (!(Var1->w[0] & 1)) JSF[lenJSF] = 0;
	else {
		JSF[lenJSF] = (Var1->w[0] & 2) ? -1 : 1;
		if ((!(Var1->w[0] & 7 ^ 3) || !(Var1->w[0] & 7 ^ 5)) && ((Var2->w[0] & 3) == 2))
			JSF[lenJSF] = -JSF[lenJSF];
	}
	if (((int)d << 1) == (JSF[lenJSF] + 1)) d = 1 - d;
	sftbit(RS, -1, RS);
}

inline void subGenJSF(big Var1, big Var2, char &JSFi, bool &d, big RS) {
	if (!(Var1->w[0] & 1)) JSFi = 0;
	else {
		JSFi = (Var1->w[0] & 2) ? -1 : 1;
		if ((!(Var1->w[0] & 7 ^ 3) || !(Var1->w[0] & 7 ^ 5)) && ((Var2->w[0] & 3) == 2))
			JSFi = -JSFi;
	}
	if (((int)d << 1) == (JSFi + 1)) d = 1 - d;
	sftbit(RS, -1, RS);
}

DWORD GenJSF(big R, big S, char *JSFr, char *JSFs)
{
	big L1 = mirvar(1), L2 = mirvar(1), R1 = mirvar(0), S1 = mirvar(0);
	bool d1 = 0, d2 = 0;
	DWORD lenJSF = 0;

	copy(R, L1); copy(S, L2); copy(R, R1); copy(S, S1);
	while (L1->len > 0 || L2->len > 0) {
		lenJSF++;
		subGenJSF(L1, L2, JSFr, lenJSF, d1, R1);
		subGenJSF(L2, L1, JSFs, lenJSF, d2, S1);
		incr(R1, d1, L1); 
		incr(S1, d2, L2);
	}
	return lenJSF;

}

DWORD GenJSF1(big R, big S, char *JSFr, char *JSFs)
{
	big L1 = mirvar(1), L2 = mirvar(1), R1 = mirvar(0), S1 = mirvar(0);
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
	return lenJSF;
}

// {P+Q, P, P-Q, Q, 0, -Q, Q-P, -P, -P-Q}
inline void PreMul_JSF(pepoint P, pepoint Q, PL *opt)
{
	epoint2_copy(P, opt->plist[0]); 
	ecurve2_padd(Q, opt->plist[0]);

	epoint2_copy(P, opt->plist[1]);
	epoint2_copy(Q, opt->plist[3]);

	epoint2_copy(Q, opt->plist[5]);
	epoint2_negate(opt->plist[5]);
	epoint2_copy(P, opt->plist[7]);
	epoint2_negate(opt->plist[7]);

	epoint2_copy(P, opt->plist[2]); 
	ecurve2_padd(opt->plist[5], opt->plist[2]);
	//ecurve2_sub(Q, opt->plist[2]);

	epoint2_copy(Q, opt->plist[6]); 
	ecurve2_padd(opt->plist[7], opt->plist[6]);
	//ecurve2_sub(P, opt->plist[6]);

	/*epoint2_set(0, 0, 1, opt->plist[8]);
	ecurve2_sub(opt->plist[0], opt->plist[8]);*/
	epoint2_copy(opt->plist[0], opt->plist[8]);
	epoint2_negate(opt->plist[8]);

	//epoint2_set(0, 0, 1, opt->plist[5]);
	//epoint2_set(0, 0, 1, opt->plist[7]);
	//ecurve2_sub(Q, opt->plist[5]);
	//ecurve2_sub(P, opt->plist[7]);
	
	for (int i = 0; i < 9; i++) {
		epoint2_norm(opt->plist[i]);
	}
}

void ShamirMul_JSF(PL *opt, big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char JSFa[300] = { 0 };
	char JSFb[300] = { 0 };
	DWORD lenJSF;
	int index;

	PreMul_JSF(P, Q, opt);
	lenJSF = GenJSF(a, b, JSFa, JSFb);
	epoint2_set(0, 0, 1, R);

	for (int i = lenJSF; i > 0; i--) {
		//index = 4 - 3 * JSFa[i] - JSFb[i];
		ecurve2_double(R);
		//if (index) ecurve2_padd(opt->plist[index], R);
		//ecurve2_add(R, R);
		ecurve2_add(opt->plist[4 - 3 * JSFa[i] - JSFb[i]], R);
	}
}