#include "option_JSF.h"
//#include "option_NAF.h"

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

inline void PreMul_JSF_origin(pepoint P, pepoint Q, pepoint *plist)
{
	for (int i = 0; i < 9; i++) {
		plist[i] = epoint_init();
	}

	epoint2_set(0, 0, 1, plist[4]);
	epoint2_copy(P, plist[0]); ecurve2_add(Q, plist[0]);
	epoint2_copy(P, plist[2]); ecurve2_sub(Q, plist[2]);
	epoint2_copy(Q, plist[6]); ecurve2_sub(P, plist[6]);
	ecurve2_sub(plist[0], plist[8]);
	epoint2_copy(P, plist[1]);
	epoint2_copy(Q, plist[3]);
	ecurve2_sub(Q, plist[5]);
	ecurve2_sub(P, plist[7]);

	/*epoint2_copy(P, plist[0]);
	ecurve2_padd(Q, plist[0]);

	epoint2_copy(P, plist[1]);
	epoint2_copy(Q, plist[3]);

	epoint2_copy(Q, plist[5]);
	epoint2_negate(plist[5]);
	epoint2_copy(P, plist[7]);
	epoint2_negate(plist[7]);

	epoint2_copy(P, plist[2]);
	ecurve2_padd(plist[5], plist[2]);

	epoint2_copy(Q, plist[6]);
	ecurve2_padd(plist[7], plist[6]);

	epoint2_copy(plist[0], plist[8]);
	epoint2_negate(plist[8]);*/
	for (int i = 0; i < 9; i++) {
		epoint2_norm(plist[i]);
	}
}

void ShamirMul_JSF_origin(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	pepoint plist[9];
	//char JSFa[MAX_M + 1] = { 0 };
	//char JSFb[MAX_M + 1] = { 0 };
	char JSFa[300] = { 0 };
	char JSFb[300] = { 0 };
	DWORD lenJSF;

	PreMul_JSF_origin(P, Q, plist);
	lenJSF = GenJSF(a, b, JSFa, JSFb);
	epoint2_set(0, 0, 1, R);

	for (int i = lenJSF; i > 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(plist[4 - 3 * JSFa[i] - JSFb[i]], R);
	}
	for (int i = 0; i < 9; i++) epoint_free(plist[i]);
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

	epoint2_copy(Q, opt->plist[6]); 
	ecurve2_padd(opt->plist[7], opt->plist[6]);

	epoint2_copy(opt->plist[0], opt->plist[8]);
	epoint2_negate(opt->plist[8]);
	
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
	mirkill(L1); mirkill(L2);
	mirkill(R1); mirkill(S1);
	return lenJSF;
}

// {P+Q, P, P-Q, Q, 0, -Q, Q-P, -P, -P-Q}
inline void PreMul_JSF1(pepoint P, pepoint Q, PL *opt)
{
	pepoint tmp = epoint_init();
	epoint2_copy(P, opt->plist[0]);
	ecurve2_padd(Q, opt->plist[0]);
	epoint2_norm(opt->plist[0]);

	epoint2_copy(P, opt->plist[1]);
	epoint2_copy(Q, opt->plist[3]);

	epoint2_copy(Q, opt->plist[5]);
	epoint2_negate(opt->plist[5]);		//plist[5] is normalized

	epoint2_copy(P, opt->plist[7]);
	epoint2_negate(opt->plist[7]);

	epoint2_copy(P, opt->plist[2]);
	ecurve2_padd(opt->plist[5], opt->plist[2]);
	epoint2_norm(opt->plist[2]);

	epoint2_copy(Q, opt->plist[6]);
	ecurve2_padd(opt->plist[7], opt->plist[6]);
	epoint2_norm(opt->plist[6]);

	epoint2_copy(opt->plist[0], opt->plist[8]);
	epoint2_negate(opt->plist[8]);
}

void ShamirMul_JSF1(PL *opt, big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char JSFa[300] = { 0 };
	char JSFb[300] = { 0 };
	DWORD lenJSF;
	int index;

	PreMul_JSF1(P, Q, opt);
	lenJSF = GenJSF1(a, b, JSFa, JSFb);
	epoint2_set(0, 0, 1, R);

	for (int i = lenJSF; i > 0; i--) {
		//index = 4 - 3 * JSFa[i] - JSFb[i];
		ecurve2_double(R);
		//if (index) ecurve2_padd(opt->plist[index], R);
		//ecurve2_add(R, R);
		ecurve2_add(opt->plist[4 - 3 * JSFa[i] - JSFb[i]], R);
	}
}