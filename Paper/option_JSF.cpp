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

/*************************************************************************************************/
/**______________________________________In Use_________________________________________________**/
/**                                                                                             **/

// {P+Q, P, P-Q, Q, 0, -Q, Q-P, -P, -P-Q}
inline void PreMul_JSF(pepoint P, pepoint Q, PL *opt)					//in use
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

void ShamirMul_JSF(PL *opt, big a, pepoint P, big b, pepoint Q, pepoint R)	//in use
{
	char JSFa[1000] = { 0 };
	char JSFb[1000] = { 0 };
	DWORD lenJSF;
	int index;

	PreMul_JSF(P, Q, opt);
	lenJSF = GenJSF(a, b, JSFa, JSFb);
	index = 4 - 3 * JSFa[lenJSF] - JSFb[lenJSF];
	epoint_copy(opt->plist[index], R);
	for (int i = lenJSF - 1; i > 0; i--) {
		index = 4 - 3 * JSFa[i] - JSFb[i];
		ecurve2_double(R);
		if (index != 4) ecurve2_padd(opt->plist[index], R);
	}
}

/**__________________________________End In Use_________________________________________________**/
/*************************************************************************************************/



/*************************************************************************************************/
/**__________________________________Another Option_____________________________________________**/
/**                                                                                             **/

// {P+Q, P, P-Q, Q, 0, -Q, Q-P, -P, -P-Q}
// 8	7	6	5	4	3	2	1	0
inline void PreMul_JSF2(pepoint P, pepoint Q, PL *opt)
{
	pepoint tmp = epoint_init();
	epoint2_copy(P, opt->plist[8]);
	ecurve2_padd(Q, opt->plist[8]);
	epoint2_norm(opt->plist[8]);	//+

	epoint2_copy(P, opt->plist[7]);
	epoint2_copy(Q, opt->plist[5]);

	epoint2_copy(Q, opt->plist[3]);
	epoint2_negate(opt->plist[3]);		//plist[5] is normalized

	epoint2_copy(P, opt->plist[1]);
	epoint2_negate(opt->plist[1]);

	epoint2_copy(P, opt->plist[6]);
	ecurve2_padd(opt->plist[3], opt->plist[6]);
	epoint2_norm(opt->plist[6]);

	epoint2_copy(Q, opt->plist[2]);
	ecurve2_padd(opt->plist[1], opt->plist[2]);
	epoint2_norm(opt->plist[2]);

	epoint2_copy(opt->plist[8], opt->plist[0]); //+
	epoint2_negate(opt->plist[0]);
}

void ShamirMul_JSF2(PL *opt, big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char JSFa[300] = { 0 };
	char JSFb[300] = { 0 };
	DWORD lenJSF;
	int index;
	int JSFai[3] = { 1, 4, 7 };

	PreMul_JSF2(P, Q, opt);
	lenJSF = GenJSF(a, b, JSFa, JSFb);
	epoint2_set(0, 0, 1, R);

	for (int i = lenJSF; i > 0; i--) {
		//index = 4 - 3 * JSFa[i] - JSFb[i];
		index = JSFai[JSFa[i] + 1] + JSFb[i];
		ecurve2_double(R);
		//if (index) ecurve2_padd(opt->plist[index], R);
		//ecurve2_add(R, R);
		//ecurve2_add(opt->plist[4 - 3 * JSFa[i] - JSFb[i]], R);
		ecurve2_add(opt->plist[index], R);
	}
}

/**_______________________________End "Another Option"__________________________________________**/
/*************************************************************************************************/



/*************************************************************************************************/
/**______________________________The Very First Option__________________________________________**/
/**                                                                                             **/

inline void subGenJSF_old(big Var1, big Var2, char *JSF, DWORD &lenJSF, bool &d, big RS) {

	if (!(Var1->w[0] & 1)) JSF[lenJSF] = 0;
	else {
		JSF[lenJSF] = (Var1->w[0] & 2) ? -1 : 1;
		if ((!(Var1->w[0] & 7 ^ 3) || !(Var1->w[0] & 7 ^ 5)) && ((Var2->w[0] & 3) == 2))
			JSF[lenJSF] = -JSF[lenJSF];
	}
	if (((int)d << 1) == (JSF[lenJSF] + 1)) d = 1 - d;
	sftbit(RS, -1, RS);
}

DWORD GenJSF_old(big R, big S, char *JSFr, char *JSFs)
{
	big L1 = mirvar(1), L2 = mirvar(1), R1 = mirvar(0), S1 = mirvar(0);
	bool d1 = 0, d2 = 0;
	DWORD lenJSF = 0;

	copy(R, L1); copy(S, L2); copy(R, R1); copy(S, S1);
	while (L1->len > 0 || L2->len > 0) {
		lenJSF++;
		subGenJSF_old(L1, L2, JSFr, lenJSF, d1, R1);
		subGenJSF_old(L2, L1, JSFs, lenJSF, d2, S1);
		incr(R1, d1, L1);
		incr(S1, d2, L2);
	}
	return lenJSF;

}

inline void PreMul_JSF_old(pepoint P, pepoint Q, pepoint *plist)
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

	for (int i = 0; i < 9; i++) {
		epoint2_norm(plist[i]);
	}
}

void ShamirMul_JSF_old(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	pepoint plist[9];
	//char JSFa[MAX_M + 1] = { 0 };
	//char JSFb[MAX_M + 1] = { 0 };
	char JSFa[300] = { 0 };
	char JSFb[300] = { 0 };
	DWORD lenJSF;

	PreMul_JSF_old(P, Q, plist);
	lenJSF = GenJSF_old(a, b, JSFa, JSFb);
	epoint2_set(0, 0, 1, R);

	for (int i = lenJSF; i > 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(plist[4 - 3 * JSFa[i] - JSFb[i]], R);
	}
	for (int i = 0; i < 9; i++) epoint_free(plist[i]);
}

/**_____________________________End "The very First Option"_____________________________________**/
/*************************************************************************************************/

