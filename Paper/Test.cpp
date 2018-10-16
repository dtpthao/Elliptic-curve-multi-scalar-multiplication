#include "Test.h"

/*
 * Compare single scalar multiplication L2R and R2L
 */
void TestSclBin(csprng &Rng, pepoint P, big n, Result &res)
{
	pepoint R = epoint_init(),
		R1 = epoint_init();
	big k = mirvar(1);
	PointList plShrBin(4);
	res.t[0] = res.t[1] = res.t[2] = 0;
	pepoint Q = epoint_init();
	big a = mirvar(1), b = mirvar(1);
	cout << "Compare single scalar multiplication L2R and R2L" << endl;
	for (int i = 0; i < TESTBIN; i++) {
		strong_bigrand(&Rng, n, k);

		if (k->len == 0 || (k->len == 1 && k->w[0] == 1)) {
			i--;
			continue;
		}
		SclDuration(k, P, R, ecurve2_mult, res.t[2]);

		SclDuration(k, P, R1, ScalarMul_Bin_R2L, res.t[1]);
		res.c[1] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_Bin_L2R, res.t[0]);
		res.c[0] += epoint2_comp(R1, R);
	}
	res.t[2] /= TESTBIN;
	res.t[0] /= TESTBIN;
	res.t[1] /= TESTBIN;
	res.p[0] = 100;
	res.p[1] = (res.t[1] / res.t[0]) * 100;
	res.p[2] = (res.t[2] / res.t[0]) * 100;
	plShrBin.Destructor();
	epoint_free(R); epoint_free(R1);
	mirkill(k);
}

/*
 * Compare single scalar multiplication via shamir method 
 * (as 100%, includes shamir decomposit and double scalar mul.)
 * and it's sub-functions.
 * NOT FOUND the problem
 */
//void TestBin1(csprng &Rng, pepoint P, big n, Result &res)
//{
//	pepoint R = epoint_init(),
//		R1 = epoint_init();
//	big k = mirvar(1);
//	PointList plShrBin(4);
//	res.t[0] = res.t[1] = res.t[2] = 0;
//	pepoint Q = epoint_init();
//	big a = mirvar(1), b = mirvar(1);
//
//	for (int i = 0; i < TESTBIN; i++) {
//		strong_bigrand(&Rng, n, k);
//
//		if (k->len == 0 || (k->len == 1 && k->w[0] == 1)) {
//			i--;
//			continue;
//		}
//
//		ShrDuration(k, P, R1, &plShrBin, res.t[0]);
//		ShrDecompDuration(k, P, a, Q, b, res.t[1]);
//		ShrDuration(&plShrBin, a, P, b, Q, R1, res.t[2]);
//	}
//	res.t[2] /= TESTBIN;
//	res.t[0] /= TESTBIN;
//	res.t[1] /= TESTBIN;
//	res.p[0] = 100;
//	res.p[1] = (res.t[1] / res.t[0]) * 100;
//	res.p[2] = (res.t[2] / res.t[0]) * 100;
//	plShrBin.Destructor();
//	epoint_free(R); epoint_free(R1);
//	mirkill(k);
//}

void Test(csprng &Rng, pepoint P, big n, Result &res, string &msg)
{
	pepoint R = epoint_init(),
		R1 = epoint_init();
	big k = mirvar(1);
	PointList plShrBin(9);
	PointList plShrJSF(9);

	for (int i = 0; i < TEST; i++) {
		strong_bigrand(&Rng, n, k);

		if (k->len == 0 || (k->len == 1 && k->w[0] == 1)) {
			i--;
			continue;
		}

		SclDuration(k, P, R, ecurve2_mult, res.t[LIB]);

		ShrDuration(k, P, R1, &plShrBin, ShamirMul_Bin_ptr, res.t[1]);
		res.c[1] += epoint2_comp(R1, R);

		//ShrDuration(k, P, R1, ShamirMul_JSF, res.t[1]);
		//res.c[1] += epoint2_comp(R1, R);

		ShrDuration(k, P, R1, &plShrJSF, ShamirMul_JSF1, res.t[2]);
		res.c[2] += epoint2_comp(R1, R);

		//ShrDuration(k, P, R1, &plShrJSF, ShamirMul_JSF, res.t[0]);
		//res.c[0] += epoint2_comp(R1, R);

		/*ShrDuration(k, P, R1, ShamirMul_JSF_origin, res.t[0]);
		res.c[0] += epoint2_comp(R1, R);*/

		SclDuration(k, P, R1, ScalarMul_Bin_L2R, res.t[0]);
		res.c[0] += epoint2_comp(R1, R);
	}
	msg = "\t\tScalarBin	ShamirBin	ShamirJSF	 Lib";
	for (int i = 0; i <= LIB; i++) {
		res.t[i] /= TESTBIN;
		res.p[i] = (res.t[i] / res.t[0]) * 100;
	}
	
	plShrJSF.Destructor();
	plShrBin.Destructor();
	epoint_free(R); epoint_free(R1);
	mirkill(k);
}

/*_______________________Old code________________________________________*/

//void TestAllScalar(big k, pepoint P, csprng &Rng, big n,
//	double tmin[LIB + 1], double per[LIB + 1], unsigned int correct[LIB])
//{
//	pepoint R = epoint_init();
//	pepoint R1 = epoint_init();
//	for (int i = 0; i <= LIB; i++) {
//		tmin[i] = 100;
//	}
//	for (int i = 0; i < TESTALL; i++) {
//		strong_bigrand(&Rng, n, k);
//		//strong_bigdig(&Rng, 7, 16, k);
//		SclDuration(k, P, R, ecurve2_mult, tmin[LIB]);
//
//		SclDuration(k, P, R1, ScalarMul_Bin, tmin[SclBin]);
//		correct[SclBin] += epoint2_comp(R1, R);
//
//		ShrDuration(k, P, R1, ShamirMul_Bin, tmin[ShrBin]);
//		correct[ShrBin] += epoint2_comp(R1, R);
//
//		ShrDuration(k, P, R1, ShamirMul_JSF, tmin[ShrJSF]);
//		correct[ShrJSF] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_NAF, tmin[SclNAF]);
//		correct[SclNAF] += epoint2_comp(R1, R);
//
//		ShrDuration(k, P, R1, ShamirMul_NAF, tmin[ShrNAF]);
//		correct[ShrNAF] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_wNAF, tmin[SclwNAF]);
//		correct[SclwNAF] += epoint2_comp(R1, R);
//
//		ShrDuration(k, P, R1, ShamirMul_wNAF, tmin[ShrwNAF]);//*********
//		correct[ShrwNAF] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_JT, tmin[SclJT]);
//		correct[SclJT] += epoint2_comp(R1, R);
//
//		ShrDuration(k, P, R1, ShamirMul_JT, tmin[ShrJT]);
//		correct[ShrJT] += epoint2_comp(R1, R);
//
//		ShrDuration(k, P, R1, ShamirMul_BNBC, tmin[ShrBNBC]);
//		correct[ShrBNBC] += epoint2_comp(R1, R);
//
//		ShrDuration(k, P, R1, ShamirMul_AKDAC, tmin[ShrAKDAC]);
//		correct[ShrAKDAC] += epoint2_comp(R1, R);
//	}
//	for (int i = 0; i <= LIB; i++) {
//		per[i] = (tmin[i] / tmin[0]) * 100;
//	}
//	epoint_free(R); epoint_free(R1);
//}
//
//void TestScalar(big k, pepoint P, csprng &Rng, big n,
//	double tmin[LIB + 1], double per[LIB + 1], unsigned int correct[LIB])
//{
//	pepoint R = epoint_init(),
//			R1 = epoint_init();
//	for (int i = 0; i <= LIB; i++) {
//		tmin[i] = 100;
//	}
//	for (int i = 0; i < TESTSCL; i++) {
//		strong_bigrand(&Rng, n, k);
//		//strong_bigdig(&Rng, 7, 16, k);
//		SclDuration(k, P, R, ecurve2_mult, tmin[LIB]);
//
//		SclDuration(k, P, R1, ScalarMul_Bin, tmin[SclBin]);
//		correct[SclBin] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_NAF, tmin[SclNAF]);
//		correct[SclNAF] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_wNAF, tmin[SclwNAF]);
//		correct[SclwNAF] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_JT, tmin[SclJT]);
//		correct[SclJT] += epoint2_comp(R1, R);
//	}
//	for (int i = 0; i <= LIB; i++) {
//		per[i] = (tmin[i] / tmin[0]) * 100;
//	}
//	epoint_free(R); epoint_free(R1);
//}
//
//void TestScalarwNAF(big k, pepoint P, csprng &Rng, big n,
//	double tmin[7], double per[7], unsigned int correct[7])
//{
//	pepoint R = epoint_init(),
//		R1 = epoint_init();
//	for (int i = 0; i <= 6; i++) tmin[i] = 100;
//	for (int i = 0; i < TESTwNAF; i++) {
//		strong_bigrand(&Rng, n, k);
//		SclDuration(k, P, R, ecurve2_mult, tmin[6]);
//
//		SclDuration(k, P, R1, ScalarMul_Bin, tmin[SclBin]);
//		correct[SclBin] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_w3NAF, tmin[1]);
//		correct[1] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_w4NAF, tmin[2]);
//		correct[2] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_w5NAF, tmin[3]);
//		correct[3] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_w6NAF, tmin[4]);
//		correct[4] += epoint2_comp(R1, R);
//
//		SclDuration(k, P, R1, ScalarMul_w7NAF, tmin[5]);
//		correct[5] += epoint2_comp(R1, R);
//	}
//	for (int i = 0; i <= 6; i++) {
//		per[i] = (tmin[i] / tmin[0]) * 100;
//	}
//	epoint_free(R); epoint_free(R1);
//}
//
//void TestShrwNAF(big k, pepoint P, csprng &Rng, big n,
//	double tmin[7], double per[7], unsigned int correct[7])
//{
//	pepoint R = epoint_init(),
//		R1 = epoint_init(),
//		Q = epoint_init();
//	big a = mirvar(1);
//	big b = mirvar(1);
//	for (int i = 0; i <= 6; i++) tmin[i] = 100;
//	for (int i = 0; i < TESTwNAF; i++) {
//		strong_bigrand(&Rng, n, k);
//		strong_bigrand(&Rng, n, a);
//		strong_bigrand(&Rng, n, b);
//		ecurve2_mult(k, P, Q);
//		ShrDuration(a, P, b, Q, R, ecurve2_mult2, tmin[6]);
//		epoint2_set(0, 0, 1, R1);
//		ShrDuration(a, P, b, Q, R1, ShamirMul_Bin, tmin[0]);
//		correct[0] += epoint2_comp(R1, R);
//		epoint2_set(0, 0, 1, R1);
//		ShrDuration(a, P, b, Q, R1, ShamirMul_w3NAF, tmin[1]);//********
//		correct[1] += epoint2_comp(R1, R);
//
//		epoint2_set(0, 0, 1, R1);
//		ShrDuration(a, P, b, Q, R1, ShamirMul_w4NAF, tmin[2]);//********
//		correct[2] += epoint2_comp(R1, R);
//
//		epoint2_set(0, 0, 1, R1);
//		ShrDuration(a, P, b, Q, R1, ShamirMul_w5NAF, tmin[3]);//********
//		correct[3] += epoint2_comp(R1, R);
//
//		epoint2_set(0, 0, 1, R1);
//		ShrDuration(a, P, b, Q, R1, ShamirMul_w6NAF, tmin[4]);//********
//		correct[4] += epoint2_comp(R1, R);
//
//		epoint2_set(0, 0, 1, R1);
//		ShrDuration(a, P, b, Q, R1, ShamirMul_w7NAF, tmin[5]);//********
//		correct[5] += epoint2_comp(R1, R);
//	}
//	for (int i = 0; i <= 6; i++) {
//		per[i] = (tmin[i] / tmin[0]) * 100;
//	}
//	epoint_free(R); epoint_free(R1);
//}
//
//void TestShamir(big k, pepoint P, csprng &Rng, big n,
//	double tmin[LIB + 1], double per[LIB + 1], unsigned int correct[LIB]) 
//{
//	pepoint R = epoint_init(),
//		R1 = epoint_init(),
//		Q = epoint_init();
//	big a = mirvar(1), b = mirvar(1);
//	for (int i = 0; i <= LIB; i++) tmin[i] = 100;
//	
//	for (int i = 0; i < TESTSHR; i++) {
//		strong_bigrand(&Rng, n, k);
//		strong_bigrand(&Rng, n, a);
//		strong_bigrand(&Rng, n, b);
//		//ShamirMul(k, P, a, Q, b);
//		ecurve2_mult(k, P, Q);
//		ShrDuration(a, P, b, Q, R, ecurve2_mult2, tmin[LIB]);
//
//		ShrDuration(a, P, b, Q, R1, ShamirMul_Bin, tmin[ShrBin]);
//		correct[ShrBin] += epoint2_comp(R1, R);
//
//		ShrDuration(a, P, b, Q, R1, ShamirMul_JSF, tmin[ShrJSF]);
//		correct[ShrJSF] += epoint2_comp(R1, R);
//
//		epoint2_set(0, 0, 1, R1);
//		ShrDuration(a, P, b, Q, R1, ShamirMul_NAF, tmin[ShrNAF]);
//		correct[ShrNAF] += epoint2_comp(R1, R);
//
//		ShrDuration(a, P, b, Q, R1, ShamirMul_wNAF, tmin[ShrwNAF]);//*******
//		correct[ShrwNAF] += epoint2_comp(R1, R);
//
//		ShrDuration(a, P, b, Q, R1, ShamirMul_JT, tmin[ShrJT]);
//		correct[ShrJT] += epoint2_comp(R1, R);
//
//		ShrDuration(a, P, b, Q, R1, ShamirMul_BNBC, tmin[ShrBNBC]);
//		correct[ShrBNBC] += epoint2_comp(R1, R);
//
//		ShrDuration(a, P, b, Q, R1, ShamirMul_AKDAC, tmin[ShrAKDAC]);
//		correct[ShrAKDAC] += epoint2_comp(R1, R);
//	}
//	for (int i = ShrBin; i <= LIB; i++) {
//		per[i] = (tmin[i] / tmin[ShrBin]) * 100;
//	}
//
//	epoint_free(Q);
//	mirkill(a); mirkill(b);
//	epoint_free(R); epoint_free(R1);
//}
//
//void TestGenForm(big k, csprng &Rng, big n, double tmin[LIB + 1], DWORD tick[LIB + 1]) {
//
//	big a = mirvar(1), b = mirvar(1);
//	double dur;
//	DWORD qPart;
//	stopWatch timer;
//	Vector* CS;
//	DiffSeq DS;
//	char *JSFa = new char[(MAX_M + 1) / 2];
//	char *JSFb = new char[(MAX_M + 1) / 2];
//	for (int i = 0; i <= LIB; i++) {
//		tmin[i] = 100; tick[i] = 1000;
//	}
//	for (int i = 0; i < TESTGF; i++) {
//		strong_bigrand(&Rng, n, k);
//		ShamirDecomposit(k, a, b);
//		GenFormSclDuration(k, GenNAF, tmin[SclNAF], tick[SclNAF]);
//		GenFormSclDuration(k, GenwNAF, tmin[SclwNAF], tick[SclwNAF]);
//		GenFormSclDuration(k, GenJT, tmin[SclJT], tick[SclJT]);
//		
//		GenFormShrDuration(a, b, GenNAF, tmin[ShrNAF], tick[ShrNAF]);
//		GenFormShrDuration(a, b, GenwNAF, tmin[ShrwNAF], tick[ShrwNAF]);
//		GenFormShrDuration(a, b, GenJT, tmin[ShrJT], tick[ShrJT]);
//
//		startTimer(&timer);
//		GenJSF(a, b,JSFa, JSFb);
//		stopTimer(&timer);
//		dur = getElapsedTime(&timer) * 1000;
//		tmin[ShrJSF] = (tmin[ShrJSF] < dur) ? tmin[ShrJSF] : dur;
//		qPart = getQuadPart(&timer);
//		tick[ShrJSF] = (tick[ShrJSF] < qPart) ? tick[ShrJSF] : qPart;
//
//		CS = new Vector[k->len <<4];
//		startTimer(&timer);
//		GenBNBC(a, b, CS, DS);
//		stopTimer(&timer);
//		dur = getElapsedTime(&timer) * 1000;
//		tmin[ShrBNBC] = (tmin[ShrBNBC] < dur) ? tmin[ShrBNBC] : dur;
//		qPart = getQuadPart(&timer);
//		tick[ShrBNBC] = (tick[ShrBNBC] < qPart) ? tick[ShrBNBC] : qPart;
//	}
//	mirkill(a); mirkill(b);
//}
//
