#include "option_JSF.h"
#include "Test.h"

#define TESTJSFs 5000
#define REPEATJSFs 10
void compare_GendJSFs(csprng &Rng, pepoint P, big n, Result &res)
{
	stopWatch timer1, timer2, timer3, timer4, timer5, timer6;
	LONGLONG dur1, min1 = LONG_MAX,
		dur2, min2 = LONG_MAX,
		dur3, min3 = LONG_MAX,
		dur5, min5 = LONG_MAX,
		dur6, min6 = LONG_MAX,
		dur4, min4 = LONG_MAX;
	const int dimens = 5;
	big *kx = new big[dimens];
	pepoint *Px = new pepoint[dimens];
	for (int i = 0; i < dimens; i++) {
		kx[i] = mirvar(0);
		Px[i] = epoint_init();
	}
	big k = mirvar(0);
	pepoint	R = epoint_init(),
		R1 = epoint_init(),
		R2 = epoint_init(),
		R3 = epoint_init(),
		R4 = epoint_init(),
		R5 = epoint_init();

	big k1 = mirvar(0),
		k2 = mirvar(0),
		k3 = mirvar(0),
		a = mirvar(0), b = mirvar(0);
	pepoint P1 = epoint_init(),
		P2 = epoint_init(),
		P3 = epoint_init(),
		Q = epoint_init();

	char **dJSF = new char*[5];
	for (int i = 0; i < 5; i++) {
		dJSF[i] = new char[1600];
	}

	int count1 = 0, count2 = 0;
	for (int i = 0; i < TESTJSFs; i++) {
		strong_bigrand(&Rng, n, k);

		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer1);
			ScalarMul_Bin_L2R(k, P, R1);
			stopTimer(&timer1);
			dur1 = getTickCount(&timer1);
			min1 = (min1 < dur1) ? min1 : dur1;
		}

		ShamirDecompose(k, a, b);
		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer2);
			GenJSF(a, b, dJSF[0], dJSF[1]);
			stopTimer(&timer2);
			dur2 = getTickCount(&timer2);
			min2 = (min2 < dur2) ? min2 : dur2;
		}

		//ShamirDecomposit_n(3, k, g, r, y, P);
		ShamirDecompose_nk(3, k, kx);
		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer3);
			GendJSF(3, kx, dJSF);
			stopTimer(&timer3);
			dur3 = getTickCount(&timer3);
			min3 = (min3 < dur3) ? min3 : dur3;
		}

		//ShamirDecomposit_n(4, k, g, r, y, P);
		ShamirDecompose_nk(4, k, kx);
		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer4);
			GendJSF(4, kx, dJSF);
			stopTimer(&timer4);
			dur4 = getTickCount(&timer4);
			min4 = (min4 < dur4) ? min4 : dur4;
		}

		ShamirDecompose_nk(5, k, kx);
		//ShamirDecomposit_n(5, k, g, r, y, P);
		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer5);
			GendJSF(5, kx, dJSF);
			stopTimer(&timer5);
			dur5 = getTickCount(&timer5);
			min5 = (min5 < dur5) ? min5 : dur5;
		}

		copy(k, kx[0]);
		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer6);
			GendJSF(1, kx, dJSF);
			stopTimer(&timer6);
			dur6 = getTickCount(&timer6);
			min6 = (min6 < dur6) ? min6 : dur6;
		}
		res.t[0] += min1;
		res.t[1] += min2;
		res.t[2] += min3;
		res.t[3] += min4;
		res.t[4] += min5;
		res.t[5] += min6;
	}
	res.t[0] /= TESTJSFs;
	res.t[1] /= TESTJSFs;
	res.t[2] /= TESTJSFs;
	res.t[3] /= TESTJSFs;
	res.t[4] /= TESTJSFs;
	res.t[5] /= TESTJSFs;
	res.p[1] = (res.t[1] / res.t[0]) * 100;
	res.p[2] = (res.t[2] / res.t[0]) * 100;
	res.p[3] = (res.t[3] / res.t[0]) * 100;
	res.p[4] = (res.t[4] / res.t[0]) * 100;
	res.p[5] = (res.t[5] / res.t[0]) * 100;

	mirkill(a); mirkill(b); mirkill(k);
	for (int i = 0; i < dimens; i++) {
		mirkill(kx[i]);
		epoint_free(Px[i]);
	}
	epoint_free(Q);
	epoint_free(R);
}

#define NUM_OF_P 5
void printcompares_GendJSFs(Result res[NUM_OF_P + 1], int *m)
{
	const int lib = 5;
	res[NUM_OF_P].p[0] = 100;

	for (int i = 0; i <= NUM_OF_P; i++) {
		res[i].p[0] = 100;
	}

	for (int i = 1; i <= lib; i++) {
		res[NUM_OF_P].p[i] = res[0].p[i] + res[1].p[i]
			+ res[2].p[i] + res[3].p[i] + res[4].p[i];
		res[NUM_OF_P].p[i] /= NUM_OF_P;
	}

	cout << setw(19) << "" << "Bin1         GenJSF       GenJSF3       GenJSF4       GenJSF5       GenJSF1\n";
	for (int i = 0; i < NUM_OF_P; i++) {
		printf("%2d | %4d: %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f\n",
			i, m[i], res[i].t[0], res[i].t[1], res[i].t[2], res[i].t[3],  res[i].t[4], res[i].t[lib]);
	}
	cout << endl;
	for (int i = 0; i < NUM_OF_P; i++) {
		printf("%2d | %4d: %12.2f%% %12.2f%% %12.2f%% %12.2f%% %12.2f%% %12.2f%%\n",
			i, m[i], res[i].p[0], res[i].p[1], res[i].p[2], res[i].p[3], res[i].p[4], res[i].p[lib]);
	}
	printf("    avg  : %12.2f%% %12.2f%% %12.2f%% %12.2f%% %12.2f%% %12.2f%%\n",
		res[NUM_OF_P].p[0], res[NUM_OF_P].p[1], res[NUM_OF_P].p[2], res[NUM_OF_P].p[3], res[NUM_OF_P].p[4], res[NUM_OF_P].p[lib]);
	cout << endl;
	cout << "Test times: " << TESTJSFs << endl;
}
