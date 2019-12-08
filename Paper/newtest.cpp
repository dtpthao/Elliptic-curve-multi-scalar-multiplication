#include "Test.h"

#define TESTS 5000
#define REPEAT 10
void compares(csprng &Rng, pepoint P, big n, Result &res)
{
	stopWatch timer1, timer2, timer3, timer4, timer5, timer6;
	LONGLONG dur1, min1 = LONG_MAX,
		dur2, min2 = LONG_MAX,
		dur3, min3 = LONG_MAX,
		dur4, min4 = LONG_MAX,
		dur5, min5 = LONG_MAX,
		dur6, min6 = LONG_MAX;
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
		k3 = mirvar(0);
	pepoint P1 = epoint_init(),
		P2 = epoint_init(),
		P3 = epoint_init(),
		Q = epoint_init();

	big a = mirvar(0), b = mirvar(0);

	for (int i = 0; i < TESTS; i++) {
		strong_bigrand(&Rng, n, k);

		for (int j = 0; j < REPEAT; j++) {
			// Bin1
			startTimer(&timer1);
			ScalarMul_Bin_L2R(k, P, R1);
			stopTimer(&timer1);
			dur1 = getTickCount(&timer1);
			min1 = (min1 < dur1) ? min1 : dur1; 

			// Bin2
			ShamirDecompose(k, P, a, Q, b);
			startTimer(&timer2);
			//ShamirDecompose(k, P, a, Q, b);
			ShamirMul_Bin(a, P, b, Q, R2);
			//ShamirMul_JSF(a, P, b, Q, R2);
			stopTimer(&timer2);
			dur2 = getTickCount(&timer2);
			min2 = (min2 < dur2) ? min2 : dur2; 

			// Bin3
			ShamirDecompose3(k, k1, k2, k3, P, P1, P2, P3);
			startTimer(&timer3);
			//ShamirDecompose3(k, k1, k2, k3, P, P1, P2, P3);
			ShamirMul_Bin3(k1, k2, k3, P1, P2, P3, R3);
			//ShamirDecompose_n(3, k, kx, P, Px);
			//ShamirMul_Bin_n(3, kx, Px, R3);
			//ShamirMul_dJSF(3, kx, Px, R3);
			stopTimer(&timer3);
			dur3 = getTickCount(&timer3);
			min3 = (min3 < dur3) ? min3 : dur3; 

			// Bin4
			ShamirDecompose_n(4, k, kx, P, Px);
			startTimer(&timer4);
			//ShamirDecompose_n(4, k, kx, P, Px);
			ShamirMul_Bin_n(4, kx, Px, R4);
			//ShamirDecompose(k, P, a, Q, b);
			//ShamirMul_JSF(&shrOpt, a, P, b, Q, R4);
			//ShamirMul_dJSF(4, kx, Px, R4);
			stopTimer(&timer4);
			dur4 = getTickCount(&timer4);
			min4 = (min4 < dur4) ? min4 : dur4; 

			// Bin5
			ShamirDecompose_n(5, k, kx, P, Px);
			startTimer(&timer5);
			//ShamirDecompose_n(5, k, kx, P, Px);
			ShamirMul_Bin_n(5, kx, Px, R5);
			//ShamirMul_dJSF(5, kx, Px, R5);
			//ShamirDecompose(k, P, a, Q, b);
			//ecurve2_mult2(a, P, b, Q, R5);
			stopTimer(&timer5);
			dur5 = getTickCount(&timer5);
			min5 = (min5 < dur5) ? min5 : dur5;

			// lib1
			startTimer(&timer6);
			//ShamirDecompose_n(3, k, kx, P, Px);
			//ecurve2_multn(3, kx, Px, R);
			//ShamirDecompose(k, P, a, Q, b);
			//ecurve2_mult2(a, P, b, Q, R);
			ecurve2_mult(k, P, R);
			stopTimer(&timer6);
			dur6 = getTickCount(&timer6);
			min6 = (min6 < dur6) ? min6 : dur6;
		}
		res.c[0] += epoint2_comp(R, R1);
		res.c[1] += epoint2_comp(R, R2);
		res.c[2] += epoint2_comp(R, R3);
		res.c[3] += epoint2_comp(R, R4);
		res.c[4] += epoint2_comp(R, R5);
		res.t[0] += min1;
		res.t[1] += min2;
		res.t[2] += min3;
		res.t[3] += min4;
		res.t[4] += min5;
		res.t[5] += min6;
	}
	res.t[0] /= TESTS;
	res.t[1] /= TESTS;
	res.t[2] /= TESTS;
	res.t[3] /= TESTS;
	res.t[4] /= TESTS;
	res.t[5] /= TESTS;
	res.p[1] = (res.t[1] / res.t[0]) * 100;
	res.p[2] = (res.t[2] / res.t[0]) * 100;
	res.p[3] = (res.t[3] / res.t[0]) * 100;
	res.p[4] = (res.t[4] / res.t[0]) * 100;
	res.p[5] = (res.t[5] / res.t[0]) * 100;

	for (int i = 0; i < dimens; i++) {
		mirkill(kx[i]);
		epoint_free(Px[i]);
	}
	mirkill(a); mirkill(b); 
	epoint_free(P1); epoint_free(P2); epoint_free(P3); epoint_free(Q);
	mirkill(k1); mirkill(k2); mirkill(k3);
	epoint_free(R); epoint_free(R1); epoint_free(R2);
	epoint_free(R3); epoint_free(R4); epoint_free(R5);
}


void printcompares(Result res[NUM_OF_EC + 1], int *m)
{
	const int ilib = 5;

	for (int i = 1; i <= ilib; i++) {
		res[NUM_OF_EC].p[i] = res[0].p[i] + res[1].p[i]
			+ res[2].p[i] + res[3].p[i] + res[4].p[i];
		res[NUM_OF_EC].p[i] /= NUM_OF_EC;
	}
	for (int i = 0; i <= NUM_OF_EC; i++) {
		res[i].p[0] = 100;
	}
	for (int i = 0; i < NUM_OF_EC; i++) {
		printf("%2d | %4d: %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f\n",
			i, m[i], res[i].t[0], res[i].t[1], res[i].t[2], res[i].t[3], res[i].t[4], res[i].t[ilib]);
	}
	cout << endl;
	for (int i = 0; i < NUM_OF_EC; i++) {
		printf("%2d | %4d: %12.1f%% %12.1f%% %12.1f%% %12.1f%% %12.1f%% %12.1f%%\n",
			i, m[i], res[i].p[0], res[i].p[1], res[i].p[2], res[i].p[3], res[i].p[4], res[i].p[ilib]);
	}
	printf("    avg  : %12.1f%% %12.1f%% %12.1f%% %12.1f%% %12.1f%% %12.1f%%\n",
		res[NUM_OF_EC].p[0], res[NUM_OF_EC].p[1], res[NUM_OF_EC].p[2], res[NUM_OF_EC].p[3], res[NUM_OF_EC].p[4], res[NUM_OF_EC].p[ilib]);
	cout << endl;
	for (int i = 0; i < NUM_OF_EC; i++) {
		printf("%2d | %4d: %13d %13d %13d %13d %13d %13d\n",
			i, m[i], res[i].c[0], res[i].c[1], res[i].c[2], res[i].c[3], res[i].c[4], res[i].c[ilib]);
	}
}

#define TESTJSFs 5000
#define REPEATJSFs 10
void compare_prepowmodJSFs(csprng &Rng, pepoint P, big n, Result &res)
{
	stopWatch timer1, timer2, timer3, timer4, timer5;
	LONGLONG dur1, min1 = LONG_MAX,
		dur2, min2 = LONG_MAX,
		dur3, min3 = LONG_MAX,
		dur5, min5 = LONG_MAX,
		dur4, min4 = LONG_MAX;
	big a = mirvar(1), b = mirvar(1);
	big	g = mirvar(1), k = mirvar(1);
	const int dimens = 5;
	big *kx = new big[dimens];
	pepoint *Px = new pepoint[dimens];
	for (int i = 0; i < dimens; i++) {
		kx[i] = mirvar(0);
		Px[i] = epoint_init();
	}
	pepoint	R = epoint_init();
	pepoint	Q = epoint_init();

	int count1 = 0, count2 = 0;
	for (int i = 0; i < TESTJSFs; i++) {
		strong_bigrand(&Rng, n, k);

		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer1);
			ScalarMul_Bin_L2R(k, P, R);
			stopTimer(&timer1);
			dur1 = getTickCount(&timer1);
			min1 = (min1 < dur1) ? min1 : dur1;
		}

		ShamirDecompose(k, P, a, Q, b);
		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer2);
			PreMul_JSF(P, Q, glob_epoints);
			stopTimer(&timer2);
			dur2 = getTickCount(&timer2);
			min2 = (min2 < dur2) ? min2 : dur2;
		}

		ShamirDecompose_n(3, k, kx, P, Px);
		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer3);
			PreMul_dJSF(3, 27, Px, glob_epoints);
			stopTimer(&timer3);
			dur3 = getTickCount(&timer3);
			min3 = (min3 < dur3) ? min3 : dur3;
		}

		ShamirDecompose_n(4, k, kx, P, Px);
		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer4);
			PreMul_dJSF(4, 81, Px, glob_epoints);
			stopTimer(&timer4);
			dur4 = getTickCount(&timer4);
			min4 = (min4 < dur4) ? min4 : dur4;
		}

		ShamirDecompose_n(5, k, kx, P, Px);
		for (int j = 0; j < REPEATJSFs; j++) {
			startTimer(&timer5);
			PreMul_dJSF(5, 243, Px, glob_epoints);
			stopTimer(&timer5);
			dur5 = getTickCount(&timer5);
			min5 = (min5 < dur5) ? min5 : dur5;
		}
		res.t[0] += min1;
		res.t[1] += min2;
		res.t[2] += min3;
		res.t[3] += min4;
		res.t[4] += min5;
	}
	res.t[0] /= TESTJSFs;
	res.t[1] /= TESTJSFs;
	res.t[2] /= TESTJSFs;
	res.t[3] /= TESTJSFs;
	res.t[4] /= TESTJSFs;
	res.p[1] = (res.t[1] / res.t[0]) * 100;
	res.p[2] = (res.t[2] / res.t[0]) * 100;
	res.p[3] = (res.t[3] / res.t[0]) * 100;
	res.p[4] = (res.t[4] / res.t[0]) * 100;

	mirkill(a); mirkill(b);
	mirkill(g); mirkill(k);
	for (int i = 0; i < dimens; i++) {
		mirkill(kx[i]);
		epoint_free(Px[i]);
	}
	epoint_free(Q);
	epoint_free(R);
}

void printcompares_JSFs(Result res[NUM_OF_EC + 1], int *m)
{
	const int lib = 4;

	for (int i = 0; i <= NUM_OF_EC; i++) {
		res[i].p[0] = 100;
	}

	for (int i = 1; i <= lib; i++) {
		res[NUM_OF_EC].p[i] = res[0].p[i] + res[1].p[i]
			+ res[2].p[i] + res[3].p[i] + res[4].p[i];
		res[NUM_OF_EC].p[i] /= NUM_OF_EC;
	}

	cout << setw(19) << "SclMul          " << "Bin1         preJSF       preJSF3       preJSF4       preJSF5\n";
	for (int i = 0; i < NUM_OF_EC; i++) {
		printf("%2d | %4d: %13.3f %13.3f %13.3f %13.3f %13.3f\n",
			i, m[i], res[i].t[0], res[i].t[1], res[i].t[2], res[i].t[3], res[i].t[lib]);
	}
	cout << endl;
	for (int i = 0; i < NUM_OF_EC; i++) {
		printf("%2d | %4d: %12.1f%% %12.1f%% %12.1f%% %12.1f%% %12.1f%%\n",
			i, m[i], res[i].p[0], res[i].p[1], res[i].p[2], res[i].p[3], res[i].p[lib]);
	}
	printf("    avg  : %12.1f%% %12.1f%% %12.1f%% %12.1f%% %12.1f%%\n",
		res[NUM_OF_EC].p[0], res[NUM_OF_EC].p[1], res[NUM_OF_EC].p[2], res[NUM_OF_EC].p[3], res[NUM_OF_EC].p[lib]);
	cout << endl;
	cout << "Test times: " << TESTJSFs << endl;
}

#define TEST_A_D 50000
void compare_Doub_Add(csprng &Rng, pepoint P, Result &res, int m)
{
	pepoint	A = epoint_init(), D = epoint_init();
	big k = mirvar(0), a = mirvar(0), b = mirvar(0);
	stopWatch timer1, timer2;
	LONGLONG dur1, dur2;
	double t1 = 0, t2 = 0;

	for (int i = 0; i < TEST_A_D; i++) {
		strong_bigdig(&Rng, m >> 1, 2, a);
		strong_bigdig(&Rng, m >> 1, 2, b);
		ecurve2_mult(a, P, A);
		ecurve2_mult(b, P, D);
		dur1 = 0; dur2 = 0;
		for (int j = 0; j < 10; j++) {
			startTimer(&timer1);
			ecurve2_double(D);
			stopTimer(&timer1);
			dur1 += getTickCount(&timer1);

			startTimer(&timer2);
			ecurve2_padd(D, A);
			stopTimer(&timer2);
			dur2 += getTickCount(&timer2);
		}
		res.t[0] += (double)dur1 / 10;
		res.t[1] += (double)dur2 / 10;
	}
	res.t[0] /= TEST_A_D;
	res.t[1] /= TEST_A_D;
	res.p[1] = (res.t[1] / res.t[0]) * 100;

	mirkill(a); mirkill(b); mirkill(k);
	epoint_free(A); epoint_free(D);
}

void print_A_D(Result res[NUM_OF_EC + 1], int *m) 
{
	res[NUM_OF_EC].p[1] = 0;
	res[NUM_OF_EC].p[0] = 100;
	for (int i = 0; i < NUM_OF_EC; i++) {
		res[i].p[0] = 100;
		res[NUM_OF_EC].p[1] += res[i].p[1];
	}
	res[NUM_OF_EC].p[1] /= NUM_OF_EC;

	cout << setw(17) << "" << "Double       Addition\n";

	for (int i = 0; i < NUM_OF_EC; i++) {
		printf("%2d | %4d: %13.3f %13.3f\n", i, m[i], res[i].t[0], res[i].t[1]);
	}
	cout << endl;
	for (int i = 0; i < NUM_OF_EC; i++) {
		printf("%2d | %4d: %12.1f%% %12.1f%%\n", i, m[i], res[i].p[0], res[i].p[1]);
	}
	printf("    avg  : %12.1f%% %12.1f%%\n", res[NUM_OF_EC].p[0], res[NUM_OF_EC].p[1]);

	cout << "Test times: " << (TEST_A_D * 10) << endl;
}


