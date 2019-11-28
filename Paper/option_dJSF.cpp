#include "option_JSF.h"

// Generate d-dimensional Simple JSF Form
// d < 32
DWORD GendJSF(int d, big *r, int **dJSF)
{
	int i, j, a0, a1, tmp;
	int A[2];
	DWORD lenJSF;
	big *x = new big[d];

	i = j = lenJSF = a0 = A[0] = 0;
	bool loop = false;
	for (; i < d; i++) {
		x[i] = mirvar(0);
		copy(r[i], x[i]);
		if (x[i]->w[0] & 1) A[a0] += 1 << i;
		loop |= !(x[i]->len == 1 && x[i]->w[0] == 0 || x[i]->len == 0);
	}

	while (loop) {
		a1 = a0 ^ 1;
		A[a1] = 0;
		for (i = 0; i < d; i++) {
			dJSF[i][lenJSF] = x[i]->w[0] & 1;
			if (x[i]->w[0] & 2) A[a1] += 1 << i;
		}

		/* check if A[a1] is a subset of A[a0] */
		if ((A[a1] | A[a0]) == A[a0]) {
			tmp = A[a1];
			A[a1] = 0;
		}
		else { /* relative complement A[a0]\A[a1] */
			tmp = (A[a0] | A[a1]) ^ A[a1];
			A[a1] |= A[a0];
		}

		i = 0;
		while (tmp) {
			if (tmp & 1) dJSF[i][lenJSF] = -dJSF[i++][lenJSF];
			tmp >>= 1;
		}

		loop = false;
		for (i = 0; i < d; i++) {
			sftbit(x[i], -1, x[i]);
			if (dJSF[i][lenJSF] == -1) {
				if (x[i]->len == 0) x[i]->len = 1;
				if (x[i]->w[0] ^ 0xffffffff) x[i]->w[0]++;
				else incr(x[i], 1, x[i]);
			}
			loop |= !(x[i]->len == 1 && x[i]->w[0] == 0 || x[i]->len == 0);
		}
		lenJSF++;
		a0 = a1;
	}
	for (i = 0; i < d; i++) mirkill(x[i]);

	return lenJSF;
}

void PreMul_dJSF(int d, int len, pepoint *P, pepoint *plist)
{
	int i, j, k, i0 = len >> 1, upi, downi;

	for (i = 0, j = 1; i < d; i++, j *= 3) {
		upi = i0 + j;
		downi = i0 - j;
		epoint_copy(P[i], plist[downi]);								// y[i]
		epoint_copy(P[i], plist[upi]);
		epoint2_negate(plist[upi]);		//negated point is normalized
		for (k = 1; k <= (j >> 1); k++) {
			epoint_copy(plist[upi], plist[upi + k]);
			ecurve2_padd(plist[i0 + k], plist[upi + k]);

			epoint_copy(plist[upi], plist[upi - k]);
			ecurve2_padd(plist[i0 - k], plist[upi - k]);

			epoint_copy(plist[downi], plist[downi + k]);
			ecurve2_padd(plist[i0 + k], plist[downi + k]);

			epoint_copy(plist[downi], plist[downi - k]);
			ecurve2_padd(plist[i0 - k], plist[downi - k]);

			// checked - no epoint2_norm() needed here
			//epoint2_norm(plist[upi + k]);
			//epoint2_norm(plist[upi - k]);
			//epoint2_norm(plist[downi + k]);
			//epoint2_norm(plist[downi - k]);
		}
	}
}

// R = k1*P1 + k2*P2 + ... + kd*Pd
void ShamirMul_dJSF(int d, PL *shrJSF, big *k, pepoint *P, pepoint R)
{
	int tmp = 1, I0, idx = 0, i, j;
	int **dJSF = new int*[d];
	for (i = 0; i < d; i++) dJSF[i] = new int[700];
	DWORD lendJSF;
	for (i = 0; i < d; i++) tmp *= 3;

	idx = I0 = tmp >> 1;
	PreMul_dJSF(d, tmp, P, shrJSF->plist);
	lendJSF = GendJSF(d, k, dJSF);
	//for (i = 0; i <= 13; i++) {
	//	cout << "pl[" << (13 - i) << "]: \n"; cotnumEp(shrJSF->plist[13 - i]);
	//	cout << "pl[" << (13 + i) << "]: \n"; cotnumEp(shrJSF->plist[13 + i]);
	//}

	for (j = 0, tmp = 1; j < d; j++, tmp *= 3) {
		idx -= tmp * dJSF[j][lendJSF - 1];
	}
	epoint2_copy(shrJSF->plist[idx], R);
	for (i = lendJSF - 2; i >= 0; i--) {
		idx = I0;
		for (j = 0, tmp = 1; j < d; j++, tmp *= 3) {
			idx -= tmp * dJSF[j][i];
		}
		ecurve2_double(R);
		if (idx != I0) ecurve2_padd(shrJSF->plist[idx], R); // no epoint2_norm() needed here
	}

	for (i = 0; i < d; i++) delete[] dJSF[i];
	delete[] dJSF;
}

#include "option_Bin.h"
void test_dJSF(int d, csprng &Rng, pepoint P, big n, string msg)
{
	big *kx = new big[d];
	pepoint *Px = new pepoint[d];
	for (int i = 0; i < d; i++) {
		kx[i] = mirvar(0);
		Px[i] = epoint_init();
	}
	big k = mirvar(0);
	pepoint	R = epoint_init(),
		R1 = epoint_init(),
		R2 = epoint_init();

	big k1 = mirvar(0),
		k2 = mirvar(0),
		k3 = mirvar(0);
	pepoint P1 = epoint_init(),
		P2 = epoint_init(),
		P3 = epoint_init();

	msg = "Test ShrMul_JSF\n";
	PL shrJSF(27);
	PL shrBin(8);
	//ecurve2_mult(a, P, Q);
	int count = 0, cmp = 0;
	for (int i = 0; i < 5000; i++) {
		strong_bigrand(&Rng, n, k);
		ShamirDecompose_n(d, k, kx, P, Px);
		//ShamirDecompose3(k, k1, k2, k3, P, P1, P2, P3);

		ShamirMul_dJSF(d, &shrJSF, kx, Px, R);
		//ShamirMul_Bin_n(d, &shrBin, kx, Px, R1);
		//ShamirMul_Bin3_ptr(&shrBin2, k1, k2, k3, P1, P2, P3, R1);
		ecurve2_mult(k, P, R2);
		//std::cout << "R : \n"; cotnumEp(R);
		//std::cout << "R1: \n"; cotnumEp(R1);
		//std::cout << "R2: \n"; cotnumEp(R2);
		cmp = epoint2_comp(R2, R);
		count += cmp;
	}
	std::cout << "Cmp: " << count << std::endl;
	shrJSF.Destructor();
	shrBin.Destructor();
	for (int i = 0; i < d; i++) {
		mirkill(kx[i]);
		epoint_free(Px[i]);
	}
	mirkill(k);
	epoint_free(R); epoint_free(R1); epoint_free(R2);
}