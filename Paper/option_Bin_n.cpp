#include "option_Bin.h"

inline void PreMul_Bin_n(int n, pepoint *P, pepoint *plist)
{
	epoint2_copy(P[1], plist[1]);

	int idx;
	for (int i = 1; i < n; i++) {
		idx = 1 << i;
		epoint_copy(P[i], plist[idx]);
		for (int j = 1; j < idx; j++) {
			epoint_copy(plist[idx], plist[idx + j]);
			ecurve2_padd(plist[j], plist[idx + j]);
			epoint2_norm(plist[idx + j]);
		}
	}
	//for (int i = 0; i < 8; i++) {
	//	std::cout << "pl[" << i << "] :\n"; cotnumEp(plist[i]);
	//}
}

void ShamirMul_Bin_n(int n, PL *shrBin, big *k, pepoint *P, pepoint R)
{
	int i, j = 0, ii = 0, index = 0, tmp;
	DWORD *w = new DWORD[n];

	i = k[n - 1]->len - 1;
	tmp = k[n - 1]->w[i];
	while (tmp >> j && j != 32) j++;

	PreMul_Bin_n(n, P, shrBin->plist);

	for (; ii < n; ii++) index += (k[ii]->w[i] >> j) & 1 << ii;
	if (j == 0) { i--; j = 32; }
	epoint2_copy(shrBin->plist[index], R);

	for (--j; i >= 0; i--, j = 31) {
		for (ii = 0; ii < n; ii++) w[ii] = k[ii]->w[i];
		while (j) {
			index = (w[0] >> j) & 1;
			for (ii = 1; ii < n; ii++) {
				index += (w[ii] >> j) & 1 << ii;
			}
			ecurve2_double(R);
			if (index) ecurve2_padd(shrBin->plist[index], R);	// it should be also normalized, but I'll see later
			//cout << "R\n"; cotnumEp(R);
			j--;
		}
		/* in case j = 0 */
		index = w[0] & 1;
		for (ii = 1; ii < n; ii++) index += (w[ii] & 1) << ii;
		ecurve2_double(R);
		if (index) ecurve2_padd(shrBin->plist[index], R);
		//cout << "R\n"; cotnumEp(R);
	}
}

void test_bin_n(int d, csprng &Rng, pepoint P, big n, string msg)
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

	msg = "Test ShrMul_Bin\n";
	PL shrBin(8);
	PL shrBin2(8);
	//ecurve2_mult(a, P, Q);
	int count = 0, cmp = 0;
	for (int i = 0; i < 1; i++) {
		strong_bigrand(&Rng, n, k);
		ShamirDecompose_n(d, k, kx, P, Px);
		ShamirDecompose3(k, k1, k2, k3, P, P1, P2, P3);
		//std::cout << "k1: "; cotnum(kx[0], stdout);
		//std::cout << "k1: "; cotnum(k1, stdout);
		//std::cout << "k2: "; cotnum(kx[1], stdout);
		//std::cout << "k2: "; cotnum(k2, stdout);
		//std::cout << "k3: "; cotnum(kx[2], stdout);
		//std::cout << "k3: "; cotnum(k3, stdout);
		//std::cout << "P1: "; cotnumEp(Px[0]);
		//std::cout << "P1: "; cotnumEp(P1);
		//std::cout << "P2: "; cotnumEp(Px[1]);
		//std::cout << "P2: "; cotnumEp(P2);
		//std::cout << "P3: "; cotnumEp(Px[2]);
		//std::cout << "P3: "; cotnumEp(P3);

		ShamirMul_Bin3_ptr(&shrBin2, k1, k2, k3, P1, P2, P3, R);
		ShamirMul_Bin_n(d, &shrBin, kx, Px, R);
		for (int i = 0; i < 8; i++) {
			std::cout << "pl[" << i << "] :\n"; cotnumEp(shrBin2.plist[i]);
			std::cout << "pl[" << i << "] :\n"; cotnumEp(shrBin.plist[i]);
		}
		ecurve2_mult(k, P, R2);
		//cmp = epoint2_comp(R2, R);
		count += cmp;
	}
	std::cout << "Cmp: " << count << std::endl;
	shrBin.Destructor();
	shrBin2.Destructor();
	for (int i = 0; i < d; i++) {
		mirkill(kx[i]);
		epoint_free(Px[i]);
	}
	mirkill(k);
	epoint_free(R); epoint_free(R1); epoint_free(R2);
}
