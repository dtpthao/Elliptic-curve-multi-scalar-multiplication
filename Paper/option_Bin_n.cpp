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
	msg = "Test ShrMul_Bin\n";
	PL shrBin(8);
	//ecurve2_mult(a, P, Q);
	int count = 0, cmp = 0;
	for (int i = 0; i < 1000; i++) {
		strong_bigrand(&Rng, n, k);
		//ShamirDecompose_n(k, kx, P, Px);		// need to write this first

		ShamirMul_Bin_n(d, &shrBin, kx, Px, R);
		//ecurve2_mult2(a, P, b, Q, R1);
		ecurve2_mult(k, P, R2);
		cmp = epoint2_comp(R2, R);
		count += cmp;
	}
	std::cout << "Cmp: " << count << std::endl;
	shrBin.Destructor();
	for (int i = 0; i < d; i++) {
		mirkill(kx[i]);
		epoint_free(Px[i]);
	}
	mirkill(k);
	epoint_free(R); epoint_free(R1); epoint_free(R2);
}

