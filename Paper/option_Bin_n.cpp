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
	int i, j = 0, ii, index = 0, tmp;
	int *bit = new int[n];
	DWORD *w = new DWORD[n];

	i = k[n - 1]->len - 1;
	tmp = k[n - 1]->w[i];
	while (tmp >> j && j != 32) j++;

	PreMul_Bin_n(n, P, shrBin->plist);
	epoint_set(0, 0, 1, R);

	for (ii = 0; ii < n; ii++) {
		bit[ii] = (k[ii]->w[i] >> j) & 1 << ii;
		index += bit[ii];
	}
	epoint2_copy(shrBin->plist[index], R);
	if (j == 0) {
		i--; j = 32;
	}
	for (--j; i >= 0; i--, j = 31) {
		for (ii = 0; ii < n; ii++) w[ii] = k[ii]->w[i];
		while (j) {
			index = (w[0] >> j) & 1;
			for (ii = 1; ii < n; ii++) {
				bit[ii] = (w[ii] >> j) & 1 << ii;
				index += bit[ii];
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

