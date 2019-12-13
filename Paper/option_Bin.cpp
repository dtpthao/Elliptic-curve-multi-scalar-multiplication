#include "option_Bin.h"

void ScalarMul_Bin_L2R(big k, pepoint P, pepoint R)
{
	unsigned int len, i = 31;

	len = k->w[k->len - 1];
	while (len && !(len & (1 << i))) i--;
	len = (k->len << 5) - (31 - i);
	epoint2_copy(P, R);
	for (int j = len - 2; j >= 0; j--) {
		ecurve2_double(R);
		if (k->w[j >> 5] & (1 << j)) {
			ecurve2_padd(P, R);
		}
	}
}

//This R2L option is much worse than L2R
void ScalarMul_Bin_R2L(big k, pepoint P, pepoint R)
{
	unsigned int len, i = 31;
	pepoint D = epoint_init();
	epoint2_copy(P, D);
	epoint_set(0, 0, 1, R);
	len = k->w[k->len - 1];
	while (!(len & (1 << i))) i--;
	len = (k->len << 5) - (31 - i);
	for (i = 0; i < len; i++) {
		if (k->w[i >> 5] & (1 << i)) {
			ecurve2_add(D, R);
		}
		ecurve2_double(D);
	}
	epoint_free(D);
}
