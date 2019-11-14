#include "option_Bin.h"

// 0,  1,  2,    3   ,  4,    5   ,    6   ,       7
// O, P1, P2, P1 + P2, P3, P1 + P3, P2 + P3, P1 + P2 + P3
inline void PreMul_Bin3(pepoint P1, pepoint P2, pepoint P3, pepoint *plist)
{
	epoint2_copy(P1, plist[1]);

	epoint2_copy(P2, plist[2]);

	epoint2_copy(P1, plist[3]);
	ecurve2_padd(P2, plist[3]);  //checked plist[3] is not normalized
	epoint2_norm(plist[3]);

	epoint2_copy(P3, plist[4]);
	
	epoint2_copy(P1, plist[5]);
	ecurve2_padd(P3, plist[5]);
	epoint2_norm(plist[5]);

	epoint2_copy(P2, plist[6]);
	ecurve2_padd(P3, plist[6]);
	epoint2_norm(plist[6]);

	epoint2_copy(plist[6], plist[7]);
	ecurve2_padd(P1, plist[7]);
	epoint2_norm(plist[7]);
}

void ShamirMul_Bin3_ptr(PL *shrBin, big k1, big k2, big k3,
	pepoint P1, pepoint P2, pepoint P3, pepoint R)
{
	int tmp1, tmp2, tmp3, i, j = 0, index;
	DWORD shift = 1, w1, w2, w3;

	i = k3->len - 1;
	tmp1 = k3->w[i];
	while (tmp1 >> j && j != 32) j++;

	PreMul_Bin3(P1, P2, P3, shrBin->plist);
	epoint_set(0, 0, 1, R);

	for (--j; i >= 0; i--, j = 31) {
		w1 = k1->w[i];
		w2 = k2->w[i];
		w3 = k3->w[i];
		while (j) {
			tmp1 = (w1 >> j) & 1;
			tmp2 = ((w2 >> j) & 1) << 1;
			tmp3 = ((w3 >> j) & 1) << 2;
			index = tmp1 + tmp2 + tmp3;
			ecurve2_double(R);
			if (index)
				ecurve2_padd(shrBin->plist[index], R);
			j--;
		}
		index = (w1 & 1) + ((w2 & 1) << 1) + ((w3 & 1) << 2);
		ecurve2_double(R);
		if (index)
			ecurve2_padd(shrBin->plist[index], R);
	}
	//if (j != 1) {
	//	bit1 = (a->w[i] >> --j) & 1;
	//	bit2 = (b->w[i] >> (j - 1)) & 2;	/** !!! j can be 0 ==>> fixed **/
	//	index = bita + bitb;
	//}
	//else {
	//	index = (a->w[i] & 1) + ((b->w[i] & 1) << 1);
	//	i--; j = 32;
	//}
	//epoint2_copy(shrBin->plist[index], R);
	//for (--j; i >= 0; i--, j = 31) {
	//	a1 = a->w[i];
	//	b1 = b->w[i];
	//	while (j) {
	//		bita = (a1 >> j) & 1;
	//		bitb = (b1 >> (j - 1)) & 2;
	//		ecurve2_double(R);
	//		index = bitb + bita;
	//		if (index)
	//			ecurve2_padd(shrBin->plist[index], R);
	//		j--;
	//	}
	//	index = (a1 & 1) + ((b1 & 1) << 1);
	//	ecurve2_double(R);
	//	if (index) ecurve2_padd(shrBin->plist[index], R);
	//}
}
