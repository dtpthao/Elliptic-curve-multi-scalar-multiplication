#include "option_Bin.h"

// 0||  1||  2||    3   ||  4||    5   ||    6   ||       7
// O|| P1|| P2|| P1 + P2|| P3|| P1 + P3|| P2 + P3|| P1 + P2 + P3
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

	//for (int i = 0; i < 8; i++) {
	//	std::cout << "pl[" << i << "] :\n"; cotnumEp(plist[i]);
	//}
}

void ShamirMul_Bin3_ptr(PL *shrBin, big k1, big k2, big k3,
	pepoint P1, pepoint P2, pepoint P3, pepoint R)
{
	int tmp, i, j = 0, index;
	DWORD w1, w2, w3;

	i = k3->len - 1;
	tmp = k3->w[i];
	while (tmp >> j && j != 32) j++;
	
	PreMul_Bin3(P1, P2, P3, shrBin->plist);

	/* *
	* ecurve2_padd doesn't work with point at infinity
	* therefore, R must be set with an initial value
	* which is not "point at infinity" before the loop
	* */
	index = (k1->w[i] >> --j) & 1;
	index += ((k2->w[i] >> j) & 1) << 1;
	index += ((k3->w[i] >> j) & 1) << 2;
	epoint2_copy(shrBin->plist[index], R);
	if (j == 0) {
		i--; j = 32;
	}
	for (--j; i >= 0; i--, j = 31) {
		w1 = k1->w[i];
		w2 = k2->w[i];
		w3 = k3->w[i];
		while (j) {
			index = (w1 >> j) & 1;
			index += ((w2 >> j) & 1) << 1;
			index += ((w3 >> j) & 1) << 2;
			ecurve2_double(R);
			if (index)
				ecurve2_padd(shrBin->plist[index], R);
			//cout << "R\n"; cotnumEp(R);
			j--;
		}
		index = (w1 & 1) + ((w2 & 1) << 1) + ((w3 & 1) << 2);
		ecurve2_double(R);
		if (index)
			ecurve2_padd(shrBin->plist[index], R);
		//cout << "R\n"; cotnumEp(R);
	}
}

void test_bin3(csprng &Rng, pepoint P, big n, string msg) 
{
	big k1 = mirvar(0),
		k2 = mirvar(0),
		k3 = mirvar(0),
		k = mirvar(0);
	pepoint P1 = epoint_init(),
		P2 = epoint_init(),
		P3 = epoint_init(),
		R = epoint_init(),
		R1 = epoint_init(),
		R2 = epoint_init();
	msg = "Test ShrMul_Bin\n";
	PL shrBin(8);
	//ecurve2_mult(a, P, Q);
	int count = 0, cmp = 0;
	for (int i = 0; i < 10000; i++) {
		//strong_bigdig(&Rng, 4, 16, k);
		//std::cout << "k: "; cotnum(k, stdout);
		strong_bigrand(&Rng, n, k);

		/*strong_bigdig(&Rng, 100, 10, a);
		strong_bigdig(&Rng, 100, 10, b);*/
		ShamirDecompose3(k, k1, k2, k3, P, P1, P2, P3);

		ShamirMul_Bin3_ptr(&shrBin, k1, k2, k3, P1, P2, P3, R);
		//ecurve2_mult2(a, P, b, Q, R1);
		ecurve2_mult(k, P, R2);
		//std::cout << "R: "; cotnumEp(R);
		//std::cout << "R1: "; cotnumEp(R1);
		//std::cout << "R2: "; cotnumEp(R2);
		cmp = epoint2_comp(R2, R);
		//if (!cmp) {
		//	std::cout << "k : "; cotnum(k, stdout);
		//	std::cout << "k1: "; cotnum(k1, stdout);
		//	std::cout << "k2: "; cotnum(k2, stdout);
		//	std::cout << "k3: "; cotnum(k3, stdout);
		//	std::cout << "R : "; cotnumEp(R);
		//	//std::cout << "R1: "; cotnumEp(R1);
		//	std::cout << "R2: "; cotnumEp(R2);
		//	//ShamirMul_Bin3_ptr(&shrBin, k1, k2, k3, P1, P2, P3, R);
		//	break;
		//}
		count += cmp;
	}
	std::cout << "Cmp: " << count << std::endl;
	shrBin.Destructor();
	mirkill(k1); mirkill(k2); mirkill(k3); mirkill(k);
	epoint_free(P1); epoint_free(P2); epoint_free(P3);
	epoint_free(R); epoint_free(R1); epoint_free(R2);
}