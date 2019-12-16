#include "option_Bin.h"

inline void PreMul_Bin(pepoint P, pepoint Q, pepoint *plist)
{
	epoint2_copy(P, plist[1]);
	epoint2_copy(Q, plist[2]);
	epoint2_copy(P, plist[3]);
	ecurve2_padd(Q, plist[3]);  //checked plist[3] != normalize
	epoint2_norm(plist[3]);
}

void ShamirMul_Bin(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	int i, j = 0, index;
	DWORD shift = 1, lastw, a1, b1;

	i = b->len - 1;
	lastw = b->w[i];
	while (lastw >> j && j != 32) j++;

	PreMul_Bin(P, Q, glob_epoints);

	/* *
	 * ecurve2_padd doesn't work with point at infinity
	 * therefore, R must be set with an initial value
	 * which is not "point at infinity" before the loop
	 * */
	index = (a->w[i] >> --j) & 1;
	index += ((b->w[i] >> j) & 1) << 1;
	if (j == 0) {
		i--; j = 32;
	}
	epoint2_copy(glob_epoints[index], R);
	for (--j; i >= 0; i--, j = 31) {
		a1 = a->w[i];
		b1 = b->w[i];
		while (j) {
			index = (a1 >> j) & 1;
			index += (b1 >> (j - 1)) & 2;
			ecurve2_double(R);
			if (index) ecurve2_padd(glob_epoints[index], R);
			j--;
		}
		index = (a1 & 1) + ((b1 & 1) << 1);
		ecurve2_double(R);
		if (index) ecurve2_padd(glob_epoints[index], R);
	}
}

void TestShrMul_Bin(csprng &Rng, pepoint P, big n, std::string &msg) {
	big a = mirvar(0),
		b = mirvar(0),
		k = mirvar(0);
	pepoint Q = epoint_init(),
		R = epoint_init(),
		R1 = epoint_init(),
		R2 = epoint_init();
	msg = "Test ShrMul_Bin\n";
	PL shrBin(4);
	//ecurve2_mult(a, P, Q);
	int count = 0, cmp = 0;
	for (int i = 0; i < 10000; i++) {
		//strong_bigdig(&Rng, 16, 16, k);
		//std::cout << "k: "; cotnum(k, stdout);
		strong_bigrand(&Rng, n, k);

		/*strong_bigdig(&Rng, 100, 10, a);
		strong_bigdig(&Rng, 100, 10, b);*/
		ShamirDecompose(k, P, a, Q, b);

		ShamirMul_Bin(a, P, b, Q, R);
		ecurve2_mult2(a, P, b, Q, R1);
		//ecurve2_mult(k, P, R2);
		//std::cout << "R: "; cotnumEp(R);
		//std::cout << "R1: "; cotnumEp(R1);
		//std::cout << "R2: "; cotnumEp(R2);
		cmp = epoint2_comp(R1, R);
		if (!cmp) {
			std::cout << "k: "; cotnum(k, stdout);
			std::cout << "a: "; cotnum(a, stdout);
			std::cout << "b: "; cotnum(b, stdout);
			std::cout << "R: "; cotnumEp(R);
			std::cout << "R1: "; cotnumEp(R1);
			//std::cout << "R2: "; cotnumEp(R2);
			ShamirMul_Bin(a, P, b, Q, R);
			break;
		}
		count += cmp;
	}
	std::cout << "Cmp: " << count << std::endl;
	mirkill(a); mirkill(b); mirkill(k);
	epoint_free(Q); epoint_free(R);
	epoint_free(R1); epoint_free(R2);
}

