#include "option_Bin.h"

inline void PreMul_Bin(pepoint P, pepoint Q, pepoint *plist)
{
	epoint2_copy(P, plist[1]);
	epoint2_copy(Q, plist[2]);
	epoint2_copy(P, plist[3]);
	ecurve2_padd(Q, plist[3]);  //checked plist[3] != normalize
	epoint2_norm(plist[3]);
}

void ShamirMul_Bin_ptr(PL *shrBin, big a,
	pepoint P, big b, pepoint Q, pepoint R)					//in use
{
	int bita, bitb, i, j = 0, index;
	DWORD shift = 1, lastw, a1, b1;

	i = b->len - 1;
	lastw = b->w[i];
	while (lastw >> j && j != 32) j++;

	PreMul_Bin(P, Q, shrBin->plist);
	epoint_set(0, 0, 1, R);

	/* *
	 * ecurve2_padd doesn't work with point at infinity
	 * therefore, R must be set with an initial value 
	 * which is not "point at infinity" before the loop
	 * */
	if (j != 1) {
		bita = (a->w[i] >> --j) & 1;
		bitb = (b->w[i] >> (j - 1)) & 2;	/** !!! j can be 0 ==>> fixed **/
		index = bita + bitb;
	}
	else {
		index = (a->w[i] & 1) + ((b->w[i] & 1) << 1);
		i--; j = 32;
	}
	epoint2_copy(shrBin->plist[index], R);
	for (--j; i >= 0; i--, j = 31) {
		a1 = a->w[i];
		b1 = b->w[i];
		while (j) {
			bita = (a1 >> j) & 1;
			bitb = (b1 >> (j - 1)) & 2;
			ecurve2_double(R);
			index = bitb + bita;
			if (index)
				ecurve2_padd(shrBin->plist[index], R);
			j--;
		}
		index = (a1 & 1) + ((b1 & 1) << 1);
		ecurve2_double(R);
		if (index) ecurve2_padd(shrBin->plist[index], R);
	}
}

void ShamirMul_Bin(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	int bita, bitb, i, j = 0, index;
	DWORD shift = 1, lastw, a1, b1;
	PL globalShrBin(4);
	i = b->len - 1;
	lastw = b->w[i];
	while (lastw >> j && j != 31) j++;

	PreMul_Bin(P, Q, globalShrBin.plist);
	epoint_set(0, 0, 1, R);
	bita = (a->w[i] >> --j) & 1;
	bitb = (b->w[i] >> (j - 1)) & 2;	/** !!!!!!! j can be 0 **/
	epoint2_copy(globalShrBin.plist[bitb + bita], R);
	for (--j; i >= 0; i--, j = 31) {
		a1 = a->w[i];
		b1 = b->w[i];
		while (j) {
			bita = (a1 >> j) & 1;
			bitb = (b1 >> (j - 1)) & 2;
			ecurve2_double(R);
			index = bitb + bita;
			if (index)
				ecurve2_padd(globalShrBin.plist[index], R);
			//ecurve2_add(R, R);
		//ecurve2_add(globalShrBin.plist[bitb + bita], R);
			j--;
		}
		index = (a1 & 1) + ((b1 & 1) << 1);
		ecurve2_double(R);
		if (index) ecurve2_padd(globalShrBin.plist[index], R);
		//ecurve2_add(globalShrBin.plist[index], R);
	}
	globalShrBin.Destructor();
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
	for (int i = 0; i < 1000; i++) {
		//strong_bigdig(&Rng, 16, 16, k);
		//std::cout << "k: "; cotnum(k, stdout);
		strong_bigrand(&Rng, n, k);

		/*strong_bigdig(&Rng, 100, 10, a);
		strong_bigdig(&Rng, 100, 10, b);*/
		ShamirDecompose(k, P, a, Q, b);

		ShamirMul_Bin_ptr(&shrBin, a, P, b, Q, R);
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
			ShamirMul_Bin_ptr(&shrBin, a, P, b, Q, R);
			break;
		}
		count += cmp;
	}
	std::cout << "Cmp: " << count << std::endl;
	mirkill(a); mirkill(b); mirkill(k);
	epoint_free(Q); epoint_free(R);
	epoint_free(R1); epoint_free(R2);
}

