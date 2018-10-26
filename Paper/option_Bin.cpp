#include "option_Bin.h"

void ScalarMul_Bin_L2R(big k, pepoint P, pepoint R)
{
	unsigned int len, i = 31;
	epoint_set(0, 0, 1, R);

	len = k->w[k->len - 1];
	while (!(len & (1 << i))) i--;
	len = (k->len << 5) - (31 - i);
	epoint2_copy(P, R);
	for (int j = len - 2; j >= 0; j--) {
		ecurve2_double(R);
		if (k->w[j >> 5] & (1 << j)) {
			ecurve2_padd(P, R);
		}
	}
}

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
		ShamirDecomposit(k, P, a, Q, b);

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


/*************************************************************************************************/
/**______________________________The Very First Option__________________________________________**/
/**                                                                                             **/
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

void PreMul_Bin_old(pepoint P, pepoint Q, pepoint *plist)
{
	for (int i = 0; i < 4; i++) {
		plist[i] = epoint_init();
	}
	epoint2_set(0, 0, 1, plist[0]);
	epoint2_copy(P, plist[1]);
	epoint2_copy(Q, plist[2]);
	epoint2_copy(P, plist[3]);
	ecurve2_add(Q, plist[3]);
}

void ShamirMul_Bin_old(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	pepoint plist[4];
	int bita, bitb, i, j = 0;
	DWORD shift = 1, lastw;
	lastw = (compare(a, b) >= 0) ? 
		a->w[a->len - 1] : b->w[b->len - 1];
	while (lastw >> j && j != 31) {
		j++;
	}
	PreMul_Bin_old(P, Q, plist);
	epoint_set(0, 0, 1, R);
	i = (lastw == a->w[a->len - 1]) ? a->len - 1 : b->len - 1;
	for (shift <<= j; i >= 0; i--, shift = 0x80000000, j = 31) {
		while (shift) {
			bita = (a->w[i] & shift) >> j;
			bitb = (b->w[i] & shift) >> j;
			ecurve2_add(R, R);
			ecurve2_add(plist[(bitb << 1) + bita], R);
			shift >>= 1;
			j--;
		}
	}
	for (int i = 0; i < 4; i++) epoint_free(plist[i]);
}

/**_____________________________End "The very First Option"_____________________________________**/
/*************************************************************************************************/


//WORD GetBin2(big k, char*Bit)
//{
//	DWORD i = 0;
//	big d = mirvar(1);
//	copy(k, d);
//	while (d->len != 1 || d->w[0] != 1) {
//		if (d->w[0] & 1) Bit[i] = 1;
//		else Bit[i] = 0;
//		sftbit(d, -1, d);
//		i++;
//	}
//	Bit[i] = 1;
//	mirkill(d);
//	return ++i;
//}


//void ShamirMul_Bin(big a, pepoint P, big b, pepoint Q, pepoint R, DWORD len)
//{
//	pepoint plist[4];
//	//char Bita[(MAX_M + 1) / 2] = { 0 };
//	//char Bitb[(MAX_M + 1) / 2] = { 0 };
//	PreMul_Bin(P, Q, plist);
//	epoint_set(0, 0, 1, R);
//	//DWORD lenBita, lenBitb;
//	//lenBita = GetBin(a, Bita);
//	//lenBitb = GetBin(b, Bitb);
//	int bita, bitb, lena = 0, lenb = 0;
//	int i, j = 0; len = 0;// (a->len - 1) << 5;
//	if (a->len == b->len) {
//		i = max(a->w[a->len - 1], b->w[b->len - 1]);
//	}
//	else if (a->len < b->len) i = b->w[b->len - 1];
//	else i = a->w[a->len - 1];
//	//i = a->w[a->len - 1];
//
//	while (i) {
//		len++; i >>= 1;
//	}
//	/*i = b->w[b->len - 1];
//	while (i) {
//	lenb++; i >>= 1;
//	}
//	len = max(lena, lenb);*/
//	//while (j < len) {
//	//	i = j / 32;
//	//	bita = (a->w[i] & (1 << j)) ? 1 : 0;
//	//	bitb = (b->w[i] & (1 << j)) ? 1 : 0;
//	//	ecurve2_add(R, R);
//	//	ecurve2_add(plist[(bitb << 1) + bita], R);
//	//	j++;
//	//	//if (j / 32) { i++; j = 0; }
//	//}
//	DWORD shift = 1 << (len - 1);
//	for (i = a->len - 1; i >= 0; i--) {
//		while (shift) {
//			bita = (a->w[i] & shift) ? 1 : 0;
//			bitb = (b->w[i] & shift) ? 1 : 0;
//			ecurve2_add(R, R);
//			ecurve2_add(plist[(bitb << 1) + bita], R);
//			shift >>= 1;
//		}
//		shift = 0x80000000;
//	}
//	/*for (int i = lenBita - 1; i >= 0; i--) {
//	ecurve2_add(R, R);
//	ecurve2_add(plist[(Bitb[i] << 1) + Bita[i]], R);
//	}*/
//	for (int i = 0; i < 4; i++) epoint_free(plist[i]);
//}



//DWORD GetBin(big k, char*Bit)
//{
//	DWORD lenbit = 0;
//	DWORD shift = 1;
//
//	for (int i = 0; i < k->len - 1; i++, shift = 1) {
//		while (shift) {
//			Bit[lenbit++] = (k->w[i] & shift) ? 1 : 0;
//			shift <<= 1;
//		}
//	}
//	shift = k->w[k->len - 1];
//	while (shift) {
//		Bit[lenbit++] = (shift & 1) ? 1 : 0;
//		shift >>= 1;
//	}
//	return lenbit;
//}
