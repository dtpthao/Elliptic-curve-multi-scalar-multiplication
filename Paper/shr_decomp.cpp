#include "shr_decomp.h"
#include <iostream>

//void ShamirMul(big k, pepoint P, pepoint R, PL *opt,
//	void(*func) (PL *, big, pepoint, big, pepoint, pepoint))
//{
//	if (k->len == 0) return;
//	pepoint Q = epoint_init();
//	big a = mirvar(1);
//	big b = mirvar(1);
//	DWORD len, i = 31;
//	big tmp2l = mirvar(1);
//	while (!(k->w[k->len - 1] & (1 << i))) i--;
//	len = (k->len << 4) - ((31 - i) >> 1);
//	sftbit(tmp2l, len, tmp2l);
//	copy(k, a);
//	divide(a, tmp2l, b);
//	
//	epoint2_copy(P, Q);
//	for (i = 0; i < len; i++) {
//		ecurve2_add(Q, Q);
//	}
//	//ecurve2_mult(tmp2l, P, Q);
//	(*func) (opt, a, P, b, Q, R);
//
//	epoint_free(Q);
//	mirkill(a); mirkill(b); mirkill(tmp2l);
//}

void ShamirMul(big k, pepoint P, pepoint R, PL *opt,
	void(*func) (PL *, big, pepoint, big, pepoint, pepoint))
{
	DWORD len, i = 31, rlen;
	pepoint Q = epoint_init();
	big a = mirvar(1);
	big b = mirvar(1);
	big tmp2l = mirvar(1);

	len = k->w[k->len - 1];
	while (!(len & (1 << i))) i--;
	len = (k->len << 4) - ((31 - i) >> 1) - !(i & 1);

	sftbit(k, 0 - len, b);
	a->len = len >> 5;
	for (i = 0; i < a->len; i++) a->w[i] = k->w[i];
	rlen = len & 0x1f;
	if (rlen) {
		a->len++;
		a->w[i] = k->w[i] & ((1 << rlen) - 1);
	}
	epoint2_copy(P, Q);
	for (i = 0; i < len; i++) ecurve2_double(Q);
	epoint2_norm(Q);
	(*func) (opt, a, P, b, Q, R);
	epoint_free(Q);
	mirkill(a); mirkill(b); mirkill(tmp2l);
}

void ShamirDecompose(big k, pepoint P, big a, pepoint Q, big b) 
{
	DWORD len, i = 31, rlen;

	len = k->w[k->len - 1];
	while (!(len & (1 << i))) i--;
	len = (k->len << 4) - ((31 - i) >> 1) - !(i & 1);

	sftbit(k, 0-len, b);
	a->len = len >> 5;
	for (i = 0; i < a->len; i++) a->w[i] = k->w[i];
	rlen = len & 0x1f;
	if (rlen) {
		a->len++;
		a->w[i] = k->w[i] & ((1 << rlen) - 1);
	}
	epoint2_copy(P, Q);
	for (i = 0; i < len; i++) ecurve2_double(Q);
	epoint2_norm(Q);
}

void ShamirDecompose(big k, big a, big b) 
{
	if (k->len == 0) return;
	DWORD len, i = 31;
	big tmp2l = mirvar(1);
	while (!(k->w[k->len - 1] & (1 << i))) i--;
	len = (k->len << 4) - ((31 - i) >> 1) - !(i & 1);
	sftbit(tmp2l, len, tmp2l);
	copy(k, a);
	divide(a, tmp2l, b);
	mirkill(tmp2l);
}

void ShamirDecompose3(big k, big k1, big k2, big k3, 
	pepoint P, pepoint P1, pepoint P2, pepoint P3)
{
	if (k->len == 0) return;
	DWORD len, i = 31;
	big tmp = mirvar(1);
	len = k->w[k->len - 1];
	while (!(len & (1 << i))) i--;
	len = (k->len << 5) - 31 + i;
	len /= 3;
	sftbit(tmp, len, tmp);
	copy(k, k1);
	divide(k1, tmp, k2);
	divide(k2, tmp, k3);

	epoint2_copy(P, P1);
	epoint2_copy(P, P2);
	for (i = 0; i < len; i++) {
		ecurve2_double(P2);
	}
	epoint2_copy(P2, P3);
	for (i = 0; i < len; i++) {
		ecurve2_double(P3);
	}

	mirkill(tmp);
}
