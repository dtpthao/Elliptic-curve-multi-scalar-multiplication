#include "option_Bin.h"


inline void PreMul_Bin_n(int n, pepoint P, pepoint *Px, pepoint *plist)
{
	epoint2_copy(Px[1], plist[1]);

	int idx;
	for (int i = 1; i < n; i++) {
		idx = 1 << i;
		epoint_copy(Px[i], plist[idx]);
		for (int j = 1; j < idx; j++) {
			epoint_copy(plist[idx], plist[idx + j]);
			ecurve2_padd(plist[j], plist[idx + j]);
			epoint2_norm(plist[idx + j]);
		}
	}
}