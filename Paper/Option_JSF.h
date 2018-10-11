#ifndef _OPTION_JSF_H
#define _OPTION_JSF_H

#include "Ellipse.h"

//inline void subGenJSF(big Var1, big Var2, char *JSF, DWORD &lenJSF, bool &d, big RS);

DWORD GenJSF(big R, big S, char *JSFr, char *JSFs);

inline void PreMul_JSF(pepoint P, pepoint Q, PL *opt);

void ShamirMul_JSF(PL *opt, big a, pepoint P, big b, pepoint Q, pepoint R);

#endif
