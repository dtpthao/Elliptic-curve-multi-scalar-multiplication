#ifndef _OPTION_JSF_H
#define _OPTION_JSF_H

#include "Ellipse.h"
#include "Duration.h"

/* option_JSF.cpp */
DWORD GenJSF(big R, big S, char *JSFr, char *JSFs);
inline void PreMul_JSF(pepoint P, pepoint Q, pepoint *plist);
void ShamirMul_JSF(big a, pepoint P, big b, pepoint Q, pepoint R);

/* option_dJSF.cpp */
DWORD GendJSF(int d, big *r, char **dJSF);
void PreMul_dJSF(int d, int len, pepoint *P, pepoint *plist);
void ShamirMul_dJSF(int d, big *k, pepoint *P, pepoint R);
void test_dJSF(int d, csprng &Rng, pepoint P, big n, string msg);


void compare_GendJSFs(csprng &Rng, pepoint P, big n, Result &res);
void printcompares_GendJSFs(Result res[NUM_OF_EC + 1], int *m);

#endif
