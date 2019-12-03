#ifndef _OPTION_JSF_H
#define _OPTION_JSF_H

#include "Ellipse.h"
#include "Duration.h"


/**************************************************************************************/
/**________________________________In Use____________________________________________**/
/**                                                                                  **/
/* option_JSF.cpp */
DWORD GenJSF(big R, big S, char *JSFr, char *JSFs);
inline void PreMul_JSF(pepoint P, pepoint Q, pepoint *plist);
void ShamirMul_JSF(big a, pepoint P, big b, pepoint Q, pepoint R);

/* option_dJSF.cpp */
void ShamirMul_dJSF(int d, big *k, pepoint *P, pepoint R);
void test_dJSF(int d, csprng &Rng, pepoint P, big n, string msg);
/**______________________________End In Use__________________________________________**/
/**************************************************************************************/



/**************************************************************************************/
/**____________________________Another Option________________________________________**/
/**                                                                                  **/
inline void PreMul_JSF2(pepoint P, pepoint Q, PL *opt);
void ShamirMul_JSF2(PL *opt, big a, pepoint P, big b, pepoint Q, pepoint R);
/**__________________________End "Another Option"____________________________________**/
/**************************************************************************************/



/**************************************************************************************/
/**_________________________The Very First Option____________________________________**/
/**                                                                                  **/
inline void subGenJSF_old(big Var1, big Var2, char *JSF, DWORD &lenJSF, bool &d, big RS);
DWORD GenJSF_old(big R, big S, char *JSFr, char *JSFs);
inline void PreMul_JSF_old(pepoint P, pepoint Q, pepoint *plist);
void ShamirMul_JSF_old(big a, pepoint P, big b, pepoint Q, pepoint R);

/**______________________End "The very First Option"_________________________________**/
/**************************************************************************************/

#endif
