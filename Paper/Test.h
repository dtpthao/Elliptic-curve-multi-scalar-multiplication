#ifndef _TEST_H
#define _TEST_H

#include <math.h>
#include <iomanip>
#include <string>

#include "option_Bin.h"
#include "option_JSF.h"
#include "shr_decomp.h"


#define LIB		5			//16

struct Result {
	double t[LIB + 1] = { 0 };
	double p[LIB + 1] = { 0 };
	unsigned int c[LIB + 1] = { 0 };
};

void compares(csprng &Rng, pepoint P, big n, Result &res);
void printcompares(Result res[NUM_OF_EC + 1], int *m);

void compare_prepowmodJSFs(csprng &Rng, pepoint P, big n, Result &res);
void printcompares_JSFs(Result res[NUM_OF_EC + 1], int *m);

void compare_Doub_Add_Neg(csprng &Rng, pepoint P, Result &res, int m);
void print_A_D_N(Result res[NUM_OF_EC + 1], int *m);

#endif
