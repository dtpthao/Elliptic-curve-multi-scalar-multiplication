#include <iostream>

//struct IUnknown;// Workaround for "combaseapi.h(229): error C2187: syntax error: 'identifier' was unexpected here" when using /permissive-

extern "C" {
#include "miracl.h"
	FILE _iob[] = { *stdin, *stdout, *stderr };
	extern "C" FILE * __cdecl __iob_func(void) { return _iob; }
}
using namespace std;

#include "Ellipse.h"
#include "Test.h"
#include "PrintTest.h"

int main()
{
	srand(time(NULL));
	miracl *M = mirsys(100, 0);
	M->IOBASE = 16;
	csprng Rng; InitStrongRNG(&Rng);
	big a = mirvar(1), x = mirvar(1), y = mirvar(1), n = mirvar(1);
	big b = mirvar(0xe);
	EC_CONSTANTS_F2m_POLY EC = {};

	pepoint P = epoint_init();
	big k = mirvar(1);
	int m[6] = { 163,233,283,409,571,0 };
	
	Result testBin[6];
	string msg;

	for (int i = 0; i < 5; i++) {
		GetConstainsEC(EC, m[i]);
		if (!GenEC(EC, a, b, P, x, y, n)) return 1;
		//TestBin(Rng, P, n, testBin[i]);
		//TestSclBin(Rng, P, n, testBin[i]);
		TestJSF(Rng, P, n, testBin[i], msg);
	}
	printBinOption(m, testBin);
	cout << endl << endl << msg;
	cout << endl << "\a\a\a\a\a\a\a\a\a\a\a" << endl;
	cout << "\a" << endl;

	epoint_free(P);
	mirkill(a); mirkill(b); mirkill(k);
	mirexit();
	system("pause");
	return 0;
}

//__________________________________________________

