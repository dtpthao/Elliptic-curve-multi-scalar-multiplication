#include <iostream>
//struct IUnknown;// Workaround for "combaseapi.h(229): error C2187: syntax error: 'identifier' was unexpected here" when using /permissive-

extern "C" {
#include "miracl.h"
	FILE _iob[] = { *stdin, *stdout, *stderr };
	extern "C" FILE * __cdecl __iob_func(void) { return _iob; }
}
using namespace std;

//#include "Ellipse.h"
#include "Test.h"

int main()
{
	srand(time(NULL));
	miracl *M = mirsys(100, 0);
	M->IOBASE = 16;
	csprng Rng; InitStrongRNG(&Rng);
	big a = mirvar(1), x = mirvar(1), y = mirvar(1), n = mirvar(1);
	big b = mirvar(0xe);
	EC_CONSTANTS_F2m_POLY EC[NUM_OF_EC] = {};

	pepoint P = epoint_init();
	big k = mirvar(1);
	//int m[NUM_OF_EC + 1] = { 0 };// { 163, 233, 283, 409, 571, 0 };
	int m[NUM_OF_EC + 1] = { 163, 233, 283, 409, 571, 0 };

	for (int i = 0; i < LEN_GLOB_EPOINTS; i++) glob_epoints[i] = epoint_init();
	for (int i = 0; i < LEN_GLOB_BIGS; i++) glob_bigs[i] = mirvar(0);

	string msg;
	int len = 0;
	//readFile("DSTU4145TablePrameters.txt", EC);
	Result res[NUM_OF_EC + 1];
	for (int i = 0; i < NUM_OF_EC; i++) {
		//m[i] = EC[i].m;
		GetConstantsEC(EC[i], m[i]);
		if (!GenEC(EC[i], a, b, P, x, y, n))
			return 1;
		//cout << m[i] << "\t";
		//TestShrMul_Bin(Rng, P, n, msg);
		//test_bin3(Rng, P, n, msg);
		//test_bin_n(3, Rng, P, n, msg);
		//test_dJSF(3, Rng, P, n, msg);
		//compares(Rng, P, n, res[i]);
		compare_GendJSFs(Rng, P, n, res[i]);
		//compare_prepowmodJSFs(Rng, P, n, res[i]);
		//compare_Doub_Add_Neg(Rng, P, res[i], m[i]);
	}

	//cout << endl << "without Shamir decomposition" << endl;
	////cout << setw(19) << "" << "n=1           n=2           n=3           n=3           JSF3        lib(n=3)\n";
	//cout << setw(19) << "SclMul          " <<   "Bin1          JSF2          JSF3          JSF4          JSF5          JSF1\n";
	//printcompares(res, m);
	//printcompares_JSFs(res, m);
	printcompares_GendJSFs(res, m);
	//print_A_D_N(res, m);
	
	//cout << endl << "\a\a\a\a\a\a\a\a\a\a\a" << endl;

	for (int i = 0; i < LEN_GLOB_EPOINTS; i++) epoint_free(glob_epoints[i]);
	for (int i = 0; i < LEN_GLOB_BIGS; i++) mirkill(glob_bigs[i]);
	epoint_free(P);
	mirkill(a); mirkill(b); mirkill(k);
	mirexit();
	system("pause");
	return 0;
}

