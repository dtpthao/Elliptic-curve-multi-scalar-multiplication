#include <iostream>
//struct IUnknown;// Workaround for "combaseapi.h(229): error C2187: syntax error: 'identifier' was unexpected here" when using /permissive-

extern "C" {
#include "miracl.h"
	FILE _iob[] = { *stdin, *stdout, *stderr };
	extern "C" FILE * __cdecl __iob_func(void) { return _iob; }
}
using namespace std;

//#include "Ellipse.h"
//#include "Test.h"
#include "PrintTest.h"

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
		compares(Rng, P, n, res[i]);
	}
	//cout << endl << "without PreMul_Bin_n norm" << endl;
	//cout << setw(19) << "" << "n=1           n=2           n=3           n=3           JSF3        lib(n=3)\n";
	cout << setw(19) << "" <<   "Bin1         JSF2         dJSF2          JSF4           JSF5           lib2\n";
	printcompares_bin(res, m);

	//cout << endl << "\a\a\a\a\a\a\a\a\a\a\a" << endl;

	epoint_free(P);
	mirkill(a); mirkill(b); mirkill(k);
	mirexit();
	system("pause");
	return 0;
}

int old_main()
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
	int m[NUM_OF_EC + 1] = { 0 };// { 163, 233, 283, 409, 571, 0 };
	
	Result testBin[NUM_OF_EC];
	string msg;
	int len = 0;
	readFile("DSTU4145TablePrameters.txt", EC);
	/*printf("%d\n\t%s\n%d %d %d %d %d\n\t%s\n\t%s\n\t%s\n\n",
		EC.a, EC.b, EC.m, EC.k1, EC.k2,
		EC.k3, EC.h, EC.n, EC.Gx, EC.Gy);*/
	//GetConstainsEC(EC, m[0]);
	//for (int i = 0; i < NUM_OF_EC; i++) {
	//	for (int i = 0; i < NUM_OF_EC; i++) {
	//		printf("%4d: %13.3f %13.3f %13.3f %13.3f\n",
	//			0, testBin[i].c[0], testBin[i].c[1], testBin[i].c[2], testBin[i].c[LIB]);
	//	}
	//}
	
	for (int i = 0; i < NUM_OF_EC; i++) {
		m[i] = EC[i].m;
		/*printf("%d\n%d\n\t%s\n%d %d %d %d %d\n%d\n\t%s\n%d\n\t%s\n%d\n\t%s\n\n",
			EC[i].a, len, EC[i].b, EC[i].m, EC[i].k1, EC[i].k2,
			EC[i].k3, EC[i].h, len, EC[i].n, len, EC[i].Gx, len, EC[i].Gy);*/
		//GetConstainsEC(EC, m[i]);
		if (!GenEC(EC[i], a, b, P, x, y, n)) 
			return 1;
		//TestBin(Rng, P, n, testBin[i]);
		//TestSclBin(Rng, P, n, testBin[i]);
		//TestJSF(Rng, P, n, testBin[i], msg);
		Test(Rng, P, n, testBin[i], msg);
		//TestShrMul_Bin(Rng, P, n, msg);
	}
	cout << endl << msg << endl;
	
	printBinOption(m, testBin);
	
	cout << endl << "\a\a\a\a\a\a\a\a\a\a\a" << endl;

	epoint_free(P);
	mirkill(a); mirkill(b); mirkill(k);
	mirexit();
	system("pause");
	return 0;
}

//__________________________________________________

