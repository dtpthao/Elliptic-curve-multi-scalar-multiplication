#ifndef _ELLIPSE_H
#define _ELLIPSE_H

#include <iostream>
#include <Windows.h>
#include <time.h>
#include <stdio.h>
extern "C" {
#include"miracl.h"
}

using namespace std;

#define NUM_OF_EC	5 //10

#define EC163 0
#define EC233 1
#define EC283 2
#define EC409 3
#define EC571 4

#define MAX_M 571
#define LEN_HEX_F2m (((MAX_M+3)/4)+1)	// 143+1
#define LEN_DWORD_F2m (((MAX_M+31)/32)*3+1)	// 18∙3+1=55

typedef struct {	//	E: y2 + x∙y = (x3 + a∙x2 + b) mod f (t)  (#E — порядок ЭК)
	DWORD m, k3, k2, k1,	// f (t) = 
		a,	// коэффициент ЭК (a  GF(2m), a = 1)
		h;	// кофактор (#E = n∙h) (h = 2)
	CHAR b[LEN_HEX_F2m],	// коэффициент ЭК b  GF(2m))
		n[LEN_HEX_F2m],	// порядок точки G (#E = n∙h) (простое число)
		Gx[LEN_HEX_F2m], Gy[LEN_HEX_F2m];	// базовая точка G = (xG , yG)
} EC_CONSTANTS_F2m_POLY;//EC_163, EC_233, EC_283, EC_409, EC_571;

typedef epoint* pepoint;

typedef struct PointList {
	int len;
	pepoint* plist;
	PointList(int len);
	void Destructor();
} PL;

#define LEN_GLOB_EPOINTS 300
//extern PointList glob_epoints;
extern pepoint glob_epoints[LEN_GLOB_EPOINTS];

void InitStrongRNG(csprng *Rng);

int GenEC(EC_CONSTANTS_F2m_POLY, big, big, pepoint, big, big, big);		//Initialize EC

void cotnumEp(pepoint P);	//display P.x, P.y

void GetConstantsEC(EC_CONSTANTS_F2m_POLY &, unsigned int);		//get const Input params EC

void readFile(const char *, EC_CONSTANTS_F2m_POLY[NUM_OF_EC]);	//get list of input params EC 

void readFile(const char *, EC_CONSTANTS_F2m_POLY &);

// doesn't work with point at infinity
BOOL ecurve2_padd(_MIPD_ epoint *p, epoint *pa); // internal add function of ms32.lib

#endif
