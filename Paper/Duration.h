#ifndef _COUNT_TIME_H
#define _COUNT_TIME_H

#include <Windows.h>
#include <time.h>
#include "Ellipse.h"
#include "ShamirMul.h"

typedef epoint* pepoint;

typedef struct {
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
} stopWatch;

inline void startTimer(stopWatch *timer);
inline void stopTimer(stopWatch *timer);
inline double LIToSecs(LARGE_INTEGER *L);
inline double getElapsedTime(stopWatch *timer);	//get time in secs/msecs
inline LONGLONG getTickCount(stopWatch *timer);	//get time in tick


void SclDuration(big k, pepoint P, pepoint R,
	void(*func) (big, pepoint, pepoint), double &min);

void ShrDuration(big k, pepoint P, pepoint R, PL *opt,
	void(*func) (PL *, big, pepoint, big, pepoint, pepoint), double &t);

void ShrDuration(big k, pepoint P, pepoint R,
	void(*func) (big, pepoint, big, pepoint, pepoint), double &min);

void Shr2Duration(big a, pepoint P, big b, pepoint Q, pepoint R, 
	void(*func) (big, pepoint, big, pepoint, pepoint), double &min);

void ShrDuration(big k, pepoint P, pepoint R, PL *opt, double &t);

void ShrDuration(PL *opt, big a, pepoint P, big b, pepoint Q, pepoint R, double &t);

//Draft
/******************************************************************************/
void ShrDecompDuration(big k, pepoint P, big a, pepoint Q, big b, double &t);

/******************************************************************************/

//
//void GenFormSclDuration(big k, DWORD (*func) (big, char*), double &min, DWORD&);
//
//void GenFormShrDuration(big a, big b, DWORD(*func) (big, char*), double &min, DWORD&);

#endif
