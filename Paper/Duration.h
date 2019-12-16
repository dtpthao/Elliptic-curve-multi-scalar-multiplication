#ifndef _COUNT_TIME_H
#define _COUNT_TIME_H

#include "Ellipse.h"
#include "shr_decomp.h"

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


#endif
