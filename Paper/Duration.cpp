#include "Duration.h"
#include <iostream>
#include "option_Bin.h"

inline void startTimer(stopWatch *timer) {
	QueryPerformanceCounter(&timer->start);
}

inline void stopTimer(stopWatch *timer) {
	QueryPerformanceCounter(&timer->stop);
}

inline double LIToSecs(LARGE_INTEGER *L) {
	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	return ((double)L->QuadPart / (double)frequency.QuadPart);
}

inline double getElapsedTime(stopWatch *timer) {
	LARGE_INTEGER time;
	time.QuadPart = timer->stop.QuadPart - timer->start.QuadPart;
	return LIToSecs(&time);
}

inline LONGLONG getTickCount(stopWatch *timer) {
	return timer->stop.QuadPart - timer->start.QuadPart;
}


#define REPEAT 10
/*
 * Counting execution time of simple scalar multiplication
 * R = kP
 */
void SclDuration(big k, pepoint P, pepoint R, 
	void(*func) (big, pepoint, pepoint), double &t)
{
	LONGLONG dur, min = LONG_MAX;
	stopWatch timer;

	for (int i = 0; i < REPEAT; i++) {
		startTimer(&timer);
		(*func)(k, P, R);
		stopTimer(&timer);

		dur = getTickCount(&timer);
		min = (min < dur) ? min : dur;
	}
	t += min;
}

/*
 * Counting execution time of single scalar multiplication via Shamir method
 * USING POINTER parameter (PL *opt)
 * R = kP = aP + bQ
 */
void ShrDuration(big k, pepoint P, pepoint R, PL *opt,
	void(*func) (PL *, big, pepoint, big, pepoint, pepoint), double &t)
{
	LONGLONG dur, min = LONG_MAX;
	big a = mirvar(1),
		b = mirvar(1);
	pepoint Q = epoint_init();
	stopWatch timer;
	for (int i = 0; i < REPEAT; i++) {
		startTimer(&timer);
		/*ShamirDecomposit(k, P, a, Q, b);
		(*func) (opt, a, P, b, Q, R);*/
		ShamirMul(k, P, R, opt, (*func));
		stopTimer(&timer);

		dur = getTickCount(&timer);
		min = (min < dur) ? min : dur;
	}
	t += min;
	epoint_free(Q);
	mirkill(a); mirkill(b);
}

/*
 * Counting execution time of single scalar multiplication via Shamir method
 * R = kP = aP + bQ
 */
void ShrDuration(big k, pepoint P, pepoint R,
	void(*func) (big, pepoint, big, pepoint, pepoint), double &t)
{
	LONGLONG dur, min = LONG_MAX;
	big a = mirvar(1),
		b = mirvar(1);
	pepoint Q = epoint_init();
	stopWatch timer;
	for (int i = 0; i < REPEAT; i++) {
		startTimer(&timer);
		ShamirDecomposit(k, P, a, Q, b);
		(*func) (a, P, b, Q, R);
		stopTimer(&timer);

		dur = getTickCount(&timer);
		min = (min < dur) ? min : dur;
	}
	t += min;
	epoint_free(Q);
	mirkill(a); mirkill(b);
}

/*
 * Counting execution time of double scalar multiplication
 */
void Shr2Duration(big a, pepoint P, big b, pepoint Q, pepoint R, 
	void(*func) (big, pepoint, big, pepoint, pepoint), double &t)
{
	LONGLONG dur, min = LONG_MAX;
	stopWatch timer;
	
	for (int i = 0; i < REPEAT; i++) {
		startTimer(&timer);
		(*func) (a, P, b, Q, R);
		stopTimer(&timer);
		dur = getTickCount(&timer);
		min = (min < dur) ? min : dur;
	}
	t += min;
}

/*
 * Counting execution time of shamir decomposition
 */
void ShrDecompDuration(big k, pepoint P, big a, pepoint Q, big b, double &t)
{
	LONGLONG dur, min = LONG_MAX;
	stopWatch timer;

	for (int i = 0; i < REPEAT; i++) {
		startTimer(&timer);
		ShamirDecomposit(k, P, a, Q, b);
		stopTimer(&timer);
		dur = getTickCount(&timer);
		min = (min < dur) ? min : dur;
	}
	t += min;
}


//void GenDuration(big k, char* r, int code,
//	void(*func) (big, char*)) {
//	double dur;
//	stopWatch timer;
//	startTimer(&timer);
//	(*func)(k, r);
//	stopTimer(&timer);
//	dur = getElapsedTime(&timer) * 1000;
//};

