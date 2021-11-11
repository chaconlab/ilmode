#include "libnma_time.h"

#if defined  _WIN32 || defined _WIN64
// Windows version of the timer
void timerReal::startTimer()
{
    QueryPerformanceCounter(&start);
}

void timerReal::stopTimer()
{
    QueryPerformanceCounter(&stop);
}

double timerReal::getElapsedTime()
{
	LARGE_INTEGER time;
	time.QuadPart = stop.QuadPart - start.QuadPart;
    return ((double)time.QuadPart /(double)frequency.QuadPart);
}





#else
// Linux version of the timer
void timerReal::startTimer()
{
    clock_gettime(CLOCK_REALTIME,&start);
}

void timerReal::stopTimer()
{
    clock_gettime(CLOCK_REALTIME,&stop);
}

double timerReal::getElapsedTime()
{
	double elapsed;
	elapsed = (double)(stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
	return elapsed;
}
#endif
