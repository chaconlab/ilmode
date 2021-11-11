#pragma once

#if defined  _WIN32 || defined _WIN64
// Windows version of timer
#include <windows.h>
class timerReal
{
private:
	LARGE_INTEGER start;
    LARGE_INTEGER stop;
	LARGE_INTEGER frequency;
public:
	void startTimer();
	void stopTimer();
	double getElapsedTime(); // Elapsed time in seconds
};
#else
// Linux version of timer
#include <time.h>
class timerReal
{
private:
	struct timespec start;
	struct timespec stop;
public:
	void startTimer();
	void stopTimer();
	double getElapsedTime(); // Elapsed time in seconds
};
#endif
