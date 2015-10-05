/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Instructor Bryan Mills, PhD
 * Timing operations.
 */

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#define NANOS 1000000000LL
struct timespec start, diff;
// Starts timer and resets the elapsed time
void timerStart(){
    int rvalue = clock_gettime(CLOCK_REALTIME, &start);
}

// Stops the timer and returns elapsed time in msec
long timerStop(){
    int rvalue = clock_gettime(CLOCK_REALTIME, &diff);
    return (diff.tv_sec*NANOS + diff.tv_nsec - start.tv_sec*NANOS - start.tv_nsec) / 1000000.0;
}
