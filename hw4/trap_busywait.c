/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Instructor Bryan Mills, PhD
 * Student: Qiming Guo
 * Implement Pthreads version of trapezoidal approximation.
 * Busy-Waiting version
 */

#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include "timer.h"

// Global variables to make coverting to pthreads easier :)
double a;
double b;
double h;
int n;
double approx;
int num_thread;
int flag;

// Actual areas under the f(x) = x^2 curves, for you to check your
// values against.
double static NEG_1_TO_POS_1 = 0.66666666666667;
double static ZERO_TO_POS_10 = 333.333;

// f function is defined a x^2
double f(double a) {
  return a * a;
}

// Serial implementation of trapezoidal approximation. You should
// refactor the loop in this function to be parallized using pthread.
void* trap(void* rank) {
    int my_rank = *((int*)rank);
    printf("My rank is %d\n", my_rank);
    double my_approx = 0.0;
    for (int i = my_rank+1; i < n-1; i += num_thread) {
        my_approx += f(a + i*h);
    }
    while (flag != my_rank);
    approx += my_approx;
    flag++;
}

void trap_busywait() {
  h = (b-a) / n;
  approx = ( f(a) + f(b) ) / 2.0;
  flag = 0;
  pthread_t ids[num_thread];
  int new_ids[num_thread];
  for (int i = 0; i < num_thread; ++i) {
      new_ids[i] = i;
      pthread_create(&ids[i], NULL, trap, &new_ids[i]);
  }
  for (int i = 0; i < num_thread; ++i) {
      pthread_join(ids[i], NULL);
  }
  approx = h*approx;
}

int main(int argc, char *argv[]) {
    if (argc == 1)
        num_thread = 1;
    else
        num_thread = atoi(argv[1]);
    printf ("Num of Thread is %d\n", num_thread);
  // Example 1 [-1,1]
  a = -1.0;
  b = 1.0;
  n = 1000000000;
  timerStart();
  trap_busywait();
  printf("Took %ld ms\n", timerStop());
  printf("a:%f\t b:%f\t n:%d\t actual:%f\t approximation:%f\n", a, b, n, NEG_1_TO_POS_1, approx);

  // Example 2 [0,10]
  a = 0.0;
  b = 10.0;
  n = 1000000000;
  timerStart();
  trap_busywait();
  printf("Took %ld ms\n", timerStop());
  printf("a:%f\t b:%f\t n:%d\t actual:%f\t approximation:%f\n", a, b, n, ZERO_TO_POS_10, approx);

  return 0;
}
