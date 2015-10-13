/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Instructor Bryan Mills, PhD (bmills@cs.pitt.edu)
 * STUDENTS: Implement OpenMP parallel shear sort.
 */
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include "timer.h"
#include "io.h"
#include <stdlib.h>

#define MAX_VALUE 10000
int N, M;
//DEBUG
int thread_count = 10;

inline int compare (const void * a, const void * b)
{
    return ( *(int*)a - *(int*)b  );
}

void LineSort(int **A, int begin, int step, int end) {
    int *temp_array = new int[M];
    for (int i = 0; i < M; ++i) {
        temp_array[i] = A[(begin+i*step)/M][(begin+i*step)%M];
    }
    qsort(temp_array, M, sizeof(int), compare);
    for (int i = 0; i < M; ++i) {
        A[(begin+i*step)/M][(begin+i*step)%M] = temp_array[i];
    }
}
void OddEvenRowSort(int **A, int M) {
    int i;
    #pragma omp parallel for \
        default(none) shared(A, M) private(i)
    for (i = 0; i < M; ++i) {
       if (i % 2 == 0) {
            LineSort(A, i*M, 1, i*M + M - 1);
        } else {
            LineSort(A, i*M + M - 1, -1, i*M);
        }
    }
}
void OddEvenColSort(int **A, int M) {
    int i;
    #pragma omp parallel for \
        default(none) shared(A, M) private(i)
    for (i = 0; i < M; ++i) {
            LineSort(A, i, M, M*(M - 1) + i);
    }
}

void shear_sort(int **A, int M) {
  // Students: Implement parallel shear sort here.
  int total_step = (int) ceil(log2(M*M)) + 1;
  for (int i = 0; i < total_step; ++i) {
      if (i % 2 == 0) { // RowSort
        OddEvenRowSort(A, M);
      } else { // ColumnSort
        OddEvenColSort(A, M);
      }
  }
}

// Allocate square matrix.
int **allocMatrix(int size) {
  int **matrix;
  matrix = (int **)malloc(size * sizeof(int *));
  for (int row = 0; row < size; row++) {
    matrix[row] = (int *)malloc(size * sizeof(int));
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      matrix[i][j] = 0;
    }
  }
  return matrix;
}

// Main method
int main(int argc, char* argv[]) {
  int **A;
  double elapsedTime;

  // checking parameters
  if (argc != 2 && argc != 3) {
    printf("Parameters: <N> [<file>]\n");
    return 1;
  }
  N = atoi(argv[1]);
  M = (int) sqrt(N);
  if(N != M*M){
    printf("N has to be a perfect square!\n");
    exit(1);
  }

  // allocating matrix A
  A = allocMatrix(M);

  // reading files (optional)
  if(argc == 3){
    readMatrixFile(A,M,argv[2]);
  } else {
    srand (time(NULL));
    // Otherwise, generate random matrix.
    for (int i=0; i<M; i++) {
      for (int j=0; j<M; j++) {
	A[i][j] = rand() % MAX_VALUE;
      }
    }
  }

  // starting timer
  timerStart();

  // calling shear sort function
  shear_sort(A,M);
  // stopping timer
  elapsedTime = timerStop();

// DEBUG: Test accuracy of algorithm
  for (int i = 0; i < M; i += 2) {
      for (int j = 0; j < M-1; ++j) {
          if (A[i][j] > A[i][j+1]) printf("Error: %d, %d\n", i, j);
  }
}
  for (int i = 1; i < M; i += 2) {
      for (int j = 0; j < M-1; ++j) {
          if (A[i][j] < A[i][j+1]) printf("Error: %d, %d\n", i, j);
  }
}

  for (int i = 0; i < M-1; ++i) {
      if (A[i][M-1] > A[i+1][M-1] || A[i][0] > A[i+1][0])
	  printf("Error: %d\n", i);
}


  // print if reasonably small
  if (M <= 10) {
    printMatrix(A,M);
  }
  char* pNumThreads;
  pNumThreads = getenv("OMP_NUM_THREADS"); 
  printf("(%s, %ld)\n", pNumThreads, timerStop());

  // releasing memory
  for (int i=0; i<M; i++) {
    delete [] A[i];
  }
  delete [] A;

  return 0;
}
