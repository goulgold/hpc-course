/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Instructor Bryan Mills, PhD
 * Student:
 * Implement Pthreads version of Strassen algorithm for matrix multiplication.
 */

#include "timer.h"
#include "io.h"
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Make these globals so threads can operate on them. You will need to
// add additional matrixes for all the M and C values in the Strassen
// algorithms.
int **A;
int **B;
int **C;
// addition matrixes
int **M[7];
int **C_sub[4];
// Reference matrix, call simpleMM to populate.
int **R;

// Calculation function
int **allocMatrix(int size);
void mat_pad(int ***mat, int N, int new_N) {
    int **new_mat = allocMatrix(new_N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            new_mat[i][j] = (*mat)[i][j];
        }
    }
    int **temp;
    temp = *mat;
    *mat = new_mat;
    free(temp);
}
static inline void mat_add(int **dest, int **source1, int loc1, int **source2, int loc2, int N) {
    int row_1 = loc1 / N;
    int col_1 = loc1 % N;
    int row_2 = loc2 / N;
    int col_2 = loc2 % N;
    for (int i = 0; i < N / 2; ++i) {
        for (int j = 0; j < N / 2; ++j) {
            dest[i][j] = source1[row_1 + i][col_1 + j] + source2[row_2 + i][col_2 + j];
        }
    }
}
static inline void mat_sub(int **dest, int **source1, int loc1, int **source2, int loc2, int N) {
    int row_1 = loc1 / N;
    int col_1 = loc1 % N;
    int row_2 = loc2 / N;
    int col_2 = loc2 % N;
    for (int i = 0; i < N / 2; ++i) {
        for (int j = 0; j < N / 2; ++j) {
            dest[i][j] = source1[row_1 + i][col_1 + j] - source2[row_2 + i][col_2 + j];
        }
    }

}
static inline void mat_mul(int **dest, int **source1, int loc1, int **source2, int loc2, int N){
    int row_1 = loc1 / N;
    int col_1 = loc1 % N;
    int row_2 = loc2 / N;
    int col_2 = loc2 % N;
    for (int i=0; i<N / 2; i++) {
        for (int j=0; j<N / 2; j++) {
            for (int k=0; k<N / 2; k++) {
                dest[i][j] += source1[row_1 + i][col_1 + k] * source2[row_2 + k][col_2 + j];
            }
        }
    }

}

//M function
void* M_func1(void* N);
void* M_func2(void* N);
void* M_func3(void* N);
void* M_func4(void* N);
void* M_func5(void* N);
void* M_func6(void* N);
void* M_func7(void* N);
void* (*M_func[7])(void*) = {M_func1, M_func2, M_func3, M_func4, M_func5, M_func6, M_func7};
void C_func1(int N);
void C_func2(int N);
void C_func3(int N);
void C_func4(int N);
void (*C_func[4])(int) = {C_func1, C_func2, C_func3, C_func4};
void combine_C(int N);

void* M_func1(void* rank) {
    int* rank_ptr = (int *) rank;
    int N = *rank_ptr;
    int **temp_M_1 = allocMatrix(N);
    int **temp_M_2 = allocMatrix(N);
    mat_add(temp_M_1, A, 0, A, N * N / 2 + N / 2, N);
    mat_add(temp_M_2, B, 0, B, N * N / 2 + N / 2, N);
    mat_mul(M[0], temp_M_1, 0, temp_M_2, 0, N);
    free(temp_M_1);
    free(temp_M_2);
}
void* M_func2(void* rank) {
    int* rank_ptr = (int *) rank;
    int N = *rank_ptr;
    int **temp_M_1 = allocMatrix(N);
    mat_add(temp_M_1, A, N * N /2, A, N * N / 2 + N / 2, N);
    mat_mul(M[1], temp_M_1, 0, B, 0, N);
    free(temp_M_1);
}
void* M_func3(void* rank) {
    int* rank_ptr = (int *) rank;
    int N = *rank_ptr;
    int **temp_M_1 = allocMatrix(N);
    mat_sub(temp_M_1, B, N / 2, B, N * N / 2 + N / 2, N);
    mat_mul(M[2], A, 0, temp_M_1, 0, N);
    free(temp_M_1);
}
void* M_func4(void* rank) {
    int* rank_ptr = (int *) rank;
    int N = *rank_ptr;
    int **temp_M_1 = allocMatrix(N);
    mat_sub(temp_M_1, B, N * N / 2, B, 0, N);
    mat_mul(M[3], A, N * N / 2 + N / 2, temp_M_1, 0, N);
    free(temp_M_1);
}
void* M_func5(void* rank) {
    int* rank_ptr = (int *) rank;
    int N = *rank_ptr;
    int **temp_M_1 = allocMatrix(N);
    mat_add(temp_M_1, A, 0, A, N / 2, N);
    mat_mul(M[4], temp_M_1, 0, B, N * N / 2 + N / 2, N);
    free(temp_M_1);
}
void* M_func6(void* rank) {
    int* rank_ptr = (int *) rank;
    int N = *rank_ptr;
    int **temp_M_1 = allocMatrix(N);
    int **temp_M_2 = allocMatrix(N);
    mat_sub(temp_M_1, A, N * N / 2, A, 0, N);
    mat_add(temp_M_2, B, 0, B, N / 2, N);
    mat_mul(M[5], temp_M_1, 0, temp_M_2, 0, N);
    free(temp_M_1);
    free(temp_M_2);
}
void* M_func7(void* rank) {
    int* rank_ptr = (int *) rank;
    int N = *rank_ptr;
    int **temp_M_1 = allocMatrix(N);
    int **temp_M_2 = allocMatrix(N);
    mat_sub(temp_M_1, A, N / 2, A, N * N / 2 + N / 2, N);
    mat_add(temp_M_2, B, N * N / 2, B, N * N / 2 + N / 2, N);
    mat_mul(M[6], temp_M_1, 0, temp_M_2, 0, N);
    free(temp_M_1);
    free(temp_M_2);
}
void C_func1(int N) {
    mat_add(C_sub[0], M[0], 0, M[3], 0, N);
    mat_add(C_sub[0], C_sub[0], 0, M[6], 0, N);
    mat_sub(C_sub[0], C_sub[0], 0, M[4], 0, N);
}
void C_func2(int N) {
    mat_add(C_sub[1], M[2], 0, M[4], 0, N);
}
void C_func3(int N) {
    mat_add(C_sub[2], M[1], 0, M[3], 0, N);
}
void C_func4(int N) {
    mat_add(C_sub[3], M[0], 0, M[2], 0, N);
    mat_add(C_sub[3], C_sub[3], 0, M[5], 0, N);
    mat_sub(C_sub[3], C_sub[3], 0, M[1], 0, N);
}
void combine_C(int N) {
    for (int i = 0; i < 4; ++i) {
        int row = (i / 2) * N / 2;
        int col = (i % 2) * N / 2;
        for (int j = 0; j < N / 2; ++j) {
            for (int k = 0; k < N / 2; ++k) {
                C[row+j][col+k] = C_sub[i][j][k];
            }
        }
    }
}
// Stupid simple Matrix Multiplication, meant as example.
void simpleMM(int N) {
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
	R[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

// WRITE YOUR CODE HERE, you will need to also add functions for each
// of the sub-matrixes you will need to calculate but you can create your
// threads in this fucntion.
void strassenMM(int N) {
    int new_N = 1;
    while (new_N < N) {
        new_N *= 2;
    }
    // pad matrix

    timerStart();
    if (N != new_N) {
        mat_pad(&A, N, new_N);
        mat_pad(&B, N, new_N);
        mat_pad(&C, N, new_N);
        for (int i = 0; i < 7; ++i) {
            mat_pad(&M[i], N, new_N);
        }
        for (int i = 0; i <4; ++i) {
            mat_pad(&C_sub[i], N, new_N);
        }

    }
    printf("matrix pad took %ld ms\n", timerStop());
    N = new_N;
    pthread_t ids[7];
    timerStart();
    for (int i = 0; i < 7; ++i) {
        pthread_create(&ids[i], NULL, *M_func[i], &N);
    }
    for (int i=0; i < 7; i++) {
        pthread_join(ids[i], NULL);
    }
    printf("seven M funcs took %ld ms\n", timerStop());
    C_func1(N);
    C_func2(N);
    C_func3(N);
    C_func4(N);
    combine_C(N);
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

// Allocate memory for all the matrixes, you will need to add code
// here to initialize any matrixes that you need.
void initMatrixes(int N) {
  A = allocMatrix(N); B = allocMatrix(N); C = allocMatrix(N); R = allocMatrix(N);
  for (int i = 0; i < 7; ++i) {
       M[i] = allocMatrix(N);
  }
  for (int i = 0; i < 4; ++i) {
      C_sub[i] = allocMatrix(N);
  }
}

// Free up matrixes.
void cleanup() {
  free(A);
  free(B);
  free(C);
  free(R);
  for (int i = 0; i < 7; ++i) {
      free(M[i]);
  }
  for (int i = 0; i < 4; ++i) {
      free(C_sub[i]);
  }
}

// Main method
int main(int argc, char* argv[]) {
  int N;
  double elapsedTime;

  // checking parameters
  if (argc != 2 && argc != 4) {
    printf("Parameters: <N> [<fileA> <fileB>]\n");
    return 1;
  }
  N = atoi(argv[1]);
  initMatrixes(N);

  // reading files (optional)
  if(argc == 4){
    readMatrixFile(A,N,argv[2]);
    readMatrixFile(B,N,argv[3]);
  } else {
    // Otherwise, generate two random matrix.
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
	A[i][j] = rand() % 5;
	B[i][j] = rand() % 5;
      }
    }
  }

  // Do simple multiplication and time it.
  timerStart();
  simpleMM(N);
  printf("Simple MM took %ld ms\n", timerStop());

  // Do strassen multiplication and time it.
  timerStart();
  strassenMM(N);
  printf("Strassen MM took %ld ms\n", timerStop());

  if (compareMatrix(C, R, N) != 0) {
    if (N < 20) {
      printf("\n\n------- MATRIX C\n");
      printMatrix(C,N);
      printf("\n------- MATRIX R\n");
      printMatrix(R,N);
    }
    printf("Matrix C doesn't match Matrix R, if N < 20 they will be printed above.\n");
  }

  // stopping timer

  cleanup();
  return 0;
}
