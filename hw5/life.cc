/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Instructor Bryan Mills, PhD (bmills@cs.pitt.edu)
 * Students:
 * Implement openmp verions of conway's game of life.
 */

#include "timer.h"
#include "omp.h"
#include "io.h"

//DEBUG
int thread_count = 10;


int **allocMatrix(int size);

// Function deciding next state of certain node
int nextState(int **World, int row, int col, int N);

// Function that counts how many neighbors a node has
int getNeighbor(int **World,int row,int col,int N);

// Function that determine whether a node lives or not.
inline bool isLive(int value) {
    return (value != 0);
}
// Function implementing Conway's Game of Life
void conway(int **World, int N, int M){
  // STUDENT: IMPLEMENT THE GAME HERE, make it parallel!
    int now_m; // what generate it is now.
    int **nextWorld = allocMatrix(N);
    int i;
    for (now_m = 0; now_m < M; ++now_m) {
        printMatrix(World, N);
        printf("\n");
        # pragma omp parallel for num_threads(thread_count) \
        default(none) shared(World, nextWorld, N) private(i)
        for (i = 0; i < N*N; ++i) {
            nextWorld[i/N][i%N] = nextState(World, i/N, i%N, N);
        }
        # pragma omp parallel for num_threads(thread_count) \
        default(none) shared(World, nextWorld, N) private(i)
        for (i = 0; i < N*N; ++i) {
            World[i/N][i%N] = nextWorld[i/N][i%N];
        }
    }

    removeMatrix(nextWorld, N);
}

int getNeighbor(int **World,int row,int col,int N) {
    int neighbor = 0;
    for (int i = row-1; i <= row+1; ++i) {
        for (int j = col-1; j <= col+1; ++j) {
            if (i == row && j == col) continue;
            else if(i < 0 || j < 0 || i > N-1 || j > N-1) continue;
            else if(isLive(World[i][j])) neighbor++;
            else continue;
        }
    }
    return neighbor;
}
int nextState(int **World, int row, int col, int N) {
    int neighbor = getNeighbor(World, row, col, N);
    if (neighbor < 2 && isLive(World[row][col])) {
        return 0;
    } else if ((neighbor == 2 || neighbor == 3) && isLive(World[row][col])) {
        return World[row][col] + 1;
    } else if (neighbor > 3 && isLive(World[row][col])) {
         return 0;
    } else if (neighbor == 3 && !isLive(World[row][col])) {
        return 1;
    } else {
        return World[row][col];
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
  int N,M;
  int **World;
  double elapsedTime;
  // checking parameters
  if (argc != 3 && argc != 4) {
    printf("Parameters: <N> <M> [<file>]\n");
    return 1;
  }
  N = atoi(argv[1]);
  M = atoi(argv[2]);

  // allocating matrices
  World = allocMatrix(N);

  // reading files (optional)
  if(argc == 4){
    readMatrixFile(World,N,argv[3]);
  } else {
    // Otherwise, generate two random matrix.
    srand (time(NULL));
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
	World[i][j] = rand() % 2;
      }
    }
  }

  // starting timer
  timerStart();

  // calling conway's game of life
  conway(World,N,M);

  // stopping timer
  elapsedTime = timerStop();

  printMatrix(World,N);

  printf("Took %ld ms\n", timerStop());

  // releasing memory
  for (int i=0; i<N; i++) {
    delete [] World[i];
  }
  delete [] World;

  return 0;
}
