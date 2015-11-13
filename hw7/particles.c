/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Student:
 * Instructor: Bryan Mills, University of Pittsburgh
 * MPI particle-interaction code.
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TAG 7
#define CONSTANT 777

// Particle-interaction constants
#define A 10250000.0
#define B 726515000.5
#define MASS 0.1
#define DELTA 1

// Random initialization constants
#define POSITION 0
#define VELOCITY 1

// Structure for shared properties of a particle (to be included in messages)
struct Particle{
  float x;
  float y;
  float mass;
  float fx;
  float fy;
};

// Headers for auxiliar functions
float random_value(int type);
void print_particles(struct Particle *particles, int n);
void interact(struct Particle *source, struct Particle *destination);
void compute_interaction(struct Particle *source, struct Particle *destination, int size1, int size2);
void compute_self_interaction(struct Particle *set, int size);
void merge(struct Particle *first, struct Particle *second, int limit);
int read_file(struct Particle *set, int size, char *file_name);

void StartSimulation();
void Initial();
void GatherResult();

  int pro_num; // Number of process should be passed.
  const int sizeofp = (sizeof(struct Particle)) / sizeof(float); // How much float a particle is
  int n;// Number of total particles
  int tag = TAG;// Tag for message
  int err; // return value of any func.
  int *displs; // Because size of particles are different. it's displacement of every groups.
  int *rcounts; // How many Particles in each group.
  struct Particle *globals;// Array of all particles in the system
  struct Particle *locals;// Array of local particles
  struct Particle *remotes;// Array of foreign particles
  int myRank;// Rank of process
  int p;// Number of processes
  int previous;// Previous rank in the ring
  int next;// Next rank in the ring
  int number;// Number of local particles
  char *file_name;// File name
  MPI_Status status;// Return status for receive
  int j, rounds, initiator, sender;
  double start_time, end_time;


// Main function
int main(int argc, char** argv){

  // checking the number of parameters
  if(argc < 2){
    printf("ERROR: Not enough parameters\n");
    printf("Usage: %s <number of particles> [<file>]\n", argv[0]);
    exit(1);
  }

  // getting number of particles
  n = atoi(argv[1]);

  // initializing MPI structures and checking p is odd
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  if(p % 2 == 0){
    p = p - 1;
    if(myRank == p){
      MPI_Finalize();
      return 0;
    }
  }
  // (p-1)/2 processes should be passed.
  pro_num = (p - 1) / 2;
  srand(myRank+myRank*CONSTANT);

  // acquiring memory for particle arrays
  number = n / p + 1;
  locals = (struct Particle *) malloc(number * sizeof(struct Particle));
  remotes = (struct Particle *) malloc(number * sizeof(struct Particle));
  globals = (struct Particle *) malloc(n * sizeof(struct Particle));
  displs = (int *) malloc(p * sizeof(int));
  rcounts = (int *) malloc(p * sizeof(int));

  // checking for file information
  if(argc == 3){
    if(myRank == 0){
      // YOUR CODE GOES HERE (reading particles from file)
    file_name = argv[2];
    read_file(globals, n, file_name);
    printf("Number of process: %d\n", p);
    //printf("Gloal: \n");
    //print_particles(globals, n);
    }

    // To send/recv (or scatter/gather) you will need to learn how to
    // transfer structs of floats, treat it as a contiguous block of
    // floats. Here is an example:
    // MPI_Send(locals,
    //          number * (sizeof (struct Particle)) / sizeof(float),
    //          MPI_FLOAT,
    //          next_rank,
    //          tag,
    //          MPI_COMM_WORLD)
    // MPI_Recv(remotes,
    //          number * (sizeof (struct Particle)) / sizeof(float),
    //          MPI_FLOAT,
    //          previous_rank,
    //          tag,
    //          MPI_COMM_WORLD,
    //          &status);
    // hint: because your nodes need to both send and receive you
    // might consider asyncronous send/recv.

    // YOUR CODE GOES HERE (distributing particles among processors)

    // Random Particles
  } else {
    // random initialization of local particle array
    for(j = 0; j < number; j++){
      locals[j].x = random_value(POSITION);
      locals[j].y = random_value(POSITION);
      locals[j].fx = 0.0;
      locals[j].fy = 0.0;
      locals[j].mass = MASS;
    }
  }
    // Copy Gloals to every locals
    Initial();
  // starting timer
  if(myRank == 0){
    start_time = MPI_Wtime();
  }

  // YOUR CODE GOES HERE (ring algorithm)
  StartSimulation();
  // stopping timer
  if(myRank == 0){
    end_time = MPI_Wtime();
    printf("Duration: %f seconds\n", (end_time-start_time));
  }

    GatherResult();
  // printing information on particles
  if(argc == 3){

    // YOUR CODE GOES HERE (collect particles at rank 0)
    if(myRank == 0) {
      printf("Result:\n");
      print_particles(globals,n);
    }
  }

  // finalizing MPI structures
  MPI_Finalize();
  return 0;
}

// Function for random value generation
float random_value(int type){
  float value;
  switch(type){
  case POSITION:
    value = (float)rand() / (float)RAND_MAX * 100.0;
    break;
  case VELOCITY:
    value = (float)rand() / (float)RAND_MAX * 10.0;
    break;
  default:
    value = 1.1;
  }
  return value;
}

// Function for printing out the particle array
void print_particles(struct Particle *particles, int n){
  int j;
  printf("Index\tx\ty\tmass\tfx\tfy\n");
  for(j = 0; j < n; j++){
    printf("%d\t%f\t%f\t%f\t%f\t%f\n",j,particles[j].x,particles[j].y,particles[j].mass,particles[j].fx,particles[j].fy);
  }
}

//Scatter Data to every process
void Initial() {
    int num_q = n / p;
    int num_r = n % p;
    for (int i = 0; i < p; ++i) {
      if (i < num_r) {
          rcounts[i] = (num_q + 1) * sizeofp;
      } else {
        rcounts[i] = num_q * sizeofp;
      }
      if (i == 0) {
           displs[i] = 0;
      } else {
           displs[i] = (displs[i-1] + rcounts[i-1]);
      }
    }
    MPI_Scatterv(globals,
                 rcounts,
                 displs,
                 MPI_FLOAT,
                 locals,
                 rcounts[myRank],
                 MPI_FLOAT,
                 0,
                 MPI_COMM_WORLD
                 );
}

void GatherResult() {
    MPI_Gatherv(locals,
                rcounts[myRank],
                MPI_FLOAT,
                globals,
                rcounts,
                displs,
                MPI_FLOAT,
                0,
                MPI_COMM_WORLD
                );
    //printf("%d Gather completed.\n", myRank);
}

// Start Simulation. Follow 1-8 steps.
void StartSimulation() {
    // next_rank;
    printf("%d in %d\n", rcounts[myRank] / sizeofp, myRank);
    print_particles(locals, rcounts[myRank] / sizeofp);
    int next_rank = (myRank + 1) % p;
    MPI_Request request;
    MPI_Status status;
    // Send local to next process
    MPI_Isend(locals,
              rcounts[myRank],
              MPI_FLOAT,
              next_rank,
              myRank+1,
              MPI_COMM_WORLD,
              &request);
    //printf("process %d send to %d\n", myRank, next_rank);
    int recv_tag;
    // receive all particles from all possible pro_num
    for (int i = 1; i <= pro_num; ++i) {
        tag = (myRank - i + p) % p + 1;
        MPI_Irecv(remotes,
                  rcounts[tag-1],
                  MPI_FLOAT,
                  (myRank - 1 + p) % p, // All particles comes from previous process
                  tag, // All possible tag
                  MPI_COMM_WORLD,
                  &request);
        MPI_Wait(&request, &status);
        // compute and store
        compute_interaction(locals, remotes, rcounts[myRank] / sizeofp,rcounts[tag-1] / sizeofp);
        // If this is last node of ring, send to original one.
        if (tag == (myRank - pro_num + p) % p + 1) {
            //printf("%d begin send tag %d back\n", myRank, tag);
            MPI_Isend(remotes,
                      rcounts[tag-1],
                      MPI_FLOAT,
                      tag-1,
                      tag,
                      MPI_COMM_WORLD,
                      &request);
        } else { // else send to next_rank
            MPI_Isend(remotes,
                      rcounts[tag-1],
                      MPI_FLOAT,
                      next_rank,
                      tag,
                      MPI_COMM_WORLD,
                      &request);
        }
    }

    MPI_Irecv(remotes,
              rcounts[myRank],
              MPI_FLOAT,
              (myRank + pro_num) % p,
              myRank + 1,
              MPI_COMM_WORLD,
              &request);
    MPI_Wait(&request, &status);
    compute_self_interaction(locals, rcounts[myRank] / sizeofp);
    merge(locals, remotes, rcounts[myRank] / sizeofp);
    //printf("%d process completed.\n", myRank);

}

// Function for computing interaction among two particles
// There is an extra test for interaction of identical particles, in which case there is no effect over the destination
void interact(struct Particle *first, struct Particle *second){
  float rx,ry,r,fx,fy,f;

  // computing base values
  rx = first->x - second->x;
  ry = first->y - second->y;
  r = sqrt(rx*rx + ry*ry);

  if(r == 0.0)
    return;

  f = A / pow(r,6) - B / pow(r,12);
  fx = f * rx / r;
  fy = f * ry / r;

  // updating sources's structure
  first->fx = first->fx + fx;
  first->fy = first->fy + fy;

  // updating destination's structure
  second->fx = second->fx - fx;
  second->fy = second->fy - fy;

}

// Function for computing interaction between two sets of particles
void compute_interaction(struct Particle *first, struct Particle *second, int size1, int size2){
  int j,k;

  for(j = 0; j < size1; j++){
      for(k = 0; k < size2; k++){
      interact(&first[j],&second[k]);
    }
  }
}

// Function for computing interaction between two sets of particles
void compute_self_interaction(struct Particle *set, int size){
  int j,k;

  for(j = 0; j < size; j++){
    for(k = j+1; k < size; k++){
      interact(&set[j],&set[k]);
    }
  }
}

// Function to merge two particle arrays
// Permanent changes reside only in first array
void merge(struct Particle *first, struct Particle *second, int limit){
  int j;

  for(j = 0; j < limit; j++){
    first[j].fx += second[j].fx;
    first[j].fy += second[j].fy;
  }
}

// Reads particle information from a text file
int read_file(struct Particle *set, int size, char *file_name){
  FILE *ifp, *ofp;
  char *mode = "r";
  ifp = fopen(file_name, mode);

  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file!\n");
    return 1;
  }

  // reading particle values
  for(int i=0; i<size; i++){
    err = fscanf(ifp, "%f\t%f\t%f", &set[i].x, &set[i].y, &set[i].mass);
    set[i].fx = 0.0;
    set[i].fy = 0.0;
  }

  // closing file
  fclose(ifp);

  return 0;
}

