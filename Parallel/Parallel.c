#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define T 3
#define DIMX 4
#define DIMY 4
#define MAX_ITER 1

void perturbate(double u[]){
  int x=0;
  int y=0;

  //Generate random number for perturbations
  double num=rand() % 1;
  if(num<0.02){
    x=rand() % ((DIMX-5) - 5 + 1) + 5;
    y=rand() % ((DIMY-5) - 5 + 1) + 5;

    for(int i = (x-2); i < (x+2); i+=1){
      for(int j = (y-2); j < (y+2); j+=1){
        u[(i * DIMY) + j] = 120;
      }//end-for
    }//end-for
  }//end-for
}

void matrixShifting(double *buffer, int disp, int sendCount, const double v[]){
  for(int i=0; i<sendCount; i+=1){
    buffer[i]=v[disp+i];
  }//end-for
}

void printMatrix(double u[]){
  for(int i = 0; i < DIMX; i+=1){
    for(int j = 0; j < DIMY; j+=1){
      printf("%3.1f ", u[(i*DIMY) + j]);
    }//end-for
    printf("\n");
  }//end-for
  printf("\n\n\n");
}


// Initializing u and alpha
void initialize(double *u[T], double *alpha, int sendCounts[], int displs[], int size){

  int h=1;
  int l=1;
  double c=0.5;

  int nElPerProc=0;
  int remainder=0;

  // Initializing sendCounts with the total n of element divide by n of processes and the last sendCounts cell with the remainder too
  nElPerProc=(DIMX*DIMY)/size;
  remainder=(DIMX*DIMY)%size;
  for(int i = 0; i < size-1; i+=1) {
    sendCounts[i]=nElPerProc;

    // printf("sendCounts[%d]=%d\n",i,sendCounts[i]);
  }//end-for
  sendCounts[size-1]=nElPerProc+remainder;
  // printf("sendCounts[%d]=%d\n",size-1,sendCounts[size-1]);

  displs[0]=0;
  for(int i = 1; i < size; i+=1) {
    displs[i]=sendCounts[i-1]*i;

    // printf("displs[%d]=%d\n",i,displs[i]);
  }//end-for

  *alpha = pow(((c*l)/h), 2);

  // printf("alpha=%lf\n",*alpha);

  for(int i = 0; i < T; i+=1) {
    u[i]= (double *)malloc(sizeof(double)*DIMX*DIMY);
  }//end-for

  for(int i = 0; i < T; i+=1){
    for(int j = 0; j < DIMX*DIMY; j++){
        u[i][j] = 0;
    }//end-for
  }//end-for
}

void main(int argc, char *argv[]){
  /* 1. Initialize MPI */
  MPI_Init(&argc, &argv);

  // MPI variables
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Domain variables
  srand(1);
  double *u[T];
  double alpha;
  int nIterations=0;

  int *sendCounts=malloc(sizeof(int) * size);
  int *displs=malloc(sizeof(int) * size);

  initialize(u, &alpha, sendCounts, displs, size);

  //SENDCOUNT NEEDS TO BE INITALIZED BEFORE!!!!!!!!!!!!!!!
  double *bufferToRecv=malloc(sizeof(double) * sendCounts[rank]);

  // boh potrebbe servire
  MPI_Barrier(MPI_COMM_WORLD);

  double tmp[]={1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};

  while (nIterations < MAX_ITER){
    perturbate(u[0]);

    // printf("rank: %d, displs: %d\n", rank, displs[rank]);
    //Scatter u[2]
    MPI_Scatterv(u[2], sendCounts, displs, MPI_DOUBLE, bufferToRecv, sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // matrixShifting(bufferToRecv, displs[rank], sendCounts[rank], tmp);
    //
    // for(int i = 0; i < sendCounts[rank]; i+=1){
    //   printf("%3.1f ", bufferToRecv[i]);
    // }//end-for

    //printMatrix(u[0]);
    nIterations+=1;
  }//end-while

  free(sendCounts);
  free(displs);
  for(int i = 0; i < T; i+=1){
    free(u[i]);
  }//end-for
  free(bufferToRecv);

  /* Terminate MPI */
  MPI_Finalize();
}
