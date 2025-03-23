#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>


#define T 3
#define DIMX 500
#define DIMY 500
#define MAX_ITER 10000


// Initializing u and alpha
void initialize(double *u[T], double *alpha, int sendCounts[], int displs[], int size) {

  int h=1;
  int l=1;
  double c=0.5;

  int nElPerProc=0;
  int remainder=0;

  // Initializing sendCounts with the total n of element divide by n of processes and the last sendCounts cell with the remainder too
  nElPerProc=(DIMX*DIMY)/size;
  remainder=(DIMX*DIMY)%size;
  for(int i = 0; i < size-1; i++) {
    sendCounts[i]=nElPerProc;

    printf("sendCounts[%d]=%d\n",i,sendCounts[i]);
  }
  sendCounts[size-1]=nElPerProc+remainder;
  printf("sendCounts[%d]=%d\n",size-1,sendCounts[size-1]);

  displs[0]=0;
  for(int i = 1; i < size; i++) {
    displs[i]=sendCounts[i-1]*i;

    printf("displs[%d]=%d\n",i,displs[i]);
  }

  *alpha = pow(((c*l)/h), 2);

  printf("alpha=%lf\n",*alpha);

  for(int i = 0; i < T; i++) {
    u[i]= (double *) malloc(sizeof(double)*DIMX*DIMY);
  }

  for(int i = 0; i < T; i++){
    for(int j = 0; j < DIMX*DIMY; j++){
        u[i][j] = 0;
    }
  }
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
  double *u[T];
  double alpha;
  int nIterations=0;

  int *sendCounts=malloc(sizeof(int) * size);
  int *displs=malloc(sizeof(int) * size);
  int *bufferToRecv=0;

  initialize(u, &alpha, sendCounts, displs, size);

  if(rank==0){
    printf("Entered master!\n");
  }else{
    printf("Worker: %d\n", rank);
  }//end-if

  // boh potrebbe servire
  MPI_Barrier(MPI_COMM_WORLD);

  while (nIterations < MAX_ITER) {


    nIterations+=1;
  }


  /*bufferToRecv=malloc(4 * sizeof(int));
  for(int i=0; i<4; i+=1){
    bufferToRecv[i]=0;
    // printf("rank %d bufferToRecv[%d] = %d\n", rank, i, bufferToRecv[i]);
  }

  // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Scatterv(mat, sendcounts, displs, MPI_INT, bufferToRecv,4, MPI_INT, 0, MPI_COMM_WORLD);
  // MPI_Scatter(mat, 4, MPI_INT, bufferToRecv, 4, MPI_INT, 0, MPI_COMM_WORLD);

  for(int i=0; i<4; i++){
    printf("%d ", bufferToRecv[i]);
  }//end-for
  printf("\n");*/


  free(sendCounts);
  free(displs);
  for(int i = 0; i < T; i+=1){
    free(u[i]);
  }
  free(bufferToRecv);

  /* Terminate MPI */
  MPI_Finalize();
}
