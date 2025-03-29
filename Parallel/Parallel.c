#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define T 3
#define DIMX 8
#define DIMY 8
#define MAX_ITER 10

void perturbate(double *u){
  int x=0;
  int y=0;

  //Generate random number for perturbations
  double num=rand() % 1;
  if(num<0.02){
    x=rand() % ((DIMX-5) - 5 + 1) + 5;
    y=rand() % ((DIMY-5) - 5 + 1) + 5;

    for(int i = (x-2); i < (x+2); i+=1){
      for(int j = (y-2); j < (y+2); j+=1){
        // printf("position %d,%d converted in %d\n", i, j, i*DIMY+j);

        u[(i * DIMY) + j] = 120.00;
      }//end-for
    }//end-for
  }//end-for
}

void matrixShifting(double *u, int sendCount, const double v[]){
  for(int i=0; i<sendCount; i+=1){
    u[i]=v[i];
  }//end-for
}

void update(double *buffer, int disp, int sendCount, const double u1[], const double u2[], double alpha){

  //true update
  int x, y=0;
  for(int i=0; i<sendCount; i+=1){
    //matrix notation
    x = (i+disp+1)/(DIMY);
    y = ((i+disp+1)%(DIMY))-1;

    if (x!=0 && y!=0 && x!=DIMX-1 && y!=DIMY-1) {
      buffer[i]= alpha*(u1[(x-1)*DIMY + y] + u1[(x+1)*DIMY + y] + u1[(x)*DIMY +(y-1)] + u1[(x)*DIMY + (y+1)] -4*u1[x*DIMY+y]);
      buffer[i]+= 2*u1[x*DIMY+y] -u2[i];
    }


  }



  // for(int i = 1; i < DIMX-1; i++){
  //   for(int j = 1; j < DIMY-1; j++){
  //     u[0][i][j] = alpha * (u[1][i-1][j] + u[1][i+1][j] + u[1][i][j-1] + u[1][i][j+1] - 4*u[1][i][j]);
  //     u[0][i][j] += 2*u[1][i][j] - u[2][i][j];
  //
  //     //double prova= (alpha[i][j] * (u[1][i-1][j] + u[1][i+1][j] + u[1][i][j-1] + u[1][i][j+1] - 4*u[1][i][j])) + 2*u[1][i][j] - u[2][i][j];
  //   }//end-for
  // }//end-for
}

void absorbingEnergy(double *u, int sendCount, int disp) {
  int x, y=0;
  for(int i=0; i<sendCount; i+=1){
    //matrix notation
    x = (i+disp+1)/(DIMY);
    y = ((i+disp+1)%(DIMY))-1;

    if (x!=0 && y!=0 && x!=DIMX-1 && y!=DIMY-1) {
      u[i] *= 0.995;
    }

  }
}

/*void absorbingEnergy(double u[DIMX][DIMY]){

  //Absorbing
  for(int i = 1; i < DIMX-1; i++){
    for(int j = 1; j < DIMY-1; j++){
      u[i][j] *= 0.995;
    }//end-for
  }//end-for
}*/

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
void initialize(double *u[], double *alpha, int sendCounts[], int displs[], int size, int rank){

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

  if(rank==0) {

    for(int i = 0; i < T; i+=1) {
      u[i]=(double *)malloc(sizeof(double)*DIMX*DIMY);
    }//end-for

    for(int i = 0; i < T; i+=1){
      for(int j = 0; j < DIMX*DIMY; j+=1){
        u[i][j]=0;
      }//end-for
    }//end-for
  }else {
    u[1]=(double *)malloc(sizeof(double)*DIMX*DIMY);

    for(int j = 0; j < DIMX*DIMY; j+=1){
      u[1][j]=0;
    }//end-for

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
  srand(1);
  double *u[T];
  double alpha;
  int nIterations=0;

  int *sendCounts=malloc(sizeof(int) * size);
  int *displs=malloc(sizeof(int) * size);


  initialize(u, &alpha, sendCounts, displs, size, rank);

  //SENDCOUNT NEEDS TO BE INITALIZED BEFORE!!!!!!!!!!!!!!!
  double *bufferToRecvU2=malloc(sizeof(double) * sendCounts[rank]);
  double *bufferToRecvU1=malloc(sizeof(double) * sendCounts[rank]);
  double *bufferToRecvU0=malloc(sizeof(double) * sendCounts[rank]);



  // boh potrebbe servire
  MPI_Barrier(MPI_COMM_WORLD);

  while (nIterations < MAX_ITER){

    if (rank==0) {
      perturbate(u[0]);
    }

    // printf("rank: %d, displs: %d\n", rank, displs[rank]);
    //Scatter u[2]
    MPI_Scatterv(u[2], sendCounts, displs, MPI_DOUBLE, bufferToRecvU2, sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(u[1], sendCounts, displs, MPI_DOUBLE, bufferToRecvU1, sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);


    matrixShifting(bufferToRecvU2, sendCounts[rank], bufferToRecvU1);

    MPI_Gatherv(bufferToRecvU2, sendCounts[rank], MPI_DOUBLE, u[2], sendCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Bcast(u[2], DIMX*DIMY, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //gather u[2]


    MPI_Scatterv(u[0], sendCounts, displs, MPI_DOUBLE, bufferToRecvU0, sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    matrixShifting(bufferToRecvU1, sendCounts[rank], bufferToRecvU0);

    MPI_Gatherv(bufferToRecvU1, sendCounts[rank], MPI_DOUBLE, u[1], sendCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(u[1], DIMX*DIMY, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    update(bufferToRecvU0, displs[rank], sendCounts[rank], u[1], bufferToRecvU2, alpha);

    absorbingEnergy(bufferToRecvU0, sendCounts[rank], displs[rank]);



    //update & absorbing energy



    MPI_Gatherv(bufferToRecvU0, sendCounts[rank], MPI_DOUBLE, u[0], sendCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //gather u[0]


    printMatrix(u[1]);
    nIterations+=1;
  }//end-while

  free(sendCounts);
  free(displs);
  if(rank==0) {
    for(int i = 0; i < T; i+=1){
      free(u[i]);
    }//end-for

  }else {
    free(u[1]);
  }

  free(bufferToRecvU0);
  free(bufferToRecvU1);
  free(bufferToRecvU2);

  /* Terminate MPI */
  MPI_Finalize();
}
