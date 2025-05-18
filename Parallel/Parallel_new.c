#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define T 3
#define DIMX 10000
#define DIMY 10000
#define RAINDROP 120.00
#define MAX_ITER 40

/**********************************************************/
// Method to initialize u, alpha, sendCounts, displs.
/**********************************************************/
void initialize(double *u[], double *alpha, int sendCounts[], int displs[], const int size, const int rank){
/**********************************************************/
//LOCAL VARIABLES
/**********************************************************/
  const int h=1;
  const int l=1;
  const double c=0.5;

  int nElPerProc=0;
  int remainder=0;
/**********************************************************/
//METHOD BODY
/**********************************************************/

  //Initializing sendCounts with the total n. of elements divided by n. of processes and the last sendCounts cell with the remainder too
  nElPerProc=(DIMX*DIMY)/size;
  remainder=(DIMX*DIMY)%size;
  for(int i = 0; i < size-1; i+=1){
    sendCounts[i]=nElPerProc;
  }//end-for
  sendCounts[size-1]=nElPerProc+remainder;

  displs[0]=0;
  for(int i = 1; i < size; i+=1){
    displs[i]=sendCounts[i-1]*i;
  }//end-for

  //Initializing alpha
  *alpha = ((c*l)/h)*((c*l)/h);

  //Master initializes u0, u1, u2, workers initialize u1 only
  if(rank==0){
    for(int i = 0; i < T; i+=2){
      u[i]=(double *)malloc(sizeof(double)*DIMX*DIMY);
    }//end-for

    for(int i = 0; i < T; i+=2){
      for(int j = 0; j < DIMX*DIMY; j+=1){
        u[i][j]=0;
      }//end-for
    }//end-for
  }//end-if
  u[1]=(double *)malloc(sizeof(double)*DIMX*DIMY);

  for(int j = 0; j < DIMX*DIMY; j+=1){
    u[1][j]=0;
  }//end-for
/**********************************************************/
}

/**********************************************************/
// Method to virtually place a raindrop
/**********************************************************/
void perturbate(double *u, const int sendCount, const int disp){
/**********************************************************/
//LOCAL VARIABLES
/**********************************************************/
  const double num=rand() % 1;
/**********************************************************/
//METHOD BODY
/**********************************************************/
  if(num<0.02){
    const int pos=rand() % ((sendCount-2)+1);

    for(int i=pos; i<pos+2; i+=1){
      const int x = (i+disp)/(DIMY);
      const int y = ((i+disp)%(DIMY));

      if(x!=0 && y!=0 && x!=DIMX-1 && y!=DIMY-1){
        u[i]=RAINDROP;
      }//end-if
    }//end-for
  }//end-if
/**********************************************************/
}

/**********************************************************/
// Method to shift matrices
/**********************************************************/
void matrixShifting(double *u, const int sendCount, const double v[]){
/**********************************************************/
//METHOD BODY
/**********************************************************/
  for(int i=0; i<sendCount; i+=1){
    u[i]=v[i];
  }//end-for
/**********************************************************/
}

/**********************************************************/
// Method implementing the actual 2D wave equation
/**********************************************************/
void update(double *buffer, const int disp, const int sendCount, const double u1[], const double u2[], const double alpha){
/**********************************************************/
//LOCAL VARIABLES
/**********************************************************/
  int x=0, y=0;
/**********************************************************/
//METHOD BODY
/**********************************************************/
  for(int i=0; i<sendCount; i+=1){
    //matrix notation
    x = (i+disp)/(DIMY);
    y = ((i+disp)%(DIMY));

    if(x!=0 && y!=0 && x!=DIMX-1 && y!=DIMY-1){
      buffer[i]= alpha*(u1[(x-1)*DIMY + y] + u1[(x+1)*DIMY + y] + u1[(x)*DIMY +(y-1)] + u1[(x)*DIMY + (y+1)] -4*u1[x*DIMY+y]);
      buffer[i]+= 2*u1[x*DIMY+y] - u2[i];
    }//end-if
  }//end-for
/**********************************************************/
}

/**********************************************************/
// Method used to remove energy from the system
/**********************************************************/
void absorbingEnergy(double *u, const int sendCount, const int disp){
/**********************************************************/
//LOCAL VARIABLES
/**********************************************************/
  int y=0;
/**********************************************************/
//METHOD BODY
/**********************************************************/
  for(int i=0; i<sendCount; i+=1){
    //matrix notation
    const int x = (i + disp)/(DIMY);
    y=((i+disp)%(DIMY));

    if(x!=0 && y!=0 && x!=DIMX-1 && y!=DIMY-1){
      u[i] *= 0.995;
    }//end-if
  }//end-for
/**********************************************************/
}

/**********************************************************/
// Method used to visualize the matrices
/**********************************************************/
void printMatrix(double u[]){
/**********************************************************/
//METHOD BODY
/**********************************************************/
  for(int i = 0; i < DIMX; i+=1){
    for(int j = 0; j < DIMY; j+=1){
      printf("%3.1f ", u[(i*DIMY) + j]);
    }//end-for
    printf("\n");
  }//end-for
  printf("\n\n\n");
/**********************************************************/
}

/**********************************************************/
int main(int argc, char *argv[]){
/**********************************************************/
  /* 1. Initialize MPI */
  MPI_Init(&argc, &argv);
/**********************************************************/
// MPI variables
/**********************************************************/
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
/**********************************************************/
// Domain variables
/**********************************************************/
  srand(1);
  double *u[T];
  double alpha;
  int nIterations=0;

  int *sendCounts=malloc(sizeof(int) * size);
  int *displs=malloc(sizeof(int) * size);
  if(!sendCounts || !displs){
    perror("Error when allocating sendCounts or displs!\n");
    exit(EXIT_FAILURE);
  }//end-if

  //SENDCOUNT NEEDS TO BE INITALIZED BEFORE!!!!!!!!!!!!!!!
  initialize(u, &alpha, sendCounts, displs, size, rank);
  double *bufferToRecvU2=malloc(sizeof(double) * sendCounts[rank]);
  double *bufferToRecvU1=malloc(sizeof(double) * sendCounts[rank]);
  double *bufferToRecvU0=malloc(sizeof(double) * sendCounts[rank]);
  if(!bufferToRecvU2 || !bufferToRecvU1 || !bufferToRecvU0){
    perror("Error when allocating one of the buffers!\n");
    exit(EXIT_FAILURE);
  }//end-if
/**********************************************************/
// Main body
/**********************************************************/
  MPI_Scatterv(u[2], sendCounts, displs, MPI_DOUBLE, bufferToRecvU2, sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(u[1], sendCounts, displs, MPI_DOUBLE, bufferToRecvU1, sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(u[0], sendCounts, displs, MPI_DOUBLE, bufferToRecvU0, sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  while(nIterations < MAX_ITER){
    perturbate(bufferToRecvU0, sendCounts[rank], displs[rank]);

    //u[2]=u[1]
    matrixShifting(bufferToRecvU2, sendCounts[rank], bufferToRecvU1);

    //u[1]=u[0]
    matrixShifting(bufferToRecvU1, sendCounts[rank], bufferToRecvU0);

    //All gather of u[1]
    MPI_Allgatherv(bufferToRecvU1, sendCounts[rank], MPI_DOUBLE, u[1], sendCounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    //update & absorbing energy
    update(bufferToRecvU0, displs[rank], sendCounts[rank], u[1], bufferToRecvU2, alpha);
    absorbingEnergy(bufferToRecvU0, sendCounts[rank], displs[rank]);

    // if (rank==0){
    //   printMatrix(u[1]);
    // }
    nIterations+=1;
  }//end-while

  if(rank==0){
    for(int i = 0; i < T; i+=1){
      free(u[i]);
    }//end-for
  }else{
    free(u[1]);
  }//end-if

  free(bufferToRecvU0);
  free(bufferToRecvU1);
  free(bufferToRecvU2);
  free(sendCounts);
  free(displs);
/**********************************************************/
  /* Terminate MPI */
  MPI_Finalize();
/**********************************************************/
  return 0;
}
