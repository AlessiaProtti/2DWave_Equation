#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

void main(int argc, char *argv[]){
  /* 1. Initialize MPI */
  MPI_Init(&argc, &argv);

  int **mat=0;
  //int mat[4][4]={{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
  int sendcounts[4]={4, 4, 4, 4};
  int displs[4]={0, 4, 8, 12};
  int *bufferToRecv=0;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0){
    printf("Entered master!\n");

    mat=malloc(4 * sizeof(int *));
    for(int i = 0; i < 4; i++) mat[i] = (int *)malloc(4 * sizeof(int));

    int con=1;
    for(int i=0; i<4; i+=1){
      for(int j=0; j<4; j+=1){
        mat[i][j] = con;
        // printf("mat[%d][%d] = %d\n", i, j, mat[i][j]);
        con+=1;
      }
    }
  }else{
    printf("Worker: %d\n", rank);
  }//end-if

  bufferToRecv=malloc(4 * sizeof(int));
  for(int i=0; i<4; i+=1){
    bufferToRecv[i]=0;
    // printf("rank %d bufferToRecv[%d] = %d\n", rank, i, bufferToRecv[i]);
  }

  // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Scatterv(mat, sendcounts, displs, MPI_INT, bufferToRecv,4, MPI_INT, 0, MPI_COMM_WORLD);

  for(int i=0; i<4; i++){
    printf("%d ", bufferToRecv[i]);
  }//end-for
  printf("\n");

  if (rank==0)
  {
    free(mat);
  }
  free(bufferToRecv);

  /* Terminate MPI */
  MPI_Finalize();
}
