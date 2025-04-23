#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define T 3
#define DIMX 10000
#define DIMY 10000
#define MAX_ITER 100


void initialize(double *u[T], double *alpha){
    int h=1;
    int l=1;
    double c=0.5;


    *alpha = pow(((c*l)/h), 2);

    for(int i=0; i<T; i+=1){
        u[i]=(double *)malloc(sizeof(double)*DIMX*DIMY);
    }//end-for

    for(int k=0; k<T; k+=1){
        for(int i=0; i<DIMX*DIMY; i+=1){
            u[k][i] = 0;
        }//end-for
    }//end-for
}

void matrixShifting(double u[], double v[]){
    for(int i = 0; i < DIMX*DIMY; i+=1){
        u[i] = v[i];
    }//end-for
}

void absorbingEnergy(double u[]){

    //Absorbing
    int x, y=0;
    for(int i=0; i<DIMX*DIMY; i+=1){
        //matrix notation
        x = (i+1)/(DIMY);
        y = ((i+1)%(DIMY))-1;

        if (x!=0 && y!=0 && x!=DIMX-1 && y!=DIMY-1) {
            u[i] *= 0.995;
        }//end-if
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



void update(double u0[], double u1[], double u2[], double alpha){

    //true update
     // for(int i = 1; i < DIMX-1; i++){
     //     for(int j = 1; j < DIMY-1; j++){
     //         u[0][i][j] = alpha * (u[1][i-1][j] + u[1][i+1][j] + u[1][i][j-1] + u[1][i][j+1] - 4*u[1][i][j]);
     //         u[0][i][j] += 2*u[1][i][j] - u[2][i][j];
     //
     //         //double prova= (alpha[i][j] * (u[1][i-1][j] + u[1][i+1][j] + u[1][i][j-1] + u[1][i][j+1] - 4*u[1][i][j])) + 2*u[1][i][j] - u[2][i][j];
     //     }//end-for
     // }//end-for

    int x, y=0;
    for(int i=0; i<DIMX*DIMY; i+=1){
        //matrix notation
        x = (i+1)/(DIMY);
        y = ((i+1)%(DIMY))-1;

        if (x!=0 && y!=0 && x!=DIMX-1 && y!=DIMY-1){
            u0[i]= alpha*(u1[(x-1)*DIMY + y] + u1[(x+1)*DIMY + y] + u1[(x)*DIMY +(y-1)] + u1[(x)*DIMY + (y+1)] -4*u1[x*DIMY+y]);
            u0[i]+= 2*u1[x*DIMY+y] - u2[i];
        }//end-if
    }//end-for
}

void perturbate(double u[]){
    int x=0;
    int y=0;

    //Generate random number for perturbations
    double num=rand() % 1;
    if(num<0.02){
        x=rand() % ((DIMX-5) - 5 + 1) + 5;
        y=rand() % ((DIMY-5) - 5 + 1) + 5;

        for(int i = (x-2); i < (x+2); i++){
            for(int j = (y-2); j < (y+2); j++){
                u[(i * DIMY) + j] = 120.00;
                // u[0][i][j] = 120;
            }
        }
    }
}


void createFile(FILE *fp, double u[], char filename[]) {

    fp=fopen(filename,"a+");
    for(int i=0; i<DIMX; i+=1) {
        for(int j=0; j<DIMY; j+=1) {
            fprintf(fp,"%f ", u[(i * DIMY) + j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}


int main(int argc, char **argv){
    //Initializing variables
    srand(1);
    double *u[T];  //3 bc time has 3 parameters (u0, u1, u2)
    double alpha;
    //
    // FILE *fp0;
    // FILE *fp1;
    // FILE *fp2;

    initialize(u, &alpha);

    int count=0;
    while(count<MAX_ITER){
        perturbate(u[0]);

        //u[2]=u[1]
        matrixShifting(u[2], u[1]);

        //u[1]=u[0]
        matrixShifting(u[1], u[0]);

        update(u[0], u[1], u[2], alpha);

        absorbingEnergy(u[0]);

        // createFile(fp0, u[0], "u0.txt");
        // createFile(fp1, u[1], "u1.txt");
        // createFile(fp2, u[2], "u2.txt");

        //Visualizing
        //printMatrix(u[1]);

        count+=1;
    }//end-while

    for(int i=0; i<T; i+=1){
        free(u[i]);
    }//end-for

    return 0;
}