#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define T 3
#define DIMX 300
#define DIMY 300
#define MAX_ITER 1000


void initialize(double u[T][DIMX][DIMY], double alpha[DIMX][DIMY]){
    int h=1;
    int l=1;

    for(int i = 0; i < DIMX; i++){
        for(int j = 0; j < DIMY; j++){
            for(int k = 0; k < T; k++){
                u[k][i][j] = 0;
            }
        }
    }
    double c=0.5;

    for(int i=0;i<DIMX;i++){
        for(int j=0;j<DIMY;j++){
            alpha[i][j]=pow(((c*l)/h), 2);
        }
    }
}

void update(double u[T][DIMX][DIMY], double alpha[DIMX][DIMY]){

    //u[2]=u[1]
    for(int i = 0; i < DIMX; i++){
        for(int j = 0; j < DIMY; j++){
            u[2][i][j] = u[1][i][j];
        }
    }

    //u[1]=u[0]
    for(int i = 0; i < DIMX; i++){
        for(int j = 0; j < DIMY; j++){
            u[1][i][j] = u[0][i][j];
        }
    }

    //true update
     for(int i = 1; i < DIMX-1; i++){
         for(int j = 1; j < DIMY-1; j++){
             u[0][i][j] = alpha[i][j] * (u[1][i-1][j] + u[1][i+1][j] + u[1][i][j-1] + u[1][i][j+1] - 4*u[1][i][j]);
             u[0][i][j] += 2*u[1][i][j] - u[2][i][j];
         }
     }

    //visualizing
    for(int i = 1; i < DIMX-1; i++){
        for(int j = 1; j < DIMY-1; j++){
            u[0][i][j] *= 0.995;            //togliere energia per cagate fisiche
        }
    }

    for(int i = 0; i < DIMX; i++){
        for(int j = 0; j < DIMY; j++){
            printf("%-6f ", u[1][i][j]);
        }
        printf("\n");
    }
    printf("\n\n\n");


}

void perturbate(double u[T][DIMX][DIMY]){
    int x=0;
    int y=0;

    //Generate random number for perturbations
    double num=rand() % 1;
    if(num<0.02){
        x=rand() % ((DIMX-5) - 5 + 1) + 5;
        y=rand() % ((DIMY-5) - 5 + 1) + 5;

        for(int i = (x-2); i < (x+2); i++){
            for(int j = (y-2); j < (y+2); j++){
                u[0][i][j] = 120;
            }
        }
    }

}


int main(void){
    //Initializing variables
    srand(time(0));
    double u[T][DIMX][DIMY];  //3 bc time has 3 parameters (u0, u1, u2)
    double alpha[DIMX][DIMY];

    initialize(u, alpha);

    int count=0;
    while(count<MAX_ITER){
        perturbate(u);
        update(u, alpha);
        count+=1;
    }//end-while

    return 0;
}