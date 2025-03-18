#include <math.h>
#include <stdlib.h>
#include <time.h>


#define T 3
#define DIMX 300
#define DIMY 300


void update(double u[T][DIMX][DIMY], double alpha[DIMX][DIMY]){

    for(int i = 0; i < DIMX; i++){
        for(int j = 0; j < DIMY; j++){
            u[2][i][j] = u[1][i][j];
        }
    }

    for(int i = 0; i < DIMX; i++){
        for(int j = 0; j < DIMY; j++){
            u[1][i][j] = u[0][i][j];
        }
    }

    for(int i = 1; i < DIMX-1; i++){
        for(int j = 1; j < DIMY-1; j++){
            u[0][i][j] = alpha[i][j] * (u[1][i-1][j] + u[1][i+1][j] + u[1][i][j-1] + u[1][i][j+1] - 4*u[1][i][j]);
            u[0][i][j] += 2*u[1][i][j] - u[2][i][j];
        }
    }

    for(int i = 1; i < DIMX-1; i++){
        for(int j = 1; j < DIMY-1; j++){
            u[0][i][j] *= 0.995;            //togliere energia per cagate fisiche
        }
    }

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
    int h=1;
    int l=1;


    double u[T][DIMX][DIMY];  //3 bc time has 3 parameters (u0, u1, u2)
    for(int i = 0; i < DIMX; i++){
        for(int j = 0; j < DIMY; j++){
            for(int k = 0; k < T; k++){
                u[k][i][j] = 0;
            }
        }
    }
    double c=0.5;

    double alpha[DIMX][DIMY];
    for(int i=0;i<DIMX;i++){
        for(int j=0;j<DIMY;j++){
            alpha[i][j]=pow(((c*l)/h), 2);
        }
    }

    //Random perturbations
    srand(time(0));

    while(1){
        perturbate(u);
        update(u, alpha);
    }



    return 0;
}