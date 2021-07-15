//--------------------------------------------------------//
//--------------------------------------------------------//
//--------------------------------------------------------//
//  3D FINITE DIFFERENCE ANALYSIS OF ARC WELDING          //
//  National Technical University of Athens               //
//  School of Mechanical Engineering                      //
//  Manufacturing Technology Lab                          //
//                                                        //
//  created by Dimitrios T. Pallas,                       //
//  under the supervision of George C. Vosniakos (2021)   //
//                                                        //
//                                                        //
//                                                        //
//--------------------------------------------------------//
//--------------------------------------------------------//
//--------------------------------------------------------//

#include <stdio.h>
#include "source.h"
#include <stdlib.h>
#include <time.h>

int main() {

    clock_t TIME_START, TIME_STOP;
    TIME_START = clock();

    //Input Variables , user defined
    double dt = 0.01;
    double sigma = 2e-3;
    double h = 30;
    double Ta = 20;
    double Q = 990;
    double electrode_vel = 40e-3 / 12;
    double length = 40e-3;
    double width  = 50e-3;
    double thickness = 3e-3;
    double total_time = 20;
    double dx =  1e-3;
    double dy =  1e-3;
    double dz =  1e-3;
    int Nx = (int) (length / dx) + 1;
    int Ny = (int) (width / dy) + 1;
    int Nz = (int) (thickness / dz) + 1;
    int Nt = (int) (total_time/dt) +1;
    //END of input variables

    double grid[Nz][Nx][Ny];
    double *TEMP_A;
    TEMP_A = malloc(sizeof(double) * Nt);
    for (int k = 0; k < Nz; k++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                grid[k][i][j] = Ta;
            }
        }
    }

    diffusion3d(Nz , Nx ,  Ny ,  Nt ,  grid,  dz ,   dx ,  dy ,  TEMP_A , electrode_vel , dt ,  sigma , Q ,  Ta ,  h);

    FILE *TEMP_POINT_A;

    TEMP_POINT_A = fopen("C:/Numerical_data/welding_sim/TEMP_POINT_A.dat","wb");

    for(int i = 0 ; i < Nt ; i++) {
        fprintf(TEMP_POINT_A,"%lf\n",TEMP_A[i]);
    }

    fclose(TEMP_POINT_A);
    free(grid);
    TIME_STOP = clock();
    printf("%lf [s] program time \n", (double) (TIME_STOP-TIME_START)/CLOCKS_PER_SEC);
    return 0;
}
