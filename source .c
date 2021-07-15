//--------------------------------------------------------//
//--------------------------------------------------------//
//--------------------------------------------------------//
//  3D FINITE DIFFERENCE ANALYSIS OF ARC WELDING          //
//  National Technical University of Athens               //
//  School of Mechanical Engineering                      //
//  Manufacturing Technology Lab                          //
//                                                        //
// created by Dimitrios T. Pallas,                        //
// under the supervision of George C. Vosniakos (2021)    //
//                                                        //
//                                                        //
//                                                        //
//--------------------------------------------------------//
//--------------------------------------------------------//
//--------------------------------------------------------//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define rho 8050 //steel density , considered constant


void copy_matrix_3d(int blocks , int rows , int cols , double original[][rows][cols] , double copy[][rows][cols]){
    register int i ; register int j;register int k;
    for(k = 0 ; k < blocks ; k++) {
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                copy[k][i][j] = original[k][i][j];
            }
        }
    }
}

//The four following functions describe the material's thermal conductivity and specific heat , as well as their derivatives with respect to temperature.
//The functions are obtained throught analytical relations , created by experimental data.
double lambda( double theta) {
    if (theta >= 20 && theta < 800) {
        return 54 - 3.33 * (pow(10,-2)) * theta;
    }
    if (theta >= 800 ) {
        return 27.3;
    }
    printf("ERROR:TERMINATE PROGRAMM\n");
    exit(-1);
}

double lambda_prime(double theta){
    if(theta >= 20 && theta <800){
        return -3.33 * (pow(10,-2));
    }
    if(theta >= 800){
        return 0;
    }
    exit(-1);
}

double c( double theta) {

    if (theta >= 20 && theta < 600) {
        return 425 + 7.73 * (pow(10, -1)) * theta - 1.69 * (pow(10, -3)) * pow(theta, 2) +
               2.22 * (pow(10, -6)) * pow(theta, 3);
    }
    if (theta >= 600 && theta < 735) {
        return 666 + 13002 / (738 - theta);
    }

    if (theta >= 735 && theta < 900) {
        return 545 + 17820 / (theta - 731);
    }
    if (theta >= 900) {
        return 650;
    }
    printf("ERROR:TERMINATE PROGRAMM\n");
    exit(-1);
}

double c_prime( double theta) {

    if (theta >= 20 && theta < 600) {
        return  7.73 * (pow(10, -1))  - 1.69 * (pow(10, -3)) * theta +
                2.22 * (pow(10, -6)) * 3 * pow(theta,2);
    }
    if (theta >= 600 && theta < 735) {
        return 13002 / pow(738 - theta,2);
    }

    if (theta >= 735 && theta < 900) {
        return  17820 / pow(theta - 731,2);
    }
    if (theta >= 900) {
        return 0;
    }
    exit(-1);
}

double newton_boundary(double dz , double h , double Temp , double Ta , double k){
    //Since transfer through radiation is considered , the equation describing the ghost boundaries cannot be analytically solved (for the ghost boundary term)
    double Sigma   = 5.670374 * 1e-8; // Latest measurements
    double epsilon = 1;
    double dT = 1e-2;
    double GB_old = 1;
    double GB_new = 0;
    int count = 0;
    //newton iterations ...
    while(count <= 100){
        GB_new = GB_old - (pow(GB_old-Ta,4)*epsilon*Sigma + GB_old*(h + lambda(Temp)/dz)  - Ta*h - lambda(Temp)/dz*Temp)/(4*epsilon*Sigma*pow(GB_old,3) + h + lambda(Temp)/dz);
        if(fabs(GB_new - GB_old)/fabs(GB_new) <= 0.1){
            GB_old = GB_new;
            break ;
        }
        GB_old = GB_new;
        count++;
    }
    return GB_old;
}


double  timestep( double weld_time , double v , double time, double sigma ,double h,double Ta ,double Q ,int Nz, int Nx , int Ny , double grid[Nz][Nx][Ny]  ,double dz,double dx , double dy , double dt  , double grid_new[Nz][Nx][Ny]) {

    double GBZ1 , GBZ2 , GBX1 , GBX2 , GBY1 , GBY2;
    double TEMP_A;

    for (int k = 0; k < Nz ; k++) {

        for (int i = 0; i < Nx ; i++) {

            for (int j = 0; j < Ny ; j++) {

                //Setting up B.C.s , using ghost boundaries. Heat transfer through convection and radiation is considered.
                if (k == 0) {
                    //GBZ2 = (Ta * h * dz + lambda(grid[k][i][j]) * grid[k][i][j]) / (dz * h + lambda(grid[k][i][j]));
                    GBZ2 = newton_boundary(dz,h,grid[k][i][j],Ta,lambda(grid[k][i][j]));
                    GBZ1 = grid[k + 1][i][j];

                }
                else if  (k == Nz - 1) {
                    //GBZ1 = (Ta * h * dz + lambda(grid[k][i][j]) * grid[k][i][j]) / (dz * h + lambda(grid[k][i][j]));
                    GBZ1 = newton_boundary(dz,h,grid[k][i][j],Ta,lambda(grid[k][i][j]));
                    GBZ2 = grid[k - 1][i][j];
                } else {
                    GBZ1 = grid[k + 1][i][j];
                    GBZ2 = grid[k - 1][i][j];
                }
                if (i == 0) {
                    //GBX2 = (Ta * h * dx + lambda(grid[k][i][j]) * grid[k][i][j]) / (dx * h + lambda(grid[k][i][j]));
                    GBX2 = newton_boundary(dx,h,grid[k][i][j],Ta,lambda(grid[k][i][j]));
                    GBX1 = grid[k][i+1][j];
                }
                else if(i == Nx-1) {
                    //GBX1 = (Ta * h * dx + lambda(grid[k][i][j]) * grid[k][i][j]) / (dx * h + lambda(grid[k][i][j]));
                    GBX1 = newton_boundary(dx,h,grid[k][i][j],Ta,lambda(grid[k][i][j]));
                    GBX2 = grid[k][i-1][j];
                } else {
                    GBX1 = grid[k][i+1][j];
                    GBX2 = grid[k][i-1][j];
                }
                if (j == 0) {
                    //GBY2 = (Ta * h * dy + lambda(grid[k][i][j]) * grid[k][i][j]) / (dy * h + lambda(grid[k][i][j]));
                    GBY2 = newton_boundary(dy,h,grid[k][i][j],Ta,lambda(grid[k][i][j]));
                    GBY1 = grid[k][i][j+1];
                }
                else if(j == Ny-1) {
                    //GBY1 = (Ta * h * dy + lambda(grid[k][i][j]) * grid[k][i][j]) / (dy * h + lambda(grid[k][i][j]));
                    GBY1 = newton_boundary(dy,h,grid[k][i][j],Ta,lambda(grid[k][i][j]));
                    GBY2 = grid[k][i][j+1];
                } else {
                    GBY1 = grid[k][i][j+1];
                    GBY2 = grid[k][i][j-1];
                }
                //Discretised Heat Transfer Eq. , material properties are not considered constant
                grid_new[k][i][j] = dt / (rho * (c(grid[k][i][j]) + c_prime(grid[k][i][j]))) *
                                    (lambda(grid[k][i][j]) *
                                     ((GBX1 - 2 * grid[k][i][j] + GBX2) / (dx * dx) +
                                      (GBY1 - 2 * grid[k][i][j] + GBY2) / (dy * dy) +
                                      (GBZ1 - 2 * grid[k][i][j] + GBZ2) / (dz * dz)) +
                                     lambda_prime(grid[k][i][j]) *
                                     ((GBX1 - GBX2) / (2 * dx) *(GBX1 - GBX2) /(2 * dx) +
                                      (GBY1 - GBY2) / (2 * dy) *(GBY1 - GBY2) /(2 * dy) +
                                      (GBZ1 - GBZ2) / (2 * dz) *(GBZ1 - GBZ2) / (2 * dz))) + grid[k][i][j];

                //Source term is only required when the electrodes add energy to the system
                if (time <= weld_time) {
                    grid_new[k][i][j] += dt / (rho * (c(grid[0][i][j]) + c_prime(grid[0][i][j]))) * (Q / v *exp(-(pow(v * time -i * dx,2)
                            +pow(25e-3 -j * dy,2) + pow(k*dz,2)) /(2 * sigma *sigma)) /(2 * sigma *sigma *M_PI));
                }
                //Newmann (second type) Boundary Condition
                //ADIABATIC POINT ( Q = 0 , no heat transfer)
                //Required for numerical stability
                //Only used at point (Nx-1,Ny-1,Nz-1)
                if(i== Nx-1 && j ==Ny-1 && k == Nz-1 ){
                    grid_new[k][i][j] = grid[k][i][j];

                }
            }
        }
    }
    for (int k = 0; k < Nz; k++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                //This constraint is required , since we only know the material properties up to 1200 Celsius.
                if (grid_new[k][i][j] >= 1500) {
                    grid_new[k][i][j] = 1500;
                }
                if (grid_new[k][i][j] < 20) {
                    printf("ERROR : BAD TEMP @ %d,%d,%d\n", i, j, k);
                    exit(-1);
                }
            }
        }
    }

    TEMP_A = grid_new[0][20][25];

    return TEMP_A;
}

void diffusion3d(int Nz , int Nx , int Ny , int Nt , double grid[Nz][Nx][Ny], double dz ,  double dx , double dy , double *TEMP_A,double electrode_vel,double dt , double sigma ,double Q , double Ta , double h){


    double weld_time = 40e-3/(electrode_vel);
    double current_temp ;
    current_temp = Ta;

    double grid_new[Nz][Nx][Ny];

    FILE *grid_to_print , *paths;
    paths = fopen("C:/numerical_data/welding_sim/paths.dat","wb");
    char path1[] = "C:/Numerical_data/welding_sim/field_01";
    char path2[] = "C:/Numerical_data/welding_sim/field_02";
    char path3[] = "C:/Numerical_data/welding_sim/field_03";
    char path4[] = "C:/Numerical_data/welding_sim/field_04";
    char loc[] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','J','K','L','Z','X','V','N','M','Q'};
    int count = 0;
    int snap = 10 ;

    for(int t = 0 ; t <= Nt  ; t++) {


        TEMP_A[t] = current_temp;
        current_temp = timestep(weld_time,electrode_vel,dt * t, sigma, h , Ta ,Q, Nz , Nx, Ny, grid, dz , dx, dy, dt, grid_new);
        copy_matrix_3d(Nz,Nx, Ny, grid_new, grid);

        if(t%snap==0) {
            if(count <= 25){
                path1[36] = loc[count];
                fprintf(paths,"%s\n",path1);
                grid_to_print = fopen(path1, "wb");
            }
            if(count <= 50 && count > 25){
                path2[36] = loc[count-25];
                fprintf(paths,"%s\n",path2);
                grid_to_print = fopen(path2, "wb");
            }
            if(count <= 75 && count > 50){
                path3[36] = loc[count-50];
                fprintf(paths,"%s\n",path3);
                grid_to_print = fopen(path3, "wb");
            }
            if(count <= 100 && count > 75){
                path4[36] = loc[count-75];
                fprintf(paths,"%s\n",path4);
                grid_to_print = fopen(path4, "wb");
            }
            count++;
//            for (int i = 0; i < Nx; i++) {
//                for (int j = 0; j < Ny; j++) {
//                    fprintf(grid_to_print, "%lf ", grid[0][i][j]);
//                }
//                fprintf(grid_to_print, "\n");
//            }
//            fclose(grid_to_print);
            printf("Time %lf -- T=%lf\n",t*dt,current_temp);

        }

    }
    fclose(paths);
    //free(grid_new);
}