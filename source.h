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

#ifndef HEAT_DIFFUSION_SOURCE_H
#define HEAT_DIFFUSION_SOURCE_H

void diffusion3d(int Nz , int Nx , int Ny , int Nt , double grid[Nz][Nx][Ny], double dz ,  double dx , double dy , double *TEMP_A,double electrode_vel,double dt , double sigma ,double Q , double Ta , double h);
double  timestep( double weld_time , double v , double time, double sigma , double Q ,int Nz, int Nx , int Ny , double grid[Nz][Nx][Ny] ,double dz, double dx , double dy , double dt  , double grid_new[Nz][Nx][Ny]);
void copy_matrix(int rows , int cols , double original[][cols] , double copy[][cols]);

#endif //HEAT_DIFFUSION_SOURCE_H
