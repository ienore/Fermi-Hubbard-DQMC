#pragma once
#include <iostream>
#include <armadillo>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
using namespace std;
using namespace arma;


int Bin(double* bin_mean, double* bin_stdvar, int bin_num, double* sample, long len_sample);
int IndexToCoor_3d(int* x, int* y, int* z, int index, int NumInEdge);
mat Get_K_3d(int NumInEdge);
int CoorToIndex_3d(int x, int y, int z, int NumInEdge);
mat Get_B_L(int alpha, int L, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene);
mat Get_B_L_inv(int alpha, int L, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene);
mat Get_B_L2(int alpha, int L2, int L1, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene);
mat Get_B_L2_TEST(int alpha, int L2, int L1, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene);
mat Get_B_L2_inv(int alpha, int L2, int L1, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene);
int Get_B_L2_svd(mat* U_out, mat* S_out, mat* V_out, int alpha, int L2, int L1, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene);
mat Get_G_L(int alpha, int L, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene);
mat WarmUp(int zjy_index, int N_wrap, mat Sigma, mat id_mat, int NumInEdge, int NumOfWarm, int NumOfEpoch, mat K, int TempSlice, int NumOfVertexs, double Miu, double Uene, double D_Tau, double lambda, double T_hop);


