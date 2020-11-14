#include "DQMC.h"
mat Get_B_L2_TEST(int alpha, int L2, int L1, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene)
{
	mat B_L2(NumOfVertexs, NumOfVertexs, fill::eye);
	for (int i = L1; i < L2; i++)
	{
		mat B;
		B = Get_B_L(alpha, i, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
		B_L2 = B * B_L2;
	};
	return B_L2;
}