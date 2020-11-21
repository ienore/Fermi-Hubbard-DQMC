#include "DQMC.h"
mat Get_B_L(int alpha, int L, int NumOfVertexs, mat Sigma,double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene)
{
	mat B_L(NumOfVertexs, NumOfVertexs);
	mat Sigma_L = Sigma.row(L);
	mat V_L(NumOfVertexs, NumOfVertexs, fill::zeros);
	for (int i = 0; i < NumOfVertexs; i++)
	{
		V_L(i, i) = exp(lambda * alpha * Sigma_L[i]);
	};
	B_L = V_L * expmat(D_Tau * T_hop * K);
	return B_L;
}