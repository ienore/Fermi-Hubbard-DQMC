#include "DQMC.h"

mat Get_G_L(int alpha, int L, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene)
{
	mat U_R, D_R, V_R;
	mat V_L, D_L, U_L;
	Get_B_L2_svd(&U_R, &D_R, &V_R, alpha, L+1, 0, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
	Get_B_L2_svd(&V_L, &D_L, &U_L, alpha, TempSlice, L+1, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
	V_R = trans(V_R);
	U_L = trans(U_L);
	mat D_R_min(NumOfVertexs, NumOfVertexs, fill::eye);
	mat D_R_max(NumOfVertexs, NumOfVertexs, fill::eye);
	mat D_R_max_inv(NumOfVertexs, NumOfVertexs, fill::eye);
	mat D_L_min(NumOfVertexs, NumOfVertexs, fill::eye);
	mat D_L_max(NumOfVertexs, NumOfVertexs, fill::eye);
	mat D_L_max_inv(NumOfVertexs, NumOfVertexs, fill::eye);
	for (int i = 0; i < NumOfVertexs; i++)
	{
		if (D_R(i, i) > 1)
		{
			D_R_max_inv(i, i) = 1.0 / D_R(i, i);
			D_R_max(i, i) = D_R(i, i);
		}
		else
		{
			D_R_min(i, i) = D_R(i, i);
		}
		if (D_L(i, i) > 1)
		{
			D_L_max_inv(i, i) = 1.0 / D_L(i, i);
			D_L_max(i, i) = D_L(i, i);
		}
		else
		{
			D_L_min(i, i) = D_L(i, i);
		}
	};

	mat G_L;
	mat temp;
	temp = D_R_max_inv * inv(U_L * U_R) * D_L_max_inv + D_R_min * V_R * V_L * D_L_min;
	G_L = inv(U_L) * D_L_max_inv * inv(temp) * D_R_max_inv * inv(U_R);

	return G_L;
}