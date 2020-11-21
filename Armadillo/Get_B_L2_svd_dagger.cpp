#include "DQMC.h"
int Get_B_L2_svd_dagger(mat* U_out, mat* S_out, mat* V_out, int alpha, int L2, int L1, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene)
{
	int Bin_Size = 4;// larger the value, faster but less stable the calculation
	mat B_L2(NumOfVertexs, NumOfVertexs, fill::eye);
	mat U;
	vec s;
	mat V;
	svd(U, s, V, B_L2);
	mat S;
	S = diagmat(s);
	mat Binned_B_L(NumOfVertexs, NumOfVertexs, fill::eye);
	for (int Bin_index = 0; Bin_index < (int)((L2 - L1) / Bin_Size) + 1; Bin_index++)
	{
		int L_start = L1 + Bin_index * Bin_Size;
		if (Bin_index < (int)((L2 - L1) / Bin_Size))
		{
			for (int index = 0; index < Bin_Size; index++)
			{
				Binned_B_L = Get_B_L(alpha, L_start + index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene) * Binned_B_L;
			}
		}
		else
		{
			for (int index = 0; index < (L2 - L_start); index++)
			{
				Binned_B_L = Get_B_L(alpha, L_start + index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene) * Binned_B_L;
			}
		};
		mat U_new;
		vec s_new;
		mat V_new;
		svd(U_new, s_new, V_new, Binned_B_L);
		mat S_new;
		S_new = diagmat(s_new);

		mat U_L;
		vec s_L;
		mat V_L;
		svd(U_L, s_L, V_L, S_new * (trans(V_new) * U) * S);
		U = U_new * U_L;
		S = diagmat(s_L);
		V = V * V_L;
	};
	*U_out = U;
	*S_out = S;
	*V_out = V;
	return 0;
}