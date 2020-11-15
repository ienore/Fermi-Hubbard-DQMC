#include "DQMC.h"
mat Get_B_L2_inv(int alpha, int L2, int L1, int NumOfVertexs, mat Sigma, double D_Tau, double lambda, int TempSlice, mat K, double T_hop, double Miu, double Uene)
{
	int Bin_Size = 4;// larger the value, faster but less stable the calculation
	mat B_L2(NumOfVertexs, NumOfVertexs, fill::eye);
	mat U;
	vec S_vec;
	mat V;
	svd(U, S_vec, V, B_L2);
	mat S;
	S = diagmat(S_vec);
	mat Binned_B_L(NumOfVertexs, NumOfVertexs, fill::eye);
	int i;
	for (int Bin_index = 0; Bin_index < (L2 - L1) / Bin_Size; Bin_index++)
	{
		Binned_B_L = eye(NumOfVertexs, NumOfVertexs);
		for (i = L1 + Bin_index * Bin_Size; i < L1 + (Bin_index + 1) * Bin_Size && i < L2; i++)
		{
			Binned_B_L = Binned_B_L*Get_B_L(alpha, i, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
		};
		if (Bin_index == (L2 - L1) / Bin_Size - 1 && i < L2)
		{
			for (int j = i; j < L2; j++)
			{
				Binned_B_L = Binned_B_L*Get_B_L(alpha, j, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
			}
		};
		mat U_new;
		vec S_new_vec;
		mat V_new;
		svd(U_new, S_new_vec, V_new, Binned_B_L);
		mat S_new;
		S_new = diagmat(S_new_vec);

		mat u;
		vec s_vec;
		mat v;
		svd(u, s_vec, v,S*(trans(V)*U_new)*S_new);
		U = U * u;
		S = diagmat(s_vec);
		V = V_new * v;

	};
	B_L2 = U * S * trans(V);
	return B_L2;
}