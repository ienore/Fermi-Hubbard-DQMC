#include "DQMC.h"

int main()
{
	int NumInEdge = 4;
	int NumOfVertexs = NumInEdge * NumInEdge * NumInEdge;
	double Beta = 4.0;
	double D_Tau = 0.1;
	double Uene = 0.0;
	double T_hop = 1.0;
	int alpha = 1;
	int L = 10;

	double Miu = Uene / 2.0;
	double lambda = 2.0 * atanh(sqrt(tanh(D_Tau * Uene / 4.0)));
	int TempSlice = Beta / (double)D_Tau;
	mat Sigma(TempSlice, NumOfVertexs, fill::zeros);
	Sigma = Sigma + 1;
	mat K(NumOfVertexs, TempSlice, fill::zeros);
	mat B_L(NumOfVertexs, TempSlice, fill::zeros);
	K = Get_K_3d(NumInEdge);
	mat G_L;
	G_L = Get_G_L(alpha,L,NumOfVertexs,Sigma, D_Tau, lambda,TempSlice,K, T_hop, Miu, Uene);
	mat B;
	mat B_old;
	B_old = Get_B_L2_inv(alpha, L+1, 0, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
	B = Get_B_L2_inv(alpha, 10, 0, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
	double sum = 0.0;
	mat sum_mat;
	sum_mat = G_L;




	for (int i = 0; i < NumOfVertexs; i++)
	{
		for (int j = 0; j < NumOfVertexs; j++)
		{
			sum = sum + sum_mat(i, j) * sum_mat(i, j);
		};
	};
	cout << sum << endl;
	return 0;
}