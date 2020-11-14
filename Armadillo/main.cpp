#include "DQMC.h"

int main()
{
	//Set up Physical parameters
	int NumInEdge = 4;
	double Beta = 4.0;
	double D_Tau = 0.1;
	double Uene = 0.0;
	double T_hop = 1.0;

	//Set up Computational parameters
	int zjy_index = 0;
	int N_wrap = 5;
	int NumOfWarm = 100;
	int NumOfEpoch = 100;

	//Generates related parametres
	int NumOfVertexs = NumInEdge * NumInEdge * NumInEdge;
	double Miu = Uene / 2.0;
	double lambda = 2.0 * atanh(sqrt(tanh(D_Tau * Uene / 4.0)));
	int TempSlice = Beta / (double)D_Tau;
	mat Sigma(TempSlice, NumOfVertexs, fill::randu);
	for (int i = 0; i < TempSlice; i++)
	{
		for (int j = 0; j < NumOfVertexs; j++)
		{
			if (Sigma(i, j) < 0.5)
			{
				Sigma(i, j) = -1.0;
			}
			else
			{
				Sigma(i, j) = 1.0;
			};
		};
	};
	mat K(NumOfVertexs, NumOfVertexs, fill::zeros);
	K = Get_K_3d(NumInEdge);
	mat id_mat(NumOfVertexs, NumOfVertexs, fill::eye);
	Sigma = WarmUp(zjy_index, N_wrap, Sigma, id_mat, NumInEdge, NumOfWarm, NumOfEpoch, K, TempSlice, NumOfVertexs, Miu, Uene, D_Tau, lambda, T_hop);

	return 0;
}