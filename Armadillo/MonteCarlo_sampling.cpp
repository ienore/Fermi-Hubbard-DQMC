#include "DQMC.h"
mat MonteCarlo_sampling(int zjy_index, int N_wrap, mat Sigma_old, mat id_mat, int NumInEdge, int NumOfWarm, int NumOfEpoch, mat K, int TempSlice, int NumOfVertexs, double Miu, double Uene, double D_Tau, double lambda, double T_hop)
{
	FILE* fp;
	errno_t err;
	err = fopen_s(&fp, "DQMC.dat", "w");

	mat Sigma;
	Sigma = Sigma_old;
	mat green_L_up, green_L_down;
	mat green_up_c, green_down_c;
	mat B_trans_up, B_trans_inv_up, B_trans_down, B_trans_inv_down;
	mat delta_up, delta_down;
	double R_up, R_down, PosToChange;
	for (int warm_index = 0; warm_index < NumOfWarm; warm_index++)
	{
		if (zjy_index % 8 == 0)
		{
			printf("D_Tau=%f, MC Sampling Ratio = %f\n", D_Tau, warm_index / (double)NumOfWarm);
		};
		for (int reverse_sign = 0; reverse_sign < 1; reverse_sign++)
		{
			for (int time_index_mother = 0; time_index_mother < TempSlice; time_index_mother++)
			{
				int time_index;
				if (reverse_sign == 0)
				{
					time_index = time_index_mother;
					if (1)//(time_index_mother % N_wrap == 0 || time_index_mother == 0)
					{
						green_L_up = Get_G_L(1.0, time_index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						green_L_down = Get_G_L(-1.0, time_index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
					}
					else
					{
						//example: B_76 * G_5 * inv(B_76) = G_6, since B_76 = B_5
						B_trans_up = Get_B_L2(1.0, time_index + 1, time_index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						B_trans_inv_up = Get_B_L2_inv(1.0, time_index + 1, time_index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						B_trans_down = Get_B_L2(-1.0, time_index + 1, time_index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						B_trans_inv_down = Get_B_L2(-1.0, time_index + 1, time_index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						green_L_up = B_trans_up * green_L_up * B_trans_inv_up;
						green_L_down = B_trans_down * green_L_down * B_trans_inv_down;
					}
				}
				else
				{
					time_index = TempSlice - time_index_mother - 1;
					if (1)//(time_index_mother % N_wrap == 0 || time_index == TempSlice - 1)
					{
						green_L_up = Get_G_L(1.0, time_index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						green_L_down = Get_G_L(-1.0, time_index, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
					}
					else
					{
						//example: inv(B_76) * G_6 * B_76 = G_5
						B_trans_up = Get_B_L2(1.0, time_index + 2, time_index + 1, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						B_trans_inv_up = Get_B_L2_inv(1.0, time_index + 2, time_index + 1, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						B_trans_down = Get_B_L2(-1.0, time_index + 2, time_index + 1, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						B_trans_inv_down = Get_B_L2(-1.0, time_index + 2, time_index + 1, NumOfVertexs, Sigma, D_Tau, lambda, TempSlice, K, T_hop, Miu, Uene);
						green_L_up = B_trans_inv_up * green_L_up * B_trans_up;
						green_L_down = B_trans_inv_down * green_L_down * B_trans_down;
					}
				}//green_L_up/down is calculated

				green_up_c = id_mat - trans(green_L_up);
				green_down_c = id_mat - trans(green_L_down);
				//Measurement
				for (int site_index_print = 0; site_index_print < NumOfVertexs; site_index_print++)
				{
					double output_value;
					output_value = green_up_c(site_index_print, site_index_print) * green_down_c(site_index_print, site_index_print);
					fprintf(fp, "%6.5e\n", output_value);
				};


				int site_index;
				for (site_index = 0; site_index < NumOfVertexs; site_index++)
				{
					delta_up = zeros(NumOfVertexs, NumOfVertexs);
					delta_down = zeros(NumOfVertexs, NumOfVertexs);;
					delta_up(site_index, site_index) = exp(2.0 * lambda * Sigma(time_index, site_index)) - 1.0;
					delta_down(site_index, site_index) = exp(-2.0 * lambda * Sigma(time_index, site_index)) - 1.0;
					R_up = 1.0 + delta_up(site_index, site_index) * (1.0 - green_L_up(site_index, site_index));
					R_down = 1.0 + delta_down(site_index, site_index) * (1.0 - green_L_down(site_index, site_index));
					PosToChange = R_up * R_down;
					if (rand() / (double)RAND_MAX < PosToChange)
					{
						Sigma(time_index, site_index) = -Sigma(time_index, site_index);
						green_L_up = green_L_up - 1.0 / R_up * green_L_up * delta_up * (id_mat - green_L_up);
						green_L_down = green_L_down - 1.0 / R_down * green_L_down * delta_down * (id_mat - green_L_down);
					}
				}
			}
		}
	}
	return Sigma;
}