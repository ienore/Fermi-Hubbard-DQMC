#include "DQMC.h"
/// <summary>
/// 
/// Run this function with:
/// Bin(&bin_mean, &bin_stdvar, bin_num, sample, len_sample);
/// 
/// </summary>
/// <param name="bin_mean">output:mean value of the sample</param>
/// <param name="bin_stdvar">output: standard variance of the sample</param>
/// <param name="bin_num">slices of the sample, 10 is recommended</param>
/// <param name="sample">data</param>
/// <param name="len_sample">the length of data</param>
/// <returns></returns>
int Bin(double* bin_mean, double* bin_stdvar, int bin_num, double* sample, long len_sample)
{	
	int max_bin_num = bin_num;
	int bin_size = int(len_sample / bin_num);
	double* bin_data = (double*)malloc(sizeof(double) * bin_num);
	double total_mean = 0.0;
	for (int data_i = 0; data_i < bin_num; data_i++)
	{
		double sum = 0.0;
		for (int j = 0; j < bin_size; j++)
		{
			sum = sum + sample[data_i * bin_size + j];
		};
		bin_data[data_i] = sum / double(bin_size);
		total_mean = total_mean + bin_data[data_i] /(double) bin_num;
	};
	double var_sum = 0.0;
	for (int i = 0; i < bin_num; i++)
	{
		var_sum = var_sum + (bin_data[i] - total_mean) * (bin_data[i] - total_mean);
	};
	*bin_mean = total_mean;
	*bin_stdvar = sqrt(var_sum / (double)bin_num);
	return 0;
}