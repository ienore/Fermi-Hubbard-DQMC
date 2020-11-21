#include "DQMC.h"

mat Get_K_3d(int NumInEdge)
{
	int NumOfVertexs = NumInEdge * NumInEdge * NumInEdge;
	mat K(NumOfVertexs,NumOfVertexs, fill::zeros);
	for (int i = 0; i < NumOfVertexs; i++)
	{
		for (int j = 0; j < NumOfVertexs; j++)
		{
			int xi, yi, zi;
			int xj, yj, zj;
			IndexToCoor_3d(&xi, &yi, &zi, i, NumInEdge);
			IndexToCoor_3d(&xj, &yj, &zj, j, NumInEdge);
			int delta_x, delta_y, delta_z;
			delta_x = min(abs(xi - xj), NumInEdge - abs(xi - xj));
			delta_y = min(abs(yi - yj), NumInEdge - abs(yi - yj));
			delta_z = min(abs(zi - zj), NumInEdge - abs(zi - zj));
			if (delta_x + delta_y + delta_z == 1)
			{
				K(i, j) = 1;
			};
		};
	};
	return K;
}