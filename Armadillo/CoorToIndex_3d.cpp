#include "DQMC.h"
int CoorToIndex_3d(int x, int y, int z, int NumInEdge)
{
	int index = NumInEdge * NumInEdge * z + NumInEdge * y + x;
	return index;
}