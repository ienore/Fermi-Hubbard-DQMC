#include "DQMC.h"
int IndexToCoor_3d(int* x, int* y, int* z, int index, int NumInEdge)
{
	*x = index % NumInEdge;
	*y = ((index - *x) / NumInEdge) % NumInEdge;
	*z = (index - *x - *y * NumInEdge) / (NumInEdge * NumInEdge);
	return 0;
}