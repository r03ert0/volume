#include "swap.h"

void swap_short(short *v)
{
	unsigned char	b[2];
	
	b[0]=((unsigned char*)v)[1];
	b[1]=((unsigned char*)v)[0];
	*v=*(short*)b;
}
void swap_int(int *v)
{
	unsigned char	b[4];
	
	b[0]=((unsigned char*)v)[3];
	b[1]=((unsigned char*)v)[2];
	b[2]=((unsigned char*)v)[1];
	b[3]=((unsigned char*)v)[0];
	*v=*(int*)b;
}
void swap_float(float *v)
{
	unsigned char	b[4];
	
	b[0]=((unsigned char*)v)[3];
	b[1]=((unsigned char*)v)[2];
	b[2]=((unsigned char*)v)[1];
	b[3]=((unsigned char*)v)[0];
	*v=*(float*)b;
}
void swap_rgbfloat(float *v)
{
	unsigned char	b[12];
	int	i;
	
	b[0]=((unsigned char*)v)[3];
	b[1]=((unsigned char*)v)[2];
	b[2]=((unsigned char*)v)[1];
	b[3]=((unsigned char*)v)[0];
	
	b[4]=((unsigned char*)v)[7];
	b[5]=((unsigned char*)v)[6];
	b[6]=((unsigned char*)v)[5];
	b[7]=((unsigned char*)v)[4];
	
	b[8]=((unsigned char*)v)[11];
	b[9]=((unsigned char*)v)[10];
	b[10]=((unsigned char*)v)[9];
	b[11]=((unsigned char*)v)[8];
	for(i=0;i<12;i++)
		((unsigned char*)v)[i]=b[i];
}
