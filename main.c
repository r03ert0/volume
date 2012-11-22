char	version[]="volume, v2, roberto toro, 10 December 2010";	// added -decompose

#include <stdio.h>
#include "Analyze.h"
#include "MGH.h"
#include "math.h"
#include <unistd.h>

#define MAX(x,y) (((x)>(y))?(x):(y))

#define kAnalyzeVolume	1
#define kMGZVolume		2

AnalyzeHeader	*hdr;
char			*img;
int				dim[4];
int				verbose=1;

void threshold(float value, int direction);

float getValue(int x, int y, int z)
{
	float		val;
	RGBValue	rgb;

	if(hdr->datatype==RGB)
	{
		rgb=((RGBValue*)img)[z*dim[1]*dim[0]+y*dim[0]+x];
		val=((int)rgb.r)>>16|((int)rgb.g)>>8|((int)rgb.b);
	}
	else
	{
		switch(hdr->datatype)
		{	case UCHAR: val=((unsigned char*)img)[z*dim[1]*dim[0]+y*dim[0]+x];	break;
			case SHORT: val=        ((short*)img)[z*dim[1]*dim[0]+y*dim[0]+x];	break;
			case INT:	val=          ((int*)img)[z*dim[1]*dim[0]+y*dim[0]+x];	break;
			case FLOAT:	val=        ((float*)img)[z*dim[1]*dim[0]+y*dim[0]+x];	break;
		}
	}
	return val;
}

void setValue(float val, int x, int y, int z)
{
	switch(hdr->datatype)
	{	case UCHAR: ((unsigned char*)img)[z*dim[1]*dim[0]+y*dim[0]+x]=val;			break;
		case SHORT: ((short*)img)[z*dim[1]*dim[0]+y*dim[0]+x]=val;					break;
		case INT:	((int*)img)[z*dim[1]*dim[0]+y*dim[0]+x]=val;					break;
		case FLOAT:	((float*)img)[z*dim[1]*dim[0]+y*dim[0]+x]=val;					break;
	}
}
float getValue2(int i, AnalyzeHeader *theHdr, char *theImg)
{
	float		val;
	RGBValue	rgb;
	
	if(theHdr->datatype==RGB)
	{
		rgb=((RGBValue*)theImg)[i];
		val=((int)rgb.r)>>16|((int)rgb.g)>>8|((int)rgb.b);
	}
	else
	{
		switch(theHdr->datatype)
		{	case UCHAR: val=((unsigned char*)theImg)[i];	break;
			case SHORT: val=        ((short*)theImg)[i];	break;
			case INT:	val=          ((int*)theImg)[i];	break;
			case FLOAT:	val=        ((float*)theImg)[i];	break;
		}
	}
	return val;
}

void setValue2(float val, int i, AnalyzeHeader *theHdr, char *theImg)
{
	switch(theHdr->datatype)
	{	case UCHAR: ((unsigned char*)theImg)[i]=val;	break;
		case SHORT: ((short*)theImg)[i]=val;			break;
		case INT:	((int*)theImg)[i]=val;				break;
		case FLOAT:	((float*)theImg)[i]=val;			break;
	}
}
#pragma mark -
#pragma mark [ Utilities ]
int     endianness;
#define kMOTOROLA   1
#define kINTEL      2
void checkEndianness(void)
{
    char    b[]={1,0,0,0};
    int     num=*(int*)b;
    
    if(num==16777216)
        endianness=kMOTOROLA;
    else
        endianness=kINTEL;
}
#pragma mark -
#define STACK 600000
int connected(int x, int y, int z, int label)
{
	int		n,i,j,k;
	int		dx,dy,dz;
	int		dire[6]={0,0,0,0,0,0};
	int		nstack=0;
	short	stack[STACK][3];
	char	dstack[STACK];
	int		empty,mark;
	char	*tmp;

	tmp=(char*)calloc(dim[0]*dim[1]*dim[2],sizeof(char));

	n=0;
	for(;;)
	{
		tmp[z*dim[1]*dim[0]+y*dim[0]+x]=1;
		n++;
		
		for(k=0;k<6;k++)
		{
			dx = (-1)*(	k==3 ) + ( 1)*( k==1 );
			dy = (-1)*( k==2 ) + ( 1)*( k==0 );
			dz = (-1)*( k==5 ) + ( 1)*( k==4 );

			if(	z+dz>=0 && z+dz<dim[2] &&
				y+dy>=0 && y+dy<dim[1] &&
				x+dx>=0 && x+dx<dim[0] )
			{
				empty=(getValue(x+dx,y+dy,z+dz)==1);
				mark = tmp[(z+dz)*dim[1]*dim[0]+(y+dy)*dim[0]+(x+dx)];
				if( empty && !mark && !dire[k])
				{
					dire[k] = 1;
					dstack[nstack] = k+1;
					stack[nstack][0]=x+dx;
					stack[nstack][1]=y+dy;
					stack[nstack][2]=z+dz;
					nstack++;
				}
				else
					if( dire[k] && !empty && !mark)
						dire[k] = 0;
			}
		}
		if(nstack>STACK-10)
		{
			printf("StackOverflow\n");
			break;
		}
		if(nstack)
		{
			nstack--;
			x = stack[nstack][0];
			y = stack[nstack][1];
			z = stack[nstack][2];
			for(k=0;k<6;k++)
				dire[k] = (dstack[nstack]==(k+1))?0:dire[k];
		}
		else
			break;
	}
	
	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	if(tmp[k*dim[1]*dim[0]+j*dim[0]+i]==1)
		setValue(label,i,j,k);
	free(tmp);
	
	return n;
}
void largestConnected(void)
{
	int	label;
	int	i,j,k,n;
	int	max,labelmax;
	
	threshold(1,1); // if vox>=1 then vox=1, else vox=0
	
	max=0;
	label=2;
	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	if(getValue(i,j,k)==1)
	{
		n=connected(i,j,k,label);
		if(n>max)
		{
			max=n;
			labelmax=label;
		}
		label++;
	}

	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	if(getValue(i,j,k)==labelmax)
		setValue(1,i,j,k);
	else
		setValue(0,i,j,k);
}
void dilate(int r)
{
	int		i,j,k;
	int		a,b,c;
	char	*tmp=calloc(dim[0]*dim[1]*dim[2],1);
	
	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
		tmp[k*dim[1]*dim[0]+j*dim[0]+i]=(getValue(i,j,k)>0);
	
	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	if(getValue(i,j,k))
	{
		for(a=-r;a<=r;a++)
		for(b=-r;b<=r;b++)
		for(c=-r;c<=r;c++)
		if(	i+a>=0 && i+a<dim[0] &&
			j+b>=0 && j+b<dim[1] &&
			k+c>=0 && k+c<dim[2])
		if(a*a+b*b+c*c<=r*r)
		if(getValue(i+a,j+b,k+c)<1)
			tmp[(k+c)*dim[1]*dim[0]+(j+b)*dim[0]+(i+a)]=1;
	}

	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	if(getValue(i,j,k)<1 && tmp[k*dim[1]*dim[0]+j*dim[0]+i]>0)
		setValue(1,i,j,k);
	free(tmp);
}
void erode(int r)
{
	int		i,j,k;
	int		a,b,c;
	char	*tmp=calloc(dim[0]*dim[1]*dim[2],1);
	
	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
		tmp[k*dim[1]*dim[0]+j*dim[0]+i]=(getValue(i,j,k)>0);
	
	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	if(getValue(i,j,k)<1)
	{
		for(a=-r;a<=r;a++)
		for(b=-r;b<=r;b++)
		for(c=-r;c<=r;c++)
		if(	i+a>=0 && i+a<dim[0] &&
			j+b>=0 && j+b<dim[1] &&
			k+c>=0 && k+c<dim[2])
		if(a*a+b*b+c*c<=r*r)
		if(getValue(i+a,j+b,k+c)>0)
			tmp[(k+c)*dim[1]*dim[0]+(j+b)*dim[0]+(i+a)]=0;
	}

	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	if(getValue(i,j,k)>0 && tmp[k*dim[1]*dim[0]+j*dim[0]+i]==0)
		setValue(0,i,j,k);
	free(tmp);
}
float max(void)
{
	int		i,j,k;
	float	val,max;
	
	max=getValue(0,0,0);
	for(i=0;i<dim[0];i++)
		for(j=0;j<dim[1];j++)
			for(k=0;k<dim[2];k++)
			{
				val=getValue(i,j,k);
				if(val>max)
					max=val;
			}
	return max;
}
float min(void)
{
	int		i,j,k;
	float	val,min;
	
	min=getValue(0,0,0);
	for(i=0;i<dim[0];i++)
		for(j=0;j<dim[1];j++)
			for(k=0;k<dim[2];k++)
			{
				val=getValue(i,j,k);
				if(val<min)
					min=val;
			}
	return min;
}
void hist(int nbins)
{
	int		i,j,k;
	float	mi,ma;
	float	*hist;
	
	hist=(float*)calloc(nbins,sizeof(float));
	
	mi=min();
	ma=max();
	
	for(i=0;i<dim[0];i++)
		for(j=0;j<dim[1];j++)
			for(k=0;k<dim[2];k++)
				hist[(int)((nbins-1)*(getValue(i,j,k)-mi)/(ma-mi))]++;
	
	for(i=0;i<nbins;i++)
	{
		printf("%g",hist[i]);
		if(i<nbins-1)
			printf(" ");
	}
	printf("\n");
	free(hist);
}
float mean(void)
{
	int		i,j,k;
	float	sum=0;

	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
		sum+=getValue(i,j,k);
	return sum/(float)(dim[0]*dim[1]*dim[2]);
}
float std(void)
{
	int		i,j,k;
	float	val;
	float	s0=dim[0]*dim[1]*dim[2];
	float	s1=0;
	float	s2=0;

	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	{
		val=getValue(i,j,k);
		s1+=val;
		s2+=val*val;
	}
	return sqrt((s0*s2-s1*s1)/(s0*(s0-1)));
}
void threshold(float value, int direction)
{
	int		i,j,k;
	float	val;

	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	{
		val=getValue(i,j,k);
		if((val>=value && direction==1)||(val<=value && direction==0))
			setValue(1,i,j,k);
	}
}
void tiff(char *path, float *m, int W, int H, char *cmapindex)
{
	FILE	*f;
	int		i,n,S;
	unsigned char	hdr[]={	0x4d,0x4d,
		0x00,0x2a,
		0x00,0x00,0x00,0x08,
		0x00,0x0d,														// declare 13 fields
		0x01,0xfe, 0x00,0x04, 0x00,0x00,0x00,0x01, 0x00,0x00,0x00,0x00,
		0x01,0x00, 0x00,0x03, 0x00,0x00,0x00,0x01, 0x00,0x40,0x00,0x00,	// image width
		0x01,0x01, 0x00,0x03, 0x00,0x00,0x00,0x01, 0x00,0x20,0x00,0x00,	// image length
		0x01,0x02, 0x00,0x03, 0x00,0x00,0x00,0x03, 0x00,0x00,0x00,0xaa,	// bits per sample [addr: aa]
		0x01,0x03, 0x00,0x03, 0x00,0x00,0x00,0x01, 0x00,0x01,0x00,0x00,	// compression
		0x01,0x06, 0x00,0x03, 0x00,0x00,0x00,0x01, 0x00,0x02,0x00,0x00,	// photometric interpretation
		0x01,0x11, 0x00,0x04, 0x00,0x00,0x00,0x01, 0x00,0x00,0x00,0xc0,	// strip offsets [addr:0xc0]
		0x01,0x15, 0x00,0x03, 0x00,0x00,0x00,0x01, 0x00,0x03,0x00,0x00,	// samples per pixel
		0x01,0x16, 0x00,0x03, 0x00,0x00,0x00,0x01, 0x00,0x20,0x00,0x00,	// rows per strip
		0x01,0x17, 0x00,0x04, 0x00,0x00,0x00,0x01, 0x00,0x00,0x18,0x00,	// strip byte counts
		0x01,0x1a, 0x00,0x05, 0x00,0x00,0x00,0x01, 0x00,0x00,0x00,0xb0,	// x resolution [addr: b0]
		0x01,0x1b, 0x00,0x05, 0x00,0x00,0x00,0x01, 0x00,0x00,0x00,0xb8,	// y resolution [addr: b8]
		0x01,0x28, 0x00,0x03, 0x00,0x00,0x00,0x01, 0x00,0x02,0x00,0x00,	// resolution unit
		0x00,0x00,0x00,0x00,
		0x00,0x08, 0x00,0x08, 0x00,0x08,								// [addr 0xAA] bits per sample: 8,8,8
		0x00,0x0a,0xfc,0x80, 0x00,0x00,0x27,0x10,						// [addr 0xB0] x resolution 72
		0x00,0x0a,0xfc,0x80, 0x00,0x00,0x27,0x10,						// [addr 0xB8] y resolution 72
	};
	// image width
	hdr[31]=((unsigned char *)&W)[0];
	hdr[30]=((unsigned char *)&W)[1];
	hdr[33]=((unsigned char *)&W)[2];
	hdr[32]=((unsigned char *)&W)[3];
	
	// image length
	hdr[43]=((unsigned char *)&H)[0];
	hdr[42]=((unsigned char *)&H)[1];
	hdr[45]=((unsigned char *)&H)[2];
	hdr[44]=((unsigned char *)&H)[3];
	
	// rows per strip = image length
	hdr[115]=((unsigned char *)&H)[0];
	hdr[114]=((unsigned char *)&H)[1];
	hdr[117]=((unsigned char *)&H)[2];
	hdr[116]=((unsigned char *)&H)[3];
	
	// strip byte counts = image length * image width * 3
	S=W*H*3;
	hdr[129]=((unsigned char *)&S)[0];
	hdr[128]=((unsigned char *)&S)[1];
	hdr[127]=((unsigned char *)&S)[2];
	hdr[126]=((unsigned char *)&S)[3];
	
	
	// make a colour map
	unsigned char cmap[256][3];
	cmap[0][0]=cmap[0][1]=cmap[0][2]=0;
	srand(0);
	for(i=1;i<256;i++)
	{
		cmap[i][0]=rand()%256;
		cmap[i][1]=rand()%256;
		cmap[i][2]=rand()%256;
	}
	
	f=fopen(path,"w");
	n=sizeof(hdr);
	for(i=0;i<n;i++)
		fputc(hdr[i],f);
	for(i=0;i<W*H;i++)
	{
		if(strcmp(cmapindex,"grey")==0)
		{
			fputc((int)(m[i]),f);
			fputc((int)(m[i]),f);
			fputc((int)(m[i]),f);
		}
		else
		if(strcmp(cmapindex,"lut")==0)
		{
			fputc(cmap[(int)(m[i])][0],f);
			fputc(cmap[(int)(m[i])][1],f);
			fputc(cmap[(int)(m[i])][2],f);
		}

	}
	fclose(f);
	
}
float volume(void)
{
	int		i,j,k;
	int		nv;
	float	vvol,vol;
	
	vvol=hdr->pixdim[1]*hdr->pixdim[2]*hdr->pixdim[3];
	nv=0;
	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	if(getValue(i,j,k))
		nv++;
	
	vol=fabs(nv*vvol);
	
	return vol;
}
void decompose(char *basename)
{
	int				i,j,k;
	int				x,y,z;
	float			val;
	AnalyzeHeader	*hdr;
	char			*addr,*img;
	char			path[1024];
	char			sizeof_hdr[4]={92,1,0,0};
	
	// create empty UCHAR volume
	addr=calloc(dim[0]*dim[1]*dim[2]+sizeof(AnalyzeHeader),sizeof(char));
	hdr=(AnalyzeHeader*)addr;
	img=addr+sizeof(AnalyzeHeader);
	hdr->sizeof_hdr=*(int*)sizeof_hdr;
	hdr->datatype=UCHAR;
	hdr->dim[0]=3;
	hdr->dim[1]=dim[0];
	hdr->dim[2]=dim[1];
	hdr->dim[3]=dim[2];
	hdr->pixdim[1]=4;
	hdr->pixdim[2]=4;
	hdr->pixdim[3]=4;
	
	// decompose volume
	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	{
		val=getValue(i,j,k);
		if(val!=0)
		{
			for(x=0;x<dim[0];x++)
			for(y=0;y<dim[1];y++)
			for(z=0;z<dim[2];z++)
			if(getValue(x,y,z)==val)
			{
				img[z*dim[1]*dim[0]+y*dim[0]+x]=1;
				setValue(0,x,y,z);
			}
			else
				img[z*dim[1]*dim[0]+y*dim[0]+x]=0;

			sprintf(path,"%s_%g.hdr",basename,val);
			Analyze_save_hdr(path,*hdr);
			sprintf(path,"%s_%g.img",basename,val);
			Analyze_save_img(path,*hdr,img);
		}
	}
}
void matchHistogram(char *volpath)
{
	int				i,j,n,n0;
	int				sz;
	int				swapped;
	char			*addr;
	AnalyzeHeader	*hdr0;
	char			*img0;
	int				*hist,*hist0;
	float			*cumsum0,cumsum,prev;
	float			val;
	float			min,min0;
	float			max,max0;
	float			d,t,*trans;
	
	// load target volume
	Analyze_load(volpath,&addr,&sz,&swapped);
	hdr0=(AnalyzeHeader*)addr;
	img0=(char*)(addr+sizeof(AnalyzeHeader));
	// find max grey value in target volume
	n0=hdr0->dim[1]*hdr0->dim[2]*hdr0->dim[3];
	min0=max0=getValue2(0,hdr0,img0);
	for(i=0;i<n0;i++)
	{
		val=getValue2(i,hdr0,img0);
		if(val>max0)
			max0=val;
		if(val<min0)
			min0=val;
	}
	// compute histogram for target volume
	hist0=(int*)calloc(max0-min0+1,sizeof(int));
	for(i=0;i<n0;i++)
		hist0[(int)(getValue2(i,hdr0,img0)-min0)]++;
	// compute cummulative histogram for target colume
	cumsum0=(float*)calloc(max0-min0+1,sizeof(float));
	cumsum0[0]=hist0[0];///(float)n0;
	for(i=1;i<max0-min0+1;i++)
		cumsum0[i]=cumsum0[i-1]+hist0[i];
	for(i=0;i<max0-min0+1;i++)
		cumsum0[i]/=(float)n0;
	
	// find max grey value in source volume
	n=hdr->dim[1]*hdr->dim[2]*hdr->dim[3];
	min=max=getValue2(0,hdr,img);
	for(i=0;i<n;i++)
	{
		val=getValue2(i,hdr,img);
		if(val>max)
			max=val;
		if(val<min)
			min=val;
	}
	// compute histogram for source volume
	hist=(int*)calloc(max-min+1,sizeof(int));
	for(i=0;i<n;i++)
		hist[(int)(getValue2(i,hdr,img)-min)]++;
	
	// compute histogram matching transformation
	trans=(float*)calloc(max,sizeof(float));
	cumsum=0;
	prev=0;
	j=0;
	for(i=0;i<max-min+1;i++)
	{
		cumsum+=hist[i];
		
		// find grey level in target volume where cumsum0==cumsum
		while(cumsum0[j]<cumsum/(float)n)
		{
			prev=cumsum0[j];
			j++;
		}
		d=cumsum0[j]-prev;
		if(d)
			t=(cumsum/(float)n-prev)/d;
		else
			t=0;
		
		val=MAX(j-1,0)*(1-t)+j*t+min0;
		
		if(val>max0)
			printf("ça va pas, non?\n");
		
		trans[i+(int)min]=val;
		
	}
	
	// apply histogram matching transformation
	for(i=0;i<n;i++)
		setValue2(trans[(int)getValue2(i,hdr,img)],i,hdr,img);
}
#pragma mark -
#pragma mark [ Format conversion ]
int getformatindex(char *path)
{
    char    *formats[]={"hdr","img","mgz"};
    int     i,j,n=3; // number of recognised formats
    int     found,index;
    char    *extension;
    
    for(i=strlen(path);i>=0;i--)
        if(path[i]=='.')
            break;
    if(i==0)
    {
        printf("ERROR: Unable to find the format extension\n");
        return 0;
    }
    extension=path+i+1;
    
    for(i=0;i<n;i++)
    {
        j=strlen(formats[i]);
        found=(strcmp(formats[i],extension)==0);
        if(found)
            break;
    }
    
    index=-1;
    if(i==0 || i==1)
    {
        index=kAnalyzeVolume;
        if(verbose)
            printf("Format: Analyze volume\n");
    }
    else
	if(i==2)
	{
		index=kMGZVolume;
		if(verbose)
			printf("Format: MGZ volume\n");
	}
	else
	{
		printf("ERROR: unrecognized format \"%s\"\n",extension);
		return 0;
	}
	
    return index;
}
int loadVolume_Analyze(char *path)
{
	int		i,sz;
	int		swapped;
	char	*addr;
	char	base[512];
	
	strcpy(base,path);
	for(i=strlen(base);i>=0;i--)
		if(base[i]=='.')
			base[i]=(char)0;
	sprintf(path,"%s.hdr",base);
	Analyze_load(path,&addr,&sz,&swapped);
	hdr=(AnalyzeHeader*)addr;
	img=(char*)(addr+sizeof(AnalyzeHeader));
	dim[0]=hdr->dim[1];
	dim[1]=hdr->dim[2];
	dim[2]=hdr->dim[3];
	
	return 0;
}
int saveVolume_Analyze(char *path)
{
	int		i;
	char	base[512];
	
	strcpy(base,path);
	for(i=strlen(base);i>=0;i--)
		if(base[i]=='.')
			base[i]=(char)0;
	sprintf(path,"%s.hdr",base);
	Analyze_save_hdr(path,*hdr);
	sprintf(path,"%s.img",base);
	Analyze_save_img(path,*hdr,img);
	
	return 0;
}
int loadVolume_MGZ(char *path)
{
	char	*addr;
	int		sz;
	
	MGH_load_GZ(path,&addr,&sz);
	hdr=(AnalyzeHeader*)addr;
	img=(char*)(addr+sizeof(AnalyzeHeader));
	dim[0]=hdr->dim[1];
	dim[1]=hdr->dim[2];
	dim[2]=hdr->dim[3];

	return 0;
}
int saveVolume_MGZ(char *path)
{
	return 0;
}
int loadVolume(char *path)
{
	int	err,format;
	
	format=getformatindex(path);
	
	switch(format)
	{
		case kAnalyzeVolume:
			err=loadVolume_Analyze(path);
			break;
		case kMGZVolume:
			err=loadVolume_MGZ(path);
			break;
		default:
			printf("ERROR: Input mesh format not recognised\n");
			return 1;			
	}
	if(err!=0)
	{
		printf("ERROR: cannot read file: %s\n",path);
		return 1;
	}
	return 0;
}
int saveVolume(char *path)
{
	int	format;
	
	format=getformatindex(path);
	
	switch(format)
	{
		case kAnalyzeVolume:
			saveVolume_Analyze(path);
			break;
		case kMGZVolume:
			saveVolume_MGZ(path);
			break;
		default:
			printf("ERROR: Output data format not recognised\n");
			break;
	}
	return 0;
}
#pragma mark -
int main (int argc, const char * argv[])
{
	printf("%s\n",version);
	
	checkEndianness();
	
	// This command performs different computations and
	// analyses in volume data.
	// The operations are carried out in the order they
	// appear in the arguments passed to the command
	
	int		i;
	

	for(i=1;i<argc;i++)
	{		
		if(strcmp(argv[i],"-i")==0)				// input volume
		{
			loadVolume((char*)argv[++i]);
		}
		else
		if(strcmp(argv[i],"-o")==0)				// output volume
		{
			saveVolume((char*)argv[++i]);
		}
		else
		if(strcmp(argv[i],"-largest6con")==0)	// largest 6 connected component
		{
			largestConnected();
		}
		else
		if(strcmp(argv[i],"-dilate")==0)		// dilate(size)
		{
			int		size;
			sscanf(argv[++i]," %i ",&size);
			dilate(size);
		}
		else
		if(strcmp(argv[i],"-erode")==0)			// erode(size)
		{
			int		size;
			sscanf(argv[++i]," %i ",&size);
			erode(size);
		}
		else
		if(strcmp(argv[i],"-hist")==0)				// hist(#bins)
		{
			int		nbins;
			sscanf(argv[++i]," %i ",&nbins);
			hist(nbins);
		}
		else
		if(strcmp(argv[i],"-matchHist")==0)				// matchHistogram(another_mri)
		{
			matchHistogram((char*)argv[++i]);
		}
		else
		if(strcmp(argv[i],"-stats")==0)				// stats, returns mean, std, min, max
		{
			float	x=mean();
			float	s=std();
			float	mi=min();
			float	ma=max();
			
			printf("mean %g\nstd %g\nmin %g\nmax %g\n",x,s,mi,ma);
		}
		else
		if(strcmp(argv[i],"-tiff")==0)				// write slice as tiff file
		{
			char	*path=(char*)argv[i+1];
			char	*cmap=(char*)argv[i+2];
			float	*m=(float*)calloc(dim[0]*dim[1],sizeof(float));
			int		x=dim[0]/2,y,z;
			for(y=0;y<dim[1];y++)
				for(z=0;z<dim[2];z++)
					m[z+y*dim[2]]=getValue(x,y,z);
			tiff(path,m,dim[1],dim[2],cmap);
			i+=2;
		}
		else
		if(strcmp(argv[i],"-info")==0)				// information: dimensions, data type, pixel size
		{
			printf("dim %i,%i,%i\n",hdr->dim[1],hdr->dim[2],hdr->dim[3]);
			printf("dataType ");
			switch(hdr->datatype)
			{	case UCHAR:		printf("uchar\n"); break;
				case SHORT:		printf("short\n"); break;
				case FLOAT:		printf("float\n"); break;
				case INT:		printf("int\n"); break;
				case RGB:		printf("rgb\n"); break;
				case RGBFLOAT:	printf("rgbfloat\n"); break;
			}
			printf("voxelSize %g,%g,%g\n",hdr->pixdim[1],hdr->pixdim[2],hdr->pixdim[3]);
		}
		else
		if(strcmp(argv[i],"-threshold")==0)				// threshold(level,direction)
		{
			float	level;
			int		direction;
			sscanf(argv[++i]," %f , %i ",&level,&direction);
			threshold(level,direction);
		}
		else
		if(strcmp(argv[i],"-volume")==0)				// calculate volume
		{
			float	vol=volume();
			printf("volume=%g\n",vol);
		}
		else
		if(strcmp(argv[i],"-decompose")==0)		// decompose(basename) a volume with many values into volumes with one single value
		{
			decompose((char*)argv[++i]);
		}
	}
    return 0;
}
