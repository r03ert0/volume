#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Nifti.h"
#include "swap.h"

int NiftiBytesPerVoxel(nifti_1_header hdr)
{
	int	bpv=0;
	
	switch(hdr.datatype)
	{	case DT_UINT8:		bpv=1; break;
		case DT_INT16:		bpv=2; break;
		case DT_INT32:		bpv=4; break;
		case DT_FLOAT32:	bpv=4; break;
		case DT_COMPLEX64:	bpv=8; break;
		case DT_FLOAT64:	bpv=8; break;
		case DT_RGB24:		bpv=3; break;
		case DT_INT8:		bpv=1; break;
		case DT_UINT16:		bpv=2; break;
		case DT_UINT32:		bpv=4; break;
		case DT_INT64:		bpv=8; break;
		case DT_UINT64:		bpv=8; break;
		case DT_FLOAT128:	bpv=16; break;
		case DT_COMPLEX128:	bpv=16; break;
		case DT_COMPLEX256:	bpv=32; break;
		case DT_RGBA32:		bpv=4; break;
	}
	
	return bpv;
}
void Nifti_load(char *path, int svol, char **addr,int *sz, int *swapped)
{
	nifti_1_header	hdr;
	FILE	*f;
	int		size;
	char    *img;

	f=fopen(path,"r");
	if(f==NULL){printf("ERROR: Cannot open Nifti file for loading\n");return;}
	
	fread(&hdr,sizeof(hdr),sizeof(char),f);
	*swapped=0;
	if(hdr.sizeof_hdr!=348) // different endianness
	{
		swap_nifti_hdr(&hdr);
		*swapped=1;
//		printf("ERROR: Cannot open Nifti files with Motorola endianess\n");
	}
	if(hdr.dim[1]==3 && hdr.dim[4]>1)
		size=hdr.dim[2]*hdr.dim[3]*hdr.dim[4];
	else
		size=hdr.dim[1]*hdr.dim[2]*hdr.dim[3];

	*sz=size*NiftiBytesPerVoxel(hdr)+sizeof(hdr);
	*addr=calloc(*sz,1);
	memcpy(*addr,&hdr,sizeof(hdr));
	
	fseek(f, hdr.vox_offset-sizeof(hdr)
		+ svol*size*NiftiBytesPerVoxel(hdr)
	,SEEK_CUR);
	
	fread(*addr+sizeof(hdr),size,NiftiBytesPerVoxel(hdr),f);
	fclose(f);

	img = *addr + (int)hdr.vox_offset;
	if(*swapped)	// different endianness
	{
		swap_nifti_img(img,hdr);
	}

}
void Nifti_save(char *path, char *addr)
{
	nifti_1_header	*hdr;
	FILE	*f;
	int		size;
	char	extension[4]={0,0,0,0};
	
	hdr=(nifti_1_header*)addr;

	if(hdr->dim[1]==3 && hdr->dim[4]>1)
		size=hdr->dim[2]*hdr->dim[3]*hdr->dim[4];
	else
		size=hdr->dim[1]*hdr->dim[2]*hdr->dim[3];

	hdr->vox_offset=352;
    if(hdr->scl_slope==0)
    {
        hdr->scl_inter=0;
        hdr->scl_slope=1;
    }
	f=fopen(path,"w");
	if(f==NULL){printf("ERROR: Cannot open Nifti file for saving\n");return;}
	fwrite(addr,sizeof(*hdr),1,f);
	fwrite(extension,4,1,f);
	fwrite(addr+sizeof(*hdr),size,NiftiBytesPerVoxel(*hdr),f);
	fclose(f);
}

#pragma mark -

void swap_nifti_hdr(nifti_1_header *hdr)
{
	int i;
	swap_int(&(*hdr).sizeof_hdr);
	swap_int(&(*hdr).extents);
	swap_short(&(*hdr).session_error);
	for(i=0;i<8;i++) swap_short(&(*hdr).dim[i]);
	swap_float(&(*hdr).intent_p1);
	swap_float(&(*hdr).intent_p2);
	swap_float(&(*hdr).intent_p3);
	swap_short(&(*hdr).intent_code);
	swap_short(&(*hdr).datatype);
	swap_short(&(*hdr).bitpix);
	swap_short(&(*hdr).slice_start);
	for(i=0;i<8;i++) swap_float(&(*hdr).pixdim[i]);
	swap_float(&(*hdr).vox_offset);
	swap_float(&(*hdr).scl_slope);
	swap_float(&(*hdr).scl_inter);
	swap_short(&(*hdr).slice_end);
	swap_float(&(*hdr).cal_max);
	swap_float(&(*hdr).cal_min);
	swap_float(&(*hdr).slice_duration);
	swap_float(&(*hdr).toffset);
	swap_int(&(*hdr).glmax);
	swap_int(&(*hdr).glmin);
	swap_short(&(*hdr).qform_code);
	swap_short(&(*hdr).sform_code);
	swap_float(&(*hdr).quatern_b);
	swap_float(&(*hdr).quatern_c);
	swap_float(&(*hdr).quatern_d);
	swap_float(&(*hdr).qoffset_x);
	swap_float(&(*hdr).qoffset_y);
	swap_float(&(*hdr).qoffset_z);
	for(i=0;i<4;i++) swap_float(&(*hdr).srow_x[i]);
	for(i=0;i<4;i++) swap_float(&(*hdr).srow_y[i]);
	for(i=0;i<4;i++) swap_float(&(*hdr).srow_z[i]);
}
void swap_nifti_img(char *img, nifti_1_header hdr)
{
	int		i,sz=hdr.dim[1]*hdr.dim[2]*hdr.dim[3];

	switch(hdr.datatype)
	{
		case DT_UINT8:
		case DT_INT8:
			break;
		case DT_UINT16:
		case DT_INT16:
			for(i=0;i<sz;i++) swap_short(&((short*)img)[i]);
			break;
		case DT_UINT32:
		case DT_INT32:
			for(i=0;i<sz;i++) swap_int(&((int*)img)[i]);
			break;
		case DT_FLOAT32:
			for(i=0;i<sz;i++) swap_float(&((float*)img)[i]);
			break;
		case DT_COMPLEX64:
		case DT_FLOAT64:
		case DT_RGB24:
		case DT_INT64:
		case DT_UINT64:
		case DT_FLOAT128:
		case DT_COMPLEX128:
		case DT_COMPLEX256:
		case DT_RGBA32:
			printf("ERROR: Unsupported hdr.datatype %i\n", (int)hdr.datatype);
			break;
	}
}

