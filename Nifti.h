#ifndef __Nifti__
#define __Nifti__

#include "nifti1.h"

void Nifti_load(char *path, int svol, char **addr,int *sz, int *swapped);
void Nifti_save(char *path, char *addr);

void swap_nifti_hdr(nifti_1_header *hdr);
void swap_nifti_img(char *img, nifti_1_header hdr);
#endif