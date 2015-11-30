/*
 * Perceptual image hash calculation library based on the phash algorithm
 *
 * Copyright (c) 2015 Commons Machinery http://commonsmachinery.se/
 * Distributed under an MIT license, please see LICENSE in the top dir.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <errno.h>
#include "phs.h"
#include "phs_endian.h"
#include "phs_swap.h"


static float* rgb_to_y(const unsigned char *data, int width, int height)
{
    float* result = NULL;
    int pixel_count = width * height;
    if(pixel_count > 0) {
        result = (float*) malloc(pixel_count * sizeof(float));
        if(result) {
            float* p = result;
            const unsigned char *data_end = data + 3 * pixel_count;
            while(data < data_end) {
                const float
                    R = (float)*data++,
                    G = (float)*data++,
                    B = (float)*data++,
                    Y = (66*R + 129*G + 25*B + 128)/256 + 16;
                *p++ = (float)(Y<0?0:(Y>255?255:Y));
            }
        }
    }
    return result;
}


static float* convolve(
    const float* image, int image_width, int image_height,
    const float* mask, int mask_width, int mask_height,
    int boundary_conditions
)
{
    float* result = malloc(image_width * image_height * sizeof(float));
    
    if(result) {
        /* generic version for other masks and boundary conditions */
        const int
            mx2 = mask_width/2, 
            my2 = mask_height/2, 
            mx1 = mx2 - 1 + (mask_width%2), 
            my1 = my2 - 1 + (mask_height%2), 
            mxe = image_width - mx2, 
            mye = image_height - my2; 
            
        int x, y, xm, ym;
                    
        /* classical correlation */
        {
            float* p = result + my1*image_width;
            const float* pimg = image;
#ifdef PHASH_SIMPLE_USE_OPENMP
#pragma omp parallel for collapse(2) if (image_width * image_height >= 32768)
#endif
            for (y = my1; y<mye; ++y, p += image_width, pimg += image_width) {
                for (x = mx1; x<mxe; ++x) {
                    const float* pmasky = mask;
                    const float* pimgy = pimg - mx1 + x;
                    float val = 0;
                    for (ym = -my1; ym<=my2; ++ym, pimgy += image_width, pmasky += mask_width) {
                        const float* pmaskx = pmasky;
                        const float* pimgx = pimgy - mx1 + x;
                        for (xm = -mx1; xm<=mx2; ++xm, ++pimgx, ++pmaskx) {
                            val+= (*pimgx) * (*pmaskx);
                        }
                    }
                    p[x] = val;
                }
            }
        }
    
        if (boundary_conditions) {
            float* p = result;
#ifdef PHASH_SIMPLE_USE_OPENMP
#pragma omp parallel for collapse(1) if (image_width >= 256 && image_height >= 128)
#endif
            for(y = 0; y < image_height; ++y, p += image_width) {
                for (x = 0; x < image_width; (y<my1 || y>=mye )?++x:((x<mx1 - 1 || x>=mxe)?++x:(x=mxe))) {
                    const float* pmasky = mask;
                    float val = 0;
                    for (ym = -my1; ym<=my2; ++ym, pmasky += mask_width) {
                        const float* pmaskx = pmasky;
                        for (xm = -mx1; xm<=mx2; ++xm, ++pmaskx) {
                            /* uses [0...last_coord] if OOR */
                            int xx = x + xm; xx = xx <= 0 ? 0 : (xx >= image_width ? image_width - 1 : xx);
                            int yy = y + ym; yy = yy <= 0 ? 0 : (yy >= image_height ? image_height - 1 : yy);
                            val+= *(image + yy*image_width + xx) * (*pmaskx);
                        }
                    }
                    p[x] = val;
                }
            }
        } else {
            float* p = result;
#ifdef PHASH_SIMPLE_USE_OPENMP
#pragma omp parallel for collapse(1) if (image_width >= 256 && image_height >= 128)
#endif
            for(y = 0; y < image_height; ++y, p += image_width) {
                for (x = 0; x < image_width; (y<my1 || y>=mye)?++x:((x<mx1 - 1 || x>=mxe)?++x:(x=mxe))) {
                    const float* pmasky = mask;
                    float val = 0;
                    for (ym = -my1; ym<=my2; ++ym) {
                        const float* pmaskx = pmasky;
                        for (xm = -mx1; xm<=mx2; ++xm, ++pmaskx) {
                            /* uses zero value if OOR */
                            int xx = x + xm; 
                            int yy = y + ym; 
                            if (xx >= 0 && xx < image_width && yy >= 0 && yy < image_height)
                                val+= *(image + yy*image_width + xx) * (*pmaskx);
                        }
                    }
                    p[x] = val;
                }
            }
        }
    }
    
    return result;
}


static int kth_smallest(const float* data, size_t size, size_t k, float* result)
{
    size_t nbytes = sizeof(float) * size;
    float *arr = malloc(nbytes);
    if(!arr) return -1;
    memcpy(arr, data, nbytes);
    {
	size_t l = 0, ir = size - 1;
	while(1) {
	    if (ir<=l + 1) {
		    if (ir==l + 1 && arr[ir]<arr[l]) 
			phs_swap_floats(arr + l, arr + ir);
		    *result = arr[k];
		    break;
	    } else {
		const size_t mid = (l + ir)>>1;
		phs_swap_floats(arr + mid,arr + l + 1);
		if (arr[l]>arr[ir]) phs_swap_floats(arr + l,arr + ir);
		if (arr[l + 1]>arr[ir]) phs_swap_floats(arr + l + 1, arr + ir);
		if (arr[l]>arr[l + 1]) phs_swap_floats(arr + l, arr + l + 1);
		{
		    size_t i = l + 1, j = ir;
		    const float pivot = arr[l + 1];
		    while (1) {
			do ++i; while (arr[i]<pivot);
			do --j; while (arr[j]>pivot);
			if (j<i) break;
			phs_swap_floats(arr + i, arr + j);
		    }
		    arr[l + 1] = arr[j];
		    arr[j] = pivot;
		    if (j>=k) ir = j - 1;
		    if (j<=k) l = i;
		}
	    }
	}
    }
    free(arr);
    return 0;
}


static int median(const float* data, size_t size, float* result)
{
    float r1 = 0.0;
    if(kth_smallest(data, size, size>>1, &r1))
	return -1;
    else if(size % 2) {
	*result = r1;
	return 0;
    } else {
	float r2 = 0.0;
	if(kth_smallest(data, size, (size>>1) - 1, &r2))
	    return -1;
	else {
	    *result = (r1 + r2) / 2;
	    return 0;
	}
    }
}


float* mul_matrix (const float* m1, int w1, int h1, const float* m2, int w2, int h2) {
    if(w1 == h2) {
	float* r = malloc(sizeof(float) * w2 * h1);
	if(r) {
	    int i, j, k;
	    const float* p1 = m1;
	    const float* p2 = m2;
	    float* pr = r;
	    for(i = 0; i < w2; ++i, ++p2) {
		const float* p2c = p2;
		for(j = 0; j < h1; ++j, ++pr) {
		    float v = 0.0;
		    for(k = 0; k < w1; ++k, ++p1, p2c += w2) {
			v += (*p1) * (*p2c);
		    }
		    *pr = v;
		}
	    }
	}
	return r;
    } else {
	errno = EINVAL;
	return NULL;
    }
}


/* nearest-neighbor interpolation resize */
static float* nni_resize(const float* img, int width, int height, int new_width, int new_height)
{
    const unsigned int
	_sx = (unsigned int)(new_width < 0 ? -new_width*width/100 : new_width),
	_sy = (unsigned int)(new_height < 0 ? -new_height*height/100 : new_height),
	sx = _sx ? _sx : 1, 
	sy = _sy ? _sy : 1;

    float* result;
    if (sx==width && sy == height)
	result = (float*)img;
   /* Nearest neighbor interpolation */
   else {
       result = malloc(sizeof(float)*sx*sy);
       if(result) {
           unsigned long* off_x = malloc(sizeof(unsigned long) * sx);
           if(off_x) {
               unsigned long* off_y = malloc(sizeof(unsigned long) * sy);
               if(off_y) {
                   unsigned int i;
                   
                   if (sx == width) {
                       for(i = 0; i < sx; ++i) 
                           off_x[i] = 1;
                   } else {
                       unsigned long *poff_x = off_x, curr = 0;
                       for(i = 0; i < sx; ++i) {
                           const unsigned long old = curr;
                           curr = (unsigned long)((i + 1.0)*width/sx);
                           *(poff_x++) = curr - old;
                       }
                   }
                   
                   if (sy == height) {
                       for(i = 0; i < sy; ++i) 
                           off_y[i] = width;
                   } else {
                       unsigned long *poff_y = off_y, curr = 0;
                       for(i = 0; i < sx; ++i) {
                           const unsigned long old = curr;
                           curr = (unsigned long)((i + 1.0)*height/sy);
                           *(poff_y++) = width*(curr - old);
                       }
                       *poff_y = 0;
                   }
                   
                   {
                        float *p = result;
                        const float *py = img;
                        const unsigned long* poff_y = off_y;
                        for (i = 0; i < sy; ) {
                            const float *px = py;
                            const unsigned long *poff_x = off_x;
                            unsigned long dy;
                            unsigned int j;
                            for(j = 0; j < sx; ++j) {
                                *p++ = *px;
                                px += *poff_x++;
                            }
                            ++i;
                            dy = *poff_y++;
                            for ( ; !dy && i < dy; ++i, p += sx, dy = *poff_y++)
                                memcpy(p, p - sx, sizeof(float)*sx);
                            py += dy;
                        }
                   }
                   
                   free(off_y);    
               } else {
                    free(result);
                    result = NULL;
               }
               free(off_x);
           } else {
               free(result);
               result = NULL;
           }
       }
   }
   return result;
}


#define DCT_MATRIX_SIDE_SIZE 32

static const float* get_dct_matrix()
{
    static int dct_matrix_prepared = 0;
    static float dct_matrix[DCT_MATRIX_SIDE_SIZE * DCT_MATRIX_SIDE_SIZE];
    
    if(dct_matrix_prepared == 0)
    {
	int x, y;
	float *p1 = dct_matrix;
	
	/* compute initial DCT matrix */
	
	float c1 = 1.0/sqrt(DCT_MATRIX_SIDE_SIZE);
	for (x = 0; x < DCT_MATRIX_SIDE_SIZE; ++x, ++p1)
	    *p1 = c1;
	
	c1 = 1.0/sqrt(2.0/DCT_MATRIX_SIDE_SIZE);
	for (y = 1; y < DCT_MATRIX_SIDE_SIZE; ++y) {
	    for (x = 0; x < DCT_MATRIX_SIDE_SIZE; ++x, ++p1) {
		*p1 = c1 * cos((M_PI/2/DCT_MATRIX_SIDE_SIZE) * y * (2*x + 1));
	    }
	}
	
	dct_matrix_prepared = 1;
    }
    
    return dct_matrix;
}


static const float* get_transposed_dct_matrix()
{
    static int transposed_dct_matrix_prepared = 0;
    static float transposed_dct_matrix[DCT_MATRIX_SIDE_SIZE * DCT_MATRIX_SIDE_SIZE];
    
    if(transposed_dct_matrix_prepared == 0)
    {
	int x, y;
	float *p1 = transposed_dct_matrix, *p2 = transposed_dct_matrix;
	
	memcpy(transposed_dct_matrix, get_dct_matrix(), 
	       sizeof(float) * DCT_MATRIX_SIDE_SIZE * DCT_MATRIX_SIDE_SIZE);
	
	/* transpose initial DCT matrix */
	for (x = 0; x < DCT_MATRIX_SIDE_SIZE; ++x, ++p2) {
	    float* p = p2;
	    for (y = 0; y < DCT_MATRIX_SIDE_SIZE; ++y, ++p1, p += DCT_MATRIX_SIDE_SIZE) {
		phs_swap_floats(p1, p);
	    }
	}
	
	transposed_dct_matrix_prepared = 1;
    }
    
    return transposed_dct_matrix;
}


int phs_dct_image_hash(unsigned char *data, int width, int height, int **hash)
{    
    #define MEAN_FILTER_WIDTH 7
    #define MEAN_FILTER_HEIGHT 7

    static const float mean_filter[MEAN_FILTER_WIDTH * MEAN_FILTER_HEIGHT] = { 
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 
    };
    
    int save_errno;
    float median_value;
    float *img1, *img2;
    const float* dct_matrix = get_dct_matrix();
    const float* transposed_dct_matrix = get_transposed_dct_matrix();
    
    img1 = rgb_to_y(data, width, height);
    if(!img1)
        return -1;
    
    img2 = convolve(img1, width, height, mean_filter, MEAN_FILTER_WIDTH, MEAN_FILTER_HEIGHT, 1);
    save_errno = errno;
    free(img1);
    if(!img2) {
        errno = save_errno;
        return -1;
    }
    
    img1 = nni_resize(img2, width, height, DCT_MATRIX_SIDE_SIZE, DCT_MATRIX_SIDE_SIZE);
    save_errno = errno;
    if(img1 != img2)
	free(img2);
    if(!img1) {
        errno = save_errno;
        return -1;	
    }
    
    img2 = mul_matrix(dct_matrix, DCT_MATRIX_SIDE_SIZE, DCT_MATRIX_SIDE_SIZE,
	img1, DCT_MATRIX_SIDE_SIZE, DCT_MATRIX_SIDE_SIZE);
    save_errno = errno;
    free(img1);
    if(!img2) {
        errno = save_errno;
        return -1;	
    }
    
    img1 = mul_matrix(img2, DCT_MATRIX_SIDE_SIZE, DCT_MATRIX_SIDE_SIZE,
	transposed_dct_matrix, DCT_MATRIX_SIDE_SIZE, DCT_MATRIX_SIDE_SIZE);
    save_errno = errno;
    free(img2);
    if(!img1) {
        errno = save_errno;
        return -1;	
    }
    
    /* crop(1,1,8,8) */
    {
        int y;
        float* p = img1;
        const float* s = img1 + DCT_MATRIX_SIDE_SIZE + 1;
        for(y = 1; y < 9; ++y, p += 8, s+= DCT_MATRIX_SIDE_SIZE)
            memcpy(p, s, sizeof(float)*8);
    }
    
    if(median(img1, 64, &median_value)) {
	save_errno = errno;
	free(img1);
        errno = save_errno;
        return -1; 	
    } else {
	int i;
	int* result;
	uint64_t one = 1;
	uint64_t hash_value = 0;
	
	for (i = 0; i < 64; i++) {
	    if (img1[i] > median_value)
		hash_value |= one;
	    one = one << 1;
	}
	free(img1);
	
	result = (int*) malloc(8);
	if(!result)
	    return -1;
	else {
	    *((uint64_t*)result) = hash_value;
#if !defined(BIG_ENDIAN)
	    phs_swap_ints(result, result + 1);
#endif
	    *hash = result;
	    return 0;
	}
    }
}
