/*
 * Perceptual image hash calculation library based on the phash algorithm
 *
 * Copyright (c) 2015 Commons Machinery http://commonsmachinery.se/
 * Distributed under an MIT license, please see LICENSE in the top dir.
 */

#ifndef __LIBPHASH_SIMPLE_PHS_H__
#define __LIBPHASH_SIMPLE_PHS_H__

#ifdef __cplusplus
extern "C" {
#endif

/** Calculate perceptual hash for an RGB image using robust phash DCT method with fixed 64-bit hash.
*
* Parameters:
*
* data - RGB image data.
* width - image width.
* height - image height.
* hash - the resulting hash will be allocated and stored in the given array as bits.
* 
* Returns: 0 on success, -1 on failure.
*/
int phs_dct_image_hash(unsigned char *data, int width, int height, int **hash);

#ifdef __cplusplus
extern "C" {
#endif

#endif
