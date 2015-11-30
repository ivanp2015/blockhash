/*
 * Perceptual image hash calculation tool based on algorithm descibed in
 * Block Mean Value Based Image Perceptual Hashing by Bian Yang, Fan Gu and Xiamu Niu
 *
 * Copyright (c) 2014-2015 Commons Machinery http://commonsmachinery.se/
 * Distributed under an MIT license, please see LICENSE in the top dir.
 */

#ifndef __BLOCKHASH_PROCESSING_H__
#define __BLOCKHASH_PROCESSING_H__

#define DEFAULT_BITS 16


enum _HashingMethod {
    HM_BLOCKHASH = 0,
    HM_BLOCKHASH_QUICK,
    HM_PHASH_DCT64
};

typedef enum _HashingMethod HashingMethod;


struct _hash_computation_task {
    const char* src_file_name;
    int bits;
    HashingMethod hashing_method;
    int debug;
    int video;
};

typedef struct _hash_computation_task hash_computation_task;


struct _MagickWand;
typedef struct _MagickWand MagickWand;


int process_image_file(const hash_computation_task* task);

int process_video_file(const hash_computation_task* task);

MagickWand* new_magic_wand();

MagickWand* load_image_from_file(const char* src_file_name);

int compute_image_hash(MagickWand* magick_wand, int bits, HashingMethod hashing_method, int** hash, int* hash_size);

void debug_print_hash(const int* hash, int bits);

void print_hash(const char* opt_filename, const int* hash, int size);

#endif
