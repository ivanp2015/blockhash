/*
 * Perceptual image hash calculation tool based on algorithm descibed in
 * Block Mean Value Based Image Perceptual Hashing by Bian Yang, Fan Gu and Xiamu Niu
 *
 * Copyright (c) 2014-2015 Commons Machinery http://commonsmachinery.se/
 * Distributed under an MIT license, please see LICENSE in the top dir.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include "processing.h"

void debug_print_hash(const int* hash, int bits) 
{
    int i;
    for (i = 0; i < bits * bits; i++) {
        if (i != 0 && i % bits == 0)
            printf("\n");
        printf("%d", hash[i]);
    }
    printf("\n");
}

void print_hash(const char* opt_filename, const int* hash, int size) 
{
    int i;
    for (i = 0; i < size; i++) {
        printf("%08x", hash[i]);
    }
    if(opt_filename)
        printf("  %s", opt_filename);
    printf("\n");
}
