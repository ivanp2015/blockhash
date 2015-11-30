/*
 * Perceptual image hash calculation library based on the phash algorithm
 *
 * Copyright (c) 2015 Commons Machinery http://commonsmachinery.se/
 * Distributed under an MIT license, please see LICENSE in the top dir.
 */

#ifndef __LIBPHASH_SIMPLE_PHS_SWAP_H__
#define __LIBPHASH_SIMPLE_PHS_SWAP_H__

#define PHS_SWAP_FUNCTION(_type, _name) \
static inline void _name(_type* a, _type* b) { \
    register _type tmp = *a; *a = *b; *b = tmp; \
}

PHS_SWAP_FUNCTION(int, phs_swap_ints)

PHS_SWAP_FUNCTION(float, phs_swap_floats)

#endif
