/*
 * Perceptual image hash calculation library based on the phash algorithm
 *
 * Copyright (c) 2015 Commons Machinery http://commonsmachinery.se/
 * Distributed under an MIT license, please see LICENSE in the top dir.
 */

#ifndef __LIBPHASH_SIMPLE_PHS_ENDIAN_H__
#define __LIBPHASH_SIMPLE_PHS_ENDIAN_H__

#ifdef __BIG_ENDIAN__
# define BIG_ENDIAN
#elif defined __LITTLE_ENDIAN__
/* override */
#elif defined __BYTE_ORDER__
# if __BYTE_ORDER__ ==  __ORDER_BIG_ENDIAN__
# define BIG_ENDIAN
# endif
#else // ! defined __LITTLE_ENDIAN__
#ifdef HAVE_MACHINE_ENDIAN_H
# include <machine/endian.h>
#else
# include <endian.h>
#endif
# if __BYTE_ORDER__ ==  __ORDER_BIG_ENDIAN__
#  define BIG_ENDIAN
# endif
#endif

#define WRITE_REORDERED_16( p, n ) \
    ( \
                ( (uint8_t*) (p) )[ 0 ] = (uint8_t) (n), \
                ( (uint8_t*) (p) )[ 1 ] = (uint8_t) ( (n) >> 8) \
    )
        
#define WRITE_REORDERED_32( p, n ) \
    ( \
                ( (uint8_t*) (p) )[ 0 ] = (uint8_t) (n), \
                ( (uint8_t*) (p) )[ 1 ] = (uint8_t) ( (n) >> 8), \
                ( (uint8_t*) (p) )[ 2 ] = (uint8_t) ( (n) >> 16), \
                ( (uint8_t*) (p) )[ 3 ] = (uint8_t) ( (n) >> 24) \
        )

#ifdef BIG_ENDIAN

#define SET_LE_16(x, v) WRITE_REORDERED_16(&x, (v))
#define SET_LE_32(x, v) WRITE_REORDERED_32(&x, (v))

#else

#define SET_LE_16(x, v) x = (v)
#define SET_LE_32(x, v) x = (v)

#endif

#endif
