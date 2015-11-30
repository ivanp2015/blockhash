lib_LTLIBRARIES = libblockhash.la libphash_simple.la

libblockhash_la_SOURCES = blockhash.c
libblockhash_la_LDFLAGS = -no-undefined

libphash_simple_la_SOURCES = phs_dct_image_hash.c
libphash_simple_la_LDFLAGS = -no-undefined

include_HEADERS = blockhash.h phs.h phs_endian.h phs_swap.h

bin_PROGRAMS = blockhash
blockhash_SOURCES = main.c misc.c process_image.c process_video.c processing.h
blockhash_LDADD = libblockhash.la libphash_simple.la
