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
#include <stdint.h>
#include <malloc.h>
#include <wand/MagickWand.h>
#include <libavformat/avformat.h>
#include <libavutil/imgutils.h>
#include <libswscale/swscale.h>
#include "blockhash.h"
#include "bmp_image.h"
#include "processing.h"
#include "phs_endian.h"

static MagickWand* load_image_from_frame(const char* src_file_name, size_t frame_number, void* frame_data, size_t frame_data_size) 
{
    // Create MagickWand object
    MagickWand* magick_wand = new_magic_wand();
    if(magick_wand) {
        // Read an image file into memory
        MagickBooleanType status = MagickReadImageBlob(magick_wand, frame_data, frame_data_size);
        if (status == MagickFalse) {
            fprintf(stderr, "Error reading converted frame #%llu of video file '%s'.\n", 
                    (unsigned long long)frame_number,
                    src_file_name);
            DestroyMagickWand(magick_wand);
            magick_wand = NULL;
        }
    }
    return magick_wand;    
}

int* process_video_frame(const hash_computation_task* task, size_t frame_number, void* frame_data, size_t frame_data_size, int* hash_size)
{
    int result = 0;
    int *hash;
    MagickWand *magick_wand;
    
    // Load Image
    magick_wand = load_image_from_frame(task->src_file_name, frame_number, frame_data, frame_data_size);
    if(!magick_wand) return NULL;
    
    // Compute Image Hash
    result = compute_image_hash(magick_wand, task->bits, task->hashing_method, &hash, hash_size);
    
    switch(result)
    {
        case 0: {
            // Show debug output
            if (task->debug && task->hashing_method < HM_PHASH_DCT64) {
                printf("Dump of the frame#%llu blockhash:\n", (unsigned long long)frame_number);
                debug_print_hash(hash, task->bits);
            } else {
                printf("Dump of the frame#%llu hash:\n", (unsigned long long)frame_number);
                print_hash(NULL, hash, *hash_size);                
            }
            break;
        }
        
        case 1: {
            fprintf(stderr, "Error computing blockhash for the zero-sized frame #%llu of the video file '%s'.\n", 
                    (unsigned long long)frame_number, task->src_file_name);
            break;
        }
        
        case 2: {
            fprintf(stderr, "Error converting image data to RGBA for the frame #%llu of the video file '%s'.\n", 
                    (unsigned long long)frame_number, task->src_file_name);
            break;
        }
        default: {
            fprintf(stderr, "Error computing blockhash for the frame #%llu of the video file '%s'.\n", 
                    (unsigned long long)frame_number, task->src_file_name);
            break;
        }
    }
    
    // Cleanup temporary data
    DestroyMagickWand(magick_wand);
    
    // Report the result
    return hash;
}


struct _video_decoding_state {
    const hash_computation_task* task;
    AVFormatContext *fmt_ctx;
    int video_stream_idx;
    AVCodecContext *original_video_codec_ctx;
    AVCodecContext *video_codec_ctx;
    AVCodec* video_codec;
    AVFrame* frame;
    AVPicture picture;
    AVPacket packet;
    struct SwsContext *sws_ctx;
    int got_frame;
    int cached;
};

typedef struct _video_decoding_state video_decoding_state;


typedef struct _video_frame_info {
    size_t frame_number;
    int* hash;
    int hash_size;
} video_frame_info;


#define HASH_PART_COUNT 4 

#define BMP_PIXEL_FORMAT AV_PIX_FMT_BGR24

int process_video_file(const hash_computation_task* task)
{
    size_t i;
    int restart_flag;
    int result = 0;
    int frames_counted = 0;
    uint64_t frame_count = 0;
    uint64_t next_hash_frame;
    uint64_t current_frame;
    int* hash = NULL;
    video_decoding_state state;
    video_frame_info hash_frames[HASH_PART_COUNT];
    int total_hash_size;

    if(task->debug)
        printf("Processing video file '%s'...\n", task->src_file_name);
    
restart:
    // Initialize data
    restart_flag = 0;
    memset(&hash_frames[0], 0, sizeof(hash_frames));
    memset(&state, 0, sizeof(state));
    state.task = task;
    
    // Coded using examples from http://dranger.com/ffmpeg/tutorial01.html
    
    // Open video file
    if (avformat_open_input(&state.fmt_ctx, task->src_file_name, NULL, NULL) < 0) {
        fprintf(stderr, "Error opening video file '%s'.\n", task->src_file_name);
        return -1;
    }
    
    // Populate streams info
    if(avformat_find_stream_info(state.fmt_ctx, NULL)<0) {
        result = -1;
        fprintf(stderr, "Error getting data streams information for the video file '%s'.\n", task->src_file_name);
        goto cleanup_close_video_file;
    }
    
    // Dump information for debug purposes
    if(task->debug) {
        av_dump_format(state.fmt_ctx, 0, task->src_file_name, 0);
    }
    
    // Find the first video stream
    state.video_stream_idx = -1;
    for(i = 0; i < state.fmt_ctx->nb_streams; i++) {
        if(state.fmt_ctx->streams[i]->codec->codec_type == AVMEDIA_TYPE_VIDEO) {
            state.video_stream_idx = i;
            break;                
        }
    }
    
    // Check if we have found a video stream
    if( state.video_stream_idx == -1) {
        result = -1;
        fprintf(stderr, "Error: Couldn't find video data stream in the video file '%s'.\n", task->src_file_name);
        goto cleanup_close_video_file;
    }
    
    // Get a pointer to the codec context for the video stream
    state.original_video_codec_ctx = state.fmt_ctx->streams[state.video_stream_idx]->codec;
    
    // Find the decoder for the video stream
    state.video_codec = avcodec_find_decoder(state.original_video_codec_ctx->codec_id);
    if(state.video_codec == NULL) {
        result = -1;
        fprintf(stderr, "Error: Unsupported codec in the video file '%s'.\n", task->src_file_name);
        goto cleanup_close_video_file;
    }
    
    // Allocate new context
    state.video_codec_ctx = avcodec_alloc_context3(state.video_codec);
    if(!state.video_codec_ctx) {
        result = -1;
        fprintf(stderr, "Error: Couldn't allocate new video codec context for the video file '%s'.\n", task->src_file_name);
        goto cleanup_close_video_file;
    }
    
    // Copy context
    if(avcodec_copy_context(state.video_codec_ctx, state.original_video_codec_ctx) != 0) {
        result = -1;
        fprintf(stderr, "Error: Couldn't copy video codec context for the video file '%s'.\n", task->src_file_name);
        goto cleanup_free_video_codec_context;
    }
    
    // Open codec
    if(avcodec_open2(state.video_codec_ctx, state.video_codec, NULL) < 0) {
        result = -1;
        fprintf(stderr, "Error: Couldn't open video codec for the video file '%s'.\n", task->src_file_name);
        goto cleanup_free_video_codec_context;
    }
    
    // Allocate storage for the original frame
    state.frame = avcodec_alloc_frame();
    if(!state.frame) {
        result = -1;
        fprintf(stderr, "Error: Couldn't allocate input frame for the video file '%s'.\n", task->src_file_name);
        goto cleanup_free_video_codec_context;        
    }

    // Get frame count
    if(!frames_counted) {     
        frame_count = state.fmt_ctx->streams[state.video_stream_idx]->nb_frames;
        if(frame_count == 0) {
            if(task->debug)
                printf("Frame count not available from stream, going to count frames manually...\n");
                
            // Frame count not available, need to count frames manually
            while(av_read_frame(state.fmt_ctx, &state.packet) >=0) {
                // Is this a packet from the video stream?
                if(state.packet.stream_index == state.video_stream_idx) {
                    // Decode video frame
                    state.got_frame = 0;
                    avcodec_decode_video2(state.video_codec_ctx, state.frame, &state.got_frame, &state.packet);
                    
                    // Did we get a video frame?
                    if(state.got_frame) {
                        ++frame_count;
                    }
                }
                
                // Free the packet that was allocated by av_read_frame
                av_free_packet(&state.packet);
            }
    
            if(task->debug)
                printf("Counted frames manually.\n");
            
            // [IP] Have to do all this trick with "restart" due to ffmpeg error
            // [H264] illegal short term buffer state detected
            frames_counted = 1;
            restart_flag = 1;
            goto cleanup_free_frame;
        }
    }
    
    // Print frame count for debug
    if(task->debug) {
        printf("Video file '%s' has %llu frames.\n", task->src_file_name, (unsigned long long)frame_count);
    }

    // Determine first and last hash frames
    if(frame_count < 11) {
        if(frame_count > 0) {
            hash_frames[0].frame_number = 0;
            hash_frames[3].frame_number = frame_count - 1;        
        } else {
            // Handle special zero frames case
            size_t hash_size = task->bits * task->bits * sizeof(int); 
            hash = (int*) malloc(hash_size);
            if(!hash) {
                result = -1;
                fprintf(stderr, "Error computing blockhash for the zero-frame video file '%s'.\n", 
                        task->src_file_name);
                goto cleanup;
            }
            memset(hash, 0, hash_size);
            goto hash_computed;
        }
    } else {
        hash_frames[0].frame_number = 10;
        hash_frames[3].frame_number = frame_count - 11;
    }
    
    // Determine middle hash frames
    hash_frames[1].frame_number = (size_t)floor(frame_count * 0.35);
    hash_frames[2].frame_number = (size_t)floor(frame_count * 0.7);

    if(av_image_alloc(state.picture.data, state.picture.linesize, state.video_codec_ctx->width, 
                   state.video_codec_ctx->height, BMP_PIXEL_FORMAT, 1) < 0) {
        result = -1;
        fprintf(stderr, "Error: Couldn't allocate picture buffer for the video file '%s'.\n", task->src_file_name);
        goto cleanup_free_frame;        
    }
    
    // initialize SWS context for software scaling
    state.sws_ctx = sws_getContext(state.video_codec_ctx->width,
        state.video_codec_ctx->height,
        state.video_codec_ctx->pix_fmt,
        state.video_codec_ctx->width,
        state.video_codec_ctx->height,
        BMP_PIXEL_FORMAT,
        SWS_BILINEAR,
        NULL,
        NULL,
        NULL
        );
    
    if(!state.sws_ctx) {
        result = -1;
        fprintf(stderr, "Error: Couldn't create image transformation context for the video file '%s'.\n", task->src_file_name);
        goto cleanup_free_picture_buffer;                                            
    }

    next_hash_frame = 0;
    current_frame = 0;

    // Read frames one by one and compute hash for the selected frames
    while(av_read_frame(state.fmt_ctx, &state.packet) >=0) {
        // Is this a packet from the video stream?
        if(state.packet.stream_index == state.video_stream_idx) {
            // Decode video frame
            state.got_frame = 0;
            avcodec_decode_video2(state.video_codec_ctx, state.frame, &state.got_frame, &state.packet);
            
            // Did we get a video frame?
            if(state.got_frame) {
                size_t bmp_size = 0;
                bitmap_image* bmp = NULL;
                for(i = next_hash_frame; i < HASH_PART_COUNT; ++i) {
                    if(hash_frames[i].frame_number == current_frame) {
                        if(!bmp) {
                            uint32_t tmp;
                            
                            // Convert the image from its native format to RGB
                            int height = sws_scale(state.sws_ctx, 
                                     (uint8_t const * const *)state.frame->data,
                                    state.frame->linesize, 
                                    0, 
                                    state.video_codec_ctx->height,
                                    state.picture.data, 
                                    state.picture.linesize);

                            // Create BMP
                            size_t image_line_size = state.video_codec_ctx->width * 3;
                            size_t image_line_size_with_padding =  image_line_size + ((4-(image_line_size)%4)%4);
                            size_t image_data_size = height * image_line_size_with_padding;
                            bmp_size = sizeof(bitmap_image) - 1 + image_data_size;
                            bmp = (bitmap_image*) malloc(bmp_size);
                            if(!bmp) {
                                result = -1;
                                fprintf(stderr, "Error: Couldn't allocate frame image buffer for the video file '%s'.\n", task->src_file_name);
                                av_free_packet(&state.packet);
                                goto cleanup;                                                                
                            }

                            SET_LE_16(bmp->bfh.bfType, 0x4D42); // 'B' 'M'
                            SET_LE_32(bmp->bfh.bfSize, bmp_size);
                            bmp->bfh.bfReserved1 = 0; // not using SET_LE_xx due to zero value
                            bmp->bfh.bfReserved2 = 0; // not using SET_LE_xx due to zero value
                            tmp = sizeof(bitmap_image) - 1;
                            SET_LE_32(bmp->bfh.bfOffBits, tmp);
                            SET_LE_32(bmp->bih.biSize, sizeof(bmp->bih));
                            SET_LE_32(bmp->bih.biWidth, state.video_codec_ctx->width);
                            SET_LE_32(bmp->bih.biHeight, height);
                            SET_LE_16(bmp->bih.biPlanes, 1);
                            SET_LE_16(bmp->bih.biBitCount, 24);
                            bmp->bih.biCompression = 0; // not using SET_LE_xx due to zero value
                            SET_LE_32(bmp->bih.biSizeImage, image_data_size);
                            bmp->bih.biXPelsPerMeter = 0; // not using SET_LE_xx due to zero value
                            bmp->bih.biYPelsPerMeter = 0; // not using SET_LE_xx due to zero value
                            bmp->bih.biClrUsed = 0; // not using SET_LE_xx due to zero value
                            bmp->bih.biClrImportant = 0; // not using SET_LE_xx due to zero value
                           
                            // Extract pixels
                            {
                                int y = 0;
                                size_t frame_line_size = state.picture.linesize[0];
                                uint8_t* bmp_data = &bmp->data[0] + (height-1)*image_line_size_with_padding;
                                uint8_t* picture_data = state.picture.data[0];
                                while(y < height)
                                {
                                    memcpy(bmp_data, picture_data, image_line_size);
                                    bmp_data -= image_line_size_with_padding;
                                    picture_data += frame_line_size;
                                    ++y;
                                }
                            }

                            if(task->debug) {
                                size_t file_name_buffer_size = strlen(task->src_file_name) + 128;
                                char* file_name_buffer = (char*) malloc(file_name_buffer_size);
                                strcpy(file_name_buffer, task->src_file_name);
                                char buffer[48];
                                snprintf(buffer, 48, "-frm-%"PRIu64".bmp", current_frame);
                                strcat(file_name_buffer, buffer);
                                printf("Saving frame #%"PRIu64" into file '%s'...\n", current_frame, file_name_buffer);
                                FILE* f = fopen(file_name_buffer, "wb");
                                if(f) {
                                    if(fwrite(bmp, bmp_size, 1, f) != 1) {
                                        fprintf(stderr, "Warning: Failed to write data into file '%s' for writing.\n", file_name_buffer);
                                    } else {
                                        fclose(f);
                                        printf("Saved frame #%"PRIu64" into file '%s'.\n", current_frame, file_name_buffer);
                                    }
                                }
                                else fprintf(stderr, "Warning: Failed to open file '%s' for writing.\n", file_name_buffer);
                            }
                        }
                        
                        {
                            int* hash = process_video_frame(task, current_frame, bmp, bmp_size, &hash_frames[i].hash_size);
                            if(!hash) {
                                free(bmp);
                                av_free_packet(&state.packet);
                                goto cleanup;
                            }
                            hash_frames[i].hash = hash;
                        }
                        
                        ++next_hash_frame;
                        // Do not break here, becasue in certain cases 
                        // a single frame may happen multipl times as hash_frame
                    }
                }
                
                // Free frame image
                if(bmp)
                    free(bmp);
                
                // Advance frame number
                ++current_frame;
            }
        }
            
        // Free the packet that was allocated by av_read_frame
        av_free_packet(&state.packet);
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    total_hash_size = 0;
    for(i = 0; i < HASH_PART_COUNT; ++i)
        total_hash_size += hash_frames[i].hash_size;
    
    hash = (int*) malloc(total_hash_size * sizeof(int) * HASH_PART_COUNT);
    if(!hash) {
        fprintf(stderr, "Error creating hash for video file '%s'.\n", task->src_file_name);
        goto cleanup;
    } else {
        int* dest = hash;
        for(i = 0; i < HASH_PART_COUNT; ++i) {
            size_t block_size = hash_frames[i].hash_size * sizeof(int);
            memcpy(dest, hash_frames[i].hash, block_size);
            dest += hash_frames[i].hash_size;
        }
    }
    
hash_computed:
    // Show debug output
    if (task->debug && task->hashing_method < HM_PHASH_DCT64)
        debug_print_hash(hash, task->bits);

    if(task->hashing_method < HM_PHASH_DCT64) {
        // Print blockhash string
        char* hex = blockhash_to_str(hash, HASH_PART_COUNT * task->bits * task->bits);
        if(hex) {
        result = 0;
        printf("%s  %s\n", hex, task->src_file_name);
        free(hex);
        } else {
        result = -1;
        fprintf(stderr, "Error converting blockhash value to string for the image file '%s'.\n", 
                task->src_file_name);
        }
    } else {
        print_hash(task->src_file_name, hash, total_hash_size);
    }

    
    // Free hash buffer
    free(hash);
    
cleanup:
    // Free partial hash buffers
    for(i = 0; i < HASH_PART_COUNT; ++i) {
        int* h = hash_frames[i].hash; 
        if(h) free(h);
    }

//cleanup_free_sws_context:
    // Free image transformation context
    if(state.sws_ctx)
        sws_freeContext(state.sws_ctx);
    
cleanup_free_picture_buffer:
    // Free picture buffer
    if(state.picture.data[0])
        av_freep(&state.picture.data[0]);
    
cleanup_free_frame:
    // Free original frame structure
    avcodec_free_frame(&state.frame);
    
cleanup_free_video_codec_context:
    // Free video codec context
    avcodec_close(state.video_codec_ctx);
    free(state.video_codec_ctx);

cleanup_close_video_file:
    // Close video file
    avformat_close_input(&state.fmt_ctx);
    
    // Restart after counting frames
    if(restart_flag) 
        goto restart;

    // report the result
    return result;
}
