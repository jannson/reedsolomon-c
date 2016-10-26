#ifndef __RS_H_

/* use small value to save memory */
#define DATA_SHARDS_MAX (255)

/* use other memory allocator */
#define RS_MALLOC(x)    malloc(x)
#define RS_FREE(x)      free(x)
#define RS_CALLOC(n, x) calloc(n, x)

typedef struct _reed_solomon {
    int data_shards;
    int parity_shards;
    int shards;
    gf* m;
    gf* parity;
} reed_solomon;

reed_solomon* reed_solomon_new(int data_shards, int parity_shards);
void reed_solomon_release(reed_solomon* rs);

/**
 * input:
 * rs
 * data_blocks[rs->data_shards][block_size]
 * fec_blocks[rs->data_shards][block_size]
 * */
int reed_solomon_encode(reed_solomon* rs,
        unsigned char** data_blocks,
        unsigned char** fec_blocks,
        int block_size);


/** input:
 * rs
 * original data_blocks[rs->data_shards][block_size]
 * dec_fec_blocks[nr_fec_blocks][block_size]
 * fec_block_nos: fec pos number in original fec_blocks
 * erased_blocks: erased blocks in original data_blocks
 * nr_fec_blocks: the number of erased blocks
 * */
int reed_solomon_decode(reed_solomon* rs,
        unsigned char **data_blocks,
        int block_size,
        unsigned char **dec_fec_blocks,
        unsigned int *fec_block_nos,
        unsigned int *erased_blocks,
        int nr_fec_blocks);

/** 
 * input:
 * rs
 * shards[rs->shards][block_size]
 * */
int reed_solomon_encode2(reed_solomon* rs, unsigned char** shards, int block_size);

/** 
 * input:
 * rs
 * shards[rs->shards][block_size]
 * */
int reed_solomon_reconstruct(reed_solomon* rs, unsigned char** shards, int block_size);
#endif

