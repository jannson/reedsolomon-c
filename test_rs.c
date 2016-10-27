#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

#define PROFILE
#include "rs.h"
#include "rs.c"

void print_matrix1(gf* matrix, int nrows, int ncols);
void print_matrix2(gf** matrix, int nrows, int ncols);

void print_buf(gf* buf, char *fmt, size_t len) {
    size_t i = 0;
    while(i < len) {
        printf(fmt, buf[i]);
        i++;
        if((i % 16) == 0) {
            printf("\n");
        }
    }
    printf("\n");
}

void test_001(void) {
    reed_solomon* rs = reed_solomon_new(11, 6);
    print_matrix1(rs->m, rs->data_shards, rs->data_shards);
    print_matrix1(rs->parity, rs->parity_shards, rs->data_shards);
    reed_solomon_release(rs);
}

void test_002(void) {
    char text[] = "hello world", output[256];
    int block_size = 1;
    int nrDataBlocks = sizeof(text)/sizeof(char) - 1;
    unsigned char* data_blocks[128];
    unsigned char* fec_blocks[128];
    int nrFecBlocks = 6;

    //decode
    unsigned int fec_block_nos[128], erased_blocks[128];
    unsigned char* dec_fec_blocks[128];
    int nr_fec_blocks;

    int i;
    reed_solomon* rs = reed_solomon_new(nrDataBlocks, nrFecBlocks);

    printf("%s:\n", __FUNCTION__);

    for(i = 0; i < nrDataBlocks; i++) {
        data_blocks[i] = (unsigned char*)&text[i];
    }

    memset(output, 0, sizeof(output));
    memcpy(output, text, nrDataBlocks);
    for(i = 0; i < nrFecBlocks; i++) {
        fec_blocks[i] = (unsigned char*)&output[i + nrDataBlocks];
    }

    reed_solomon_encode(rs, data_blocks, fec_blocks, block_size);
    print_buf((gf*)output, "%d ", nrFecBlocks+nrDataBlocks);

    text[1] = 'x';
    text[3] = 'y';
    text[4] = 'z';
    erased_blocks[0] = 4;
    erased_blocks[1] = 1;
    erased_blocks[2] = 3;

    fec_block_nos[0] = 1;
    fec_block_nos[1] = 3;
    fec_block_nos[2] = 5;
    dec_fec_blocks[0] = fec_blocks[1];
    dec_fec_blocks[1] = fec_blocks[3];
    dec_fec_blocks[2] = fec_blocks[5];
    nr_fec_blocks = 3;

    printf("erased:%s\n", text);

    reed_solomon_decode(rs, data_blocks, block_size, dec_fec_blocks,
            fec_block_nos, erased_blocks, nr_fec_blocks);

    printf("fixed:%s\n", text);

    reed_solomon_release(rs);
}

void test_003(void) {
    char text[] = "hello world hello world ", output[256];
    int block_size = 2;
    int nrDataBlocks = (sizeof(text)/sizeof(char) - 1) / block_size;
    unsigned char* data_blocks[128];
    unsigned char* fec_blocks[128];
    int nrFecBlocks = 6;

    //decode
    unsigned int fec_block_nos[128], erased_blocks[128];
    unsigned char* dec_fec_blocks[128];
    int nr_fec_blocks;

    int i;
    reed_solomon* rs = reed_solomon_new(nrDataBlocks, nrFecBlocks);

    printf("%s:\n", __FUNCTION__);
    //printf("text size=%d\n", (int)(sizeof(text)/sizeof(char) - 1) );

    for(i = 0; i < nrDataBlocks; i++) {
        data_blocks[i] = (unsigned char*)&text[i*block_size];
    }

    memset(output, 0, sizeof(output));
    memcpy(output, text, nrDataBlocks*block_size);
    //print_matrix1((gf*)output, nrDataBlocks + nrFecBlocks, block_size);
    for(i = 0; i < nrFecBlocks; i++) {
        fec_blocks[i] = (unsigned char*)&output[i*block_size + nrDataBlocks*block_size];
    }
    reed_solomon_encode(rs, data_blocks, fec_blocks, block_size);
    //print_buf((gf*)output, "%d ", nrFecBlocks*block_size + nrDataBlocks*block_size);
    //print_matrix1((gf*)output, nrDataBlocks + nrFecBlocks, block_size);

    //decode
    text[1*block_size] = 'x';
    text[10*block_size+1] = 'y';
    text[4*block_size] = 'z';
    erased_blocks[0] = 4;
    erased_blocks[1] = 1;
    erased_blocks[2] = 10;

    fec_block_nos[0] = 1;
    fec_block_nos[1] = 3;
    fec_block_nos[2] = 5;
    dec_fec_blocks[0] = fec_blocks[1];
    dec_fec_blocks[1] = fec_blocks[3];
    dec_fec_blocks[2] = fec_blocks[5];
    nr_fec_blocks = 3;

    printf("erased:%s\n", text);

    reed_solomon_decode(rs, data_blocks, block_size, dec_fec_blocks,
            fec_block_nos, erased_blocks, nr_fec_blocks);

    printf("fixed:%s\n", text);

    reed_solomon_release(rs);
}

void test_004(void) {
    char text[] = "hello world hello world ";
    int dataShards = 6;
    int parityShards = 4;
    int blockSize = 2;
    struct timeval tv;
    int i, j, n, seed, size, nrShards, nrBlocks, nrFecBlocks;
    unsigned char *origin, *data;
    unsigned char **data_blocks;
    unsigned char *zilch;
    reed_solomon *rs;

    gettimeofday(&tv, 0);
    seed = tv.tv_sec ^ tv.tv_usec;
    srandom(seed);

    fec_init();

    size = sizeof(text)/sizeof(char)-1;
    //size = 1024*1024;
    origin = malloc(size);
    memcpy(origin, text, size);
    /* for(i = 0; i < size; i++) {
        origin[i] = (unsigned char)(random() % 255);
    } */

    nrBlocks = (size+blockSize-1) / blockSize;
    nrBlocks = ((nrBlocks+dataShards-1)/dataShards) * dataShards;
    n = nrBlocks / dataShards;
    nrFecBlocks = n*parityShards;
    nrShards = nrBlocks + nrFecBlocks;
    data = malloc(nrShards * blockSize);
    memcpy(data, origin, size);
    memset(data + size, 0, nrShards*blockSize - size);
    printf("nrBlocks=%d nrFecBlocks=%d nrShards=%d n=%d left=%d\n", nrBlocks, nrFecBlocks, nrShards, n, nrShards*blockSize - size);
    print_buf(origin, "%d ", size);
    print_buf(data, "%d ", nrShards*blockSize);

    data_blocks = (unsigned char**)malloc( nrShards * sizeof(unsigned char**) );
    for(i = 0; i < nrShards; i++) {
        data_blocks[i] = data + i*blockSize;
    }

    rs = reed_solomon_new(dataShards, parityShards);
    reed_solomon_encode2(rs, data_blocks, nrShards, blockSize);
    i = memcmp(origin, data, size);
    assert(0 == i);
    //print_matrix2(data_blocks, nrShards, blockSize);

    zilch = (unsigned char*)calloc(1, nrShards);
    n = parityShards;

    int es[100];
    es[0] = 3;
    es[1] = 3;
    es[2] = 2;
    es[3] = 8;

    for(i = 0; i < n; i++) {
        //j = random() % (nrBlocks-1);
        j = es[i];
        memset(data + j*blockSize, 137, blockSize);
        zilch[j] = 1; //erased!
        printf("erased %d\n", j);
    }
    /*if(nrFecBlocks > 2) {
        for(i = 0; i < 2; i++) {
            j = nrBlocks + (random() % nrFecBlocks);
            memset(data + j*blockSize, 139, blockSize);
            zilch[j] = 1;
        }
    }*/

    reed_solomon_reconstruct(rs, data_blocks, zilch, nrShards, blockSize);
    i = memcmp(origin, data, size);
    //print_buf(origin, "%d ", nrBlocks);
    //print_buf(data, "%d ", nrBlocks);
    printf("rlt=%d\n", i);

    free(origin);
    free(data);
    free(data_blocks);
    free(zilch);
    reed_solomon_release(rs);
}

int main(void) {
    fec_init();
    //test_001();
    //test_002();
    //test_003();
    test_004();

    return 0;
}
