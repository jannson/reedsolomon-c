#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <fcntl.h>
#include "../rs.h"

#ifndef PROFILE
#ifdef __x86_64__
static long long rdtsc(void)
{
    unsigned long low, hi;
    asm volatile ("rdtsc" : "=d" (hi), "=a" (low));
    return ( (((long long)hi) << 32) | ((long long) low));
}
#elif __arm__
static long long rdtsc(void)
{
    u64 val;
    asm volatile("mrs %0, cntvct_el0" : "=r" (val));
    return val;
}
#endif
#endif

/*
 * ./fec_test dataShards parityShards blockSize input.txt
 * example: ./fec_test 10 3 256 input.txt > output.txt
 * */
int main(int argc, char **argv) {
    struct stat st;
    int size;
    struct timeval tv;
    int seed;
    int fd;
    int dataShards, parityShards, nrShards, blockSize, nrBlocks, nrFecBlocks;
    int i, n;
    int corrupted;
    reed_solomon* rs;
    unsigned long long begin, end;
    unsigned char *data;
    unsigned char **data_blocks;
    unsigned char *zilch;

    if(argc != 5) {
        fprintf(stderr, "example: ./fec_test 10 3 256 input.txt > output.txt\n");
        exit(1);
    }

    gettimeofday(&tv, 0);
    seed = tv.tv_sec ^ tv.tv_usec;
    srandom(seed);

    begin = rdtsc();
    fec_init();
    end = rdtsc();
    fprintf(stderr, "%d cycles to create FEC buffer\n", (unsigned int) (end-begin));

    fd = open(argv[4], O_RDONLY);
    if(fd < 0) {
        fprintf(stderr, "input file not found\n");
        exit(1);
    }

    dataShards = atoi(argv[1]);
    parityShards = atoi(argv[2]);
    blockSize = atoi(argv[3]);
    rs = reed_solomon_new(dataShards, parityShards);
    nrShards = rs->shards;

    fstat(fd, &st);
    size = st.st_size;
    nrBlocks = (size+blockSize-1) / blockSize;
    nrBlocks = ((nrBlocks+dataShards-1)/dataShards) * dataShards;
    i = nrBlocks / dataShards;
    nrFecBlocks = i*parityShards;
    nrShards = nrBlocks + nrFecBlocks;
    fprintf(stderr, "size=%d nr=%d\n", size, nrBlocks);

    data = (unsigned char*)malloc(nrShards * blockSize);
    n = read(fd, data, size);
    if(n < size) {
        fprintf(stderr, "Short read\n");
        close(fd);
        exit(1);
    }
    close(fd);
    memset(data+size, 0, nrShards*blockSize - size);

    data_blocks = (unsigned char**)malloc( nrShards * sizeof(unsigned char**) );
    for(i = 0; i < nrShards; i++) {
        data_blocks[i] = data + i*blockSize;
    }

    begin = rdtsc();
    reed_solomon_encode2(rs, data_blocks, nrBlocks+nrFecBlocks, blockSize);
    end = rdtsc();
    fprintf(stderr, "times %ld\n", (unsigned long) (end-begin));

    zilch = (unsigned char*)calloc(1, nrShards);
    memset(zilch, 0, sizeof(zilch));
    corrupted = parityShards;
    for(i = 0; i < corrupted; i++) {
        int corr = random() % nrBlocks;
        memset(data + blockSize * corr, 137, blockSize);
        fprintf(stderr, "Corrupting %d\n", corr);
        zilch[corr] = 1;
    }

    begin = rdtsc();
    reed_solomon_reconstruct(rs, data_blocks, zilch, nrShards, blockSize);
    end = rdtsc();
    fprintf(stderr, "times %ld\n", (unsigned long) (end-begin));

    write(1, data, size);

    free(data);
    free(data_blocks);
    free(zilch);
    reed_solomon_release(rs);

    return 0;
}

