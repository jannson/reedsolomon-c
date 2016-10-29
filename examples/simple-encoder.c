#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <fcntl.h>
#include "../rs.h"

int main(int argc, char *argv[]) {
    struct stat st;
    int parity_shards = 0;
    int data_shards = 0;
    int nr_shards;
    int block_size;
    char* out = NULL;
    char* filename = NULL, f[256];
    int i, n;
    int fd, size;
    unsigned char* data;
    reed_solomon* rs = NULL;
    unsigned char **data_blocks = NULL;
    char output[256];

    while(-1 != (i = getopt(argc, argv, "d:p:o:f:"))) {
        switch(i) {
            case 'd':
                data_shards = atoi(optarg);
                break;
            case 'p':
                parity_shards = atoi(optarg);
                break;
            case 'o':
                out = optarg;
                break;
            case 'f':
                strcpy(f, optarg);
                filename = f;
                break;
            default:
                fprintf(stderr, "simple-encoder -d 10 -p 3 -o output -f filename.ext\n");
                exit(1);
        }
    }

    if(0 == parity_shards || 0 == data_shards || NULL == filename) {
        fprintf(stderr, "error input, example:\nsimple-encoder -d 10 -p 3 -o output -f filename.ext\n");
        exit(1);
    }
    if(out == NULL) {
        out = "./";
    }

    fec_init();

    fd = open(filename, O_RDONLY);
    if(fd < 0) {
        fprintf(stderr, "input file: %s not found\n", filename);
        exit(1);
    }

    fstat(fd, &st);
    size = st.st_size;
    block_size = (size+data_shards-1) / data_shards;
    nr_shards = data_shards + parity_shards;
    printf("filename=%s size=%d block_size=%d nr_shards=%d\n", filename, size, block_size, nr_shards);

    data = malloc(nr_shards*block_size);
    n = read(fd, data, size);
    if(n < size) {
        fprintf(stderr, "file read error!\n");
        close(fd);
        exit(1);
    }
    close(fd);

    filename = basename(filename);

    memset(data+size, 0, nr_shards*block_size - size);
    //printf("data=%s\n", data);
    data_blocks = (unsigned char**)malloc(nr_shards * sizeof(unsigned char*));
    for(i = 0; i < nr_shards; i++) {
        data_blocks[i] = data + i*block_size;
    }

    rs = reed_solomon_new(data_shards, parity_shards);
    reed_solomon_encode2(rs, data_blocks, nr_shards, block_size);

    for(i = 0; i < nr_shards; i++) {
        sprintf(output, "%s/%s.%d", out, filename, i);
        fd = open(output, O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
        write(fd, data_blocks[i], block_size);
        close(fd);
    }

    free(data_blocks);
    free(data);
    reed_solomon_release(rs);
    return 0;
}

