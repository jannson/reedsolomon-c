# reedsolomon-c
C compatible version for [golang-klauspost-reedsolomon](https://github.com/klauspost/reedsolomon).

It can run in platform of MIPS/ARM/i386/x64 of linux, especially for some embedded projects.

It's really simple, with a little change, it can run in other platform.

# Installation
Copy rs.c and rs.h to your projects and make it.

# Usage

This section assumes you know the basics of Reed-Solomon encoding. A good start is this [Backblaze blog post](https://www.backblaze.com/blog/reed-solomon/).

To create an encoder/decoder with 10 data shards (where your data goes) and 3 parity shards (calculated):
```C
    reed_solomon* rs = reed_solomon_new(10, 3);
```

Destroy:
```C
    reed_solomon_relelase(rs);
```

Encode:

Look at [test_create_encoding](https://github.com/jannson/reedsolomon-c/blob/master/test_rs.c#L163)

```C
int test_create_encoding(
        reed_solomon *rs,
        unsigned char *data,
        int data_size,
        int block_size
        ) {
    unsigned char **data_blocks;
    int data_shards, parity_shards;
    int i, n, nr_shards, nr_blocks, nr_fec_blocks;

    data_shards = rs->data_shards;
    parity_shards = rs->parity_shards;
    nr_blocks = (data_size+block_size-1)/block_size;
    nr_blocks = ((nr_blocks+data_shards-1)/data_shards) * data_shards;
    n = nr_blocks / data_shards;
    nr_fec_blocks = n * parity_shards;
    nr_shards = nr_blocks + nr_fec_blocks;

    data_blocks = (unsigned char**)malloc(nr_shards * sizeof(unsigned char*));
    for(i = 0; i < nr_shards; i++) {
        data_blocks[i] = data + i*block_size;
    }

    n = reed_solomon_encode2(rs, data_blocks, nr_shards, block_size);
    free(data_blocks);

    return n;
}
```

Decode:

Look at [test_data_decode](https://github.com/jannson/reedsolomon-c/blob/master/test_rs.c#L192)

```C
int test_data_decode(
        reed_solomon *rs,
        unsigned char *data,
        int data_size,
        int block_size,
        int *erases,
        int erase_count) {
    unsigned char **data_blocks;
    unsigned char *zilch;
    int data_shards, parity_shards;
    int i, j, n, nr_shards, nr_blocks, nr_fec_blocks;

    data_shards = rs->data_shards;
    parity_shards = rs->parity_shards;
    nr_blocks = (data_size+block_size-1)/block_size;
    nr_blocks = ((nr_blocks+data_shards-1)/data_shards) * data_shards;
    n = nr_blocks / data_shards;
    nr_fec_blocks = n * parity_shards;
    nr_shards = nr_blocks + nr_fec_blocks;

    data_blocks = (unsigned char**)malloc(nr_shards * sizeof(unsigned char*));
    for(i = 0; i < nr_shards; i++) {
        data_blocks[i] = data + i*block_size;
    }

    zilch = (unsigned char*)calloc(1, nr_shards);
    for(i = 0; i < erase_count; i++) {
        j = erases[i];
        memset(data + j*block_size, 137, block_size);
        zilch[j] = 1; //mark as erased
    }

    n = reed_solomon_reconstruct(rs, data_blocks, zilch, nr_shards, block_size);
    free(data_blocks);
    free(zilch);

    return n;
}
```

# Performance

I have implemented a benchmarkEncode in test_rs.c. but is this implement right ? it's so quick, I think there maybe some mistake in benchmark test.

# TODO

Implement a tool to test C and golang version.

# License

published under an MIT license. See LICENSE file for more information.

