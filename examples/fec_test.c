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
static inline long long rdtsc(void)
{
    unsigned long low, hi;
    asm volatile ("rdtsc" : "=d" (hi), "=a" (low));
    return ( (((long long)hi) << 32) | ((long long) low));
}
#endif

/* example: ./fec_test 10 256 10 input.txt > output.txt */
int main(int argc, char **argv) {
    struct stat st;
    int size;
    int redundancy;
    int blocksize;
    int corrupted;
    unsigned char *data, *fec_data, *fec_data2;
    unsigned char *data_blocks[128], *fec_blocks[128], *dec_fec_blocks[128];
    unsigned int fec_block_nos[128], erased_blocks[128];
    int fec_size;
    int n;
    int nrBlocks;
    int i,j;
    int zilch[1024];
    struct timeval tv;
    int seed;
    long long begin, end;
    int fd;
    reed_solomon* rs;

    if(argc != 5) {
        fprintf(stderr, "example: ./fec_test 10 256 10 input.txt > output.txt\n");
        exit(1);
    }

    redundancy = atoi(argv[1]);
    blocksize = atoi(argv[2]);
    corrupted = atoi(argv[3]);
    fd = open(argv[4], O_RDONLY);
    if(fd < 0) {
        fprintf(stderr, "input file not found\n");
        exit(1);
    }

    rs = reed_solomon_new(redundancy, corrupted);
    
    gettimeofday(&tv, 0);
    seed = tv.tv_sec ^ tv.tv_usec;
    seed = 996035588;
    srandom(seed);
    fprintf(stderr, "%d\n", seed);

    fstat(fd, &st);
    size = st.st_size;
    if(size > blocksize * 128)
	size = blocksize * 128;
    nrBlocks = (size+blocksize-1) / blocksize;
    fprintf(stderr, "Size=%d nr=%d\n", size, nrBlocks);

    data = malloc(size + 4096);
    data += 4096 - ((unsigned long) data) % 4096; 
    fec_size = blocksize * redundancy;
    fec_data = malloc(fec_size + 4096);
    fec_data += 4096 - ((unsigned long) fec_data) % 4096; 
    fec_data2 = malloc(fec_size + 4096);
    fec_data2 += 4096 - ((unsigned long) fec_data2) % 4096; 
    n = read(fd, data, size);
    if(n < size) {
	fprintf(stderr, "Short read\n");
        close(fd);
	exit(1);
    }
    close(fd);

    begin = rdtsc();
    fec_init();
    end = rdtsc();
    fprintf(stderr, "%d cycles to create FEC buffer\n",
	    (unsigned int) (end-begin));
    for(i = 0, j = 0; i < size; i += blocksize, j++) {
	data_blocks[j] = data + i;
    }
    for(i = 0, j = 0; j < redundancy; i += blocksize, j++) {
	fec_blocks[j] = fec_data + i;
    }

    begin = rdtsc();
    reed_solomon_encode(rs, data_blocks, fec_blocks, blocksize);
    end = rdtsc();

    fprintf(stderr,"times %ld %f %f\n", 
	    (unsigned long) (end-begin),
	    ((double) (end-begin)) / 
	    blocksize / nrBlocks / redundancy,
	    ((double) (end-begin)) / blocksize / nrBlocks);

    fprintf(stderr, "old:\n");

    memset(zilch, 0, sizeof(zilch));
    for(i=0; i < corrupted; i++) {	
	int corr = random() % nrBlocks;
	memset(data + blocksize * corr, 137, blocksize);
	fprintf(stderr, "Corrupting %d\n", corr);
	zilch[corr] = 1;
    }

    {
	int i;
	int fec_pos = 0;
	int nr_fec_blocks = 0;

	for(i = 0; i < nrBlocks; i++) {
	    if(zilch[i]) {
		fec_pos++;
		//if(random() % 2)
		//    fec_pos++;
		erased_blocks[nr_fec_blocks] = i;
		fec_block_nos[nr_fec_blocks] = fec_pos;
		fprintf(stderr, "Fixing %d with %d (%p)\n",
			i, fec_pos, fec_blocks[fec_pos]);		
		dec_fec_blocks[nr_fec_blocks] = fec_blocks[fec_pos];
		nr_fec_blocks++;
		assert(fec_pos <= redundancy);
		assert(nr_fec_blocks <= redundancy);
	    }
	}

	begin = rdtsc();
	reed_solomon_decode(rs, data_blocks, blocksize,
		   dec_fec_blocks, fec_block_nos, erased_blocks, 
		   nr_fec_blocks);
	end = rdtsc();
	fprintf(stderr,"times %ld %f %f\n", 
		(unsigned long) (end-begin),
		((double) (end-begin)) / 
		blocksize / nrBlocks / redundancy,
		((double) (end-begin)) / blocksize / nrBlocks);	
    }

    write(1, data, size);
    exit(0);
}

