#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <errno.h>
#include "htslib/sam.h"
#include "common.h"

int main(int argc,char** argv){
    
    //check args
    if(argc!=3){
        fprintf(stderr, "Usage %s file.bam chr:start-stop\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    
    //these come from htslib/sam.h
	hts_itr_t *iter=NULL;
	hts_idx_t *idx=NULL;
	samFile *in = NULL;
	bam1_t *b= NULL;
    bam_hdr_t *header = NULL;

    //open the BAM file for reading (though called sam_open it opens bam files too :P)
    in = sam_open(argv[1], "r");
    errorCheckNULL(in);
    
    //get the sam header. 
    if ((header = sam_hdr_read(in)) == 0){
        fprintf(stderr,"No sam header?\n");
        exit(EXIT_FAILURE);
    }
    //print the chromosome names in the header
    //see the bam_hdr_t struct in htslib/sam.h for parsing the rest of the stuff in header
    int i;
    for(i=0; i< (header->n_targets); i++){
        printf("Chromosome ID %d = %s\n",i,(header->target_name[i]));
    } 
        
    //load the index file for BAM
	idx = sam_index_load(in, argv[1]);
	errorCheckNULL(idx);
    
    //the iterator for the BAM random access is probably initialised here. Note that we pass the region string to this function 
	iter  = sam_itr_querys(idx, header, argv[2]); 
	errorCheckNULL(iter);
    
    //this must be the initialisation for the structure that stores a read (need to verify)
	b = bam_init1();
    
    //my structure for a readbuffer (see common.h)
    struct bamReadBuffer readBuffer;
    readBuffer.reads.__capacity=1000;
    readBuffer.reads.__size=0;
    readBuffer.reads.array = (struct alignedRead*)malloc(sizeof(struct alignedRead)*1000);
    readBuffer.reads.windowStart=readBuffer.reads.array;
    readBuffer.reads.windowEnd=readBuffer.reads.array;
    
    //repeat until all reads in the region are retrieved
	while ( sam_itr_next(in, iter, b) >= 0){
        getRead(readBuffer.reads.windowEnd, b); //copy the current read to the myread structure. See common.c for information
        printRead(readBuffer.reads.windowEnd,header);  //print data in structure. See common.c for information
        readBuffer.reads.__size++;
        readBuffer.reads.windowEnd++;
        if(readBuffer.reads.__size==readBuffer.reads.__capacity){
            fprintf(stderr,"Buffer fill\n");
            exit(EXIT_FAILURE);
        }
        
	}
    
    
    //wrap up
    free(readBuffer.reads.array);
	hts_idx_destroy(idx);
	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);
    
    return 0;
}