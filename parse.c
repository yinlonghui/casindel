#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "parse.h"

static int usage()
{
	fprintf(stderr,"casindel [option]  in.normal.bam   in.tumor.bam  out.feature\n");
	fprintf(stderr,"         -l <INT>   head/tail length");
	fprintf(stderr,"         -d <INT>   distant");
	return 0 ;
}



opt_t *parse_main( int argc  , char *argv[])
{
	opt_t  *p =  calloc(1,sizeof(opt_t));
	p->len = 10 ;
	p->dist = 5 ;
	int c  ; 
	while( (c = getopt(argc ,argv ,"l:d:")) >= 0 ){
		switch(c){
			case 'l': p->len = atoi(optarg);  break;
			case 'd': p->dist = atoi(optarg); break;
			default:  break; 
		}

	}

	if(argc  != 4 + optind) {
		free(p);
		usage();
		return  NULL;
	}
	p->fn[NORMAL] = argv[optind]; 
	p->fn[TUMOR] = argv[optind+1];
	p->fn[FEATURE] = argv[optind+2];
	p->fn[REF] = argv[optind+3];
	return  p ;
}


inline void free_opt(opt_t *p){


	free(p);
}
