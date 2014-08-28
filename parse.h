#ifndef __PARSE_H
#define __PARSE_H

#define NORMAL 0
#define TUMOR  1 
#define FEATURE 2 
#define REF	3


typedef struct {
	char *fn[4] ;
	int  len ;
} opt_t ;

opt_t *parse_main( int argc  , char *argv[]);

inline void free_opt(opt_t *p);


#endif
