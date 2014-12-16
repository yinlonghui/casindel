#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>

#include <bam.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>



int	usage()
{
	fprintf(stderr,"filter [option] <input.pup> <input.normal.bam> <input.tumor.bam> <out.pup> \n");
	fprintf(stderr,"    -1  [chr-col] <default:1> \n");
	fprintf(stderr,"    -2  [pos-col] <default:2> \n");
	fprintf(stderr,"    -3  [no_ref-col] <default:42> \n");
	fprintf(stderr,"    -4  [n_col]	<defualt:63>\n");
	return -1 ;
}
typedef struct {
	int	chr_col ;
	int	pos_col ;
	int	no_ref_col ;
	int	n_col ;
} opt_t ;

opt_t	*init_opt()
{
	opt_t *p = malloc(sizeof(opt_t));
	p->chr_col =  1;
	p->pos_col =  2;
	p->no_ref_col = 42 ;
	p->n_col = 63 ;
	return p ;
}

int	is_support( bam1_t *b[2] , char c ,FILE *out , int pos , int *normal_support) ;



typedef struct{
	bam1_t **b ;
	int m ;
} bam_v ;

int	main(int argc ,  char *argv[])
{
	opt_t *p =  init_opt();
	samFile    *fp[2];
	bam_hdr_t  *h[2] ; 
	hts_idx_t *idx[2] ;
	bam1_t    *b[2] ;
	FILE *in_pup , *out_pup;
	int parse_c ;

	while(( parse_c = getopt(argc,argv,"1:2:3:4:"))>= 0){
		switch(parse_c){
			case	'1':  p->chr_col = atoi(optarg); break;
			case	'2':  p->pos_col = atoi(optarg); break;
			case	'3':  p->no_ref_col = atoi(optarg); break ;
			case	'4':  p->n_col  =  atoi(optarg); break ;
		}
	}


/* 	parse cmd */

	if(argc != 4 +optind)	return usage();


	in_pup = fopen(argv[optind],"r");
	out_pup = fopen(argv[optind+3],"w");

	if(!in_pup){
		fprintf(stderr,"can't open %s\n",argv[optind]);
		return -1 ;
	}

	if(!out_pup){
		fprintf(stderr,"can't open %s\n",argv[optind+3]);
		return -1 ;
	}
	bam_v	normal;
	normal.m = 1024 ;
	normal.b = malloc(sizeof(bam1_t *)*normal.m);

	int k ;
	for(  k = 0 ; k < normal.m ; k++){
		normal.b[k] = bam_init1();
	}

/* open normal and tumor bam */

	for( k = 0 ; k < 2 ; k++){
		fp[k] =  sam_open(argv[optind + k+ 1]);
		if(!fp[k]){
			fprintf(stderr,"can't open bam %s\n",argv[1+k+optind]);
			return -1 ;
		} 
		h[k]  =  sam_hdr_read(fp[k]);
		if(!h[k]){
			fprintf(stderr,"can't open bam header %s\n",argv[1+k+optind]);
			return -1 ;

		}
		idx[k] = sam_index_load(fp[k],argv[optind+k+1]);
		if(!idx[k]){
			fprintf(stderr,"can't open bam index %s\n",argv[1+k+optind]);
			return -1 ;
		}

		b[k] = bam_init1();
	}



	while(!feof(in_pup)){
		char	buffer[1024] ,pup_c  ;
		int	i , pup_tid , pup_pos ,sp  ,support , nor_sp , nor_support ,pre_bam_pos , result  ;
		support = sp = nor_sp = nor_support = 0 ;


		for( i = 0 ; i < p->n_col ; i++){
			fscanf(in_pup,"%s\t",buffer);
			if(i+1 == p->chr_col)
				pup_tid =  bam_name2id(h[0],buffer);
			if(i+1 == p->pos_col)
				pup_pos =  atoi(buffer);
			if(i+1 == p->no_ref_col){
				pup_c = buffer[0];
			}
			fprintf(out_pup,"%s\t",buffer);
		}

		hts_itr_t *iter[2] ; 
		pup_pos--;
		iter[0] = sam_itr_queryi( idx[0], pup_tid, pup_pos, pup_pos+1 );
		iter[1] = sam_itr_queryi( idx[1], pup_tid, pup_pos, pup_pos+1 );
		k = 0 ;
		do {

			if( k ==  normal.m )	{
				normal.m =  normal.m << 1 ;
				normal.b = realloc(normal.b , sizeof(bam1_t*)*normal.m);
				for( i = k ; i < normal.m ; i++)
					normal.b[i] = bam_init1();
			}
//			fprintf(stderr,"%d\n",k);
			
			result = sam_itr_next( fp[0] , iter[0] ,normal.b[k]) ;
			k++ ;
		}while(result >= 0);

		pre_bam_pos = b[0]->core.pos = b[1]->core.pos = -1 ;
	
		if(!iter[0] || !iter[1]){
			fprintf(stderr,"can't open find index tid:%d,pos:%d\n",pup_tid,pup_pos);
			return -1 ;
		}
		int offset = 0 ;
		while((result = sam_itr_next(fp[1],iter[1],b[1]))>=0){
			int  ns , tmp_n , tmp_ns;
			tmp_n = ns = tmp_ns= 0;
			for( k = offset ;  k < normal.m ; k++)
				if(normal.b[k]->core.pos == b[0]->core.pos) break ;
			offset =  k - 1 ;
			
		
			
//			fprintf(out_pup,"%s\t",bam_get_qname(b[1]));
//
			if( k == normal.m )  b[0]->core.pos = -1 ;
			else *b[0]  =  *normal.b[k];

			do {
				tmp_n |= is_support(b,pup_c,out_pup,pup_pos,&ns);
				tmp_ns |= ns ;
				k++ ;
				if(k >= normal.m) break ;
		
			}while(normal.b[k]->core.pos > b[1]->core.pos);

			if(pre_bam_pos != b[1]->core.pos ){
				nor_sp +=  tmp_ns;
				sp +=  tmp_n ;
			}
			support += tmp_n ;
			nor_support += tmp_ns;

			if(tmp_n) pre_bam_pos  =  b[1]->core.pos;
#if 0
			fprintf(out_pup,"pup_pos: %d \t",pup_pos);
			fprintf(out_pup,"tid:%d pos: %d \t",b[1]->core.tid , b[1]->core.pos );
			fprintf(out_pup,"\n");
#endif
		}
		fprintf(out_pup,"%d\t%d\t%d\t%d\n",support,sp,nor_support,nor_sp);
		hts_itr_destroy(iter[0]);
		hts_itr_destroy(iter[1]);
	}
	for( k = 0 ; k < normal.m ;  k++)
		free(normal.b[k]);
	free(normal.b);
	for( k = 0 ; k < 2; k++){
		bam_destroy1(b[k]);
		sam_close(fp[k]);
		bam_hdr_destroy(h[k]);
		hts_idx_destroy(idx[k]);
	}
	fclose(in_pup);
	fclose(out_pup);
	free(p);
	return 0 ;
}


typedef struct{
	int	*coor;
	int	n,m ;
} mis_v ;

int	is_support( bam1_t *b[2] , char c ,FILE *out ,int pos , int *normal_support) 
{
	char MD[2] = "MD" ;

	int	k , l ,MD_num = 0 ,count =  0 ;
	int	is_SNP = 0 ;
	uint8_t  *s = bam_get_seq(b[1]);

	mis_v tumor ;
	*normal_support  =  0 ;
	tumor.m = 64 ;
	tumor.n = 0 ;
	tumor.coor = malloc(tumor.m*sizeof(int));

	uint32_t  *cigar =  bam_get_cigar(b[1]);
	for( k = 0 ; k < b[1]->core.n_cigar ; k++)
		if(bam_cigar_op(cigar[k]) == BAM_CINS || bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP || bam_cigar_op(cigar[k]) == BAM_CSOFT_CLIP){
			return 0 ;
		}
	

//
	for( k = 0 , l = 0  ; k < bam_get_l_aux(b[1]) && l < 2 ; k++){
		if(MD[l] == bam_get_aux(b[1])[k]){
			l++ ;
		}else{
			l = 0 ;
		}
	}
	k++;

	//fprintf(out,"MD:");
	


	// parse MD ..
	for(  ; k < bam_get_l_aux(b[1]) ; k++){
		if(!bam_get_aux(b[1])[k]) break; 
		char ch = bam_get_aux(b[1])[k];
		if(isalpha(ch)) {
			count+= MD_num ;
			if(!is_SNP && b[1]->core.pos + count > pos ){
				free(tumor.coor); 
				return 0 ;
			}else if( b[1]->core.pos + count == pos  && "=ACMGRSVTWYHKDBN"[bam_seqi(s,count)] == c ) {
				is_SNP =  1 ;
				continue ;
			}
			if(tumor.n == tumor.m){
				tumor.m =  tumor.m << 1 ;
				tumor.coor =  realloc(tumor.coor,sizeof(int));
			}
			tumor.coor[tumor.n] = count ;
			count++;
			tumor.n++;
			MD_num = 0;
		}else if(isalnum(ch)) {
			MD_num = MD_num*10 + ch - '0';
		}else if(ch == '^' ){
			free(tumor.coor); 
			return 0 ;
		}
	}

	if( (is_SNP && tumor.n ==0) || !is_SNP ){
		free(tumor.coor);
		*normal_support  =  is_SNP ;
		return  is_SNP;
	}

	if( b[0]->core.pos == b[1]->core.pos){
		count = 0 , MD_num = 0 ;
		for( k = 0 , l = 0  ; k < bam_get_l_aux(b[0]) && l < 2 ; k++){
			if(MD[l] == bam_get_aux(b[0])[k]){
				l++ ;
			}else{
				l = 0 ;
			}
		}
		k++;
		int  mis_addr  = 0 ;
		for(  ; k < bam_get_l_aux(b[0]) ; k++){
			if(!bam_get_aux(b[0])[k]) break; 
			char ch = bam_get_aux(b[0])[k];
			if(isalpha(ch)) {
				count += MD_num ;
				if(count != tumor.coor[mis_addr++]){
					free(tumor.coor); 
					return 0 ;
				}
				count++ ;
				MD_num = 0 ;
			}else if(isalnum(ch)) {
				MD_num = MD_num*10 + ch - '0';
			}else if(ch == '^' ){
				free(tumor.coor); 
				return 0 ;
			}
		
		}
		if( mis_addr == tumor.n )  {
			free(tumor.coor);
			*normal_support  =  1 ;
			return  0;
		}
		free(tumor.coor);
		return 0 ;
		
	}else {
		free(tumor.coor);
		return 0 ;
	}
}
