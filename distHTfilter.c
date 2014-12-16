#include <stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<math.h>

#include <bam.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>


int	usage()
{
	fprintf(stderr,"distHTfilter [option] <input.pup> <input.tumor.bam> <out.pup> \n");
	fprintf(stderr,"    -1  [chr-col] <default:1> \n");
	fprintf(stderr,"    -2  [pos-col] <default:2> \n");
	fprintf(stderr,"    -3  [no_ref-col] <default:42> \n");
	fprintf(stderr,"    -4  [n_col]	<defualt:63>\n");
	return  1 ;
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

typedef struct {
	int	*val, m , n ;
}  diff_v ;

typedef struct{
	int	beg , len ;
} mut_t ;

typedef struct{
	int	n,m;
	mut_t   *value ;
} mut_v ;

void extract_htlen(bam1_t *b , int mut_pos , char mut_c , diff_v  *d[2] , int all_repeat , int M_repeat , int  *all_isSNP , int  *M_isSNP )
{
	uint32_t *cigar = bam_get_cigar(b);
	char  MD[2] = "MD"; 
	int   i , k , l ,read_offset ,ref_offset , n_mis  ; 
	int   other_mut ;
	int   l_seq ;
	int   count  ;
	uint8_t  *s = bam_get_seq(b);

	count =  read_offset = ref_offset = other_mut = *all_isSNP =  *M_isSNP = 0  ;


	mut_v    *mi =  malloc(sizeof(mut_v));
	mi->n = 0  , mi->m = 64 ;
	mi->value =  malloc(sizeof(mut_t)*mi->m);

	mut_v    *md = malloc(sizeof(mut_v));
	md->n = 0  , md->m = 64 ;
	md->value =  malloc(sizeof(mut_t)*md->m);
			
	if(bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP){
		count = read_offset = bam_cigar_oplen(cigar[0]);
		other_mut = 1 ;
	} else  if(bam_cigar_op(cigar[0]) == BAM_CMATCH){
		count =  bam_cigar_oplen(cigar[0]);
	} else {
		other_mut = 1 ;
	}


	for( i = 1 ;  i < b->core.n_cigar ; i++){
		if(bam_cigar_op(cigar[i])== BAM_CINS){
			if(mi->n == mi->m ){
				mi->m = mi->m >> 1 ;
				mi->value = realloc(mi->value,sizeof(mut_t)*mi->m);
			}
			mi->value[mi->n].beg = count ;
			mi->value[mi->n].len = bam_cigar_oplen(cigar[i]);
			count += bam_cigar_oplen(cigar[i]);
			mi->n++;
			other_mut =  1 ;
		}else if(bam_cigar_op(cigar[i]) == BAM_CMATCH){
			count += bam_cigar_oplen(cigar[i]);
		}else if(bam_cigar_op(cigar[i]) == BAM_CDEL){
			other_mut =  1 ;
			if(md->n == md->m ){
				md->m = md->m >> 1 ;
				md->value = realloc(md->value,sizeof(mut_t)*md->m);
			}
			md->value[md->n].beg = count ;
			md->value[md->n].len = bam_cigar_oplen(cigar[i]);
			md->n++;
		}else {
			other_mut =  1 ;
			count +=  bam_cigar_oplen(cigar[i]);
		}
	}
	
	l_seq = count ;


	for( i = 0 , k = 0; i <  bam_get_l_aux(b) && k < 2 ; i++)
		if(MD[k] == bam_get_aux(b)[i]) k++ ;
		else  k = 0 ;	
	i++;
	
	count = 0 ;
	n_mis = k = l = 0 ;
		
	char ch ;
	int  is_SNPs = 0 ;
	while((ch = bam_get_aux(b)[i])!= 0 && i < bam_get_l_aux(b)){
		if(isalpha(ch)){
			read_offset += count ;
			ref_offset += count ;
			count = 0 ;
			for(  ;  k < mi->n && mi->value[k].beg > read_offset  ; k++){
				read_offset +=  mi->value[k].len ;
			}
			if( mut_pos == b->core.pos + ref_offset && "=ACMGRSVTWYHKDBN"[bam_seqi(s,read_offset)] == mut_c ){  //  is SNPs

				is_SNPs =  1 ;
				if(!all_repeat){
					if(d[0]->n == d[0]->m ){
						d[0]->m =  d[0]->m << 1 ;
						d[0]->val =  realloc(d[0]->val , d[0]->m*sizeof(int));
					}
					d[0]->val[d[0]->n] =  l_seq - read_offset <  read_offset ?  l_seq - read_offset : read_offset ;
					d[0]->n++ ;
					*all_isSNP = 1 ;
				}
			}
			read_offset++ ;
			ref_offset++  ;
			n_mis++;
			i++;
		}else if(isalnum(ch)){
			count = 10 * count +  ch - '0' ;
			i++;
		}else if( ch == '^'){
			read_offset += count ;
			ref_offset += count ;
			count = 0 ;
			ref_offset += md->value[l].len ;
			i += md->value[l].len + 1 ;
			l++;
		}
	}
	if( is_SNPs && !other_mut && n_mis == 1 && !M_repeat ){
		*M_isSNP = 1 ;
		if(d[1]->n == d[1]->m ){
			d[1]->m =  d[1]->m << 1 ;
			d[1]->val =  realloc(d[1]->val , d[1]->m*sizeof(int));
		}
		d[1]->val[d[1]->n] = d[0]->val[d[0]->n-1]  ;
		d[1]->n++ ;
	}
	
	free(mi->value);
	free(mi);
	free(md->value);
	free(md);

}

int	main( int argc , char *argv[])
{
	int	parse_c , i ;

	FILE	*in_pup , *out_pup ;
	samFile *fp ;
	bam_hdr_t *h ;
	hts_idx_t *idx ;
	bam1_t    *b ;
	opt_t	*p = init_opt() ;
	diff_v   *d[2] ;

	for( i = 0 ; i < 2 ; i++){
		d[i] = malloc(sizeof(diff_v));
		d[i]->n = 0 ;
		d[i]->m = 1024 ;
		d[i]->val= malloc(sizeof(int)*d[i]->m);
	}

	while(( parse_c = getopt(argc,argv,"1:2:3:4:"))>= 0){
		switch(parse_c){
			case	'1':  p->chr_col = atoi(optarg); break;
			case	'2':  p->pos_col = atoi(optarg); break;
			case	'3':  p->no_ref_col = atoi(optarg); break ;
			case	'4':  p->n_col  =  atoi(optarg); break ;
		}
	}


/* 	parse cmd */

	if(argc != 3 + optind)	return usage();

	in_pup = fopen(argv[optind] , "r");

	if(!in_pup){
		fprintf(stderr,"can't open input pup %s\n",argv[optind]);
		return 1 ;
	}

	out_pup = fopen(argv[optind+2],"w");
	if(!out_pup){
		fprintf(stderr,"can't open out pup %s\n",argv[optind+2]);
		return 1;
	}

	fp = sam_open(argv[optind+1]);
	if(!fp){
		fprintf(stderr,"can't open bam %s\n", argv[optind + 1]);
		return 1;
	}
	
	h = sam_hdr_read(fp);
	if(!h){
		fprintf(stderr,"can't open bam header %s\n ",argv[optind+1]);
		return 1;
	}

	idx =  sam_index_load(fp,argv[optind+1]);
	if(!idx){
		fprintf(stderr,"can't open bam index %s\n ",argv[optind+1]);
		return 1;
	}
	b = bam_init1();

	while(!feof(in_pup)){
		char	buffer[1024] , pup_c ;
		int	pup_tid , pup_pos , result , pre_pos  ;
		double  avg ,  std ;
		int	is_SNP  , M_is_SNP ;
		avg =  std = 0  ;
		is_SNP = M_is_SNP = 0 ;
		
		pup_tid = pup_pos = pre_pos =  -1 ;
		pup_c =  -1 ;
		d[0]->n = d[1]->n = 0 ;

		for( i = 0 ;  i < p->n_col ; i++){
			fscanf(in_pup,"%s\t",buffer);
			if(i+1 == p->chr_col) pup_tid = bam_name2id(h,buffer);
			if(i+1 == p->pos_col) pup_pos = atoi(buffer); 
			if(i+1 == p->no_ref_col)  pup_c = buffer[0] ;
			fprintf(out_pup,"%s\t",buffer);
		}
	
		hts_itr_t *iter ;
		pup_pos-- ;
		iter = sam_itr_queryi(idx , pup_tid , pup_pos , pup_pos+1);
		while((result = sam_itr_next(fp , iter, b))>= 0){
			int  all_repeat = ( pre_pos  ==  b->core.pos ) && is_SNP ;
			int  M_repeat  =  ( pre_pos  ==  b->core.pos ) && M_is_SNP ;
			extract_htlen( b, pup_pos, pup_c, d, all_repeat,M_repeat,&is_SNP ,&M_is_SNP);
			
		}
		
		for( i = 0 ; i < d[0]->n ; i++)
			avg +=  d[0]->val[i];	
		
		if(d[0]->n)	avg = avg/d[0]->n ;

		for( i = 0 ; i < d[0]->n ; i++)
			std +=  (d[0]->val[i] - avg ) * (d[0]->val[i] -avg) ; 
		if(d[0]->n)	std = sqrt(std/d[0]->n);
		fprintf(out_pup,"%f\t%f\t",avg ,std);
		
		avg =  std = 0 ;
		for( i = 0 ; i < d[1]->n ; i++)
			avg +=  d[1]->val[i];	
		
		if(d[1]->n)	avg = avg/d[1]->n ;

		for( i = 0 ; i < d[1]->n ; i++)
			std +=  (d[1]->val[i] - avg ) * (d[1]->val[i] -avg) ; 
		if(d[1]->n)	std = sqrt(std/d[1]->n);

		fprintf(out_pup,"%f\t%f\n",avg ,std );
		hts_itr_destroy(iter);
	}

	for( i = 0 ; i < 2 ; i++ ){
		free(d[i]->val);
		free(d[i]);
	}
	bam_destroy1(b);
	sam_close(fp);
	bam_hdr_destroy(h);
	hts_idx_destroy(idx);

	fclose(in_pup);
	fclose(out_pup);
	free(p);

	return 0 ;
}
