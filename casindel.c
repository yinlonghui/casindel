#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <bam.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include  <time.h>

#include "parse.h"

/*
 *   This struct link_list store perfect Mapping.
 */
typedef struct {
	int  len:16,type:16;
	int  base_q:22,map_q:10;
	int  n_mis: 16 , strand:8 , HT:8 ;
	int  start_pos ;
	int  end_pos ;
	int  indel_pos ;
	int  HT_len   ;
	int  tid ;
	// 
	int  flag_rep ;
	char *ref ;
	char *reads ;
} cov_t ;

typedef struct {
	cov_t   *t ;
	int  n , m ;
	int  only_S ;
	int  tid ;
} cov_info ;  //  have sorted 

typedef struct {
	int  n , m ;
	cov_info  *c ;
} cov_info_ar ;


struct link_node {
	struct link_node *next ;
	struct link_node *pre  ;
	cov_info   p;
};

typedef struct {
	struct link_node  *head , *tail ;
} link_list ;

struct link_node_D {
	struct link_node_D  *next ;
	struct link_node_D  *pre  ;
	cov_t   p ;
};

typedef struct {
	struct link_node_D  *head , *tail ;
} link_list_D ;

#define  TS   1 
#define  TI   2
#define  TD   3


void  insert_sort(link_list  *list , cov_t em  , int tid )
{
	if(list->head == NULL  &&  list->tail == NULL){
		struct link_node *node  = malloc(sizeof(struct link_node));
		node->next = NULL ;
		node->pre  = NULL ;
		list->head = list->tail = node ;

		node->p.m = 16 ;
		node->p.t = malloc(sizeof(cov_t)*node->p.m);
		node->p.t[0] = em ;
		node->p.n =  1 ;
		node->p.only_S = (em.type == TS ) ;
		node->p.tid =  tid ;
	} else {
		struct link_node *tmp , *node ;
		node = NULL ;
		tmp = list->tail ;
		while((( tmp->p.tid > tid ) ||(tmp->p.tid == tid && em.indel_pos <  tmp->p.t[0].indel_pos ))&& tmp->pre)  tmp = tmp->pre ;

		if(tmp->p.t[0].indel_pos == em.indel_pos && tmp->p.tid == tid ){
			if(tmp->p.n == tmp->p.m){
				tmp->p.m  =  tmp->p.m << 1 ;
				tmp->p.t  =  realloc(tmp->p.t ,sizeof(cov_t)*tmp->p.m);
			}
			cov_t  *c  =  tmp->p.t + tmp->p.n ;
			em.tid = tid ;
			*c  = em  ;
			if(tmp->p.only_S)  tmp->p.only_S  =  (em.type == TS) ;
			tmp->p.n++;
		}else{
			node = malloc(sizeof(struct link_node)) ;
			if(tmp->pre == NULL && ((em.indel_pos < tmp->p.t[0].indel_pos && tid == tmp->p.tid)||  tid < tmp->p.tid)){
				node->next = tmp ;
				node->pre  = NULL ;
				tmp->pre  =  node ;
				list->head = node ;

				node->p.m = 16 ;
				node->p.t = malloc(sizeof(cov_t)*node->p.m);
				node->p.t[0] = em ;
				node->p.n =  1 ;
				node->p.only_S = (em.type == TS ) ;
				node->p.tid =  tid ;
			}else{
				if(tmp == list->tail)   {
					list->tail = node ;
				}else{
					tmp->next->pre = node ;
				}
				node->next =  tmp->next;
				tmp->next  =  node ;
				node->pre  =  tmp  ;

				node->p.m = 16 ;
				node->p.t = malloc(sizeof(cov_t)*node->p.m);
				node->p.t[0] = em ;
				node->p.n =  1 ;
				node->p.only_S = (em.type == TS ) ;
				node->p.tid =  tid ;
			}
		}
	}
}



void   insert_sort_D(link_list_D *list , cov_t  em , int tid )
{
	if(list->head == NULL && list->tail == NULL){
		struct link_node_D  *node =  malloc(sizeof(struct link_node_D));
		node->next = NULL ;
		node->pre  = NULL ;
		em.tid = tid ;
		node->p = em ;
		list->head =  list->tail  = node ;
	}else{
		struct link_node_D  *tmp , *node ;
		node =  malloc(sizeof(struct link_node_D)) ;
		tmp  =  list->tail ;
		while (  ((em.indel_pos < tmp->p.indel_pos && tmp->p.tid == tid )||(tmp->p.tid>tid)) && tmp->pre)  tmp =  tmp->pre ;
		
		if(tmp->pre == NULL && ((em.indel_pos < tmp->p.indel_pos && tmp->p.tid == tid ) ||  tid < tmp->p.tid )){
			node->next = tmp ;
			node->pre  = NULL ;
			tmp->pre  =  node ;
			em.tid = tid ;
			node->p = em ;
			list->head = node ;

		}else{
			if(tmp == list->tail)   {
				list->tail = node ;
			}else{
				tmp->next->pre = node ;
			}
			node->next =  tmp->next;
			tmp->next  =  node ;
			node->pre  =  tmp  ;
			em.tid = tid ;
			node->p = em ;
		
			
		}
	}
}
void  free_list(link_list  *list)
{
	struct link_node  *tmp , *next ;
	int i ;
	tmp  = list->head ;
	while( tmp ){
		next = tmp->next ;
		for( i = 0 ;  i  < tmp->p.n ; i++){
			free(tmp->p.t[i].ref);
			free(tmp->p.t[i].reads);
		}
		free(tmp->p.t);
		free(tmp);
		tmp = next ;
	}
	list->head  =list->tail = NULL ;
}

void  free_node(link_list *list , struct link_node *node)
{
	int  i = 0 ;
	if(list->head == list->tail){
		for( i = 0 ;  i  < node->p.n ; i++){
			free(node->p.t[i].ref);
			free(node->p.t[i].reads);
		}
		free(node->p.t);
		free(node);
		list->head = list->tail = NULL ;
	}else  if(list->head == node) {
		list->head = node->next ;
		list->head->pre = NULL ;
		for( i = 0 ;  i  < node->p.n ; i++){
			free(node->p.t[i].ref);
			free(node->p.t[i].reads);
		}
		free(node->p.t);
		free(node);
	} else if(list->tail== node ) {
		list->tail =   node->pre ;
		list->tail->next = NULL ;
		for( i = 0 ;  i  < node->p.n ; i++){
			free(node->p.t[i].ref);
			free(node->p.t[i].reads);
		}
		free(node->p.t);
		free(node);
	}else{
		node->pre->next = node->next ;
		node->next->pre = node->pre ;
		for( i = 0 ;  i  < node->p.n ; i++){
			free(node->p.t[i].ref);
			free(node->p.t[i].reads);
		}
		free(node->p.t);
		free(node);
	}

}

void free_node_D(link_list_D *list , struct link_node_D *node)
{
	if(list->head == list->tail){
		free(node);
		list->head = list->tail = NULL ;
	}else if(list->head == node){
		list->head = node->next ;
		list->head->pre = NULL ;
		free(node);
	}else if(list->tail== node ){
		list->tail =   node->pre ;
		list->tail->next = NULL ;
		free(node);
	}else{
		node->pre->next = node->next ;
		node->next->pre = node->pre ;
		free(node);
	}

}


void  free_list_D(link_list_D *list)
{
	struct link_node_D  *tmp , *next ;
	tmp =  list->head ;
	while( tmp ){
		next = tmp->next ;
		free(tmp);
		tmp = next ;
	}
	list->head = list->tail = NULL ;
}

void print_list_D(link_list_D *list)
{
	struct link_node_D  *tmp  = list->head ;
	while(tmp){
		printf("%d:%d-%d\t",tmp->p.tid ,tmp->p.indel_pos,tmp->p.indel_pos + tmp->p.len);
		tmp = tmp->next ;
	}
	printf("\n");

}

void print_list_D1(link_list_D *list)
{
	struct link_node_D  *tmp  = list->tail;
	while(tmp){
		printf("%d\t",tmp->p.indel_pos);
		tmp = tmp->pre;
	}
	printf("\n");

}


#define  TM   0
#define  TH   1 
#define  TT   2
void    print_cov(cov_t c , bam1_t *b)
{
	int  i ;
	uint32_t *cigar = bam_get_cigar(b);
	if( c.ref && c.reads )
	printf("%s\t%d\t%d\t%d\t%c\t%d\t%c\t%d\t%d\t%s\t%s\t%d\t%d\t",bam_get_qname(b),c.indel_pos,c.start_pos, c.end_pos ," SID"[c.type],c.len ,"MHT"[c.HT],c.strand,c.n_mis,c.reads,c.ref,c.base_q,c.map_q);
	else  printf("%s\t%d\t%d\t%d\t%c\t%d\t%c\t%d\t%d\t%d\t%d\t",bam_get_qname(b),c.indel_pos,c.start_pos, c.end_pos ," SID"[c.type],c.len ,"MHT"[c.HT],c.strand,c.n_mis,c.base_q,c.map_q);

	for( i = 0 ;  i < b->core.n_cigar ; i++){
		printf("%d",bam_cigar_oplen(cigar[i]));
		printf("%c",bam_cigar_opchr(cigar[i]));
	}
	printf("\n");

}


cov_t   bam2cov( bam1_t *b ,int n_cigar , int HTlen , faidx_t *fai , bam_hdr_t *h )
{
	cov_t   tmp ;
	int  i , l_qseq , pos , y_pos ,t_pos ;
	uint32_t *cigar = bam_get_cigar(b);

	l_qseq = pos = y_pos = t_pos =  0 ;

	for( i = 0 ; i < n_cigar ; i++){
		if(bam_cigar_op(cigar[i])&(BAM_CINS|BAM_CSOFT_CLIP |BAM_CHARD_CLIP) || bam_cigar_op(cigar[i]) == BAM_CMATCH )  l_qseq +=  bam_cigar_oplen(cigar[i]) ;
		if(bam_cigar_op(cigar[i])&(BAM_CDEL) || bam_cigar_op(cigar[i]) == BAM_CMATCH ) pos +=  bam_cigar_oplen(cigar[i]) ;
		if(bam_cigar_op(cigar[i])&(BAM_CINS|BAM_CSOFT_CLIP |BAM_CHARD_CLIP) || bam_cigar_op(cigar[i]) == BAM_CMATCH )  y_pos +=  bam_cigar_oplen(cigar[i]) ;
	}
	t_pos = pos ;
	for(  ; i <  b->core.n_cigar  ; i++ ){
		if(bam_cigar_op(cigar[i])&(BAM_CINS|BAM_CSOFT_CLIP |BAM_CHARD_CLIP) || bam_cigar_op(cigar[i]) == BAM_CMATCH )  l_qseq +=  bam_cigar_oplen(cigar[i]) ;
		if(bam_cigar_op(cigar[i])&(BAM_CDEL) || bam_cigar_op(cigar[i]) == BAM_CMATCH ) t_pos +=  bam_cigar_oplen(cigar[i]) ;
	}
	
	tmp.map_q  =  b->core.qual ;
	tmp.base_q = 0 ;

	/* avg base quality */
	uint8_t   *bq =   bam_get_qual(b);

	for( i = 0  ; i < b->core.l_qseq ; i++)
		tmp.base_q += bq[i];

	tmp.base_q /=  b->core.l_qseq ;
	tmp.strand = bam_is_rev(b);
	tmp.indel_pos = b->core.pos + pos ;
	tmp.start_pos = b->core.pos ;
	tmp.end_pos = b->core.pos + t_pos ;


	

	char  NM[2] = "NM";
	  
	int  j = 0 ;
	
	for( i = 0  , j = 0; i < bam_get_l_aux(b) && j < 2  ; i++){
		if(NM[j] == bam_get_aux(b)[i]){
			j++;
		}
	}
	tmp.n_mis = bam_get_aux(b)[i+1];

	if(n_cigar == 0){
		tmp.len = l_qseq ;
		tmp.ref = NULL;
		tmp.reads = NULL;
		return tmp ;
	}
	
	tmp.flag_rep = 1 ;
	tmp.len =  bam_cigar_oplen(cigar[n_cigar-1]) ;

	if(  y_pos < HTlen)  tmp.HT = TH ;
	else if(y_pos > l_qseq - HTlen) tmp.HT = TT ; 
	else  tmp.HT = TM ; 

	if(  y_pos < l_qseq/2)   tmp.HT_len = y_pos ;
	else			 tmp.HT_len = l_qseq - y_pos ;

	if(n_cigar == b->core.n_cigar )  tmp.indel_pos--;
	
	
	if( n_cigar == 1 || n_cigar == b->core.n_cigar){
		tmp.type = TS ;
	}else if( bam_cigar_op(cigar[n_cigar-1])&BAM_CDEL) {
		tmp.type = TD ;
	}else if( bam_cigar_op(cigar[n_cigar-1]) & BAM_CINS){
		tmp.type = TI ;
	}

	// char  to char   
	if( n_cigar == 1 || n_cigar == b->core.n_cigar){
		tmp.ref = NULL;
		tmp.reads = NULL ;
	}else if(bam_cigar_op(cigar[n_cigar-1])&BAM_CINS) {
		tmp.reads = malloc(sizeof(char)*(tmp.len+2));
		uint8_t  *s =  bam_get_seq(b);
		int len ;
		for( j = 0 ; j  < tmp.len + 1 ; j++)
			tmp.reads[j] = "=ACMGRSVTWYHKDBN"[bam_seqi(s, j + y_pos - tmp.len - 1)];
		tmp.reads[j] = 0 ; 
		tmp.ref = faidx_fetch_seq(fai,h->target_name[b->core.tid],tmp.indel_pos-1,tmp.indel_pos-1,&len);
	}else if( bam_cigar_op(cigar[n_cigar-1])&BAM_CDEL){

		tmp.reads = malloc(sizeof(char)*2);
		uint8_t  *s =  bam_get_seq(b);
		int len ;
		tmp.reads[0] = "=ACMGRSVTWYHKDBN"[bam_seqi(s,  y_pos - 1)];
		tmp.reads[1] = 0 ;
		tmp.ref = faidx_fetch_seq(fai,h->target_name[b->core.tid],tmp.indel_pos-tmp.len -1,tmp.indel_pos-1,&len);
		tmp.indel_pos -= tmp.len ;
	}
//	print_cov(tmp,b);
	return  tmp ;
}

int  search_list_D( link_list_D  *list , cov_t  em ,  cov_info *p)
{
	struct link_node_D  *node = list->head  ;
	p->n =  0 ;
	p->m =  16 ;
	p->t = malloc(sizeof(cov_t)*p->m);

	while(node){
		if(((node->p.indel_pos < em.indel_pos && node->p.indel_pos + node->p.len > em.indel_pos)||node->p.indel_pos == em.indel_pos) && em.tid == node->p.tid \
			&&!(node->p.flag_rep &(1<<(em.indel_pos-node->p.indel_pos)))){

			if(p->m == p->n){
				p->m =  p->m << 1 ;
				p->t = realloc(p->t , sizeof(cov_t)*p->m);
			}
			node->p.flag_rep  |=  (1<< (em.indel_pos - node->p.indel_pos));
			cov_t  *emd =  p->t + p->n ;
			*emd = node->p ;
			emd->ref = emd->reads = NULL ;
			emd->indel_pos =  em.indel_pos ;
			p->n++;
		}
		node = node->next ;
	}
	if(!p->n) free(p->t);
	
	return p->n ;
}

int search_out_D(link_list_D *list , cov_t em ,  cov_info *p )
{
	struct link_node_D  *node = list->head  ;
	p->n =  0 ;
	p->m =  16 ;
	p->t = malloc(sizeof(cov_t)*p->m);

	while(node){
		if(em.tid ==  node->p.tid && em.indel_pos < node->p.end_pos  && em.indel_pos > node->p.start_pos){

			if(p->m == p->n){
				p->m =  p->m << 1 ;
				p->t = realloc(p->t , sizeof(cov_t)*p->m);
			}
			cov_t  *emd =  p->t + p->n ;
			*emd = node->p ;
			emd->ref = emd->reads = NULL ;
			emd->indel_pos =  em.indel_pos ;
			p->n++;
		}
		node = node->next ;
	}
	if(!p->n) free(p->t);
	
	return p->n ;

}

struct   link_node *search_node_ ( struct link_node *node , cov_t em , int dist)
{
	struct link_node *n = node , *rt = NULL ;
	int	min =  dist ;
	while(n){
		if(em.tid == n->p.t[0].tid  && abs(em.indel_pos- n->p.t[0].indel_pos) < min){
			rt =  n  ;
			min  = abs(em.indel_pos- n->p.t[0].indel_pos) ;
		}
		if(em.tid < node->p.t[0].tid || em.indel_pos + min < node->p.t[0].indel_pos)  return rt ;
		n = n->next ;
	}
	return rt ;
}

struct   link_node *search_node(link_list  *list , cov_t em , int dist )
{
	struct link_node  *node = list->head;
	while(node){
		if(em.tid == node->p.t[0].tid  && abs(em.indel_pos-node->p.t[0].indel_pos) < dist) return  node ;
		if(em.tid < node->p.t[0].tid || em.indel_pos + dist < node->p.t[0].indel_pos)  return NULL ;
		node = node->next ;
	}
	return  NULL;
}

void out_data_link_D( link_list_D *list , cov_t em , int sel)
{
	struct link_node_D  *node  = list->head ;
	while(node){
		struct link_node_D *next =  node->next ;
		int len =  sel ? node->p.len + 101 : 0 ;
		if((node->p.indel_pos + len  < em.start_pos && em.tid == node->p.tid) ||  em.tid > node->p.tid){
			//if(sel) printf("De:%d:%d\n",node->p.tid,node->p.start_pos);
			free_node_D(list,node);
		}else{
			return  ;
		}
		node = next ;
	}
}

struct  link_node  *out_data_link(link_list  *list , cov_t em)
{
	struct link_node  *node = list->head;
	while(node){
		struct link_node *next =  node->next ;
		if((node->p.t[0].indel_pos + 101  < em.start_pos && em.tid == node->p.t[0].tid) ||  em.tid > node->p.t[0].tid){
			//printf("De:%d:%d\n",node->p.tid,node->p.start_pos);
			return node;
		}else{
			return NULL;
		}
		node = next ;
	}
	return NULL ;
}


typedef struct {
	int  type ;
	char	*bg ;
	int  n ;
}  t_cmp ;

int   cmp_type( int *n , t_cmp *s , char *b , int type)
{
	if( *n == 0) {
		s->bg  = strdup(b);
		s->n = 1;
		s->type = type;
		(*n)++;
		return  0 ;
	}
	int  i = 0; 

	for( i = 0 ; i < *n ; i++){
		t_cmp *p = s + i ;
		if(strcmp(p->bg,b) == 0 && p->type == type){
			p->n++;
			return 1 ;
		}
	}
	if( i == *n ){
		s[i].bg = strdup(b);
		s[i].n = 1 ;
		s[i].type = type;
		(*n)++;
	}
	return 0 ;
}


t_cmp  cmp_indel(cov_info p)
{
	t_cmp  tmp , *s ;
	s =  malloc(p.n*sizeof(t_cmp));
	int  i , n = 0 ;

	for( i = 0 ; i < p.n ; i++){
		cov_t  *ptr =  p.t + i ;
		if(ptr->ref == NULL){
			continue ;
		}if(ptr->type == TD) {
			cmp_type(&n,s,ptr->ref,TD);
		}else if( ptr->type == TI) {
			cmp_type(&n,s,ptr->reads,TI);
		}
	}

	t_cmp  *high  =  s ;
	for( i = 1 ;  i < n ; i++){
		t_cmp *ptr = s + i ;
		if(ptr->n > high->n)  high = ptr ;
	}

	tmp.n =  high->n ;
	tmp.type = high->type ;
	tmp.bg  =  strdup(high->bg);

	for( i = 0 ;  i < n ; i++){
		t_cmp *ptr = s + i ;
		free(ptr->bg);
	}
	free(s);
//	printf("%s \t %d\n",tmp.bg, tmp.n);

	return  tmp ;
}

typedef struct {
	int	tid ;
	int	pos ;
	int	type ;
	char	*ref ,*reads ;
	int	support ;
	int	len ;
	int	rev_num ;
	int	pos_num ;
	int	n_mismatch;
	int	H_num ;
	int	T_num ;
	int	map_q ;
	int	base_q ;
	int	start_point;
	int     start_pos ;
	int	HT_len ;
//   noise ...
	int     noise ;
	int     N_mis ;
	int  	N_map_q ;
	int	N_base_q ;
//   ext
	int	M_num ;
	int	M_map_q ;
	int 	M_base_q ;
	int	M_mis ;
} feature_t ;


void  cov2feature( int sel , feature_t *t, cov_t *c )
{
	if(sel == 1){// support 
		if(t->start_pos == -1){
			t->ref =  strdup(c->ref);
			t->reads = strdup(c->reads);
			t->tid = c->tid ;
			t->len = c->len ;
			t->pos = c->indel_pos ;
		}
		if(c->start_pos != t->start_pos){
			t->start_pos = c->start_pos ;
			t->start_point++ ;
		}
		if(c->strand)  t->rev_num++ ;
		else   t->pos_num++ ;

		if(c->HT == TH) t->H_num++ ; 
		else  if(c->HT == TT) t->T_num++ ;
		t->HT_len += c->HT_len ;
		
		t->n_mismatch += c->n_mis ;
		t->map_q += c->map_q ;
		t->base_q+= c->base_q;
		t->support++ ;
	
	}else if(sel == 2){
		t->noise++;
		t->N_mis += c->n_mis;
		t->N_map_q += c->map_q ;
		t->N_base_q +=  c->base_q;
	}else{
		t->M_num++ ; 
		t->M_mis += c->n_mis;
		t->M_map_q += c->map_q ;
		t->M_base_q +=  c->base_q;
	}
}



void classify_indel( struct  link_node *node , cov_info s,bam_hdr_t *h)
{
	if(!node){
		feature_t  *p = calloc(1,sizeof(feature_t));
		printf("-\t0\t-\t0\t-\t-\t0\t0\t0\t0\t");
		printf("0\t0\t0\t0\t0\t");
		int  i = 0 ;
		for( i = 0 ; i < s.n ; i++){
			cov_t  *c =  s.t + i ;
			cov2feature(0,p,c);
		}
		printf("0\t0\t");
		printf("0\t0\t0\t0\t");
		if(p->M_num) printf("%d\t%f\t%f\t%f",p->M_num,(double)((double)p->M_mis/p->M_num), (double)(p->M_base_q/p->M_num) , (double)(p->M_map_q/p->M_num));
		else  printf("0\t0\t0\t0");
		free(p);
	}else if(node->p.only_S){
		feature_t  *p = calloc(1,sizeof(feature_t));
		printf("-\t0\t-\t0\t-\t-\t0\t0\t0\t0\t");
		printf("0\t0\t");
		printf("0\t0\t0\t0\t0\t");
		int i = 0 ;
		for( i = 0 ; i < node->p.n ; i++){
			cov_t  *c =  node->p.t + i ;
			cov2feature(2,p,c);
		}
		if(p->noise) printf("%d\t%f\t%f\t%f\t",p->noise ,(double)((double)p->N_mis/p->noise), (double)((double)p->N_base_q/p->noise) , (double)((double)p->N_map_q/p->noise) );
		else  printf("0\t0\t0\t0\t");
		if(p->M_num) printf("%d\t%f\t%f\t%f",p->M_num,(double)((double)p->M_mis/p->M_num), (double)(p->M_base_q/p->M_num) , (double)(p->M_map_q/p->M_num));
		else  printf("0\t0\t0\t0");
		free(p);
	}else{
		t_cmp  tmp =  cmp_indel(node->p);
		feature_t  *p = calloc(1,sizeof(feature_t));

		p->start_pos = -1 ;
		p->type = tmp.type ;

		int i = 0 ;
		for( i = 0 ; i < node->p.n ; i++){
			cov_t  *c =  node->p.t + i ;
			if(c->ref == NULL){
				cov2feature(2,p,c);
			}else if(c->type == TI && tmp.type == c->type && strcmp(c->reads,tmp.bg) == 0){
				cov2feature(1,p,c);
			}else if(c->type == TD && tmp.type == c->type && strcmp(c->ref , tmp.bg) == 0){
				cov2feature(1,p,c);
			}else {
				cov2feature(2,p,c);
			}

		}
		for( i = 0 ; i < s.n ; i++){
			cov_t  *c =  s.t + i ;
			cov2feature(0,p,c);
		}
		printf("%s\t%d\t%c\t%d\t%s\t%s\t%d\t",h->target_name[p->tid], p->pos ," SID"[p->type],p->len,p->ref,p->reads,p->support);
		if(p->support) {
			printf("%f\t%f\t",(double)((double)(abs(p->rev_num - p->pos_num))/p->support),(double)((double)p->HT_len/p->support));
			printf("%f\t%f\t%f\t",(double)((double)p->n_mismatch/p->support),(double)((double)p->base_q/p->support),(double)((double)p->map_q/p->support));
		}
		else  {
			printf("0\t0\t");
			printf("0\t0\t0\t");
		}

		printf("%d\t%d\t%d\t%d\t%d\t",p->pos_num,p->rev_num,p->H_num,p->T_num,p->start_point);
		if(p->noise) printf("%d\t%f\t%f\t%f\t",p->noise ,(double)((double)p->N_mis/p->noise), (double)((double)p->N_base_q/p->noise) , (double)((double)p->N_map_q/p->noise) );
		else  printf("0\t0\t0\t0\t");
		if(p->M_num) printf("%d\t%f\t%f\t%f",p->M_num,(double)((double)p->M_mis/p->M_num), (double)(p->M_base_q/p->M_num) , (double)(p->M_map_q/p->M_num));
		else  printf("0\t0\t0\t0");
		free(p->ref);
		free(p->reads);
		free(p);
		free(tmp.bg);
	}
}


/*
 *   This funtion generate  support features.
 */

int  casindel_core( bam_hdr_t *h[2] , samFile *fp[2], faidx_t *fai , opt_t *opt )
{
	bam1_t *b[2];
	b[0] = bam_init1();
	b[1] = bam_init1();
	int ret[2] , i ;
	ret[0]  =  ret[1]  = 0 ;

	if(ret[0] >= 0) ret[0]  =  sam_read1(fp[0],h[0],b[0]);
	if(ret[1] >= 0) ret[1]  =  sam_read1(fp[1],h[1],b[1]);


	if(ret[0] < 0 || ret[1] < 0)  return -1 ;

	link_list   lls[2] ;
	link_list_D  lld[4] ;


	/* init  list */

	for( i = 0 ; i < 2 ; i++){
		lls[i].head =  lls[i].tail = NULL ;
		lld[i].head =  lld[i].tail = NULL ;
	}
	for(  ;i < 4 ;i++) lld[i].head =  lld[i].tail = NULL ;


	do{
		int sel  = 0 , indel = 0 ,is_indel = 0; 
		cov_t    em ;
		uint32_t *cigar = NULL;

		if( b[0]->core.tid < b[1]->core.tid  || (b[0]->core.tid == b[1]->core.tid && b[0]->core.pos <= b[1]->core.pos))	sel = 0 ;
		else sel = 1 ;
		
		cigar = bam_get_cigar(b[sel]);

		for( i = 0 ; i <  b[sel]->core.n_cigar ; i++){
			cov_info  p ;
			indel = bam_cigar_op(cigar[i]) & ( BAM_CINS |BAM_CDEL |BAM_CSOFT_CLIP |BAM_CHARD_CLIP) ;

			if(indel){

				is_indel = 1 ;
				em = bam2cov(b[sel],i+1,opt->len,fai, h[sel]);
				em.tid  = b[sel]->core.tid ;
				insert_sort(&lls[sel], em ,b[sel]->core.tid);

				if(bam_cigar_op(cigar[i]) & BAM_CDEL && em.len > 1) {
					insert_sort_D(&lld[sel] ,em ,b[sel]->core.tid);  //  Insert DEL list  , and  sort 

				}


				//  Search  DEL list ,If DEL list's element coordinate region overlap the new element , both insert INDEL list. Otherwise insert  new element
				if(em.type !=TS){
					if(search_list_D( &lld[sel] , em , &p)){  //  must change D coordinate .
						int  j ; 
						for( j =  0 ; j < p.n ; j++)
							insert_sort(&lls[sel],p.t[j],p.t[j].tid);
						free(p.t);
					}
				}

				//  Release "out of date" element  
				if(!(bam_cigar_op(cigar[i]) & BAM_CDEL)) {
					out_data_link_D(&lld[sel],em,0);
				}
				struct  link_node  *node , *tmp , *n = lls[0].head ;
				tmp = NULL ;
				if(sel == 1){
					while((node =out_data_link(&lls[sel],em))){
						//printf("*********:%d\n",node->p.t[0].indel_pos);
						if(node->p.only_S) {
							free_node(&lls[sel],node);
							continue ;
						}
						tmp = search_node_(n,node->p.t[0],opt->dist);
	//					tmp = search_node(&lls[1-sel],node->p.t[0] ,opt->dist);
						search_out_D(&lld[sel+2],node->p.t[0],&p);
#if 0
						int j = 0 ;
						for( j = 0 ; j < p.n ; j++ )
							printf("%d:%d-%d %d\n",sel,p.t[j].start_pos,p.t[j].end_pos,p.t[j].n_mis);
#endif
						classify_indel(node,p,h[1]);
						if(p.n) free(p.t);
						printf("\t");
						out_data_link_D(&lld[sel+2],node->p.t[0],1);
						
						search_out_D(&lld[3-sel],node->p.t[0],&p);
						classify_indel(tmp,p,h[0]);
						printf("\n");
						out_data_link_D(&lld[3-sel],node->p.t[0],1);
						if(p.n) free(p.t);
						free_node(&lls[1],node);
				//		n = tmp ;
						//if(tmp) free_node(&lls[0],tmp);
						//  output support 
					}
					while((node = out_data_link(&lls[0],em)))  free_node(&lls[0],node);
//					while((node = out_data_link(&lls[1],em)))  free_node(&lls[1],node);

				}

			}
		}

		if(!is_indel){
			em = bam2cov(b[sel],0,0,fai,h[sel]);
			insert_sort_D(&lld[sel+2],em,b[sel]->core.tid);
		}

		if(ret[0]  >= 0  &&  ret[1] >= 0){
			if(!sel){
				ret[0]  =  sam_read1(fp[0],h[0],b[0]);
			}else{
				ret[1]  =  sam_read1(fp[1],h[1],b[1]);
			}
		}


		if( ret[0] < 0)  ret[1]  =  sam_read1(fp[1],h[1],b[1]);
		else  if(ret[1] < 0)		 ret[0]  =  sam_read1(fp[0],h[0],b[0]);

	}while( ret[0] >= 0 || ret[1] >= 0);
/*
 * 	have  none  
 */
	struct  link_node  *node =  lls[1].head , *tmp , *n = lls[0].head ;
	tmp = NULL ;
	while(node){
		if(node->p.only_S){
			node = node->next ;
			continue ;
		}else{
			cov_info  p ;
			tmp = search_node_(n,node->p.t[0],opt->dist);
			search_out_D(&lld[3],node->p.t[0],&p);
			classify_indel(node,p,h[1]);
			printf("\t");
			if(p.n)  free(p.t);
			search_out_D(&lld[2],node->p.t[0],&p);
			classify_indel(tmp,p,h[0]);
			printf("\n");
			if(p.n)  free(p.t);
//			n = tmp ;
		}
		node = node->next ;
	}



	free_list(&lls[0]);
	free_list(&lls[1]);
	free_list_D(&lld[0]);
	free_list_D(&lld[1]);
	free_list_D(&lld[2]);
	free_list_D(&lld[3]);

	bam_destroy1(b[0]);
	bam_destroy1(b[1]);
	return 0;
}


int  casindel(opt_t *p){
	bam_hdr_t *h[2] ;
	samFile *fp[2] ;
	faidx_t  *fai ;
//  open  bam..
	fp[0] = sam_open(p->fn[0]);
	fp[1] = sam_open(p->fn[1]);
	h[0] = sam_hdr_read(fp[0]);
	h[1] = sam_hdr_read(fp[1]);


	fai  = fai_load(p->fn[3]);
	casindel_core(h,fp,fai,p);
	bam_hdr_destroy(h[0]);
	bam_hdr_destroy(h[1]);
	sam_close(fp[0]);
	sam_close(fp[1]);
	fai_destroy(fai);
	
	return 0; 
}

int main(int argc ,char *argv[])
{

	opt_t  *p  =  parse_main(argc , argv);
	if(p == NULL) return -1 ;
	
	casindel(p);
	free_opt(p);
	return 0;
}

