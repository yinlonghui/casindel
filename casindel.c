#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
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
	int  tid ;
	char *ref ;
	char *reads ;
} cov_t ;

typedef struct {
	cov_t   *t ;
	int  n , m ;
	int  only_S ;
	int  tid ;
} cov_info ;  //  have sorted 


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
			*c  = em  ;
			tmp->p.only_S  =  (em.type == TS) ;
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
	tmp  = list->head ;
#ifdef  DEBUG
	int  x = -1  ,  y = 0  ,h  = -1  ;
#endif
	while( tmp ){
#ifdef  DEBUG
		int  i ,z = tmp->p.t[0].indel_pos ;
		int  tid =  tmp->p.tid ;
//		printf("tid:%d %d\t",tid ,z);
		for(  i  =  1 ;  i <  tmp->p.n ; i++)
			if( z != tmp->p.t[i].indel_pos){
				printf("1:e:%d\t",tmp->p.t[i].indel_pos);
				printf("1:error\n");
			}

		if( tid < h  || ( (tid == h &&  x > z ) )){
			y =  1 ;
		}
		/* 
		if(y){
			printf("2:e:%d\t%de",x,z);

		} */
		x  =  z  ;
		h  =  tid ;
#endif
		next = tmp->next ;
		free(tmp->p.t);
		free(tmp);
		tmp = next ;
	}
	list->head  =list->tail = NULL ;
#ifdef DEBUG
	//printf("\n");
	if(y) printf("error");
#endif


}

void  free_node(link_list *list , struct link_node *node)
{
	if(list->head == list->tail){
		free(node->p.t);
		free(node);
	}else  if(list->head == node) {
		list->head = node->next ;
		list->head->pre = NULL ;
		free(node->p.t);
		free(node);
	} else if(list->tail== node ) {
		list->tail =   node->pre ;
		list->tail->next = NULL ;
		free(node->p.t);
		free(node);
	}else{
		node->pre->next = node->next ;
		node->next->pre = node->pre ;
		free(node->p.t);
		free(node);
	}

}

void free_node_D(link_list_D *list , struct link_node_D *node)
{
	if(list->head == list->tail){
		free(node);
	}if(list->head == node){
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
#ifdef  DEBUG
	int  x  = -1 , y = 0 , h = -1 ;
#endif
	while( tmp ){
#ifdef  DEBUG
		int tid = tmp->p.tid ;
		if( tid < h  || (tid == h &&  tmp->p.indel_pos < x)){
			y = 1  ;
		}
		h =  tid ;
//		y += (tmp->p.indel_pos < x ) ;
//		printf("%d:%d\t",tid,tmp->p.indel_pos);
		x = tmp->p.indel_pos ;
#endif
		next = tmp->next ;
		free(tmp);
		tmp = next ;
	}
	list->head = list->tail = NULL ;
#ifdef  DEBUG
//	printf("\n");
	if(y) printf("error\n");
#endif
}
void print_list_D(link_list_D *list)
{
	struct link_node_D  *tmp  = list->head ;
	while(tmp){
		printf("%d\t",tmp->p.indel_pos);
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
	printf("%d\t%d\t%d\t%c\t%d\t%c\t%d\t",c.indel_pos,c.start_pos, c.end_pos ," SID"[c.type],c.len ,"MHT"[c.HT],c.strand);
	for( i = 0 ;  i < b->core.n_cigar ; i++){
		printf("%d",bam_cigar_oplen(cigar[i]));
		printf("%c",bam_cigar_opchr(cigar[i]));
	}
	printf("\n");

}

cov_t   bam2cov( bam1_t *b ,int n_cigar , int HTlen)
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
	
	tmp.len =  bam_cigar_oplen(cigar[n_cigar-1]) ;
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


	if(  y_pos < HTlen)  tmp.HT = TH ;
	else if(y_pos > l_qseq - HTlen) tmp.HT = TT ; 
	else  tmp.HT = TM ; 
	

	if( n_cigar == 1 || n_cigar == b->core.n_cigar){
		tmp.type = TS ;
	}else if( bam_cigar_op(cigar[n_cigar-1])&BAM_CDEL) {
		tmp.type = TD ;
	}else if( bam_cigar_op(cigar[n_cigar-1]) & BAM_CINS){
		tmp.type = TI ;
	}
	print_cov(tmp,b);
	
	return  tmp ;
}

/*
 *   This funtion generate  support features.
 */



int  casindel_core( bam_hdr_t *h[2] , samFile *fp[2],faidx_t *fai , opt_t *opt)
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
	link_list_D  lld[2] ;


	/* init  list */

	for( i = 0 ; i < 2 ; i++){
		lls[i].head =  lls[i].tail = NULL ;
		lld[i].head =  lld[i].tail = NULL ;
	}


	do{
		int sel  = 0 , indel = 0; 
		cov_t    em ;
		uint32_t *cigar = NULL;

		if( b[0]->core.tid < b[1]->core.tid  || (b[0]->core.tid == b[1]->core.tid && b[0]->core.pos <= b[1]->core.pos))	sel = 0 ;
		else sel = 1 ;
		
		cigar = bam_get_cigar(b[sel]);

		for( i = 0 ; i <  b[sel]->core.n_cigar ; i++){
			indel = bam_cigar_op(cigar[i]) & ( BAM_CINS |BAM_CDEL |BAM_CSOFT_CLIP |BAM_CHARD_CLIP) ;

			if(indel){

				em = bam2cov(b[sel],i+1,opt->len);
				
				if(bam_cigar_op(cigar[i]) & BAM_CDEL && em.len > 1) {
					insert_sort_D(&lld[sel] ,em  ,b[sel]->core.tid);  //  Insert DEL list  , and  sort 
				}



				//  Search  DEL list ,If DEL list's element coordinate region overlap the new element , both insert INDEL list. Otherwise insert  new element
				//  and  release "out of date" element  



				//  all support , to support feature 
			}
		}

		if(ret[0]  >= 0  &&  ret[1] >= 0){
			if(!sel){
				ret[0]  =  sam_read1(fp[0],h[0],b[0]);
			}else{
				ret[1]  =  sam_read1(fp[1],h[1],b[1]);
			}
		}


		if( ret[0] < 0)  ret[1]  =  sam_read1(fp[1],h[1],b[1]);
		else 		 ret[0]  =  sam_read1(fp[0],h[0],b[0]);


	}while( ret[0] >= 0 || ret[1] >= 0);

	bam_destroy1(b[0]);
	bam_destroy1(b[1]);
	return 0;
}


int  casindel(opt_t *p){
	samFile *fp[2] ;
	faidx_t  *fai ;
	bam_hdr_t *h[2] ;

//  read  bam 1 

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
	

#if 0
	printf("%s\n",faidx_fetch_seq(fai,h->target_name[b->core.tid],b->core.pos,b->core.pos+99,&len));
#endif
	return 0; 
}

void  test()
{
	link_list_D  l ;
	cov_t  em;
	int  j = 0;
	l.head = l.tail = NULL ;
	srand((unsigned long)time(0));
	while(j < 1000){
		int  i = 0  ;
		while(i< 1000){
			em.indel_pos = rand()% 100 ;
			int  tid =  rand()%5;
			//scanf("%d",&em.indel_pos);
			//printf("tid:%d ,%d:\t",tid ,em.indel_pos);
			insert_sort_D(&l,em,tid);
		//	print_list_D(&l);
		//	print_list_D1(&l);
			i++ ;
		}
//		printf("\n");
		free_list_D(&l);
		j++ ;
	}

}
void test1()
{
	link_list  l ; 
	cov_t     em ;
	int   j =  0 ;
	l.head = l.tail = NULL ;
	srand((unsigned long)time(0));
	while( j < 1000){
		int i = 0 ;
		while(i< 10000){
			em.indel_pos = rand()%100;
			//scanf("%d",&em.indel_pos);
			//scanf("%d",&tid);
			int  tid =  rand()%5;
		//	printf("tid:%d ,%d:\t",tid ,em.indel_pos);
			insert_sort(&l,em,tid);
			i++;
		}
	//	printf("\n");
		j++ ;
		free_list(&l);
	}

}

void test3()
{
	link_list  l ; 
	cov_t     em ;
	l.head = l.tail = NULL ;
	em.indel_pos = 1;
	insert_sort(&l,em,0);
	em.indel_pos = 2;
	insert_sort(&l,em,0);
	em.indel_pos = 3;
	insert_sort(&l,em,0);
	free_node(&l,l.head->next);
	free_node(&l,l.tail);
	free_node(&l,l.head);
	

}

int main(int argc ,char *argv[])
{

//	test3();
//	test();
//	test1();
#if 1
	opt_t  *p  =  parse_main(argc , argv);
	if(p == NULL) return -1 ;
	casindel(p);
	free_opt(p);
#endif 
	return 0;
}

