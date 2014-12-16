#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>

#include <bam.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>


static int usage()
{
	fprintf(stderr," readsPairCheck [option]  <in.input>  <file.bam> [file2.bam] <out.put> \n");
	fprintf(stderr,"                -d <INT>   [default:500] \n");

	return 1;
}

typedef struct {
	int	dist ;
}  opt_t ;

typedef struct{
	int	tid ;
	int	pos ;
	int	direct ;
} sv_inf;


int main(int argc , char *argv[])
{

	
	int parse_c ; 
	
	FILE	*in,*out;
	
	samFile  **fp ; 
	bam_hdr_t  **header ;
	hts_idx_t  **idx ;
	bam1_t	  **b ;
	sv_inf   *sv ;
	hts_itr_t **iter ;
	
	int	i , n_platform  ;

	opt_t  p ;
	p.dist = 500 ;

	while(( parse_c = getopt(argc,argv,"d:"))>= 0){
		switch(parse_c){
			case	'd':  p.dist= atoi(optarg); break;
		}
	}

	if(argc != optind + 4  &&  argc !=  optind + 3)  return   usage();
	
	if( argc == optind + 3 )  n_platform =  1 ; 
	else if( argc == optind + 4) n_platform = 2 ;


	in = fopen( argv[optind] , "r");

	if(!in){
		fprintf(stderr,"can't open input file: %s \n",argv[optind]);
		return  1;
	}
	
	out = fopen( argv[optind+n_platform+1],"w") ;

	if(!out){
		fprintf(stderr,"can't open output file:%s \n",argv[optind]);
		return  1;
	}


	fp =  malloc(sizeof(samFile *)*n_platform);
	header = malloc(sizeof(bam_hdr_t *)*n_platform);
	idx =  malloc(sizeof(hts_idx_t *)*n_platform);
	b  =  malloc(sizeof(bam1_t *)*n_platform);
	iter = malloc(sizeof(hts_itr_t *)*n_platform);

	sv = malloc(sizeof(sv_inf)*2); 
	for( i = 0 ; i < n_platform ; i++){
		fp[i] =  sam_open(argv[optind+1+i]);
		if(!fp[i]){
			fprintf(stderr,"can't open input bam file: %s \n",argv[optind+1+i]);
			return  1;
		}
		header[i] = sam_hdr_read(fp[i]);
		if(!header[i]){
			fprintf(stderr,"can't open input bam's header: %s \n",argv[optind+1+i]);
			return  1;
		} idx[i] = sam_index_load(fp[i] , argv[optind+1+i]);
		if(!idx[i]){
			fprintf(stderr,"can't open input bam's index: %s \n",argv[optind+1+i]);
			return  1;
		}
		b[i] = bam_init1();
	}

	while(!feof(in)){
		char	buffer[1024] , c ;  // temporary various  
		int	result ;
	

		fscanf(in,"%s\t%d\t%c\t",buffer,&sv[0].pos,&c);
		sv[0].tid = bam_name2id(header[0],buffer);
		sv[0].direct  =  c == 'R' ? 1 : 0 ;
		fscanf(in,"%s\t%d\t%c\n",buffer,&sv[1].pos,&c);
		sv[1].tid = bam_name2id(header[0],buffer);
		sv[1].direct  =  c == 'R' ? 1 : 0 ;

		for( i = 0 ;  i < n_platform ; i++){
			int	rev = 0 , pos = 0 ;
		
			if(sv[0].direct == 1 )	 iter[i] = sam_itr_queryi( idx[i], sv[0].tid , sv[0].pos-1, sv[0].pos + p.dist -1);
			else	iter[i] = sam_itr_queryi( idx[i], sv[0].tid , sv[0].pos - p.dist -1, sv[0].pos - 1);
			while((result = sam_itr_next(fp[i],iter[i],b[i]))>=0){
				uint32_t   *cigar = bam_get_cigar(b[i]);
				int k = 0 ;
				int len = 0 ;
			
				if((sv[0].direct == 0 && b[i]->core.pos < sv[0].pos - p.dist - 1) || (sv[0].direct == 1 && b[i]->core.pos < sv[0].pos -1))  continue ;
			//     <=====================|	
				if( sv[0].direct == 0 ){
					for( k = 0 ; k < b[i]->core.n_cigar - 1 ; k++)
						if(bam_cigar_op(cigar[k])&BAM_CDEL || bam_cigar_op(cigar[k]) == BAM_CMATCH)  len+= bam_cigar_oplen(cigar[k]);
					len +=  bam_cigar_oplen(cigar[k]);
					if( b[i]->core.pos + len >  sv[0].pos - 1 ) continue ;
				}else {  // |===================>
					if( bam_cigar_op(cigar[0])&(BAM_CDEL & BAM_CINS & BAM_CSOFT_CLIP &BAM_CHARD_CLIP) &&  b[i]->core.pos - bam_cigar_op(cigar[0]) < sv[0].pos -1) continue ;

				}

				if( sv[1].tid ==  b[i]->core.mtid )
					if((sv[1].direct == 0 && b[i]->core.mpos > sv[1].pos - p.dist - 1 && b[i]->core.mpos < sv[1].pos -1 ) || (sv[1].direct == 1 && b[i]->core.mpos > sv[1].pos -1 && b[i]->core.mpos <  sv[1].pos +  p.dist -1 )){
						//  flags 
						int	a[2];
						a[0] = b[i]->core.flag & BAM_FREVERSE ;
						a[1] = b[i]->core.flag & BAM_FMREVERSE ; 
						if(a[0] == a[1])  pos++ ;
						else	rev++ ;
					}

			
			}
			fprintf(out,"%d\t%d",pos,rev);
			if( i == 0 )  fprintf(out,"\t");
			else	fprintf(out,"\n");
			hts_itr_destroy(iter[i]);
		}
	}

	for( i = 0 ; i < n_platform ; i++){
		bam_destroy1(b[i]);
		bam_hdr_destroy(header[i]);
		sam_close(fp[i]);
		hts_idx_destroy(idx[i]);
	}
	free(iter);
	free(sv);
	free(fp);
	free(header);
	free(idx);
	free(b);
	fclose(in);
	fclose(out);
	return 0 ;
}

