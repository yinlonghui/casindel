#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>

typedef struct{
	int      type ;
	int	pos[8] ;
	char	mut[8] ;
} vcf_info ;

int main(int argc , char *argv[])
{
	
	FILE *fin , *fout ;

	vcf_info    vcf ;

	faidx_t   *fai ;
	
	char	buffer[32];
	
	char *s = malloc(1024*sizeof(char));

	if( argc != 4 ){
		fprintf(stderr,"sv2vcf  [ref.fa]  [input] [output] \n");
		return 1 ;
	}

	fai =  fai_load(argv[1]);
	
	if(!fai){
		fprintf(stderr,"can't open the file:%s.fai\n",argv[1]);
		return  1;
	}


	fin =  fopen(argv[2],"r");

	if(!fin){
		fprintf(stderr,"can't open the file:%s\n", argv[2]);
		return 1 ;
	}

	fout = fopen(argv[3],"w");
	if(!fout){
		fprintf(stderr,"can't open the file:%s\n", argv[3]);
		return 1 ;
	}
	fprintf(fout,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILT\tINFO\n");
	while(!feof(fin)){
		char *tmp ;
		int  len ;
		fscanf(fin,"%s\t",buffer);
		if(!strcmp(buffer,"TRSdL")){
			fscanf(fin,"%s\t",buffer);
			fscanf(fin,"%d\t%d\n",&vcf.pos[0],&vcf.pos[1]);
			int	k  ;
			for(  k = 0 ;  k < 2 ; k++){
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[k],vcf.pos[k]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[k] = *tmp ;
				free(tmp);
			}
			fprintf(fout,"%s\t%d\tbnd_V\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],buffer,vcf.pos[1],vcf.mut[0]);
			fprintf(fout,"%s\t%d\tbnd_W\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],buffer,vcf.pos[0],vcf.mut[1]);

		}else if(!strcmp(buffer,"TRSdR")){
			fscanf(fin,"%s\t",buffer);
			fscanf(fin,"%d\t%d\n",&vcf.pos[0],&vcf.pos[1]);
			int	k  ;
			for(  k = 0 ;  k < 2 ; k++){
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[k],vcf.pos[k]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[k] = *tmp ;
				free(tmp);
			}
			fprintf(fout,"%s\t%d\tbnd_V\t%c\t%s:%d]%c]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],buffer,vcf.pos[1],vcf.mut[0]);
			fprintf(fout,"%s\t%d\tbnd_W\t%c\t%s:%d]%c]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],buffer,vcf.pos[0],vcf.mut[1]);

		}else if((tmp = strstr(buffer,"TandemCNV"))){
			tmp = tmp + 9 ;
			vcf.type = atoi(tmp);
			fscanf(fin,"%s\t",buffer);
			if(vcf.type == 1 ) {
				fscanf(fin,"%*[^\n]\n");
				continue ;
			}
			fscanf(fin,"%d\t%d\n",&vcf.pos[0],&vcf.pos[1]);
			int	k  ;
			for(  k = 0 ;  k < 2 ; k++){
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[k],vcf.pos[k]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[k] = *tmp ;
				free(tmp);
			}
			if(vcf.type == 2){
				fprintf(fout,"%s\t%d\tbnd_V\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],vcf.mut[1],buffer,vcf.pos[1]);
				fprintf(fout,"%s\t%d\tbnd_W\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],vcf.mut[1],buffer,vcf.pos[1]);
			}else {
				fprintf(fout,"%s\t%d\tbnd_V\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],buffer,vcf.pos[0],vcf.mut[0]);
				fprintf(fout,"%s\t%d\tbnd_W\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],buffer,vcf.pos[0],vcf.mut[0]);
			}


		}else if((tmp = strstr( buffer,"CNV"))){
			
			tmp = tmp + 3;
			vcf.type = atoi(tmp);
			fscanf(fin,"%s\t",buffer);
			fscanf(fin,"%d\t%d\t%d\n",&vcf.pos[0],&vcf.pos[1],&vcf.pos[2]);
			if( vcf.type == 1 ){

				sprintf(s,"%s:%d-%d",buffer,vcf.pos[0],vcf.pos[0]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[0] = *tmp ;
				free(tmp);
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[2],vcf.pos[2]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[2] = *tmp ;
				free(tmp);

				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[1],vcf.pos[1]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[1] = *tmp ;
				free(tmp);
				
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[2]+1,vcf.pos[2]+1);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[3] = *tmp ;
				free(tmp);

				fprintf(fout,"%s\t%d\tbnd_V\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],buffer,vcf.pos[2],vcf.mut[0]);
				fprintf(fout,"%s\t%d\tbnd_W\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V;EVENT=PRO\n",buffer,vcf.pos[2],vcf.mut[2],vcf.mut[2],buffer,vcf.pos[0]);
				fprintf(fout,"%s\t%d\tbnd_U\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_X;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],vcf.mut[1],buffer,vcf.pos[2]+1);
				fprintf(fout,"%s\t%d\tbnd_X\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_U;EVENT=PRO\n",buffer,vcf.pos[2]+1,vcf.mut[3],buffer,vcf.pos[1],vcf.mut[3]);
			}else if( vcf.type == 2 ){
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[0],vcf.pos[0]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[0] = *tmp ;
				free(tmp);

				sprintf(s,"%s:%d-%d",buffer,vcf.pos[2]+1,vcf.pos[2]+1);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[3] = *tmp ;
				free(tmp);
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[1],vcf.pos[1]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[1] = *tmp ;
				free(tmp);
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[2],vcf.pos[2]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[2] = *tmp ;
				free(tmp);

				fprintf(fout,"%s\t%d\tbnd_V\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],buffer,vcf.pos[2]+1,vcf.mut[0]);
				fprintf(fout,"%s\t%d\tbnd_W\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_U;EVENT=PRO\n",buffer,vcf.pos[2]+1,vcf.mut[3],buffer,vcf.pos[0],vcf.mut[3]);
				fprintf(fout,"%s\t%d\tbnd_U\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_X;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],vcf.mut[1],buffer,vcf.pos[2]);
				fprintf(fout,"%s\t%d\tbnd_X\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_U;EVENT=PRO\n",buffer,vcf.pos[2],vcf.mut[2],vcf.mut[2],buffer,vcf.pos[1]);


			}else if(vcf.type == 3){
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[0],vcf.pos[0]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[0] = *tmp ;
				free(tmp);
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[1],vcf.pos[1]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[1] = *tmp ;
				free(tmp);
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[0]+1,vcf.pos[0]+1);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[2] = *tmp ;
				free(tmp);
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[2],vcf.pos[2]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[3] = *tmp ;
				free(tmp);

				fprintf(fout,"%s\t%d\tbnd_V\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],vcf.mut[0],buffer,vcf.pos[1]);
				fprintf(fout,"%s\t%d\tbnd_W\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],buffer,vcf.pos[0],vcf.mut[1]);
				fprintf(fout,"%s\t%d\tbnd_U\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_X;EVENT=PRO\n",buffer,vcf.pos[0]+1,vcf.mut[2],buffer,vcf.pos[2],vcf.mut[2]);
				fprintf(fout,"%s\t%d\tbnd_X\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[2],vcf.mut[3],vcf.mut[3],buffer,vcf.pos[0]+1);

			}else if(vcf.type == 4){
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[0],vcf.pos[0]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[0] = *tmp ;
				free(tmp);
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[1],vcf.pos[1]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[1] = *tmp ;
				free(tmp);
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[0]+1,vcf.pos[0]+1);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[2] = *tmp ;
				free(tmp);
				
				sprintf(s,"%s:%d-%d",buffer,vcf.pos[2],vcf.pos[2]);
				tmp = fai_fetch(fai,s,&len);
				vcf.mut[3] = *tmp ;
				free(tmp);

				fprintf(fout,"%s\t%d\tbnd_V\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],vcf.mut[0],buffer,vcf.pos[2]);
				fprintf(fout,"%s\t%d\tbnd_W\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V;EVENT=PRO\n",buffer,vcf.pos[2],vcf.mut[3],vcf.mut[3],buffer,vcf.pos[0]);
				fprintf(fout,"%s\t%d\tbnd_U\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_X;EVENT=PRO\n",buffer,vcf.pos[0]+1,vcf.mut[2],buffer,vcf.pos[1],vcf.mut[2]);
				fprintf(fout,"%s\t%d\tbnd_X\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_U;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[3],buffer,vcf.pos[0]+1,vcf.mut[3]);
			}
		}else if((tmp =strstr(buffer,"TRS"))){
			
			tmp = tmp + 3;
			vcf.type = atoi(tmp);
			fscanf(fin,"%s\t",buffer);
			if( vcf.type == 2 || vcf.type == 3 ){
				fscanf(fin,"%d\t%d\t%d\t%d\n",&vcf.pos[1],&vcf.pos[2],&vcf.pos[5],&vcf.pos[6]);

				vcf.pos[0] = vcf.pos[1] - 1 ;
				vcf.pos[3] = vcf.pos[2] + 1 ;
				vcf.pos[4] = vcf.pos[5] - 1 ;
				vcf.pos[7] = vcf.pos[6] + 1 ;

				int  k  = 0 ;
				for( k = 0 ;  k < 8 ; k++){
					sprintf(s,"%s:%d-%d",buffer,vcf.pos[k],vcf.pos[k]);
					tmp = fai_fetch(fai,s,&len);
					vcf.mut[k] = *tmp ;
					free(tmp);
				}
				if(vcf.type == 2 ){
					fprintf(fout,"%s\t%d\tbnd_V\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],vcf.mut[0],buffer,vcf.pos[6]);
					fprintf(fout,"%s\t%d\tbnd_W\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V;EVENT=PRO\n",buffer,vcf.pos[6],vcf.mut[6],vcf.mut[6],buffer,vcf.pos[0]);
					fprintf(fout,"%s\t%d\tbnd_U\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_X;EVENT=PRO\n",buffer,vcf.pos[3],vcf.mut[3],buffer,vcf.pos[5],vcf.mut[3]);
					fprintf(fout,"%s\t%d\tbnd_X\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_U;EVENT=PRO\n",buffer,vcf.pos[5],vcf.mut[5],buffer,vcf.pos[3],vcf.mut[5]);
					
					fprintf(fout,"%s\t%d\tbnd_E\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_D;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],buffer,vcf.pos[4],vcf.mut[1]);
					fprintf(fout,"%s\t%d\tbnd_D\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_E;EVENT=PRO\n",buffer,vcf.pos[4],vcf.mut[4],vcf.mut[4],buffer,vcf.pos[1]);
					fprintf(fout,"%s\t%d\tbnd_Z\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_H;EVENT=PRO\n",buffer,vcf.pos[2],vcf.mut[2],vcf.mut[2],buffer,vcf.pos[7]);
					fprintf(fout,"%s\t%d\tbnd_H\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_Z;EVENT=PRO\n",buffer,vcf.pos[7],vcf.mut[7],buffer,vcf.pos[2],vcf.mut[7]);
					

				}else {
					fprintf(fout,"%s\t%d\tbnd_V\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],vcf.mut[0],buffer,vcf.pos[5]);
					fprintf(fout,"%s\t%d\tbnd_W\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V;EVENT=PRO\n",buffer,vcf.pos[5],vcf.mut[5],buffer,vcf.pos[0],vcf.mut[5]);
					fprintf(fout,"%s\t%d\tbnd_U\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_X;EVENT=PRO\n",buffer,vcf.pos[3],vcf.mut[3],buffer,vcf.pos[6],vcf.mut[3]);
					fprintf(fout,"%s\t%d\tbnd_X\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_U;EVENT=PRO\n",buffer,vcf.pos[6],vcf.mut[6],vcf.mut[6],buffer,vcf.pos[3]);
					
					fprintf(fout,"%s\t%d\tbnd_E\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_D;EVENT=PRO\n",buffer,vcf.pos[2],vcf.mut[2],vcf.mut[2],buffer,vcf.pos[4]);
					fprintf(fout,"%s\t%d\tbnd_D\t%c\t%c]%s:%d]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_E;EVENT=PRO\n",buffer,vcf.pos[4],vcf.mut[4],vcf.mut[4],buffer,vcf.pos[2]);
					fprintf(fout,"%s\t%d\tbnd_Z\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_H;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],buffer,vcf.pos[7],vcf.mut[1]);
					fprintf(fout,"%s\t%d\tbnd_H\t%c\t[%s:%d[%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_Z;EVENT=PRO\n",buffer,vcf.pos[7],vcf.mut[7],buffer,vcf.pos[1],vcf.mut[7]);

				}
			}else if(vcf.type == 1){
				fscanf(fin,"%d\t%d\n",&vcf.pos[1],&vcf.pos[2]);
				vcf.pos[0] = vcf.pos[1] - 1 ;
				vcf.pos[3] = vcf.pos[2] + 1 ;
				int	k ;
				for( k = 0 ;  k < 4 ; k++){
					sprintf(s,"%s:%d-%d",buffer,vcf.pos[k],vcf.pos[k]);
					tmp = fai_fetch(fai,s,&len);
					vcf.mut[k] = *tmp ;
					free(tmp);
				}
				fprintf(fout,"%s\t%d\tbnd_V\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[0],vcf.mut[0],vcf.mut[0],buffer,vcf.pos[2]);
				fprintf(fout,"%s\t%d\tbnd_W\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[2],vcf.mut[2],buffer,vcf.pos[0],vcf.mut[2]);
				fprintf(fout,"%s\t%d\tbnd_U\t%c\t]%s:%d]%c\t6\tPASS\tSVTYPE=BND;MATEID=bnd_X;EVENT=PRO\n",buffer,vcf.pos[1],vcf.mut[1],buffer,vcf.pos[3],vcf.mut[1]);
				fprintf(fout,"%s\t%d\tbnd_X\t%c\t%c[%s:%d[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_W;EVENT=PRO\n",buffer,vcf.pos[3],vcf.mut[3],vcf.mut[3],buffer,vcf.pos[1]);


			}else{
				fscanf(fin,"%*[^\n]\n");
			}
		}else{
				fscanf(fin,"%*[^\n]\n");
		}
	}
	
	free(s);
	fclose(fin);
	fclose(fout);
	return 0;
}
