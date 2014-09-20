#include<stdio.h>
#include<string.h>
#include<stdlib.h>

int main(int argc , char *argv[])
{
	FILE *fp =  fopen(argv[1],"r");
	if(!fp) return -1;
	char	name[1024]; 
	int	pos[3];
	char	c ;
	int	sel ;
	int	tmp ;
	int	chr ;
	printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILT\tINFO\n");
	while(!feof(fp)){ 
		sel = -1 ;
		fscanf(fp,"%s\t",name);
		printf("%s\t",name);
		if(!strcmp(name,"CNV1") || !strcmp(name,"CNV2")) sel = 1 ;
		if(!strcmp(name,"CNV3") || !strcmp(name,"CNV4")) sel = 2 ;
		fscanf(fp,"%s\t",name);
		chr = atoi(name);
		fscanf(fp,"%s\t",name);
		pos[0] = atoi(name);
		fscanf(fp,"%s\t",name);
		pos[1] = atoi(name);

		fscanf(fp,"%c\t",&c);
		fscanf(fp,"%s\t",name);
		pos[2] = atoi(name);
		if(set == -1) continue ;
		
		if(sel == 1){
			tmp = pos[1] ;
			pos[1] = pos [2];
			pos[2] = tmp ;
		}	
		printf("%d\t",chr);
		printf("%d\t",pos[0]);
		printf("%d\t",pos[1]);
		printf("%d\n",pos[2]);
	}
	fclose(fp);
	return 0;
}
