/*
Program for Quick short sequence search in a fasta file

Usage: program_name genome_file.fa reads.fa

Sreenu Vattipally
MRC-CVR, University of Glasgow
Glasgow, UK

16/May/2014
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "miRNA_Search.h"

#define LINE_LEN 1024

int main (int argc, char **argv){

	FILE *input;
	int j=0, k=0, found=0,w=0,x=0,y=0,z=0;
	char nuc, buffer[LINE_LEN+1],rev_seq[LINE_LEN+1],name[LINE_LEN+1],sbstr[40];

	// Define unlimited char array
	CharArray c;
	// Initialize the array with size 51200
	initCharArray(&c, 51200);

	name[0]='\0';

	// Read the first input file and fill in the array c

	if((input = fopen(argv[1], "r"))!=NULL){ 
		while((nuc = fgetc(input))!=EOF) { 
			if (feof(input)) break; 
			if(nuc!='\n') insertCharArray(&c, nuc);
			} fclose(input);
		} // End of reading the first file

		// Read the k-mer(short reads) file. Each kmer (short read) length should be less than 1024bp.

    	if((input = fopen(argv[2], "r"))!=NULL){ 
		while (!feof(input)) { 
			fgets(buffer,LINE_LEN,input); 
			if(feof(input)) break; 
			
			buffer[strlen(buffer)-1]='\0'; 

			if(buffer[0]=='>') strcpy(name,buffer);

			else{ 
				TempLate(buffer,rev_seq); 
				// Searching...
				found=KR2(buffer,strlen(buffer), c.array,c.used); 
				found+=KR2(rev_seq,strlen(rev_seq), c.array,c.used); 
				if(found!=0)  printf("%d\t%s\t%s\n",found,buffer,name);
				else  printf("0\t%s\t%s\n",buffer,name); 
			}
			found=0;
		} fclose(input);
	}
}
