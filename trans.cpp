#include<iostream.h>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<ctype.h>
#define BASE 10000
class DNA
{
 int size[BASE];
 int noDNAseq;
 int* DNA_seq[BASE],*COMP_seq[BASE];
 char firststr[BASE][100];
	public:
	 DNA(char *filename);
	 ~DNA();
	 int getsize(int k) const{return size[k];}
	 int* getseq(int k);
	 int getnodnaseq(){return noDNAseq;}
	 int* getcompseq(int k){return COMP_seq[k];}
	 char* getfirststr(int k){return firststr[k];}
	 void printseq();
	 void copyseq(char *filename);
};
//prints sequence on screen
void DNA::printseq()
{
	int i,j;
	for(int j=0;j<noDNAseq;j++)
	for(i=0;i<size[j];i++)
		switch(DNA_seq[i][j])
		{
			case 0:
				cout<<"A";
				break;
			case 1:
				cout<<"T";
				break;
			case 2:
				cout<<"G";
				break;
			case 3:
				cout<<"C";
				break;
		}
	//printf("i=%d\n",i);
}
//copies translated protein sequence into file fname
void DNA::copyseq(char *fname)
{
	FILE* fp;
	fp=fopen(fname,"w+");
	for(int j=0;j<noDNAseq;j++)
	for(int i=0;i<size[j];i++)
		switch(DNA_seq[j][i])
		{
			case '0':
				fputc('A',fp);
				break;
			case '1':
				fputc('T',fp);
				break;
			case '2':
				fputc('G',fp);
				break;
			case '3':
				fputc('C',fp);
				break;
		}
	fclose(fp);
}
class protein
{
 char* protein_seq[BASE][6];
 int size[BASE][6];
 int noprotseq;
 char firststr[BASE][100];
 char table[5][5][5];
 public:
  protein(char*,DNA);
 void printseq();
 void copyseq(char *filename);
};
void protein::printseq()
{
	for(int a=0;a<noprotseq;a++)
	{
		for(int l=0;l<6;l++)
		for(int i=0;i<size[a][l];i++)
			cout<<protein_seq[a][l][i];
	}
}
//copies sequence into file , whose name is received as argument
void protein::copyseq(char *fname)
{
	FILE* fp;
	fp=fopen(fname,"w+");
	for(int y=0;y<noprotseq;y++)
	{
		
		fputc('\n',fp);
		fputc('>',fp);
		fprintf(fp,"%s",firststr[y]);
		for(int j=0;j<6;j++)
		{
			fprintf(fp,"\n\tPROTEIN SEQUENCE READING FRAME %d\t\n",j+1);
			for(int i=0;i<size[y][j];i++)
				fputc(protein_seq[y][j][i],fp);
	
		}
	}
	fclose(fp);
}
//maps charcter n returns as 0/1/2/3 
int map(char x)
{
	if(x=='A')
		return 0;
	else if(x=='T')
		return 1;
	else if(x=='G')
		return 2;
	else if(x=='U')
		return 1;
	else if(x=='C')
		return 3;
}
//constructor of protein class, which reads genetic code and  translates DNA sequence in all 3 frames
 protein::protein(char* fname,DNA DNAclass)
 {
	 int i,j;
	 noprotseq=DNAclass.getnodnaseq();
	 FILE* fp;
	 fp=fopen(fname,"r");
	 char str1[100],codon[4],aacid;
	 while(fp)
	 {
		 if(feof(fp))
			 break;
		 fgets(str1,100,fp);
		 if(str1[11]!='i')
		 	sscanf(str1,"%s %c %*s",codon,&aacid);
		 else
		 	 {
		 	 //printf("in else\n");
			 sscanf(str1,"%s %*c %*s %c",codon,&aacid);
			 }
		 //printf("before allocating to table %d %d %d %c %s\n",map(codon[0]),map(codon[1]),map(codon[2]),aacid,codon);
		 table[map(codon[0])][map(codon[1])][map(codon[2])]=aacid;	 
	 }
	 fclose(fp); 
	 //for(int a=0;a<4;a++)
	 //for(int b=0;b<4;b++)
	 //for(int c=0;c<4;c++)
	//	 printf("%c ",table[a][b][c]);
	 for(int k=0;k<noprotseq;k++)
	 {
		strcpy(firststr[k],DNAclass.getfirststr(k));
	 	for(int l=0;l<6;l++)
	 	{
		         protein_seq[k][l]=(char *)malloc(sizeof(char)*DNAclass.getsize(k));
			 size[k][l]=DNAclass.getsize(k)/3;
	 	}
	 int *sequence=DNAclass.getseq(k);
	 int tempsize=DNAclass.getsize(k); 
	 for(j=0;j<3;j++)
	 {
	  size[k][j]=0;
	  for(i=j;(i+3)<tempsize;i+=3)
  	  {
	  	protein_seq[k][j][size[k][j]++]=table[sequence[i]][sequence[i+1]][sequence[i+2]];
	  }
	 //cout<<"After generating protein sequence for the forward strand in "<<j<<" reading frame\n";
	 }
	 //now read s complementary sequence and tranlates that in all three frames
	 sequence=DNAclass.getcompseq(k);
	 for(j=3;j<6;j++)
	 {
	  size[k][j]=0;
	  for(int i=DNAclass.getsize(k)-j+2;i-3>0;i-=3)
  	  {
		   protein_seq[k][j][size[k][j]++]=table[sequence[i]][sequence[i-1]][sequence[i-2]];
	  }
	 //cout<<"After generating protein sequence for the reverse strand in reading frame "<<j<<endl;
	 }
 }
 }
//destructor
 DNA::~DNA()
 {
 }
//accessor function
int* DNA::getseq(int k)
{
	return DNA_seq[k];
}
//constructor of class DNA
//takea filename as argument 
DNA::DNA(char* filename)
 {
 	  FILE* fp;
	  char ch;
	  int temp=0;
	  noDNAseq=-1;
	  if(!(fp=fopen(filename ,"r")))
	  {
			  cout<<"File "<<filename<<"does not exist\n";
			  exit(1);
	  }
	  //reads every character from the file and stores it as 0,1, 2 or 3
	 while(!feof(fp))
	  {	
	     ch=fgetc(fp);
	     if(feof(fp))
		     break;
	     if(ch=='>')
	     {
		if(noDNAseq!=-1)
		size[noDNAseq]=temp;
			noDNAseq++;
		temp=0;	
		if(noDNAseq>9999)
		{
			printf("Can process only 9999 sequences at a given time. Exiting......\n");
			exit(1);
		}
   		fgets(firststr[noDNAseq],100,fp);
		DNA_seq[noDNAseq]=(int *)malloc(sizeof(int)*2);
	     }
	     else if(isalpha(ch))
	     {
		     DNA_seq[noDNAseq]=(int *)realloc(DNA_seq[noDNAseq],sizeof(int)*(temp+3));
		     if(toupper(ch)=='A')
	     	    {
   	                DNA_seq[noDNAseq][temp++]=0;
		    }
	     	    else if(toupper(ch)=='T')
	    	    {
	            	DNA_seq[noDNAseq][temp++]=1;
		    }
	     	    else if(toupper(ch)=='G')
	    	    {
	            	DNA_seq[noDNAseq][temp++]=2;
		    }
	     	    else if(toupper(ch)=='U')
	    	    {
	            	DNA_seq[noDNAseq][temp++]=1;
		    }
		     else if(toupper(ch)=='C')
	    	    {
	                DNA_seq[noDNAseq][temp++]=3;
		    }
	    }
  	  }
		size[noDNAseq]=temp;
          fclose(fp);
	  //iof no of sequences in file is less that 10000 then free the rest
	  //else if number if sequences greater than 09999, then print message 
	  if(noDNAseq<BASE)
		  for(int l=noDNAseq+1;l<BASE;l++)
		  {
			  free((void *)DNA_seq[l]); 
			  free((void *)size[l]);
		  }
	  //get complementary sequence
	  for(int j=0;j<noDNAseq;j++)
	  {
       	  COMP_seq[j]=(int *)malloc(sizeof(int)*size[j]);				 
	  for(int i=size[j]-1,k=0;i>=0;i--)
	  {
		  switch(DNA_seq[j][i])
		  {
			  case 0:
		  		COMP_seq[j][k++]=1;
				break;
			  case 1:
		  		COMP_seq[j][k++]=0;
				break;
			  case 2:
		  		COMP_seq[j][k++]=3;
				break;
			  case 3:
		  		COMP_seq[j][k++]=2;
				break;
		  }
	  }
	 }
	  //noDNAseq++;
	  printf("%d \n",noDNAseq);
 }

main(int argc, char* argv[])
{
       	if(argc!=3)
	 {
	  printf("enter the input in the form <filename> <file containing genetic codes>\n"); 
	  exit(1);
 	 }
//calls constructor of class DNA with filename as argument
          DNA DNAsequence(argv[1]);
	 // printf("\n after reeading files in multi fasta format\n");
	  protein PROTsequence(argv[2],DNAsequence);
	 PROTsequence.copyseq("output.txt");
}
