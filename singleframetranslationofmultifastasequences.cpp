#include<iostream.h>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<ctype.h>
#define BASE 10000
// Single frame translation from DNA to Protein
class DNA
{
 int size[BASE];
 int noDNAseq;
 int* DNA_seq[BASE],*COMP_seq[BASE];
 char firstline[BASE][100];
 public:
	 DNA(char *filename);
	 ~DNA();
	 int getsize(int k) const{return size[k];}
	 int* getseq(int k);
	 char* getstr(int k){return firstline[k];}
	 int getnodnaseq(){return noDNAseq;}
	 int* getcompseq(int k){return COMP_seq[k];}
	 void printseq();
	 void copyseq(char *filename);
};

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
	printf("i=%d\n",i);
}

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
 char firstr[BASE][100];
 public:
  protein(DNA);
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

void protein::copyseq(char *fname)
{
	FILE* fp;
	fp=fopen(fname,"w+");
	for(int y=0;y<noprotseq;y++)
	{
		fputc('\n',fp);
		fputc('>',fp);
		fprintf(fp,"%s",firstr[y]);
		for(int j=0;j<1;j++)
		{
			for(int i=0;i<size[y][j];i++)
				fputc(protein_seq[y][j][i],fp);
	
	}
	}
	fclose(fp);
}
 protein::protein(DNA DNAclass)
 {
	 int i,j;
	 noprotseq=DNAclass.getnodnaseq();
	 for(int k=0;k<noprotseq;k++)
	 {
		strcpy(firstr[k],DNAclass.getstr(k));
	 	for(int l=0;l<6;l++)
	 	{
		         protein_seq[k][l]=(char *)malloc(sizeof(char)*DNAclass.getsize(k));
			 size[k][l]=DNAclass.getsize(k)/3;
	 	}
	 int *sequence=DNAclass.getseq(k);
	 int tempsize=DNAclass.getsize(k); 
	 for(j=0;j<1;j++)
	 {
	  size[k][j]=0;
	  for(i=j;(i+3)<tempsize;i+=3)
  	  {
		  int number=(sequence[i]*100)+(sequence[i+1]*10)+sequence[i+2];
		  switch(number)
		  {
			  case 111:
			  case 113:
				  protein_seq[k][j][size[k][j]++]='F';
				  break;
			  case 110:
			  case 112:
			  case 311:
			  case 313:
			  case 310:
			  case 312:
			          protein_seq[k][j][size[k][j]++]='L';
				  break;
			  case 131:
			  case 133:
			  case 130:
			  case 132:
			  case 021:		  
			  case 023:		  
			          protein_seq[k][j][size[k][j]++]='S';
				  break;
			  case 101:
			  case 103:
				  protein_seq[k][j][size[k][j]++]='Y';
				  break;
			  case 100:
			  case 102:
			  case 120:
			          protein_seq[k][j][size[k][j]++]='X';
			          break;
			  case 121:
			  case 123:
			          protein_seq[k][j][size[k][j]++]='C';
			          break;
			  case 122:
			           protein_seq[k][j][size[k][j]++]='W';
			           break;
			  case 331:
			  case 333:
			  case 330:
			  case 332:
				   protein_seq[k][j][size[k][j]++]='P';
				   break;
			  case 301:
			  case 303:
				   protein_seq[k][j][size[k][j]++]='H';
				   break;
			  case 300:
			  case 302:
				   protein_seq[k][j][size[k][j]++]='Q';
				   break;
			  case 321:
			  case 323:
			  case 320:
			  case 322:
			  case 20:		  
			  case 22:
				   protein_seq[k][j][size[k][j]++]='R';
				   break;
		          		   
			  case 11:		  
			  case 13:		  
			  case 10:
		 		   protein_seq[k][j][size[k][j]++]='I';
				   break;		   
			  case 12:
		 		   protein_seq[k][j][size[k][j]++]='M';
				   break;		   
			  case 31:		  
			  case 33:		  
			  case 30:		  
			  case 32:
		 		  protein_seq[k][j][size[k][j]++]='T';
       				  break;				  
			  case 1:		  
			  case 3:
				   protein_seq[k][j][size[k][j]++]='N';
		 		   break;		   
			  case 0:		  
			  case 2:
		 		   protein_seq[k][j][size[k][j]++]='K';
				   break;		   
		          case 211:		   
		          case 213:		   
		          case 210:		   
		          case 212:
				   protein_seq[k][j][size[k][j]++]='V';
				   break;		   
		          case 231:		   
		          case 233:		   
		          case 230:		   
		          case 232:
				   protein_seq[k][j][size[k][j]++]='A';
				   break;		   
		          case 201:		   
		          case 203:
				   protein_seq[k][j][size[k][j]++]='D';
				   break;		   
		          case 200:		   
		          case 202:
				   protein_seq[k][j][size[k][j]++]='E';
				   break;			   
		          case 221:		   
		          case 223:		   
		          case 220:		   
		          case 222:
				   protein_seq[k][j][size[k][j]++]='G';
				   break;			   
		  }
	  
	  }
	 }
	 }
 }
 DNA::~DNA()
 {
 }
int* DNA::getseq(int k)
{
	return DNA_seq[k];
}
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
	 while(!feof(fp))
	  {	
	     ch=fgetc(fp);
	     if(feof(fp))
		     break;
	     if(ch=='>')
	     {
		noDNAseq++;
		if(noDNAseq>9999)
		{
			printf("Can process only 9999 sequences at a given time. Exiting......\n");
			exit(1);
		}
   		fgets(firstline[noDNAseq],500,fp);
		DNA_seq[noDNAseq]=(int *)malloc(sizeof(int)*2);
		if(noDNAseq!=0)
		size[noDNAseq-1]=temp;
		temp=0;	
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
          fclose(fp);
 }

main(int argc, char* argv[])
{
       	if(argc!=2)
	 {
	  printf("enter the input in the form <filename>\n"); 
	  exit(1);
 	 }
          DNA DNAsequence(argv[1]);
	  protein PROTsequence(DNAsequence);
	 PROTsequence.copyseq("output.txt");
}
