#include<iostream.h>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<ctype.h>
class DNA
{
 int size;
 int* DNA_seq,*COMP_seq;
 public:
	 DNA(char *filename);
	 ~DNA();
	 int getsize() const{return size;}
	 int* getseq();
	 int* getcompseq(){return COMP_seq;}
	 void printseq();
	 void copyseq(char *filename);
};

void DNA::printseq()
{
	//cout<<"\n inside print seq";
	int i;
	for(i=0;i<size;i++)
		switch(DNA_seq[i])
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
	for(int i=0;i<size;i++)
		switch(DNA_seq[i])
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
 char* protein_seq[6];
 int size[6];
 public:
  protein(DNA);
 void printseq(int);
 void copyseq(char *filename);
};
void protein::printseq(int x)
{
	for(int i=0;i<size[x];i++)
		cout<<protein_seq[i];
}

void protein::copyseq(char *fname)
{
	FILE* fp;
	fp=fopen(fname,"w+");
	for(int j=0;j<6;j++)
	{
	for(int i=0;i<size[j];i++)
		fputc(protein_seq[j][i],fp);
	fputs("\n\tPROTEIN SEQUENCE\t\n",fp);
	
	}
	fclose(fp);
}
 protein::protein(DNA DNAclass)
 {
	 int i,j;
	 int *sequence=DNAclass.getseq();
	 int tempsize=DNAclass.getsize(); 
	 for(j=0;j<3;j++)
	 {
	  protein_seq[j]=(char*)malloc(sizeof(char)*(tempsize/3));
	  size[j]=0;
	  for(i=j;(i+3)<tempsize;i+=3)
  	  {
		  int number=(sequence[i]*100)+(sequence[i+1]*10)+sequence[i+2];
		  switch(number)
		  {
			  case 111:
			  case 113:
				  protein_seq[j][size[j]++]='F';
				  break;
			  case 110:
			  case 112:
			  case 311:
			  case 313:
			  case 310:
			  case 312:
			          protein_seq[j][size[j]++]='L';
				  break;
			  case 131:
			  case 133:
			  case 130:
			  case 132:
			  case 021:		  
			  case 023:		  
			          protein_seq[j][size[j]++]='S';
				  break;
			  case 101:
			  case 103:
				  protein_seq[j][size[j]++]='Y';
				  break;
			  case 100:
			  case 102:
			  case 120:
			          protein_seq[j][size[j]++]='X';
			          break;
			  case 121:
			  case 123:
			          protein_seq[j][size[j]++]='C';
			          break;
			  case 122:
			           protein_seq[j][size[j]++]='W';
			           break;
			  case 331:
			  case 333:
			  case 330:
			  case 332:
				   protein_seq[j][size[j]++]='P';
				   break;
			  case 301:
			  case 303:
				   protein_seq[j][size[j]++]='H';
				   break;
			  case 300:
			  case 302:
				   protein_seq[j][size[j]++]='Q';
				   break;
			  case 321:
			  case 323:
			  case 320:
			  case 322:
			  case 20:		  
			  case 22:
				   protein_seq[j][size[j]++]='R';
				   break;
		          		   
			  case 11:		  
			  case 13:		  
			  case 10:
		 		   protein_seq[j][size[j]++]='I';
				   break;		   
			  case 12:
		 		   protein_seq[j][size[j]++]='M';
				   break;		   
			  case 31:		  
			  case 33:		  
			  case 30:		  
			  case 32:
		 		  protein_seq[j][size[j]++]='T';
       				  break;				  
			  case 1:		  
			  case 3:
				   protein_seq[j][size[j]++]='N';
		 		   break;		   
			  case 0:		  
			  case 2:
		 		   protein_seq[j][size[j]++]='K';
				   break;		   
		          case 211:		   
		          case 213:		   
		          case 210:		   
		          case 212:
				   protein_seq[j][size[j]++]='V';
				   break;		   
		          case 231:		   
		          case 233:		   
		          case 230:		   
		          case 232:
				   protein_seq[j][size[j]++]='A';
				   break;		   
		          case 201:		   
		          case 203:
				   protein_seq[j][size[j]++]='D';
				   break;		   
		          case 200:		   
		          case 202:
				   protein_seq[j][size[j]++]='E';
				   break;			   
		          case 221:		   
		          case 223:		   
		          case 220:		   
		          case 222:
				   protein_seq[j][size[j]++]='G';
				   break;			   
		  }
	  
	  }
	 cout<<"After generating protein sequence for the forward strand in "<<j<<" reading frame\n";
	 }
	 
	 for(j=3;j<6;j++)
	 {
	  protein_seq[j]=(char*)malloc(sizeof(char)*(DNAclass.getsize()/3));
	  size[j]=0;
	  for(int i=DNAclass.getsize()-j+2;i-3>0;i-=3)
  	  {
		  //cout<<sequence[i]<<sequence[i-1]<<sequence[i-2]<<"	"<<DNAclass.getsize()<<"\n";
		  int number=(sequence[i]*100)+(sequence[i-1]*10)+sequence[i-2];
		  //cout<<"number is "<<number<<" and reading frame is "<<j<<endl;
		  switch(number)
		  {
			  case 0:
			  case 200:
				  protein_seq[j][size[j]++]='F';
				  break;
			  case 100:
			  case 300:
			  case 2:
			  case 202:
			  case 102:
			  case 302:
			          protein_seq[j][size[j]++]='L';
				  break;
			  case 20:
			  case 220:
			  case 120:
			  case 320:
			  case 31:		  
			  case 231:		  
			          protein_seq[j][size[j]++]='S';
				  break;
			  case 10:
			  case 210:
				  protein_seq[j][size[j]++]='Y';
				  break;
			  case 110:
			  case 310:
			  case 130:
			          protein_seq[j][size[j]++]='X';
			          break;
			  case 30:
			  case 230:
			          protein_seq[j][size[j]++]='C';
			          break;
			  case 330:
			           protein_seq[j][size[j]++]='W';
			           break;
			  case 22:
			  case 222:
			  case 122:
			  case 322:
				   protein_seq[j][size[j]++]='P';
				   break;
			  case 12:
			  case 212:
				   protein_seq[j][size[j]++]='H';
			   break;
			  case 112:
			  case 312:
				   protein_seq[j][size[j]++]='Q';
				   break;
			  case 32:
			  case 232:
			  case 132:
			  case 332:
			  case 131:		  
			  case 331:
				   protein_seq[j][size[j]++]='R';
				   break;
		          		   
			  case 001:		  
			  case 201:		  
			  case 101:
		 		   protein_seq[j][size[j]++]='I';
				   break;		   
			  case 301:
		 		   protein_seq[j][size[j]++]='M';
				   break;		   
			  case 21:		  
			  case 221:		  
			  case 121:		  
			  case 321:
		 		  protein_seq[j][size[j]++]='T';
       				  break;				  
			  case 11:		  
			  case 211:
				   protein_seq[j][size[j]++]='N';
		 		   break;		   
			  case 111:		  
			  case 311:
		 		   protein_seq[j][size[j]++]='K';
				   break;		   
		          case 3:		   
		          case 203:		   
		          case 103:		   
		          case 303:
				   protein_seq[j][size[j]++]='V';
				   break;		   
		          case 23:		   
		          case 223:		   
		          case 123:		   
		          case 323:
				   protein_seq[j][size[j]++]='A';
				   break;		   
		          case 13:		   
		          case 213:
				   protein_seq[j][size[j]++]='D';
				   break;		   
		          case 113:		   
		          case 313:
				   protein_seq[j][size[j]++]='E';
				   break;			   
		          case 33:		   
		          case 233:		   
		          case 133:		   
		          case 333:
				   protein_seq[j][size[j]++]='G';
				   break;			   
		  }
	  
	  }
	 cout<<"After generating protein sequence for the reverse strand in reading frame "<<j<<endl;
	 }
 }
 DNA::~DNA()
 {
 }
int* DNA::getseq()
{
	return DNA_seq;
}
 DNA::DNA(char* filename)
 {
 	  FILE* fp;
	  char str[500],ch;
	  int temp=0;
	  if(!(fp=fopen(filename ,"r")))
	  {
			  cout<<"File "<<filename<<"does not exist\n";
			  exit(1);
	  }
	  fseek(fp,0L,SEEK_END);
	  size=ftell(fp);
	  DNA_seq=(int *)(malloc(sizeof(int)*size));
	  COMP_seq=(int *)(malloc(sizeof(int)*size));
	  fseek(fp,0L,SEEK_SET);
	  fgets(str,500,fp);
	 while(!feof(fp))
	  {	
	     ch=fgetc(fp);
	     if(toupper(ch)=='A')
	     	    {
   	                DNA_seq[temp++]=0;
		    }
	     else if(toupper(ch)=='T')
	    	    {
	            	DNA_seq[temp++]=1;
		    }
	     else if(toupper(ch)=='G')
	    	    {
	            	DNA_seq[temp++]=2;
		    }
	     else if(toupper(ch)=='C')
	    	    {
	                DNA_seq[temp++]=3;
		    }
  	  }
          fclose(fp);
	  size=temp;
//	  cout<<"am here\n"<<endl;
		  int j=0;
	  for(int i=size-1;i>=0;i--)
	  {
		  switch(DNA_seq[i])
		  {
			  case 0:
		  		COMP_seq[j++]=1;
				break;
			  case 1:
		  		COMP_seq[j++]=0;
				break;
			  case 2:
		  		COMP_seq[j++]=3;
				break;
			  case 3:
		  		COMP_seq[j++]=2;
				break;
		  }
	  }
  //        cout<<"Size of the genome ="<<size<<"bp\n";
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
