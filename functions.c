#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MIN(a,b) (((a)<(b)) ? (a) : (b) ) // the usual "min" function
#define MAX(a,b) (((a)<(b)) ? (b) : (a) ) // the usual "max" function

/*
 * This function returns the smallest common multiple
 * of two integer numbers a and b
 * */
/*
int smallest_common_multiple( int a, int b ) 
{  
	int l = MIN(a,b);
  int h = MAX(a,b);
  int i;
	for(i=h; ; i+= h ) {
		if( i%l == 0 )
			return i;
	}
}
*/
/* Written by S. Brieden, 26.07.19
 * This function pastes n_files files with names
 * filename_wo_extension__id.txt (where id is an integer
 * between [1, n_files]) into a new file with name
 * filename_wo_extension.txt.
 * After copying, the old files are removed.
 *  */
void concatenate_files(int n_files, char *filename_wo_extension)
{
int ifile, c;
int status;
FILE *fnew, *fold;
char filename_new[2000];
char filename_old[2000];

// Write new file
sprintf(filename_new,"%s.txt",filename_wo_extension);
fnew = fopen(filename_new, "w"); 
if (fnew == NULL) 
{ 
printf("Cannot open file %s \n", filename_new); 
exit(0); 
} 
fclose(fnew);

printf("Writing old files to new file\n");

// Loop over old files
fnew = fopen(filename_new, "a"); 
for (ifile=0;ifile<n_files;ifile++)
{
sprintf(filename_old,"%s__%d.txt",filename_wo_extension,ifile+1);
fold = fopen(filename_old, "r");
if (fold == NULL) 
{ 
printf("Cannot open file %s \n", filename_old); 
exit(0); 
}
// Read contents from file and append to new file
c = fgetc(fold); 
while (c != EOF) 
{ 
fputc(c, fnew); 
c = fgetc(fold); 
} 
fclose(fold); 
//remove the old file from the system

status = remove(filename_old);
if (status != 0){
printf("Warning: File %s could not be deleted \n", filename_old);
}
}

printf("Contents copied to %s\n", filename_new); 

fclose(fnew);

}

/* Written by S. Brieden, 08.03.21
 * * For an already opened file* myfile, this function parses
 * * through the header until it is over and sets the file pointer
 * * starting at the line where the data begins. The structure
 * * is very similar to the function countlines()
 * */

void skipheader(FILE* myfile)
{
int ch_ascii = -1;
int line_start = -1;
int doc_start = -1;

doc_start = fgetc(myfile);//initialise

while (doc_start == ' '){//Skip blank spaces
doc_start = fgetc(myfile);
}
if (doc_start>39){// End function, if the first non-space character is not header-like.

fseek(myfile, -1, SEEK_CUR);// set the filestream one back, so that it is one position before the first relevant character

return;// finish function
}
else {// Else, loop over the following lines until reaching end of header
do
{

ch_ascii = fgetc(myfile);// ch_ascii is the ascii number of the next character of the file


if(ch_ascii == '\n'){// When line end is reached, check how next line starts
line_start = fgetc(myfile);// skip blank spaces at the beginning of the line

while (line_start == ' '){
line_start = fgetc(myfile);
}
}

} while (ch_ascii != EOF && line_start<39);// continue loop while the first character in line is header-like. Loop ends, either when the header has finished or the end of the file is reached.

fseek(myfile, -1, SEEK_CUR);// set the filestream one back, so that it is one position before the first relevant character
}
}

/*
void skipheader(FILE* myfile)//old version
{
int ch_ascii = -1;
int line_start = -1;
int doc_start = -1;

// initialise
doc_start = fgetc(myfile);

//Skip blank spaces
while (doc_start == ' '){
doc_start = fgetc(myfile);
}
// End function, if the first non-space character is not header-like.
if (doc_start>39){
return;
}
// Else, loop over the following lines until reaching end of header
else {
do
{
// ch_ascii is the ascii number of the next character of the file
ch_ascii = fgetc(myfile);

// When line end is reached, check how next line starts
if(ch_ascii == '\n'){
line_start = fgetc(myfile);
// skip blank spaces at the beginning of the line
while (line_start == ' '){
line_start = fgetc(myfile);
}
}
// continue loop while the first character in line is header-like.
// Loop ends, either when the header has finished or the end of the file is reached.
} while (ch_ascii != EOF && line_start<39);
}
}
*/


/* Written by S. Brieden, 24.07.19
 * * For an already opened file* myfile, this function parses
 * * through the lines the user wants to skip and sets the file
 * * pointer starting at the interesting line.
 * */
void skiplines(FILE* myfile, long int nlines)
{
char line[1000];
int doc_start = -1;
long int linecounter;

// initialise
doc_start = fgetc(myfile);
linecounter=0;

//Skip blank spaces
while (doc_start == ' '){
  doc_start = fgetc(myfile);
}
//Skip lines
while (linecounter<nlines){
  fgets(line, 1000, myfile);
  linecounter++;
}

}

void determine_spacing_theo2(char *spacing,double **k0,int NeffP0)
{
     double check1,check2;
       double epsilon=1e-9;
       double check1log,check2log;
       double check1log10,check2log10;
       int mode0=0;

     check1=k0[1][0]-k0[0][0];
       check2=(k0[NeffP0-1][0]-k0[0][0])/(NeffP0*1.-1.);

       check1log=log(k0[1][0])-log(k0[0][0]);
       check2log=(log(k0[NeffP0-1][0])-log(k0[0][0]))/(NeffP0*1.-1.);

       check1log10=log10(k0[1][0])-log10(k0[0][0]);
       check2log10=(log10(k0[NeffP0-1][0])-log10(k0[0][0]))/(NeffP0*1.-1.);

       if(fabs(check1-check2)<epsilon){mode0=1;}//equally linearly spaced

       if(fabs(check1log-check2log)<epsilon){mode0=2;}//equally log-spaced

       if(fabs(check1log10-check2log10)<epsilon){mode0=3;}//equally log10-spaced

if(mode0==0){sprintf(spacing,"irregular");printf("Warning, PTcool file spaced irregularly. The mcmc performance will be significantly affected.\n");}
if(mode0==1){sprintf(spacing,"linear");printf("PTcool file spaced linearly.\n");}
if(mode0==2){sprintf(spacing,"log");printf("PTcool file spaced logly.\n");}
if(mode0==3){sprintf(spacing,"log10");printf("PTcool spaced log10ly.\n");}


}


void determine_spacing_theo(char *spacing,double *k0,int NeffP0,int mode)
{
     double check1,check2;
       double epsilon=1e-9;
       double check1log,check2log;
       double check1log10,check2log10;
       int mode0=0;
       char file[200];
       if(mode==0){sprintf(file,"Plin");}
       if(mode==1){sprintf(file,"Olin");}
       if(mode==2){sprintf(file,"Mask NGC");}
       if(mode==3){sprintf(file,"Mask SGC");}

    check1=k0[1]-k0[0];
       check2=(k0[NeffP0-1]-k0[0])/(NeffP0*1.-1.);

       check1log=log(k0[1])-log(k0[0]);
       check2log=(log(k0[NeffP0-1])-log(k0[0]))/(NeffP0*1.-1.);

       check1log10=log10(k0[1])-log10(k0[0]);
       check2log10=(log10(k0[NeffP0-1])-log10(k0[0]))/(NeffP0*1.-1.);

       if(fabs(check1-check2)<epsilon){mode0=1;}//equally linearly spaced

       if(fabs(check1log-check2log)<epsilon){mode0=2;}//equally log-spaced

       if(fabs(check1log10-check2log10)<epsilon){mode0=3;}//equally log10-spaced

if(mode0==0){sprintf(spacing,"irregular");printf("Warning, %s file spaced irregularly. The mcmc performance will be significantly affected.\n",file);}
if(mode0==1){sprintf(spacing,"linear");printf("%s file spaced linearly.\n",file);}
if(mode0==2){sprintf(spacing,"log");printf("%s file spaced logly.\n",file);}
if(mode0==3){sprintf(spacing,"log10");printf("%s file spaced log10ly.\n",file);}


}
void determine_spacing(char *spacing,double *k0, double *k2, double *k4,int NeffP0,int NeffP2,int NeffP4)
{

     double check1,check2;
       double epsilon=1e-5;//1e-9;
       double check1log,check2log;
       double check1log10,check2log10;
       int mode0=0;
       int mode2=0;
       int mode4=0;

       check1=k0[1]-k0[0];
       check2=(k0[NeffP0-1]-k0[0])/(NeffP0*1.-1.);

       check1log=log(k0[1])-log(k0[0]);
       check2log=(log(k0[NeffP0-1])-log(k0[0]))/(NeffP0*1.-1.);

       check1log10=log10(k0[1])-log10(k0[0]);
       check2log10=(log10(k0[NeffP0-1])-log10(k0[0]))/(NeffP0*1.-1.);

       if(fabs(check1-check2)<epsilon){mode0=1;}//equally linearly spaced

       if(fabs(check1log-check2log)<epsilon){mode0=2;}//equally log-spaced

       if(fabs(check1log10-check2log10)<epsilon){mode0=3;}//equally log10-spaced


       check1=k2[1]-k2[0];
       check2=(k2[NeffP2-1]-k2[0])/(NeffP2*1.-1.);

       check1log=log(k2[1])-log(k2[0]);
       check2log=(log(k2[NeffP2-1])-log(k2[0]))/(NeffP2*1.-1.);

       check1log10=log10(k2[1])-log10(k2[0]);
       check2log10=(log10(k2[NeffP2-1])-log10(k2[0]))/(NeffP2*1.-1.);

       if(fabs(check1-check2)<epsilon){mode2=1;}//equally linearly spaced

       if(fabs(check1log-check2log)<epsilon){mode2=2;}//equally log-spaced

       if(fabs(check1log10-check2log10)<epsilon){mode2=3;}//equally log10-spaced

       check1=k4[1]-k4[0];
       check2=(k4[NeffP4-1]-k4[0])/(NeffP4*1.-1.);

       check1log=log(k4[1])-log(k4[0]);
       check2log=(log(k4[NeffP4-1])-log(k4[0]))/(NeffP4*1.-1.);

       check1log10=log10(k4[1])-log10(k4[0]);
       check2log10=(log10(k4[NeffP4-1])-log10(k2[0]))/(NeffP4*1.-1.);

       if(fabs(check1-check2)<epsilon){mode4=1;}//equally linearly spaced

       if(fabs(check1log-check2log)<epsilon){mode4=2;}//equally log-spaced

       if(fabs(check1log10-check2log10)<epsilon){mode4=3;}//equally log10-spaced


if(mode0 != mode2 || mode0 != mode4 || mode2 != mode4 ){printf("Error, bins across multipoles not equally spaced!. Exiting now...\n");exit(0);}

if(mode0==0){sprintf(spacing,"irregular");printf("Warning, data spaced irregularly. The mcmc performance will be significantly affected.\n");}
if(mode0==1){sprintf(spacing,"linear");printf("Data spaced linearly.\n");}
if(mode0==2){sprintf(spacing,"log");printf("Data spaced logly.\n");}
if(mode0==3){sprintf(spacing,"log10");printf("Data spaced log10ly.\n");}

}

void string_copy(char *from, char *to) {

    while ((*to++ = *from++) != '\0')
        ;
}

int countlines(char *filename)
{
FILE* myfile = fopen(filename, "r");
if(myfile==NULL){printf("Error reading %s. Exiting now...\n",filename);exit(0);}
int ch, number_of_lines = -1;

skipheader(myfile);

do
{
            ch = fgetc(myfile);
                    if(ch == '\n')
                                        number_of_lines++;
} while (ch != EOF);

if(ch != '\n' && number_of_lines != 0)
    number_of_lines++;

        fclose(myfile);
        return number_of_lines;
}

int count_all_lines(char *filename)
{
FILE* myfile = fopen(filename, "r");
int ch, number_of_lines = -1;

do
{
            ch = fgetc(myfile);
                    if(ch == '\n')
                                        number_of_lines++;
} while (ch != EOF);

if(ch != '\n' && number_of_lines != 0)
    number_of_lines++;

        fclose(myfile);
        return number_of_lines;
}

long int countlinesLI(char *filename)
{
FILE* myfile = fopen(filename, "r");
if(myfile==NULL){printf("Error reading %s. Exiting now...\n",filename);exit(0);}
long int ch, number_of_lines = -1;

skipheader(myfile);

do
{
            ch = fgetc(myfile);
                    if(ch == '\n')
                                        number_of_lines++;
} while (ch != EOF);

if(ch != '\n' && number_of_lines != 0)
    number_of_lines++;

        fclose(myfile);
        return number_of_lines;
}

void freeTokens(double **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{    
     free(tokens[i]);
}
free(tokens);

}

void freeTokens2(double ***tokens, int N1, int N2)
{
int i,j;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2;++j)
     {
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}



void freeTokensInt(int **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{
     free(tokens[i]);
}
free(tokens);
}


double P_interpolLOG(double k0, double *k, double *P, int N)
{
double P0,m,n;
int i;
i=-1;
do
{
i=i+1;
}while(k0>k[i] && i<N-1);
if(i==0)
{
m=( log10(P[i]) - log10(P[i+1]) )/( log10(k[i]) - log10(k[i+1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
else
{
m=( log10(P[i]) - log10(P[i-1]) )/( log10(k[i]) - log10(k[i-1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
return P0;
}

double P_interpol(double k0, double *k, double *P, int N)
{
double P0,m,n;
int i;

if(k0<k[0] || k0>k[N-1])
{
P0=0;
}
else
{
i=-1;
do
{
i=i+1;
}while(k0>k[i] && i<N-1);
if(i==0)
{
m=( (P[i]) - (P[i+1]) )/( (k[i]) - (k[i+1]) );
n=(P[i])-m*(k[i]);
P0=m*(k0)+n;
}
else
{
m=( (P[i]) - (P[i-1]) )/( (k[i]) - (k[i-1]) );
n=(P[i])-m*(k[i]);
P0=m*(k0)+n;
}

}
//printf("%d -> %d: %lf+(%lf-%lf)*(%lf-%lf)/(%lf-%lf), w=%lf\n",i-1,i,P[i],k0,k[i],P[i],P[i-1],k[i],k[i-1],(k0-k[i])/( (k[i]) - (k[i-1])));
return P0;
}

double P_interpol_extrapol(double k0, double *k, double *P, int N)
{
double P0,m,n;
int i;

if(k0<k[0] || k0>k[N-1])//extrapol
{
//P0=0;

if(k0<k[0])
{
m=( (P[0]) - (P[1]) )/( (k[0]) - (k[1]) );
n=(P[0])-m*(k[0]);
P0=m*(k0)+n;
}
if(k0>k[N-1])
{
m=( (P[N-1]) - (P[N-2]) )/( (k[N-1]) - (k[N-2]) );
n=(P[N-1])-m*(k[N-1]);
P0=m*(k0)+n;
}


}
else
{
i=-1;
do
{
i=i+1;
}while(k0>k[i] && i<N-1);
if(i==0)
{
m=( (P[i]) - (P[i+1]) )/( (k[i]) - (k[i+1]) );
n=(P[i])-m*(k[i]);
P0=m*(k0)+n;
}
else
{
m=( (P[i]) - (P[i-1]) )/( (k[i]) - (k[i-1]) );
n=(P[i])-m*(k[i]);
P0=m*(k0)+n;
}

}
return P0;
}

double P_interpol_doublearray2(double k0, double **k,double *P, int N)
{
double P0,m,n;
int i;

if(k0<k[0][0] || k0>k[N-1][0])
{
P0=0;
}
else
{
i=-1;
do
{
i=i+1;
}while(k0>k[i][0] && i<N-1);
if(i==0)
{
m=( (P[i]) - (P[i+1]) )/( (k[i][0]) - (k[i+1][0]) );
n=(P[i])-m*(k[i][0]);
P0=m*(k0)+n;
}
else
{
m=( (P[i]) - (P[i-1]) )/( (k[i][0]) - (k[i-1][0]) );
n=(P[i])-m*(k[i][0]);
P0=m*(k0)+n;
}

}
return P0;

}

double P_interpol_doublearray(double k0, int index, double **P, int N)
{
double P0,m,n;
int i;

if(k0<P[0][0] || k0>P[N-1][0])
{
P0=0;
}
else
{
i=-1;
do
{
i=i+1;
}while(k0>P[i][0] && i<N-1);
if(i==0)
{
m=( (P[i][index]) - (P[i+1][index]) )/( (P[i][0]) - (P[i+1][0]) );
n=(P[i][index])-m*(P[i][0]);
P0=m*(k0)+n;
}
else
{
m=( (P[i][index]) - (P[i-1][index]) )/( (P[i][0]) - (P[i-1][0]) );
n=(P[i][index])-m*(P[i][0]);
P0=m*(k0)+n;
}

}
return P0;
}


    double determine_w1_doublearray(double **Theory,double kinput,int N1,char *spacing)
   {
       double w1;
if(strcmp(spacing,"linear") == 0 ){w1=(Theory[N1+1][0]-kinput)/(Theory[N1+1][0]-Theory[N1][0]);}
if(strcmp(spacing,"log") == 0 ){w1=(log(Theory[N1+1][0])-log(kinput))/(log(Theory[N1+1][0])-log(Theory[N1][0]));}
if(strcmp(spacing,"log10") == 0 ){w1=(log10(Theory[N1+1][0])-log10(kinput))/(log10(Theory[N1+1][0])-log10(Theory[N1][0]));}

       return w1;
   }

    double determine_w1_singlearray(double *k,double kinput,int N1, char *spacing)
   {
       double w1;
if(strcmp(spacing,"linear") == 0 ){w1=(k[N1+1]-kinput)/(k[N1+1]-k[N1]);}
if(strcmp(spacing,"log") == 0 ){w1=(log(k[N1+1])-log(kinput))/(log(k[N1+1])-log(k[N1]));}
if(strcmp(spacing,"log10") == 0 ){w1=(log10(k[N1+1])-log10(kinput))/(log10(k[N1+1])-log10(k[N1]));}

       return w1;
   }

double P_interpol_w1_doublearray(double **Theory,int index, int N1,double w1,char *spacing)
{
double a;
if(strcmp(spacing,"linear") == 0 ){a=w1*Theory[N1][index]+(1.-w1)*Theory[N1+1][index];}

if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0  ){a=w1*log(Theory[N1][index])+(1.-w1)*log(Theory[N1+1][index]);a=exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0  ){a=w1*log(-Theory[N1][index])+(1.-w1)*log(-Theory[N1+1][index]);a=-exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0  ){a=w1*(Theory[N1][index])+(1.-w1)*(Theory[N1+1][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0  ){a=w1*(Theory[N1][index])+(1.-w1)*(Theory[N1+1][index]);}

if(strcmp(spacing,"log10") == 0 ){a=w1*log10(Theory[N1][index])+(1.-w1)*log10(Theory[N1+1][index]);a=pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0  ){a=w1*log10(-Theory[N1][index])+(1.-w1)*log10(-Theory[N1+1][index]);a=-pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0  ){a=w1*(Theory[N1][index])+(1.-w1)*(Theory[N1+1][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0  ){a=w1*(Theory[N1][index])+(1.-w1)*(Theory[N1+1][index]);}

return a;
}

double P_interpol_w1_singlearray(double *Theory, int N1,double w1, char *spacing)
{
double a;
if(strcmp(spacing,"linear") == 0 ){a=w1*Theory[N1]+(1.-w1)*Theory[N1+1];}

if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]>0){a=w1*log(Theory[N1])+(1.-w1)*log(Theory[N1+1]);a=exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]<0 ){a=w1*log(-Theory[N1])+(1.-w1)*log(-Theory[N1+1]);a=-exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]<0 ){a=w1*(Theory[N1])+(1.-w1)*(Theory[N1+1]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]>0 ){a=w1*(Theory[N1])+(1.-w1)*(Theory[N1+1]);}

if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]>0){a=w1*log10(Theory[N1])+(1.-w1)*log10(Theory[N1+1]);a=pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]<0){a=w1*log10(-Theory[N1])+(1.-w1)*log10(-Theory[N1+1]);a=-pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]<0){a=w1*(Theory[N1])+(1.-w1)*(Theory[N1+1]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]>0){a=w1*(Theory[N1])+(1.-w1)*(Theory[N1+1]);}

return a;
}


double P_interpol_w012_singlearray(double *Theory, int N1,double w0, double w1, double w2, char *spacing)
{

double a;
if(strcmp(spacing,"linear") == 0){a=w0*Theory[N1]+w1*Theory[N1+1]+w2*Theory[N1+2];}

if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]>0 && Theory[N1+2]>0){a=w0*log(Theory[N1])+w1*log(Theory[N1+1])+w2*log(Theory[N1+2]);a=exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]<0 && Theory[N1+2]<0){a=w0*log(-Theory[N1])+w1*log(-Theory[N1+1])+w2*log(-Theory[N1+2]);a=-exp(a);}

if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]<0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]>0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]<0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}

if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]>0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]<0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]>0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}


if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]>0 && Theory[N1+2]>0){a=w0*log10(Theory[N1])+w1*log10(Theory[N1+1])+w2*log10(Theory[N1+2]);a=pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]<0 && Theory[N1+2]<0){a=w0*log10(-Theory[N1])+w1*log10(-Theory[N1+1])+w2*log10(-Theory[N1+2]);a=-pow(10,a);}

if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]<0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]>0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]<0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}

if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]>0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]<0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]>0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}


return a;
}

double P_interpol_w012_doublearray(double **Theory,int index, int N1,double w0, double w1, double w2,char *spacing)
{
double a;
if(strcmp(spacing,"linear") == 0){a=w0*Theory[N1][index]+w1*Theory[N1+1][index]+w2*Theory[N1+2][index];}

//if(strcmp(spacing,"log") == 0){a=w0*log(Theory[N1][index])+w1*log(Theory[N1+1][index])+w2*log(Theory[N1+2][index]);a=exp(a);}
//if(strcmp(spacing,"log10") == 0){a=w0*log10(Theory[N1][index])+w1*log10(Theory[N1+1][index])+w2*log10(Theory[N1+2][index]);a=pow(10,a);}

if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0 && Theory[N1+2][index]>0){a=w0*log(Theory[N1][index])+w1*log(Theory[N1+1][index])+w2*log(Theory[N1+2][index]);a=exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0 && Theory[N1+2][index]<0){a=w0*log(-Theory[N1][index])+w1*log(-Theory[N1+1][index])+w2*log(-Theory[N1+2][index]);a=-exp(a);}

if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}

if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}


if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0 && Theory[N1+2][index]>0){a=w0*log10(Theory[N1][index])+w1*log10(Theory[N1+1][index])+w2*log10(Theory[N1+2][index]);a=pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0 && Theory[N1+2][index]<0){a=w0*log10(-Theory[N1][index])+w1*log10(-Theory[N1+1][index])+w2*log10(-Theory[N1+2][index]);a=-pow(10,a);}

if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}

if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}


return a;
}
         
  double determine_w0_2ndorder_singlearray(double *k,double kinput,int N1, char *spacing)
   {
       double w1;
           if(strcmp(spacing,"linear") == 0){w1=(kinput-k[N1+1])*(kinput-k[N1+2])/((k[N1]-k[N1+1] )*(k[N1]-k[N1+2]));}
           if(strcmp(spacing,"log") == 0){w1=(log(kinput)-log(k[N1+1]))*(log(kinput)-log(k[N1+2]))/((log(k[N1])-log(k[N1+1]) )*(log(k[N1])-log(k[N1+2])));}
           if(strcmp(spacing,"log10") == 0){w1=(log10(kinput)-log10(k[N1+1]))*(log10(kinput)-log10(k[N1+2]))/((log10(k[N1])-log10(k[N1+1]) )*(log10(k[N1])-log10(k[N1+2])));}


       return w1;
   }

    double determine_w1_2ndorder_singlearray(double *k,double kinput,int N1,char *spacing)
   {
       double w1;
      if(strcmp(spacing,"linear") == 0){w1=(kinput-k[N1])*(kinput-k[N1+2])/((k[N1+1]-k[N1] )*(k[N1+1]-k[N1+2]));}
      if(strcmp(spacing,"log") == 0){w1=(log(kinput)-log(k[N1]))*(log(kinput)-log(k[N1+2]))/((log(k[N1+1])-log(k[N1]) )*(log(k[N1+1])-log(k[N1+2])));}
      if(strcmp(spacing,"log10") == 0){w1=(log10(kinput)-log10(k[N1]))*(log10(kinput)-log10(k[N1+2]))/((log10(k[N1+1])-log10(k[N1]) )*(log10(k[N1+1])-log10(k[N1+2])));}

       return w1;
   }

    double determine_w2_2ndorder_singlearray(double *k,double kinput,int N1,char *spacing)
   {
       double w1;
       if(strcmp(spacing,"linear") == 0){w1=(kinput-k[N1])*(kinput-k[N1+1])/((k[N1+2]-k[N1])*(k[N1+2]-k[N1+1]));}
       if(strcmp(spacing,"log") == 0){w1=(log(kinput)-log(k[N1]))*(log(kinput)-log(k[N1+1]))/((log(k[N1+2])-log(k[N1]))*(log(k[N1+2])-log(k[N1+1])));}
       if(strcmp(spacing,"log10") == 0){w1=(log10(kinput)-log10(k[N1]))*(log10(kinput)-log10(k[N1+1]))/((log10(k[N1+2])-log10(k[N1]))*(log10(k[N1+2])-log10(k[N1+1])));}


       return w1;
   }
  double determine_w0_2ndorder_doublearray(double **k,double kinput,int N1,char *spacing)
   {
       double w1;
       if(strcmp(spacing,"linear") == 0){w1=(kinput-k[N1+1][0])*(kinput-k[N1+2][0])/((k[N1][0]-k[N1+1][0] )*(k[N1][0]-k[N1+2][0]));}
       if(strcmp(spacing,"log") == 0){w1=(log(kinput)-log(k[N1+1][0]))*(log(kinput)-log(k[N1+2][0]))/((log(k[N1][0])-log(k[N1+1][0]))*(log(k[N1][0])-log(k[N1+2][0])));}
       if(strcmp(spacing,"log10") == 0){w1=(log10(kinput)-log10(k[N1+1][0]))*(log10(kinput)-log10(k[N1+2][0]))/((log10(k[N1][0])-log10(k[N1+1][0]))*(log10(k[N1][0])-log10(k[N1+2][0])));}


       return w1;
   }

    double determine_w1_2ndorder_doublearray(double **k,double kinput,int N1, char *spacing)
   {
       double w1;
       if(strcmp(spacing,"linear") == 0){w1=(kinput-k[N1][0])*(kinput-k[N1+2][0])/((k[N1+1][0]-k[N1][0] )*(k[N1+1][0]-k[N1+2][0]));}
       if(strcmp(spacing,"log") == 0){w1=(log(kinput)-log(k[N1][0]))*(log(kinput)-log(k[N1+2][0]))/((log(k[N1+1][0])-log(k[N1][0]))*(log(k[N1+1][0])-log(k[N1+2][0])));}
       if(strcmp(spacing,"log10") == 0){w1=(log10(kinput)-log10(k[N1][0]))*(log10(kinput)-log10(k[N1+2][0]))/((log10(k[N1+1][0])-log10(k[N1][0]))*(log10(k[N1+1][0])-log10(k[N1+2][0])));}


       return w1;
   }

  double determine_w2_2ndorder_doublearray(double **k,double kinput,int N1, char *spacing)
   {
       double w1;
       if(strcmp(spacing,"linear") == 0){w1=(kinput-k[N1][0])*(kinput-k[N1+1][0])/((k[N1+2][0]-k[N1][0])*(k[N1+2][0]-k[N1+1][0]));}
       if(strcmp(spacing,"log") == 0){w1=(log(kinput)-log(k[N1][0]))*(log(kinput)-log(k[N1+1][0]))/((log(k[N1+2][0])-log(k[N1][0]))*(log(k[N1+2][0])-log(k[N1+1][0])));}
       if(strcmp(spacing,"log10") == 0){w1=(log10(kinput)-log10(k[N1][0]))*(log10(kinput)-log10(k[N1+1][0]))/((log10(k[N1+2][0])-log10(k[N1][0]))*(log10(k[N1+2][0])-log10(k[N1+1][0])));}


       return w1;
   }


double P_interpol_fast_doublearray(double k0, double **P, int index, int N, char *spacing,int interpolation_order, int Ninterpol,double w0,double w1,double w2 )
{
int shiftN;
double Pout;
interpolation_order=1;//1 or 2
shiftN=interpolation_order;
if(Ninterpol>=N-shiftN || Ninterpol<0 || k0<=0){Pout=0;}
else{
if(interpolation_order==1){Pout=P_interpol_w1_doublearray(P,index,Ninterpol,w1,spacing);}
if(interpolation_order==2){Pout=P_interpol_w012_doublearray(P,index,Ninterpol,w0,w1,w2,spacing);}
}

return Pout;
}

double P_interpol_fast(double k0, double *P, int N, char *spacing,int interpolation_order, int Ninterpol,double w0,double w1,double w2 )
{
int shiftN;
double Pout;

shiftN=interpolation_order;


if(Ninterpol>=N-shiftN || Ninterpol<0 || k0<=0){Pout=0;}
else{
if(interpolation_order==1){Pout=P_interpol_w1_singlearray(P,Ninterpol,w1,spacing);}
if(interpolation_order==2){Pout=P_interpol_w012_singlearray(P,Ninterpol,w0,w1,w2,spacing);}
}

return Pout;
}

double get_error2(double *alpha, double *chi2, double chi2min, int imin,int N, double Deltachi2)
{
int l;
double m,n;
int lpivot=-1;
double error=0;
for(l=N-1;l>1;l--)
{
if(chi2[l]-chi2min-Deltachi2>0 && chi2[l-1]-chi2min-Deltachi2<=0){
lpivot=l;//printf("lpivot=%d\n",l);
}
if(chi2[l]-chi2min-Deltachi2<0 && chi2[l-1]-chi2min-Deltachi2>=0){
lpivot=l;//printf("lpivot=%d\n",l);
}

}
if(lpivot==-1){/*printf("Error in get error1. Set error to 0\n");*/error=0;}
else{
m=(alpha[lpivot]-alpha[lpivot-1])/(chi2[lpivot]-chi2[lpivot-1]);
n=alpha[lpivot]-m*chi2[lpivot];
error=-( (chi2min+Deltachi2)*m+n   )+alpha[imin];
}
//printf("%lf %d\n",error,lpivot);

if(error<0){error=0;}
return error;
}

double get_error1(double *alpha, double *chi2, double chi2min, int imin, int N, double Deltachi2)
{
int l;
double m,n;
int lpivot=-1;
double error=0;

lpivot=-1;
error=0;
for(l=0;l<N-1;l++)
{
if(chi2[l]-chi2min-Deltachi2<0 && chi2[l+1]-chi2min-Deltachi2>=0){
lpivot=l;
}
if(chi2[l]-chi2min-Deltachi2>0 && chi2[l+1]-chi2min-Deltachi2<=0){
lpivot=l;
}
}
if(lpivot==-1){/*printf("Error in get error2. Set error to 0\n");*/error=0;}
else{
m=(alpha[lpivot]-alpha[lpivot+1])/(chi2[lpivot]-chi2[lpivot+1]);
n=alpha[lpivot]-m*chi2[lpivot];
error=( (chi2min+Deltachi2)*m+n   )-alpha[imin];
}

return error;
}


 int determine_N_doublearray(double **Theory,double kinput,int Nlin, char *spacing)
{
       int N1;
       int i;

  double check1;
       double check1log;
       double check1log10;
     
       check1=(Theory[Nlin-1][0]-Theory[0][0])/(Nlin*1.-1.);
       check1log=(log(Theory[Nlin-1][0])-log(Theory[0][0]))/(Nlin*1.-1.);
       check1log10=(log10(Theory[Nlin-1][0])-log10(Theory[0][0]))/(Nlin*1.-1.);

if(kinput<Theory[0][0] || kinput>Theory[Nlin-1][0] || kinput<=0){N1=-1;}
else{
       if(strcmp(spacing,"linear") == 0){N1=(int)((kinput-Theory[0][0])/check1);}
       if(strcmp(spacing,"log") == 0){N1=(int)((log(kinput)-log(Theory[0][0]))/check1log);}
       if(strcmp(spacing,"log10") == 0){N1=(int)((log10(kinput)-log10(Theory[0][0]))/check1log10);}
       if(strcmp(spacing,"irregular") == 0){
       i=-1;
            do{

               i++;
               }while(kinput>Theory[i][0]);
               N1=i-1;

       }

}
       return N1;

}

 int determine_N_singlearray(double *Theory,double kinput,int Nlin, char *spacing)
{

       int N1;
       int i;

  double check1;
       double check1log;
       double check1log10;

       check1=(Theory[Nlin-1]-Theory[0])/(Nlin*1.-1.);
       check1log=(log(Theory[Nlin-1])-log(Theory[0]))/(Nlin*1.-1.);
       check1log10=(log10(Theory[Nlin-1])-log10(Theory[0]))/(Nlin*1.-1.);

if(kinput<Theory[0] || kinput>Theory[Nlin-1] || kinput<=0){N1=-1;}
else{
       if(strcmp(spacing,"linear") == 0){N1=(int)((kinput-Theory[0])/check1);}
       if(strcmp(spacing,"log") == 0){N1=(int)((log(kinput)-log(Theory[0]))/check1log);}
       if(strcmp(spacing,"log10") == 0){N1=(int)((log10(kinput)-log10(Theory[0]))/check1log10);}
       if(strcmp(spacing,"irregular") == 0){     
       i=-1;
            do{

               i++;
               }while(kinput>Theory[i]);
               N1=i-1;

       }

}
       return N1;


}

long int getMin(long int arr[], int n) 
{   
    int i;
    int res = arr[0]; 
    for (i = 1; i < n; i++) 
        res = MIN(res, arr[i]); 
    return res; 
} 
  
long int getMax(long int arr[], int n) 
{ 
    int i;
    int res = arr[0]; 
    for (i = 1; i < n; i++) 
        res = MAX(res, arr[i]); 
    return res; 
}


int re_sampling(char *spacing, int Neff,double kmin,double kmax)
{
int N;

if( strcmp(spacing,"linear") == 0  ){
N=Neff-2+25+(int)(kmin/(kmax-kmin)*(Neff-1));
}

if( strcmp(spacing,"log") == 0  ){
N=Neff-2+60+(int)(log(kmin)/(log(kmax)-log(kmin))*(Neff-1));
}

if( strcmp(spacing,"log10") == 0  ){
N=Neff-2+60+(int)(log10(kmin)/(log10(kmax)-log10(kmin))*(Neff-1));
}
if( strcmp(spacing,"irregular") == 0  ){
N=0;
}

return N;
}

int  get_Neffmax(char *spacing_data, int modeP0, int modeP2, int modeP4, int NeffP0, int NeffP2, int NeffP4, double k0min, double k0max, double k2min, double k2max, double k4min, double k4max,int Ntheo )
{
int Neffmax;

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=re_sampling(spacing_data,NeffP0,k0min,k0max);}//    NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=re_sampling(spacing_data,NeffP0,k0min,k0max);}//    NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=re_sampling(spacing_data,NeffP2,k2min,k2max);}//    NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=re_sampling(spacing_data,NeffP2,k2min,k2max);}//    NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=re_sampling(spacing_data,NeffP4,k4min,k4max);}//    NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=re_sampling(spacing_data,NeffP4,k4min,k4max);}//    NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}

/*
if( strcmp(spacing_data,"linear") == 0  ){

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}

}

if( strcmp(spacing_data,"log") == 0  ){

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(log(k0[0])/(log(k0[NeffP0-1])-log(k0[0]))*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(log(k0[0])/(log(k0[NeffP0-1])-log(k0[0]))*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(log(k2[0])/(log(k2[NeffP2-1])-log(k2[0]))*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(log(k2[0])/(log(k2[NeffP2-1])-log(k2[0]))*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(log(k4[0])/(log(k4[NeffP4-1])-log(k4[0]))*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(log(k4[0])/(log(k4[NeffP4-1])-log(k4[0]))*(NeffP4-1));}

}

if( strcmp(spacing_data,"log10") == 0  ){

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(log10(k0[0])/(log10(k0[NeffP0-1])-log10(k0[0]))*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(log10(k0[0])/(log10(k0[NeffP0-1])-log10(k0[0]))*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(log10(k2[0])/(log10(k2[NeffP2-1])-log10(k2[0]))*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(log10(k2[0])/(log10(k2[NeffP2-1])-log10(k2[0]))*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(log10(k4[0])/(log10(k4[NeffP4-1])-log10(k4[0]))*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(log10(k4[0])/(log10(k4[NeffP4-1])-log10(k4[0]))*(NeffP4-1));}

}
*/
if( strcmp(spacing_data,"irregular") == 0  ){
Neffmax=Ntheo;
//factor_sampling_mask=1;
}



return Neffmax;
}

void get_NeffP(int params[], int NeffP0, int NeffP2, int NeffP4, char *spacing_data, double kmin, double kmax, double k0min, double k0max,double k2min,double k2max,double k4min,double k4max, int N_Plin)
{
int NeffP,NeffP_max;

NeffP=0;
//if(strcmp(spacing_data,"linear") == 0){
if(kmax==k0max && kmin==k0min){NeffP=re_sampling(spacing_data,NeffP0,k0min,k0max);/*NeffP0+25-2+(int)(k0min/(k0max-k0min)*(NeffP0-1.));*/NeffP_max=NeffP0;}
if(kmax==k2max && kmin==k2min){NeffP=re_sampling(spacing_data,NeffP2,k2min,k2max);/*NeffP2+25-2+(int)(k2min/(k2max-k2min)*(NeffP2-1.));*/NeffP_max=NeffP2;}
if(kmax==k4max && kmin==k4min){NeffP=re_sampling(spacing_data,NeffP4,k4min,k4max);/*NeffP4+25-2+(int)(k4min/(k4max-k4min)*(NeffP4-1.));*/NeffP_max=NeffP4;}
//}
/*
if(strcmp(spacing_data,"log") == 0){
if(kmax==k0max && kmin==k0min){NeffP=NeffP0+25-2+(int)(log(k0min)/(log(k0max)-log(k0min))*(NeffP0-1.));NeffP_max=NeffP0;}
if(kmax==k2max && kmin==k2min){NeffP=NeffP2+25-2+(int)(log(k2min)/(log(k2max)-log(k2min))*(NeffP2-1.));NeffP_max=NeffP2;}
if(kmax==k4max && kmin==k4min){NeffP=NeffP4+25-2+(int)(log(k4min)/(log(k4max)-log(k4min))*(NeffP4-1.));NeffP_max=NeffP4;}
}

if(strcmp(spacing_data,"log10") == 0){

if(kmax==k0max && kmin==k0min){NeffP=NeffP0+25-2+(int)(log10(k0min)/(log10(k0max)-log10(k0min))*(NeffP0-1.));NeffP_max=NeffP0;}
if(kmax==k2max && kmin==k2min){NeffP=NeffP2+25-2+(int)(log10(k2min)/(log10(k2max)-log10(k2min))*(NeffP2-1.));NeffP_max=NeffP2;}
if(kmax==k4max && kmin==k4min){NeffP=NeffP4+25-2+(int)(log10(k4min)/(log10(k4max)-log10(k4min))*(NeffP4-1.));NeffP_max=NeffP4;}

}
*/
if(strcmp(spacing_data,"irregular") == 0){

NeffP=N_Plin;
NeffP_max=0;
}

if(NeffP== 0){printf("Warning, crossing intervals among P0 (%lf<k<%lf), P2 (%lf<k<%lf), P4 (%lf<k<%lf) [%lf,%lf] make impossible to determine NeffP within mask application. Fix this\n",k0min,k0max,k2min,k2max,k4min,k4max,kmin,kmax);exit(0);}


params[0]=NeffP;
params[1]=NeffP_max;
}

void get_kmin_kmax(double params[] ,int modeP0,int modeP2,int modeP4,double k0min,double k0max,double k2min,double k2max,double k4min,double k4max)
{
double kmin,kmax;
kmin=-1;
kmax=-1;
if(k0max>=k2max && modeP0==1 && modeP2==1){kmax=k0max;}
if(k0max>=k4max && modeP0==1 && modeP4==1){kmax=k0max;}

if(k2max>=k0max && modeP0==1 && modeP2==1){kmax=k2max;}
if(k2max>=k4max && modeP4==1 && modeP2==1){kmax=k2max;}

if(k4max>=k0max && modeP0==1 && modeP4==1){kmax=k4max;}
if(k4max>=k2max && modeP4==1 && modeP2==1){kmax=k4max;}

if(k0min<=k2min && modeP0==1 && modeP2==1){kmin=k0min;}
if(k0min<=k4min && modeP0==1 && modeP4==1){kmin=k0min;}

if(k2min<=k0min && modeP0==1 && modeP2==1){kmin=k2min;}
if(k2min<=k4min && modeP4==1 && modeP2==1){kmin=k2min;}

if(k4min<=k0min && modeP0==1 && modeP4==1){kmin=k4min;}
if(k4min<=k2min && modeP4==1 && modeP2==1){kmin=k4min;}

if(modeP0==1 && modeP2==0 && modeP4==0){kmax=k0max;kmin=k0min;}
if(modeP2==1 && modeP0==0 && modeP4==0){kmax=k2max;kmin=k2min;}
if(modeP4==1 && modeP2==0 && modeP0==0){kmax=k4max;kmin=k4min;}

if(kmin<0 || kmax<0){printf("Error in function get_kmin_kmax, kmin=%lf, kmax=%lf. k0min=%lf k0max=%lf (%d); k2min=%lf k2max=%lf (%d); k4min=%lf k4max=%lf (%d). Exiting now...\n",kmin,kmax,k0min,k0max,modeP0,k2min,k2max,modeP2,k4min,k4max,modeP4);exit(0);}

params[0]=kmin;
params[1]=kmax;
}

double get_ktheo(char *spacing_data,int j_k,double kmin,double kmax,int NeffP_max,int factor_for_sampling, double *k_theory, double **Theory)
{
double k;
if(strcmp(spacing_data,"linear") == 0){k=(j_k+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);}

if(strcmp(spacing_data,"log") == 0){k=exp(  (j_k+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );}

if(strcmp(spacing_data,"log10") == 0){k=pow(10,  (j_k+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );}

if(strcmp(spacing_data,"irregular") == 0){ 

    if(Theory==NULL){k=k_theory[j_k];}

    if(k_theory==NULL){k=Theory[j_k][0];}

} 


return k;
}


void apply_mask_matrix( double *P_theo0, double *P_theo2,double *P_theo4,double **Matrix_mask, double *vector_in, int Nin, int Nout, int modeP0, int modeP2, int modeP4, int NeffP0, int NeffP2, int NeffP4, char *type_of_analysis)
{
int i,j,l;
double p_out;
int modeP0_trial,modeP2_trial,modeP4_trial;
modeP0_trial=modeP0;
modeP2_trial=modeP2;
modeP4_trial=modeP4;

l=0;
for(j=0;j<Nout;j++)
{
    p_out=0;
    for(i=0;i<Nin;i++)
    {
    p_out=p_out+Matrix_mask[i][j]*vector_in[i];
    }

//printf("P[%d]=%lf\n",l,p_out);
    if(l<NeffP4*modeP4_trial && modeP0_trial==0 && modeP2_trial==0 && modeP4_trial>0){
          P_theo4[l]=p_out;l++;
          if(l==NeffP4*modeP4_trial){l=0;modeP4_trial=0;}
    }

  if(l<NeffP2*modeP2_trial && modeP0_trial==0 && modeP2_trial>0){
          P_theo2[l]=p_out;l++;
          if(l==NeffP2*modeP2_trial){l=0;modeP2_trial=0;}
    }

    if(l<NeffP0*modeP0_trial && modeP0_trial>0){
          P_theo0[l]=p_out;l++;
          if(l==NeffP0*modeP0_trial){l=0;modeP0_trial=0;}
     }

}


//end
}


void order(double ctheta12,double ctheta13,double ctheta23,double params[])
{
double max,min,med;

if(ctheta12>=ctheta13 && ctheta12>=ctheta23)
{
   max=ctheta12;
   if(ctheta13>=ctheta23)
   {
     med=ctheta13;
     min=ctheta23;
   }
   else
   {
     min=ctheta13;
     med=ctheta23;
   }

}
else
{

   if(ctheta13>=ctheta12 && ctheta13>=ctheta23)
   {
      max=ctheta13;
      if(ctheta12>=ctheta23)
      {
         med=ctheta12;
         min=ctheta23;
      }
      else
      {
         min=ctheta12;
         med=ctheta23;
      }

   }
   else
   {
     max=ctheta23;
     if(ctheta12>=ctheta13)
     {
       med=ctheta12;
       min=ctheta13;
     }
     else
     {
       min=ctheta12;
       med=ctheta13;
     }

   }


}


params[0]=max;
params[1]=med;
params[2]=min;

}

