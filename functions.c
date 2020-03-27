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

/* Written by S. Brieden, 30.01.19
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


void determine_spacing_theo(char *spacing,double *k0,int NeffP0)
{
     double check1,check2;
       double epsilon=1e-9;
       double check1log,check2log;
       double check1log10,check2log10;
       int mode0=0;

    check1=k0[1]-k0[0];
       check2=(k0[NeffP0-1]-k0[0])/(NeffP0*1.-1.);

       check1log=log(k0[1])-log(k0[0]);
       check2log=(log(k0[NeffP0-1])-log(k0[0]))/(NeffP0*1.-1.);

       check1log10=log10(k0[1])-log10(k0[0]);
       check2log10=(log10(k0[NeffP0-1])-log10(k0[0]))/(NeffP0*1.-1.);

       if(fabs(check1-check2)<epsilon){mode0=1;}//equally linearly spaced

       if(fabs(check1log-check2log)<epsilon){mode0=2;}//equally log-spaced

       if(fabs(check1log10-check2log10)<epsilon){mode0=3;}//equally log10-spaced

if(mode0==0){sprintf(spacing,"irregular");printf("Warning, Plin/Olin file spaced irregularly. The mcmc performance will be significantly affected.\n");}
if(mode0==1){sprintf(spacing,"linear");printf("Plin/Olin file spaced linearly.\n");}
if(mode0==2){sprintf(spacing,"log");printf("Plin/Olin file spaced logly.\n");}
if(mode0==3){sprintf(spacing,"log10");printf("Plin/Olin spaced log10ly.\n");}


}
void determine_spacing(char *spacing,double *k0, double *k2, double *k4,int NeffP0,int NeffP2,int NeffP4)
{

     double check1,check2;
       double epsilon=1e-9;
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
}while(k0>k[i] && i<N);
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
}while(k0>k[i] && i<N);
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
}while(k0>P[i][0] && i<N);
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

double get_error2(double *alpha, double *chi2, double chi2min, int imin,int N)
{
int l;
double m,n;
int lpivot=-1;
double error=0;
for(l=N-1;l>1;l--)
{
if(chi2[l]-chi2min-1>0 && chi2[l-1]-chi2min-1<=0){
lpivot=l;//printf("lpivot=%d\n",l);
}
if(chi2[l]-chi2min-1<0 && chi2[l-1]-chi2min-1>=0){
lpivot=l;//printf("lpivot=%d\n",l);
}

}
if(lpivot==-1){/*printf("Error in get error1. Set error to 0\n");*/error=0;}
else{
m=(alpha[lpivot]-alpha[lpivot-1])/(chi2[lpivot]-chi2[lpivot-1]);
n=alpha[lpivot]-m*chi2[lpivot];
error=-( (chi2min+1)*m+n   )+alpha[imin];
}
//printf("%lf %d\n",error,lpivot);

if(error<0){error=0;}
return error;
}

double get_error1(double *alpha, double *chi2, double chi2min, int imin, int N)
{
int l;
double m,n;
int lpivot=-1;
double error=0;

lpivot=-1;
error=0;
for(l=0;l<N-1;l++)
{
if(chi2[l]-chi2min-1<0 && chi2[l+1]-chi2min-1>=0){
lpivot=l;
}
if(chi2[l]-chi2min-1>0 && chi2[l+1]-chi2min-1<=0){
lpivot=l;
}
}
if(lpivot==-1){/*printf("Error in get error2. Set error to 0\n");*/error=0;}
else{
m=(alpha[lpivot]-alpha[lpivot+1])/(chi2[lpivot]-chi2[lpivot+1]);
n=alpha[lpivot]-m*chi2[lpivot];
error=( (chi2min+1)*m+n   )-alpha[imin];
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

 
