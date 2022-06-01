
/*
/*% Programa para decodificar Turbo Código, con el algoritmo MAP / BCJR, con Interleaver más grande
% Autores Matlab: Damian Gustavo Levin Jorge Castiñeira   pasado a C por luciana :)
% Creado el: 5/10/2003
% Última actualización: 16/12/2003        pasado a C 16/06/2010
%--------------------------------------------------------------------------------------------
 */

#define INFINITY (1.0/0.0)
#pragma hdrstop
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
//#include <condefs.h>
//#include <fstream.h>    //para abrir y cerrar archivos
#include <time.h> // para el calculo del tiempo empleado
#include <math.h>
#include "conio.h"
//#include <iostream.h>
#include "nr.h"// la pide NUMRECN. Es NR.H modificada para long double
#include "nrutil.h"
#include "nrutil.c"
#define NR_END 1        // la usan vector, freevector, matrix y freematrix
#define FREE_ARG char*  // puntero a una variable tipo char
#define log2(val) (log((val))/log(2.0))
//------------------------------------------------------------------------------
// DECLARACION DE ARCHIVOS DE DATOS/RESULTADOS
//------------------------------------------------------------------------------
 FILE *sali = fopen("Pbe_4saltos_D30_9.txt","w");
 //FILE *salu = fopen("TWBM_4saltos.txt","w");

//-------------------------Punteros a matrices y vectores----------------------
int Co,*COEF,**Tabla_elem_bin,**Key_Coef,*ii,*GF,*q,*InputElem,*TrellisInput,**Salida,*InputElem2,*Inter,*TrellisInput2,**Salida2,**Salidabin,**In,**LukyHard;
double **RecivedBit,**Y1,**Y2,**Y,**Luk ,**Alfa,*A,**Beta,*B,**Luky,*Numer,*Denom,Pbe,**Le;
//---------declaracion de funciones----------------------------------------------

void dec2bin(long decimal, int *binary);
void duracion(void);
float gasdev(long *);
float ran1(long *);
double TWBM4(void);
double TWBM(void);
double logistico(void);
int LFSR(void);

//----------------funciones-----------------------------------------------------
// duracion
/* Esta funci¢n realiza la impresi¢n de la duraci¢n de un programa. Es la de
plares modificada */

void duracion(void)
{

static int i=1;
static time_t tinicial, tfinal;         // variables tipo tiempo
int deltah, deltam, deltas, deltac; //para calcular horas, mins y segs.
long double Interval;

switch(i)
{  case 1:
      tinicial = time(NULL); i=2; return;

   case 2:
      tfinal = time(NULL);
      Interval = difftime(tfinal,tinicial);
      deltah = Interval/3600.;
      deltam = (Interval = Interval - deltah * 3600)/60.;
      deltas = Interval - deltam * 60;
  //    sali <<"\n"<< " CORRIO HOY: "<< ctime(&tinicial)
    //      <<" TIEMPO TOTAL DE CALCULO: "<<deltah<<" hs. "<<deltam<<" ms. "
      //    <<deltas<<" segs.\n";

}   // Fin del switch nro.
}   // Fin de la funcion duracion.
//-----------------------------------------------------

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
int LFSR() {
static int a7,a6,a5,a4,a3,a2,a1,a0;
int lfsr;//0xACE1u;
  /* taps: 16 14 13 11; characteristic polynomial: x^16 + x^14 + x^13 + x^11 + 1 */
 // bit  = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 4) ) & 1; /*[8 6 5 4 0]*/
 // lfsr =  (lfsr >> 1) | (bit << 15);
if (a7==NULL && a6==NULL && a5==NULL && a4==NULL && a3==NULL && a2==NULL && a1==NULL && a0==NULL ){
   a7=1;a6=1; a5=1;a4=1;a3=1;a2=1;a1=1;a0=1;
   // printf(" nuloo primera vez \n");
}
lfsr=a7*128+a6*64+a5*32+a4*16+a3*8+a2*4+a1*2+a0;
//printf("%d; ", lfsr);

int bit;
//unsigned period = 0;
 bit  = a7+a5+a4+a3;

 bit=bit % 2;

a7=a6;
a6=a5;
a5=a4;
 a4=a3;
 a3=a2;
 a2=a1;
 a1=a0;
 a0=bit;
 //a6=a7;
 //a7=bit;
//lfsr=a7*128+a6*64+a5*32+a4*16+a3*8+a2*4+a1*2+a0;

return lfsr;
}



float ran1(long *idum) {
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    temp=(float)AM*iy;
    if (temp > RNMX) return (float)RNMX;
    else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

double logistico() {
long utime;
utime=(long)time(NULL);					/*Unix time*/
long idum2=utime;
static double semilla1;
if (semilla1==NULL) semilla1=ran1(&idum2);
if (semilla1<=0) semilla1=ran1(&idum2);
if (semilla1>=0.999) semilla1=ran1(&idum2);
//printf("semilla1= %lf   ",semilla1);
double x0=semilla1;

double k=4*x0*(1-x0);
//printf("k= %lf \n  ",k);
semilla1=k;

return k;
}


double TWBM() {
long utime;
utime=(long)time(NULL);					/*Unix time*/
long idum2=utime;
static double semilla1;
if (semilla1==NULL) semilla1=ran1(&idum2);
if (semilla1<=0) semilla1=ran1(&idum2);
if (semilla1>=0.999) semilla1=ran1(&idum2);
//printf("semilla1= %lf   ",semilla1);
double x0=semilla1;

double k,x1;

       if (x0<=0.333333) k=3*x0;
      else if (x0<=0.66667) k=3*x0-1;
      else k=3*x0-2;

semilla1=k;

return k;
}


double TWBM4() {
long utime;
utime=(long)time(NULL);					/*Unix time*/
long idum2=utime;
static double semilla1;
if (semilla1==NULL) semilla1=ran1(&idum2);
if (semilla1<=0) semilla1=ran1(&idum2);
if (semilla1>=0.999) semilla1=ran1(&idum2);
//printf("semilla1= %lf   ",semilla1);
double x0=semilla1;

double k,x1,x2,x3;

       if (x0<=0.333333) x1=3*x0;
      else if (x0<=0.66667) x1=3*x0-1;
      else x1=3*x0-2;

//       %*************

      if (x1<=0.333333) x2=3*x1;
      else if (x1<=0.66667) x2=3*x1-1;
      else x2=3*x1-2;

//        %*************

      if (x2<=0.333333) x3=3*x2;
      else if (x2<=0.66667) x3=3*x2-1;
      else x3=3*x2-2;
  //     %*************

      if (x3<=0.333333) k=3*x3;
      else if (x3<=0.66667) k=3*x3-1;
      else k=3*x3-2;

//printf("k= %lf \n  ",k);
semilla1=k;

return k;
}


#include <math.h>

float gasdev(long *idum)
{
    float ran1(long *idum);
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;

    if (iset == 0) {
        do {
            v1=2.0*ran1(idum)-1.0;
            v2=2.0*ran1(idum)-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return (float)(v2*fac);
    } else {
        iset=0;
        return (float)gset;
    }
}

//---------accepts a decimal integer and returns a binary coded string
//
void dec2bin(long decimal, int *binary)   {
int k = 0, n = 0,iii=0;
int neg_flag = 0;
int remain;
int old_decimal; // for test
char temp[80];
// take care of negative input
if (decimal < 0)
{
decimal = -decimal;
neg_flag = 1;
}
do {
iii++;
old_decimal = decimal; // for test
remain = decimal % 2;
// whittle down the decimal number
decimal = decimal / 2;
// this is a test to show the action
//printf("%d/2 = %d remainder = %d\n", old_decimal, decimal, remain);
binary[iii] = remain;
// converts digit 0 or 1 to character '0' or '1'
temp[k++] = remain + '0';
} while (decimal > 0);
if (neg_flag)
temp[k++] = '-'; // add - sign
else
temp[k++] = ' '; // space
// reverse the spelling
//while (k >= 0)

//binary[n++] = temp[--k];
//binary[n-1] = 0; // end with NULL
}

//------------------------------------------------------------------------------
//Empieza el programa
//------------------------------------------------------------------------------
int main(int argc, char **argv)
{

long nume;


double randx,iF,NF,N2F,ale2;

long utime;
utime=(long)time(NULL);					/*Unix time*/
long idum=utime;
long idum2=utime;
long idum3=utime;
long idum5=utime;
utime=(long)time(NULL);					/*Unix time*/
//printf("idum= %ld  \n idum2= %ld \n idum3= %ld  idum5=%ld \n",idum,idum2,idum3,idum5);
//duracion();       //-cuenta cuanto tarda desde aca

int N=400;     //Longitud de Palabra

float gmax=50000,dd;


float N3F=pow(N,0.5);
int  N3=N3F;
Inter=ivector(1,N);


for (int i=0;i<=N3-1;i++ ) for (int j=1;j<=N3;j++)   Inter[N3*i+j]=(i+1)+N3*(j-1);    // Con esto creamos el vector de Interleave para cualquier N

//N2 es long de ent con igual coef, L cant de saltos
float N2=N/(0.2*N3F); //CANTIDAD DE SALTOS=200  //<==De ACA SALE L!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               //LONG. DE CADA SALTO=N/N2=2      //L=4 0.2; L=1 0.05;  L=2 0.1



int N1=2; //numero de bits del campo GF


float NL=1;
int Nmax=pow(2,N1)-1;
int Ngenerador=3; 			//Este numero resulta de hacer P(alfa)=0, siendo P(X) el ploinomio que genera el campo GF
						  //		%Por ejemplo en GF(4) P(X)=X^2+X+1, luego alfa^2=alfa+1, Ngenerador=3 (11)
                  	//Se ingresan los coeficientes del polinomio de encriptamiento, 15 es coeficiente cero
						  //		%Generamos la tabla del Campo de Galois
int DosAlaN1=pow(2,N1);
int M=DosAlaN1*DosAlaN1;
COEF=ivector(1,N+1);
GF=ivector(1,DosAlaN1);
A=dvector(1,M);
Beta=dmatrix(1,N+1,1,M);
B=dvector(1,M);
Numer=dvector(1,M);
Denom=dvector(1,M);
InputElem=ivector(1,N);
InputElem2=ivector(1,N);
Salida=imatrix(1,N,1,2);
Salidabin=imatrix(1,N,1,N1*2);
Salida2=imatrix(1,N,1,2);
TrellisInput2=ivector(1,N+1);
TrellisInput=ivector(1,N+1);
RecivedBit=dmatrix(1,N,1,2*N1);
Y1=dmatrix(1,N,1,2*N1);
Y2=dmatrix(1,N,1,2*N1);
Y=dmatrix(1,N,1,2*N1);
Luk =dmatrix(1,N,1,N1);
Alfa=dmatrix(1,N,1,M);
Luky=dmatrix(1,N,1,N1);
LukyHard=imatrix(1,N,1,N1);
In=imatrix(1,N,1,N1);
Le=dmatrix(1,N,1,N1);

//parte del programa para el calculo de la matriz del trellis de codificadores GF

GF[1]=pow(2,(N1-1));
for (int i=1;i<=(Nmax-1);i++){

   if (GF[i]%2==0) GF[i+1]=GF[i]/2;
   else  GF[i+1]=((GF[i]-1)/2)^Ngenerador;


}
GF[DosAlaN1]=0; //el elemento 0 se guarda en la ultima posicion de la tabla GF


//Generacion de la Tabla de conversion de elementos del campo a su forma +-1
int tamM=DosAlaN1+5;

Tabla_elem_bin=imatrix(1,tamM+1,1,N1);   //Tabla_elem_bin contiene la conversion de elemento a palabra binaria del campo GF

for (int iu=1;iu<=tamM+1;iu++)  for (int ie=1;ie<=N1;ie++)   Tabla_elem_bin[iu][ie]=0;

long decimal;
q=ivector(1,N1);
for (int w=1;w<=N1;w++) q[w]=0;
for (int jj=1;jj<=DosAlaN1;jj++) {
dec2bin(jj-1,q);
 for (int ir=1;ir<=N1;ir++){
	   if  (q[N1+1-ir]==1)
   	   Tabla_elem_bin[jj][ir]=1;
   	else
      	Tabla_elem_bin[jj][ir]=-1;

   }

 }
free_ivector(q,1,N1);

Tabla_elem_bin[tamM+1][1]=2;
Tabla_elem_bin[tamM+1][2]=2;

//for(int g=1;g<=10;g++) {cout<<Tabla_elem_bin [g][1]<<" "<<Tabla_elem_bin [g][2]<<endl;}


                           //     16        16      3  200
int ****Mat;
int mat1=M+1,mat2=2+1,mat3=6+1;

Mat = (int ****)malloc (mat1 * sizeof(int ***));

for (int i=1; i<mat1; i++) {

                         Mat[i] = (int ***)malloc (mat1*sizeof(int**));
                                    for (int j=1; j<mat1; j++) {
                                                        Mat[i][j] = (int **)malloc (mat2*sizeof(int*));
                                                            for (int k=1; k<mat2; k++)
                                                            Mat[i][j][k] = (int *)malloc (mat3*sizeof(int)); }
                         }

int valou=DosAlaN1+5;
for (int a1=1;a1<=M;a1++) for (int b1=1;b1<=M;b1++) for (int c1=1;c1<=2;c1++) for (int d1=1;d1<=6;d1++) Mat[a1][b1][c1][d1]=valou;    //Mat(:,:,:)=2^N1+5;

Key_Coef=imatrix(1,6,1,3);

Key_Coef[1][1]=0;Key_Coef[1][2]=1;Key_Coef[1][3]=3;
Key_Coef[2][1]=0;Key_Coef[2][2]=2;Key_Coef[2][3]=3;
Key_Coef[3][1]=1;Key_Coef[3][2]=1;Key_Coef[3][3]=3;
Key_Coef[4][1]=1;Key_Coef[4][2]=2;Key_Coef[4][3]=3;
Key_Coef[5][1]=2;Key_Coef[5][2]=1;Key_Coef[5][3]=3;
Key_Coef[6][1]=2;Key_Coef[6][2]=2;Key_Coef[6][3]=3;

//for (int fil=1;fil<=6;fil++) {for (int col=1;col<=3;col++) printf("Key_Coef[%d][%d]=%d ",fil,col,Key_Coef[fil][col]);printf("\n");}
//getch();
int exp0,exp1,exp2;

for (int p=1;p<=6;p++){//------------------------------------barro los 6 coef optimos------------

exp0=Key_Coef[p][1];
exp1=Key_Coef[p][2];
exp2=Key_Coef[p][3];

int e=0;    //esta variable no necesita ser un arreglo (NL es 1)
int ed=0;
int edgf=0;
int edgfd=0;
int d=0,ExpEd,Exp1,Exp2,Exp3;

for (int ci1=0;ci1<=Nmax;ci1++){
   for (int ci2=0;ci2<=Nmax;ci2++){   // se corren todos los posibles estados, y se anota entrada-salida y proximo estado en la matriz Mat
      for (int h=1;h<=DosAlaN1;h++){ //h es la entrada, equivale a u en el esquema

        ed=ci1;
        edgfd=ci2;
  		d= ed^edgfd;
        e= (h-1)^d;

   if (ed==0)  edgf = GF[exp0+1];   //efecto del polinomio

  	else  {
     for (int m=1;m<=Nmax;m++) if (GF[m]==ed) ExpEd=m-1;

            Exp1=exp0;
            if (exp1==Nmax) Exp2=exp1;
            else   Exp2=(exp1+ExpEd)%Nmax;

            if (exp2==Nmax)  Exp3=exp2;
            else Exp3=(exp2+2*ExpEd)%Nmax;

        edgf = ((GF[Exp1+1]^GF[Exp2+1])^GF[Exp3+1]);

       }
     Mat[DosAlaN1*ci1+ci2+1][DosAlaN1*e+edgf+1][1][p]=h-1;
     Mat[DosAlaN1*ci1+ci2+1][DosAlaN1*e+edgf+1][2][p]=e;
  // printf("estado inicial:%d  entrada:%d  estado futuro:%d   salida:%d Coef:%d, \n",DosAlaN1*ci1+ci2+1,h-1,DosAlaN1*e+edgf+1,e,Coef);
   //getch();

 }}}               //16    16     3   200

 }    // for (int p=1;p<=N2;p++)  //fin parte del programa para el calculo de la matriz del trellis de codificadores GF

int num_it=1;         //<=========================================================================PASADAS!!!!!!
double Sigma,Lc,Errores;
double ***Gama;
int Nmas1=N+1,Mmas1=M+1;               //Gama(400,16,16)=Gama(N,M,M)
Gama = (double ***)malloc (Nmas1 * sizeof(double **));

for (int i=1; i<Nmas1; i++) {

                         Gama[i] = (double **)malloc (Mmas1*sizeof(double*));
                         for (int j=1; j<Mmas1; j++) {  Gama[i][j] = (double *)malloc (Mmas1*sizeof(double)); }
                         }




for (float EbNo=1.5;EbNo<=3.5;EbNo=EbNo+0.5){  //=================Este lazo recorre valores de Sigma=======================================
//printf("EbNo:%f  \n ",EbNo);
Sigma = pow(10,(-1*EbNo/20));        // Dispersión del Ruido

Lc = 2/pow(Sigma,2);               // Lc=4.a.Eb/(2.Sigma^2),   con a=1  y   Eb=1
Errores=0;


for (int g=1;g<=gmax;g++){    //Este bucle corre varias veces el programa completo*********************************************************
//printf("*** g:%d \n",g);

//for (int ree=0;ree<10000;ree++){
int pp=0,ale,ala,alu;
int Tr_tot=0;
int Tr=0;
int Nr=0;

for (int h=1;Tr_tot<N;h++){
ale=LFSR(); //sorteo la long del tramo
ala=ale % 3; //sortea cuanto tiempo mentiene el salto
Nr=ala+1;
ale=LFSR();
Tr=ala;
for(int barr=1;barr<=Nr;barr++){

if(Tr_tot+Tr>N) {Tr=N-Tr_tot;}
Tr_tot=Tr_tot+Tr;
alu=ale % 6;
Co=alu;
for (int ii=1;ii<=Tr;ii++){
pp=pp+1;
COEF[pp]=Co+1;
}
}
}

//for (int ii=1;ii<=N;ii++) fprintf("%d ",COEF[ii]);
//fprintf("\n");

//}
//for (int ii=1;ii<=N;ii++) COEF[ii]=1;
/*ale=ran1(&idum5);
for (int ii=1;ii<=100;ii++) COEF[ii]=floor(6*ale+1);
ale=ran1(&idum5);
for (int ii=101;ii<=200;ii++) COEF[ii]=floor(6*ale+1);
ale=ran1(&idum5);
for (int ii=201;ii<=300;ii++) COEF[ii]=floor(6*ale+1);
ale=ran1(&idum5);
for (int ii=301;ii<=400;ii++) COEF[ii]=floor(6*ale+1);
//printf("idum6=%dl \n",idum6);*/


//for (int ii=1;ii<=N;ii++) COEF[ii]=6;
//printf("Tr=%d ==> Tr_tot=%d \n",Tr,Tr_tot);

//getch();

 for (int ey=1;ey<=N-2;ey++){ale2=ran1(&idum3);  InputElem[ey] =floor(DosAlaN1*ale2); }//valores aleatorios de la entrada al codificador 1
//for (int ey=1;ey<=N-2;ey++){ infile3>>ale2;  InputElem[ey] =INT(DosAlaN1*ale2);  }//valores aleatorios de la entrada al codificador 1

   TrellisInput[1]=1;           // Armamos el Trellis del 1er Codificador hasta la posición N-2
   int cailDa;

for (int i=1;i<=N-2;i++){              // Los últimos dos elementos deben ser tales que terminen el Trellis en el Estado 1
    iF=i;
    NF=N;
    N2F=N2;
         cailDa=ceil(iF/(NF/N2F));

         for (int j= 1;j<=DosAlaN1*DosAlaN1;j++){
              //busca para que estado futuro (j) la entrada
   	     if (Mat[TrellisInput[i]][j][1][COEF[i]]==InputElem[i]){    //if Mat(TrellisInput(i),j,1,ceil(i/(N/N2)))==InputElem(i);
              TrellisInput[i+1]=j;
            		Salida[i][1]= Mat[TrellisInput[i]][j][1][COEF[i]];
          		Salida[i][2]= Mat[TrellisInput[i]][j][2][COEF[i]];
          		//fprintf(sali,"Salida[%d][1]=%d Salida[%d][2]=%d \n",i,Salida[i][1],i,Salida[i][2]);

        	}
   	}   // for (int j= 1;j<=DosAlaN1*DosAlaN1;j++)
    }  //for (int i=1;i<=N-2;i++)

//Ahora forzamos los últimos dos elementos de entrada al Primer Codificador

  int Daceil;
	  for (int j=1;j<=DosAlaN1*DosAlaN1;j++){
    NF=N;
    N2F=N2;
      Daceil=ceil((NF-1)/(NF/N2F));

        if ((Mat[TrellisInput[N-1]][j][1][COEF[N]]!=DosAlaN1+5) && (Mat[j][1][1][COEF[N]]!=DosAlaN1+5))  {
          TrellisInput[N]=j;
          InputElem[N-1]=Mat[TrellisInput[N-1]][j][1][COEF[N]];    //ceil((N-1)/(N/N2)));
          TrellisInput[N+1]=1;
          InputElem[N]=Mat[j][1][1][COEF[N]];   //ceil((N-1)/(N/N2)));
   	  }   // if ((Mat[TrellisInput[N-1]-1][j-1][0][Daceil]~=DosAlaN1+5) && (Mat[j-1][0][0][Daceil]~=DosAlaN1+5))
     }   //for (int j=1;j<=DosAlaN1*DosAlaN1;j++)

//   save InputElem.dat InputElem /ascii;

	for (int i= N-1;i<=N;i++){            // Hacemos el Trellis de los últimos dos Elem
		for (int j= 1;j<=DosAlaN1*DosAlaN1;j++){

   		if (Mat[TrellisInput[i]][j][1][COEF[i]]==InputElem[i]){
         	TrellisInput[i+1]=j;

         	Salida[i][1]= Mat[TrellisInput[i]][j][1][COEF[i]];
           	Salida[i][2]= Mat[TrellisInput[i]][j][2][COEF[i]];
         //  	fprintf(sali,"Salida[%d][1]=%d Salida[%d][2]=%d \n",i,Salida[i][1],i,Salida[i][2]);
      	}
   	}  //	for (int j= 1;j<=DosAlaN1*DosAlaN1;j++)
	}    //for (int i= N-1;i<=N;i++)

//   save Salida.dat Salida /ascii;


    int ceilDA;
	for (int i= 1;i<=N;i++) InputElem2[i]= InputElem[Inter[i]];    // Interleave de la Secuencia Recibida para el Segundo Decodificador

TrellisInput2[1]=1;
	for (int i= 1;i<=(N-1);i++){
    iF=i;
    NF=N;
    N2F=N2;
    ceilDA=ceil(iF/(NF/N2F));
		for (int j= 1;j<=DosAlaN1*DosAlaN1;j++){       //if Mat(TrellisInput2(i),j,1,ceil(i/(N/N2)))==InputElem2(i);
   		if (Mat[TrellisInput2[i]][j][1][COEF[i]]==InputElem2[i]){
         	TrellisInput2[i+1]=j;
          	Salida2[i][1]= Mat[TrellisInput2[i]][j][1][COEF[i]];
         	Salida2[i][2]= Mat[TrellisInput2[i]][j][2][COEF[i]];
      	}
   	}     //for (int j= 1;j<=DoaAlaN1*DosAlaN1;j++)
   } //for (int i= 1;i<=(N-1);i++)

   for (int j= 1;j<=DosAlaN1*DosAlaN1;j++){
   		if (Mat[TrellisInput2[N]][j][1][COEF[N]]==InputElem2[N]){  //ceil((N-1)/(N/N2)))    Daceil=ceil((N-1)/(N/N2))-1;
         	TrellisInput2[N+1]=j;
         	Salida2[N][1]= Mat[TrellisInput2[N]][j][1][COEF[N]];
         	Salida2[N][2]= Mat[TrellisInput2[N]][j][2][COEF[N]];
      	}
   } //for (int j= 1;j<=DosAlaN1*DosAlaN1;j++)

    for (int i= 2;i<=N;i=i+2) Salida[i][2]=Salida2[i][2];                  // Ahora ponemos los Elementos de Paridad pares, del Segundo Codificador
 //  save Salida2.dat Salida2 /ascii;

 for (int j=1;j<=N1;j++) for (int i= 1;i<=N;i++) {/*cout<<"Salida["<<i<<"][1]+1:"<<Salida[i][1]+1<<endl; getch(); */ Salidabin[i][j]=Tabla_elem_bin[Salida[i][1]+1][j];}	// Ahora se arma la salida en formato binario con los elementos de redundancia de cada codificador

for (int j=N1+1;j<=2*N1;j++) for (int i= 1;i<=N;i++) Salidabin[i][j]=Tabla_elem_bin[Salida[i][2]+1][j-N1];		// Ahora se arma la salida en formato binario con los elementos de redundancia de cada codificador

//   save Salidabin.dat Salidabin /ascii;

	//Ahora le sumamos ruido a los Elementos (bits) transmitidos
for(int i=1;i<=N;i++) for(int j=1;j<=2*N1;j++) {Y1[i][j]=0; Y2[i][j]=0;}

for(int i=1;i<=N;i++){ for(int j=1;j<=2*N1;j++) {randx=gasdev(&idum); RecivedBit[i][j] =Salidabin[i][j]+ Sigma*randx; Y1[i][j]=RecivedBit[i][j];}  }//ver distribucion normal!!!!!!!!!!!!!!!!!!

for (int i= 2;i<=N;i=i+2)  for (int j=N1+1;j<=2*N1;j++)  Y1[i][j]=0;    // Ponemos ceros en los bit de paridad pares (puncturing) para el primer Decodificador

   // Ahora vamos a generar las secuencias para cada Decodificador
   for (int i= 1;i<=N;i++) for (int j=1;j<=N1;j++)  Y2[i][j]= RecivedBit[Inter[i]][j];      // Interleave de la Secuencia Recivida para el Segundo Decodificador

   for (int i=2;i<=N;i=i+2)  for (int j=N1+1;j<=2*N1;j++)  	Y2[i][j]=RecivedBit[i][j];

	// Ahora comienza el proceso de Decodificación

for (int hh=1;hh<=N;hh++) for (int jj=1;jj<=N1;jj++) Luk[hh][jj]=0;
for (int i=1; i<Nmas1; i++) for (int j=1; j<Mmas1; j++) for (int k=1; k<Mmas1; k++) Gama[i][j][k]=-INFINITY;

//************************************************************** Cada bucle es un Decodificador ==> Cantidad de iteraciones = la mitad del h máximo
for (int pasada=1;pasada<=15;pasada++){
if (pasada%2==1)
{for (int i=1;i<=N;i++) for (int j=1;j<=N1*2;j++) Y[i][j]=Y1[i][j]; }
else
{for (int i= 1;i<=N;i++) for (int j=1;j<=N1*2;j++)   Y[i][j]=Y2[i][j]; }

 double TEB;

int ceilDDA;
// ******************************************************************************************Cálculo de los GAMAS
for (int i= 1;i<=N-1;i++){
    iF=i;
    NF=N;
    N2F=N2;
      ceilDDA=ceil(iF/(NF/N2F));

     	for (int k= 1;k<=M;k++){
           for (int j= 1;j<=M;j++){

               if (Mat[j][k][1][COEF[i]] != DosAlaN1+5){     //ceil(i/(N/N2)))~=2^N1+5;
                    Gama[i][j][k]=0;                                            // Tabla_elem_bin(Mat(j,k,2,ceil(i/(N/N2)))+1,n);
                   for (int n=1;n<=N1;n++){TEB=Tabla_elem_bin[Mat[j][k][1][COEF[i]]+1][n];  Gama[i][j][k]=Gama[i][j][k]+Y[i][n]*TEB;      }//el calculo de los Gama cambia para la generalizacion, revisar

                 for (int n=N1+1;n<=2*N1;n++){TEB=Tabla_elem_bin[Mat[j][k][2][COEF[i]]+1][n-N1]; Gama[i][j][k]=Gama[i][j][k]+Y[i][n]*TEB; }//el calculo de los Gama cambia para la generalizacion, revisar

              	Gama[i][j][k]=Lc/2*Gama[i][j][k];                //Tabla_elem_bin(Mat(j,k,1,ceil(i/(N/N2)))+1,n)
            for (int n=1;n<=N1;n++) {TEB=Tabla_elem_bin[ Mat[j][k][1][COEF[i]]+1 ][n]; Gama[i][j][k]=Gama[i][j][k]+TEB*Luk[i][n]/2; /*if (Luk[i][n]!=0){cout<<"h:"<<h<<" ,Luk["<<i<<"]["<<n<<"]:"<<Luk[i][n]<<" , Luk[i][n]/2:"<<Luk[i][n]/2<<endl;getch();}       */ }

        }   // if (Mat[j][k][1][ceilDDA] != DosAlaN1+5){


        }    // for (int j= 1;j<=DosAlaN1*DosAlaN1;j++){
        }  //for (int k= 1;k<=DosAlaN1*DosAlaN1;k++){
     }   //    for (int i= 1;i<=N-1;i++){  // Cálculo de los GAMAS


                  //M=(2^N1)*(2^N1)
    for (int k= 1;k<=M;k++){
            for (int j= 1;j<=M;j++){                                   //Daceil=ceil((NF-1)/(NF/N2F));

               if (Mat[j][k][1][COEF[N]]!=DosAlaN1+5){              //      Daceil=ceil((N-1)/(N/N2));

                  Gama[N][j][k]=0;
                  for (int n=1;n<=N1;n++){TEB=Tabla_elem_bin[Mat[j][k][1][COEF[N]]+1][n];  Gama[N][j][k]=Gama[N][j][k]+Y[N][n]*TEB;}    //el calculo de los Gama cambia para la generalizacion, revisar[

               for (int n=N1+1;n<=2*N1;n++){TEB=Tabla_elem_bin[Mat[j][k][2][COEF[N]]+1][n-N1];  Gama[N][j][k]=Gama[N][j][k]+Y[N][n]*TEB;}   //el calculo de los Gama cambia para la generalizacion, revisar

						Gama[N][j][k]=Lc/2*Gama[N][j][k];
                  for (int n=1;n<=N1;n++){TEB=Tabla_elem_bin[Mat[j][k][1][COEF[N]]+1][n]; Gama[N][j][k]=Gama[N][j][k]+TEB*Luk[N][n]/2; }


            } //          if (Mat[j-1][k-1][0][Daceil]!=DosAlaN1+5){
         } // for (int j= 1;j<=M;j++){
      }    //for (int k= 1;k<=M;k++){

                //M=(2^N1)*(2^N1)
//-------------------------------------------------------------------------------------// Cálculo de los ALFAS
for(int aa=1;aa<=N;aa++) for(int bb=1;bb<=M;bb++) Alfa[aa][bb]=0;   //Alfa(1,:)=-Inf;
for(int bb=1;bb<=M;bb++) Alfa[1][bb]=-INFINITY;   //Alfa(1,:)=-Inf;


      Alfa[1][1]=0;
double Pru;
   	for (int i=2;i<=N;i++){
      	for (int k=1;k<=M;k++){
         	for (int j=1;j<=M;j++)  A[j]= Gama[i-1][j][k] + Alfa[i-1][j];

            if (A[1]==-INFINITY && A[2]==-INFINITY)  Pru=-INFINITY;
            else  {
            if (A[1]>=A[2])   Pru= A[1] + log(1+exp(-1*fabs(A[1]-A[2])));
            else            Pru= A[2] + log(1+exp(-1*fabs(A[1]-A[2])));
            }

				for (int j=3;j<=M;j++)
                {
                if (A[j]!=-INFINITY) {
                 if (Pru>=A[j]) Pru= Pru + log(1+exp(-1*fabs(Pru-A[j])));
                 else Pru= A[j]+log(1+exp(-1*fabs(Pru-A[j])));
                }
                }

            Alfa[i][k]= Pru;


      	} //	for (int k=1;k<=M;k++){
   }  // 	for (int i=2;i<=N;i++){
//------------------------------------------------------------------------------// Cálculo de los BETAS
for(int cc=1;cc<=N+1;cc++) for(int bb=1;bb<=M;bb++) Beta[cc][bb]=-INFINITY;   //   	Beta(:,:,:)=-Inf;

   if ((pasada%2)==1)  Beta[N+1][1]=0;    // Si es el Primer Decodificador, el último Beta(1) debe ser 1 y los demás deben ser 0
   else   for(int cc=1;cc<=M;cc++) Beta[N+1][cc]=0;    // Si es el Segundo Decodificador, todos los últimos Beta debe ser 1
                         // porque el segundo Trellis no necesariamente termina en el Estado 1 (por el Interleaver)

for (int i=N;i>=2;i=i-1){
      	for (int k=1;k<=M;k++){
         	for (int j=1;j<=M;j++) B[j]= Beta[i+1][j] + Gama[i][k][j];

            if (B[1]==-INFINITY && B[2]==-INFINITY)  Pru=-INFINITY;
                            //Pru=max(B(1),B(2)) + log(1+exp(-1*abs(B(1)-B(2))));
            else { if (B[1]>=B[2]) Pru=B[1] + log(1+exp(-1*fabs(B[1]-B[2])));  else  Pru=B[2] + log(1+exp(-1*fabs(B[1]-B[2]))); }

            for (int j=3;j<=M;j++)
            if (B[j]!=-INFINITY) {
                         if (Pru>=B[j]) Pru= Pru + log(1+exp(-1*fabs(Pru-B[j])));
                         else  Pru= B[j] + log(1+exp(-1*fabs(Pru-B[j])));
                                 }

            Beta[i][k]= Pru;
      	}  //	for (int k=1;k<=M;k++){
   	}    //for (int i=N;i>=2;i=i-1){
      //save Alfa.dat Alfa /ascii;
      //%save Beta.dat Beta /ascii;
//*****************************************************************************************


for(int lif=1;lif<=N;lif++) for(int lic=1;lic<=N1;lic++)Luky[lif][lic]=0;

double  NumerAux, DenomAux;
int CCEEil;

for (int n=1;n<=N1;n++){
   for (int i=1;i<=N-1;i++){
    iF=i;
    NF=N;
    N2F=N2;
      CCEEil=ceil(iF/(NF/N2F));

      NumerAux=-INFINITY;
      DenomAux=-INFINITY;
      for (int limpia=1;limpia<=M;limpia++) { Numer[limpia]=0;   Denom[limpia]=0;}

         for (int k=1;k<=M;k++){
				for (int j=1;j<=M;j++){
                                               // Tabla_elem_bin(Mat(j,k,1,ceil(i/(N/N2)))+1,n)==1      ceil(i/(N/N2)))
               if (Tabla_elem_bin[Mat[j][k][1][COEF[i]]+1][n]==1)
               {

                  Numer[k]=Alfa[i][j] + Gama[i][j][k] + Beta[i+1][k];    //	NumerAux = max(NumerAux,  Numer(1,k)+log(1+exp(-1*abs(NumerAux-Numer(1,k)))));
                  if (Numer[k]!=-INFINITY)
                            if (NumerAux<(Numer[k]+log(1+exp(-1*fabs(NumerAux-Numer[k])))))
                            NumerAux = Numer[k]+log(1+exp(-1*fabs(NumerAux-Numer[k])));

              }
              if (Tabla_elem_bin[Mat[j][k][1][COEF[i]]+1][n]==-1) {
                  Denom[k] = Alfa[i][j] + Gama[i][j][k] + Beta[i+1][k];
                  if (Denom[k]!=-INFINITY)
                            if (DenomAux<Denom[k]+log(1+exp(-1*fabs(DenomAux-Denom[k]))))
                            DenomAux = Denom[k]+log(1+exp(-1*fabs(DenomAux-Denom[k])));
                        }  // } else if (Tabla_elem_bin[Mat[j-1][k-1][0][CCEEil]+1][n]==-1) {

               	}  //for (int j=1<=M;j++){
      	}    //for (int k=1;k<=M;k++){



         Luky[i][n]=NumerAux - DenomAux;

      } //   for (int i=1;i<=N-1;i++){


  } //for (int n=1;n<=N1;n++){

    NF=N;
    N2F=N2;
int cedil=ceil((NF-1)/(NF/N2F));
  for (int n=1;n<=N1;n++){

      NumerAux=-INFINITY;
      DenomAux=-INFINITY;
      for (int limpia=1;limpia<=M;limpia++) { Numer[limpia]=0;   Denom[limpia]=0;}

         for (int k=1;k<=M;k++){
				for (int j=1;j<=M;j++){

               if (Tabla_elem_bin[Mat[j][k][1][COEF[N]]+1][n]==1){
                  Numer[k]=Alfa[N][j] + Gama[N][j][k] + Beta[N+1][k];
                  if (Numer[k]!=-INFINITY)
                        if (NumerAux<Numer[k]+log(1+exp(-1*fabs(NumerAux-Numer[k]))))
                        NumerAux =Numer[k]+log(1+exp(-1*fabs(NumerAux-Numer[k])));
               	}
               if (Tabla_elem_bin[Mat[j][k][1][COEF[N]]+1][n]==-1){
                  Denom[k] = Alfa[N][j] + Gama[N][j][k] + Beta[N+1][k];
                  if (Denom[k]!=-INFINITY)
                        if(DenomAux<Denom[k]+log(1+exp(-1*fabs(DenomAux-Denom[k]))))
                        DenomAux = Denom[k]+log(1+exp(-1*fabs(DenomAux-Denom[k])));

            	 } //  if (Tabla_elem_bin[Mat[j-1[k-1][0][cedil]+1][n]==-1){
         }  //	for (int j=1;j<=M;j++){


      	}  //    for (int k=1;k<=M;k++){

         Luky[N][n]=NumerAux - DenomAux;
 }  // for (int n=1;n<=N1;n++){

   //save numer.dat Numer /ascii;
	//save denom.dat Denom /ascii;
    //save Luky.dat Luky /ascii;

for (int j=1;j<=N;j++) for (int f=1;f<=N1;f++) if (Luky[j][f]<=0) LukyHard[j][f]=-1; else LukyHard[j][f]=1;   //Para graficar el Trellis Decodificado
    //  save LukyHard.dat LukyHard /ascii;

for (int j=1;j<=N1;j++) for (int fil=1;fil<=N;fil++) In[fil][j]=Tabla_elem_bin[InputElem[fil]+1][j];	//In contiene los bits de entrada correpondientes a la secuencia aleatoria
         									//Sirve para comparar lo detectado con lo enviado
   //   save In.dat In /ascii  se guarda la secuencia de entrada en formato binario
   //Ahora la Información Extrínseca para pasársela al otro Decodificador

for (int i= 1;i<=N;i++) for (int j=1;j<=N1;j++) { Le[i][j] = Luky[i][j] - Lc*Y[i][j] - Luk[i][j];}

for (int j=1;j<=N1;j++)	for (int i= 1;i<=N;i++) Luk[i][j]= Le[Inter[i]][j];   // De-interleave de la Información Extrínseca que entra
                                 //como Información a Priori en el otro Decodificador

}   //for (int h=1;h<=num_it;h++) // Cada bucle es un Decodificador ==> Cantidad de iteraciones = la mitad del h máximo
//********************************** % Fin del bucle de las iteraciones de la Decodificación ****************************
float ResTEMP=0;

for (int j=1;j<=N1;j++) {

ResTEMP=0;

for (int jj=1;jj<=N;jj++)ResTEMP=ResTEMP+fabs(LukyHard[jj][j]-In[jj][j])/2;
Errores= Errores + ResTEMP;     //esta bien esto? estaba en el prog el comentario no es de luciana
}

}  //for (int g=1;g<=gmax;g++)

Pbe=Errores/(N1*N*gmax);
printf("%lf %lf %lf\n",EbNo,Pbe,Errores);

fprintf(sali,"%lf %lf %lf\n",EbNo,Pbe,Errores);
}  //for (float EbNo=1.0;EbNo<=3.5;EbNo=EbNo+0.5)

//for (int a1=1; a1<=M; a1++)  for (int b1=1; b1<=M; b1++)for (int c1=1; c1<=3; c1++) for (int d1=1; d1<=N2; d1++) sal1<<Mat[a1][b1][c1][d1]<<" ";
for (int i=1; i<Nmas1; i++)  for (int j=1; j<Mmas1; j++)    free (Gama[i][j]);
free(Gama);
for (int i=1; i<mat1; i++)  for (int j=1; j<mat1; j++) for (int k=1; k<mat2; k++)    free (Mat[i][j][k]);
free(Mat);

//duracion();

//------------------------------------------------------------------------------
free_ivector(COEF,1,N+1);
free_ivector(Inter,1,N);             //N3=pow(N,0.5)=20    TamA=N3*(N3-1)+N3=400;
free_ivector(TrellisInput,1,N+1);    //N=400
free_ivector(InputElem,1,N);
free_ivector(InputElem2,1,N);
free_dvector(A,1,M);                //M=2^N1*2^N1=4*4=16          N1=2 numeromde bits del campo de GF
free_dmatrix(Y,1,N,1,2*N1);
free_dmatrix(Y1,1,N,1,2*N1);
free_dmatrix(Y2,1,N,1,2*N1);
free_dmatrix(Luky,1,N,1,N1);
free_dmatrix(RecivedBit,1,N,1,2*N1);
free_ivector(GF,1,DosAlaN1);             //DosAlaN1=2^N1=2^2=4
free_dvector(B,1,M);
free_ivector(TrellisInput2,1,N+1);
free_imatrix(Salidabin,1,N,1,2*N1);
free_imatrix(Salida2,1,N,1,2);
free_dvector(Numer,1,M);
free_imatrix(In,1,N,1,N1);
free_dvector(Denom,1,M);
free_dmatrix(Le,1,N,1,N1);
free_imatrix(LukyHard,1,N,1,N1);
free_imatrix(Salida,1,N,1,2);
free_imatrix(Key_Coef,1,6,1,3);
free_dmatrix(Luk,1,N,1,N1);
free_dmatrix(Beta,1,N+1,1,M);
free_dmatrix(Alfa,1,N,1,M);
free_imatrix(Tabla_elem_bin,1,tamM+1,1,N1);
return 0;
}







