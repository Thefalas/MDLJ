// aleat.cpp RETURNS A RANDOMLY GENERATED AND UNIFORMILY DISTRIBUTED 
// DOUBLE PRECISION REAL NUMBER "al" between 0 and 1 (both included)

// the seed is idum (integer), and the seed is initialized each time 
// iff (integer) is set to zero


#include "params.h"
//#include "def_functions.h"

//************************************************************
//*      PARAMETROS DE LA FUNCION DE NUMERO ALEATORIO
//************************************************************ 

//	parametros internos del generador
#define MBIG 1000000000
#define	MSEED 161803398
#define MZ 0
#define FAC 1.e-9f


int iff;
double al;


double aleat(int idum)
{

	static int ia,k,ma[56],ii,mj,mk,inext,inextp;
	ma[0]=0;


do{
if(idum < 0 || iff==0)
{
	iff=1;
	mj=MSEED*-abs(idum); 
	mj%=MBIG;
	ma[55]=mj;
	mk=1;
    
	for(ia=1;ia<=54;++ia)
	{
		ii=ia*21%55;
		ma[ii]=mk;
		mk=mj-mk;
		if(mk < MZ) mk=mk+MBIG;
		mj=ma[ii];
    }
	

	for(k=1;k<=4;++k)
	{	for(ia=1;ia<=55;++ia)
	 {
			ma[ia]=ma[ia]-ma[1+(ia+30)%55];
			if(ma[ia] < MZ) ma[ia]=ma[ia]+MBIG; 
	 }
	}   


	inext=0;
	inextp=31;
//	idum=1;

} // if IDUM


++inext;
if(inext==56) inext=1 ;
++inextp;
if(inextp==56) inextp=1;
mj=ma[inext]-ma[inextp];
if(mj < MZ) mj+=MBIG;
ma[inext]=mj;
al=mj*FAC;
if(al==0.) al+=0.00001;


}while(!al); // la condicion de bucle equivale a: 'while(al=0)'
//	FIN del bucle condicional de 'al' mayor que cero 


return al;


}

