/*
 *  main.c
 *  LJMD
 *
 *  Created by Francisco Vega Reyes on 20/07/09.
 *  Copyright 2009. All rights reserved.
 *
 */
 

#include "params.h"
//#include "def_functions.h"

int i, j, k, l,m,cfuerzas;   
int it,itr;
int xi;
int MoreCycles;
int iPC;

int *lista;

double temp,rij2,ljr,ir2,ir6;
double v2,m2;
double up, energy;
double t;
double virial;

double bin, normar, *histo, rangor,*rbin;
int limr,nbin, ibin;

double cred;

ivector nvcelda;

vector Lv ;

vector r[1+npart],ro[1+npart];
vector v[1+npart], vo[1+npart];
vector a[1+npart], a1[1+npart], a2[1+npart];

vector dr,drb;

rva *mol;
rvae *mole;

vector vijt,vijn,atot,vij;
vector awtot, aw[npart+1];
vector aNe[npart+1],aNi[npart+1],aT[npart+1];
vector rot;
double ir,rij,vr;
vector w[npart+1];


time_t hora0,horaF;


int main () {

	time(&hora0); // inicia el conteo de tiempo
	
	
	iff=0;		// inicializa el generador de numero aleatorio con semilla predeterminada LR (solo una vez en cada simulacion)
	// MUY IMPORTANTE: la linea de arriba NO SE TOCA.
	
	iPC=1; // iPC=1, Predictor-Corrector; iPC=0, LeapFrog
	
	it=0; // inicia el paso temporal
	SetupJob(1.);	// inicializa en estado de equilibrio con temperatura 1. 
	ComputeForcesCS(); // calcula las aceleraciones iniciales!!
	midetemp(v);// comprueba que la temperatura inicial coincide con la esperada
	salida(it);// salva el estado inicial

	
	mensajes0();	// $ imprime en pantalla informaci—n de la evoluci—n del sistema

	

	// Bucle Predictor-Corrector
	do{
		
		it++;
		PredictorStep();
		ApplyBoundaryCond();
		ComputeForcesCS();
		rdf(); //mide la funcion de distribucion radial
		corrFCC();
		CorrectorStep();
		ApplyBoundaryCond();
		
		// funcion de salida de informacion 
		if (it/(utermo*1.) == it/utermo){
			midetemp(v); // mide la temperatura del sistema
			midemom(v);	// mide el momento lineal total
			normar=VProd(Lv)/(2.*3.14159265*Cube(bin)*Sqr(npart)*utermo);
			for (i=1;i<nbin;i++) {histo[i]*=(normar/Sqr(rbin[i]));/*printf("%g \t %g\n", bin*rbin[i], histo[i])*/;}
			cred/=(utermo*npart);
			printf("correlaciones: %g\n",cred); 
			itr=1;
			mensajes();	// $ imprime en pantalla informaci—n de la evoluci—n del sistema
			salida(it/utermo); // $ salva el estado del sistema cada 'utermo' pasos
		}
		
		
	}while(it<nt);

	
/*
	// Bucle LeapFrog
	do{
	
		// funcion de evolucion temporal 
		it++;
		LeapfrogStep(1);
		ApplyBoundaryCond();
		ComputeForces();
		LeapfrogStep(2);
		
	
		// funcion de salida de informacion 
		if (it/(utermo*1.) == it/utermo){
			midetemp(v); // mide la temperatura del sistema
			midemom(v);	// mide el momento lineal total
			mensajes();	// $ imprime en pantalla informaci—n de la evoluci—n del sistema
			salida(it/utermo); // $ salva el estado del sistema cada 'utermo' pasos
		}
		
	
	}while(it<nt);
*/
	
	//system ("pause");  // esta instrucci—n no es C estandar. no funciona en todos los SO.

	
	fprintf(stdout,"\a\a\a\n");
	time(&horaF);
	fprintf(stdout,"Tiempo de ejecucion: %g minutos\n",(float) (horaF-hora0)/60.);	
	
    return 0;
	
}



