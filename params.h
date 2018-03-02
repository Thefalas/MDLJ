/*
 *  params.h
 *  LJMD
 *
 *  Created by Juanjo Meléndez on 02/02/09.
 *  Modified by Francisco Vega Reyes (main versions) on 20/07/09, 29/07/10.
 *  Copyright 2009. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <time.h>
#include "def_functions.h"

//	Definición de constantes y parámetros 

#define LR 102133*10		/* Semilla para el generador de numeros aleatorios */
#define r2 sqrt(2.)  // raiz de 2

#define ncx 4	//numero de celdas de cristal en una direccion
#define npart 4*ncx*ncx*ncx			// Numero de particulas. TIENE QUE SER PAR (compatibilidad con initCoord)  
							// o MULTIPLO de 4 (compatibilidad con InitCristal)
#define nc ncx*ncx*ncx // numero total de celdas de cristal

#define n 0.3// la densidad del gas. m‡ximo: 0.7404 (cristal fcc)

#define vp acos(-1.0)/6. // pi/6
#define Lc r2*pow(vp*r2/n, 1./3)  //parametro reticular, determinado por la densidad
#define dx r2*Lc-2.
#define L ncx*Lc	// la longitud total del sistema es el parametro reticular x numero de celdas

//#define npart 10*9
//#define L 2.75*4

#define sigma 1.				// diametro de las particulas. si =1, unidad de distancia
#define sigma2 sigma*sigma		// diametro cuadrado de las particulas

// parametros temporales de la simulacion
#define nt 2000000	// numero de pasos temporales 
#define dt 1.e-4
#define utermo 100000  //intervalo entre medidas

//magnitudes relacionadas con la fuerza LJ
#define fact pow(2., 1./6)
#define rc fact*sigma  /* longitud de corte de la cola atractiva */
#define rc2 rc*rc	   /* idem al cuadrado */

#define N_VECIN 14
#define VECINS {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},{-1,1,0},{0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1}, {-1,-1,1}, {0,-1,1},{1,-1,1}}


// magnitudes para el gas granular de esferas blandas
#define Y 10000000
#define gamman 200
#define gammashear 100
#define K 0.1 // se define de forma que I= K m sigma^2 


typedef struct {double x,y,z;} vector;
typedef struct {int x,y,z;} ivector;

typedef struct { vector r, rv, ra;} rva; // estructura molecular -sin uso de momento-
typedef struct { vector r, rv, ra, ra1, ra2, roo, rvo;} rvae; // idem pero extendida Pred-Correc

extern ivector nvcelda;

extern vector Lv,LvR; /* Para L = 10, el numero maximo de particulas es aproximadamente 700 */
extern vector r[1+npart];
extern vector ro[1+npart];

extern vector v[1+npart];
extern vector vo[1+npart];

extern vector a[1+npart];
extern vector a1[1+npart];
extern vector a2[1+npart];

extern vector dr,drb;

extern ivector icelda;

extern rva *mol;
extern rvae *mole;

extern int *lista;

extern int i, j, k, l,m,cfuerzas;
extern int xi;
extern int it,itr;
extern int MoreCycles;
extern int iPC;

//extern double L;

extern double al;				/* el numero aleatorio entre 0 y 1 devuelto por aleat es almacenado en esta variable */
extern int iff;				/* iff inicializa, si nulo, el generador de numero aleatorio con la semilla 'idum' */
extern double temp,rij2,ljr,ir2,ir6;
extern double v2,m2;
extern double up, energy;
extern double t;
extern double virial;

// variables para la funcion de distribucion radial
extern double bin, normar, *histo, rangor,*rbin;
extern int limr,nbin, ibin;
// variable de resultado de las correlaciones de largo alcance
extern double cred;


extern vector vijt,vijn,atot,vij;
extern vector awtot, aw[npart+1];
extern vector aNe[npart+1],aNi[npart+1],aT[npart+1];
extern vector rot;
extern double ir,rij,vr;
extern vector w[npart+1];


extern double aleat(int idum);	/* aleat: generador del numero aleatorio; idum: semilla */

extern int midetemp(vector rr[npart+1]);
extern int midemom(vector rr[npart+1]);
extern int salida(int numarchiv);
extern int mensajes(void);
extern int mensajes0(void);

extern void SetupJob();
extern void InitCoords();
extern void	InitVels();
extern void	InitForce();
extern void	InitPCAccels();

extern void SingleStep();
extern void LeapfrogStep(int ileap);
extern void ApplyBoundaryCond();
extern void ComputeForces();
extern void ComputeForcesCS();
extern void EvalProps();

extern void PredictorStep();
extern void CorrectorStep();

extern void RDF();

extern void AllocArrays();

extern int rdf(void);
extern void InitFCC(void);
extern int corrFCC(void);

