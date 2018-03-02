/*
 *  SetupJob.c
 *  LJMD
 *
 *  Created by Juanjo Meléndez on 02/02/09.
 *  Modified by Francisco Vega Reyes on 20/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#include "params.h"
//#include "def_functions.h"

void AllocArrays(){
	AllocMem(lista,VProd(nvcelda)+npart,int);
	AllocMem(histo,nbin,double);
	AllocMem(rbin,nbin,double);
}


void SetupJob (double norma){
	it = 0;
	itr=1;
	nbin=500; //numero de bins para la lectura de la rdf
	
	VSetAll(Lv,L); // tamano 3D del sistema
	

	nvcelda.x=Lv.x/rc;  //numero de celdas para calculo de fuerzas
	nvcelda.y=Lv.y/rc;
	nvcelda.z=Lv.z/rc;
	
	AllocArrays();
	
	for (i=0; i<nbin; i++) {histo[i]=0;rbin[i]=i-0.5;}

	cred=0.; //inicia correlaciones FCC de largo alcance
	
	InitFCC();
	InitVels(norma);
	if (iPC) InitPCAccels();
	//InitForce();
}

void InitCoords(){		
	
	int fl,counter;
	
	fl=0;counter=0;
	for (i=1;i<=npart;i++) {

do{ // i. mientras 
		fl=0;
		VVSet(r[i],(aleat(LR)-0.5)*Lv);	// ii. la particula situada
	
	
		for(j=1; j<i; j++){ // iii. solape con las anteriores

				dr.x = r[i].x - r[j].x;dr.y = r[i].y - r[j].y;dr.z = r[i].z - r[j].z;
				if (VLengSq(dr) < sigma2)  {fl=5;} // sigma es el diametro de las particulas!!
		}		
		//printf("counter: %d\n",counter);
			
	}	while (fl > 1);	// iv. sigue generando posiciones aleatorias para esa particula
	
	//printf ("Particula no. %d\n", i);
	} // fin del primer bucle
		
}


void InitFCC(){
	
	vector cr,Lcv;
	
	VSetAll(Lcv,Lc);

	m=1;
	for (k=0; k<ncx; k++) {
		for (j=0; j<ncx; j++) {
			for (i=0; i<ncx; i++) {
				VSet(cr,(i+0.25)*Lc,(j+0.25)*Lc,(k+0.25)*Lc);
				//VMul(cr,cr,Lcv);
				VVSAdd(cr,-0.5,Lv);
				VVSet(r[m],cr);
				VSet(cr,0.5*Lc,0.5*Lc,0.);VAdd(r[m+1],cr,r[m]);
				VSet(cr,0.5*Lc,0.,0.5*Lc);VAdd(r[m+2],cr,r[m]);
				VSet(cr,0.,0.5*Lc,0.5*Lc);VAdd(r[m+3],cr,r[m]);
				m+=4;
				
			}
		}
	}
	
}


void InitVels(double norma){

	double x1, x2, w, momento[3+1], rp[3+1][npart+1];

		for (xi=1; xi<=3; xi++){

		momento[xi]=0.;	// inicializa momento de primer orden 
		
		for (i=2;i<=npart;i+=2) { 
			do {
				x1=2.0*aleat(LR)-1.0;
				x2=2.0*aleat(LR)-1.0;
				w=x1*x1+x2*x2;
			} while(w>=1.0);
			
			w = sqrt((-2.0*norma*log(w))/w);
			rp[xi][i-1] = x1*w;
			rp[xi][i] = x2*w;
				
			momento[xi] += (rp[xi][i-1]+rp[xi][i]);	// $ comienza calculo del momento de 1er orden		
		}
		    		  
		momento[xi] /= npart;	// $ acaba calculo del momento de 1er orden

		do_part rp[xi][i] -= momento[xi]; // $ corrige la distribucion de velocidad para evitar el flujo 
												 // comentar lineas marcadas $ para suprimir correccion de momento
	}	
	do_part {
		v[i].x = rp[1][i];
		v[i].y = rp[2][i];
		v[i].z = rp[3][i];
		}
	
}


void InitForce(){
	do_part {a[i].x = 0.;a[i].y = 0.;a[i].z = 0.;}
}			


void InitPCAccels(){
	do_part {a1[i].x = 0.;a1[i].y = 0.;a1[i].z = 0.;}
	do_part {a2[i].x = 0.;a2[i].y = 0.;a2[i].z = 0.;}

}			
