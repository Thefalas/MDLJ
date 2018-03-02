/*
 *  termo.c
 *  LJMD
 *
 *  Created by Francisco Vega Reyes on 20/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 *	Aqui estan todas las funciones de medida y toma de datos
 *
 */

#include "params.h"
//#include "def_functions.h"

int midetemp(vector rr[npart+1]){

	double v2x, v2y, v2z;

	v2x=0.0;v2y=0.;v2z=0.;
	v2=0.0;
	
	for (i=1;i<=npart;i++) v2+=((rr[i].x*rr[i].x)+(rr[i].y*rr[i].y)+(rr[i].z*rr[i].z));
	
	v2 /= npart;// aqui v2 es (d)(T/m)	
	v2 /= 3.;	 // aqui ya es la temperatura (T)
	
	return 0;	
}

int midemom(vector rr[npart+1]){

	double mx, my, mz;

	mx=0.0;my=0.;mz=0.;
	
	for (i=1;i<=npart;i++) {
	mx+=(rr[i].x);
	my+=(rr[i].y);
	mz+=(rr[i].z);
	}
	
	m2=mx*mx+my*my+mz*mz;
	
	//printf ("it: %d,  m: %3.3f\t %3.3f\t %3.3f \n", it, mx,my,mz);
	
	return 0;	
}


int rdf(void){
	
	double rangor2;
	
	limr=utermo; //acumula la funcion de distribucion radial en todos los 'it' entre escrituras de archivos
	
	
	rangor=Lv.x; // distancia maxima para rho(r)
	rangor2=rangor*rangor; 
	
	bin=rangor/nbin; // tamaño del bin (anchura de la division de rho(r))
	
	if (itr==1) for (i=0; i<nbin; i++) {histo[i]=0;rbin[i]=i-0.5;}
	
	for (i=1;i<npart;i++){
		for (j=i+1;j<=npart;j++){
		
			VSub(dr,r[i],r[j]);
			VWrapAll(dr);
			rij2=VLengSq(dr);
			if(rij2<rangor2) ibin=sqrt(rij2)/bin+1;
			histo[ibin]+=1.;
		
		}
		
	}
	
	++itr;

	
	return 0;
	
}


int corrFCC(void){

	vector kvec;
	double si, sr, phi;
	
	if (itr==1) cred=0.;
	
	kvec.x=2*3.14159265*ncx/Lv.x;
	kvec.y=-kvec.x;
	kvec.z=kvec.x;
	sr=0.;si=0.;
	
	do_part{
		phi=VDot(kvec,r[i]);
		sr+=cos(phi);
		si+=sin(phi);
		
	}
	
	cred+=sqrt(Sqr(sr)+Sqr(si));
		
	return 0;
	
}
