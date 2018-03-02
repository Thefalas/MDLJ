/*
 *  salida.c
 *  LJMD
 *
 *  Created by Francisco Vega Reyes on 20/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "params.h"
//#include "def_functions.h"

int salida (int numarchiv){
	
	char nombre[72];
	FILE *archivo;

	sprintf(nombre, "xyz%4.4i.dat", numarchiv);
	archivo = fopen (nombre, "w");
	for(i=1;i<=npart;i++){
	fprintf (archivo, "%d\t %5.5g \t %5.5g \t %5.5g\n", i, r[i].x, r[i].y, r[i].z);
	}
	fclose(archivo);

	sprintf(nombre, "v%4.4i.dat", numarchiv);
	archivo = fopen (nombre, "w");
	for(i=1;i<=npart;i++){
	fprintf (archivo, "%d\t %7.7g \t %7.7g \t %7.7g\n", i, v[i].x, v[i].y, v[i].z);
	}
	fclose(archivo);

	if(cfuerzas){
	sprintf(nombre,"f%4.4i.dat",numarchiv);
	archivo= fopen(nombre,"w");
	for (i=1;i<=npart;i++){
	fprintf(archivo,"%d \t %7.7g \t %7.7g \t %7.7g \n", i, a[i].x, a[i].y, a[i].z);
	}
	} // cond. if de calculo de fuerzas
	fclose(archivo);
	
	if(it){
	sprintf(nombre,"r%4.4i.dat",numarchiv);
	archivo=fopen(nombre, "w");
	for (i=1; i<nbin; i++) fprintf(archivo, "%g \t %g \n", bin*rbin[i] , histo[i]);
	}
	fclose(archivo);
	
	sprintf(nombre,"corrFCC.dat");
	archivo=fopen(nombre, "a");
	fprintf(archivo, "%d \t %g \n", numarchiv , cred);
	fclose(archivo);
	
return 0;

}

int mensajes0(void){
	
	
	printf("numero de celdas en un eje: %d\n",ncx);
	printf("npart %d\n",npart);
	printf("longitud de celda: %g\n",Lc);
	printf("longitud del sistema L: %g\n",L);
	printf("dx: %g\n",dx);
	printf("densidad: %g\n",n);
	printf("rc: %g\n",rc);
	printf("nceldas Fuerza en un eje: %d\n",nvcelda.x);
	
	printf ("it: %d\t, m2: %3.3f\t, T: %g\t, Ep: %g\t, E: %g\n", it, m2, v2, up, 1.5*v2+up);
	return 0;
}


int mensajes(void){

	printf ("it: %d\t, m2: %3.3f\t, T: %g\t, Ep: %g\t, E: %g\n", it, m2, v2, up, 1.5*v2+up);
	
	return 0;

}