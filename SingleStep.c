/*
 *  SingleStep.c
 *  LJMD
 *
 *  Created by Juanjo Meléndez on 02/02/09.
 *  Modified by Francisco Vega Reyes on 20/07/09, 19/07/10.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
 
#include "params.h"
//#include "def_functions.h"


void SingleStep(){	
	++it;
	t = it*dt;
	LeapfrogStep(1);	
	ApplyBoundaryCond();
	ComputeForces();
	LeapfrogStep(2);
	}


void LeapfrogStep(int ileap){


	if (ileap==1){
		do_part{
			v[i].x += 0.5*dt*a[i].x;
			v[i].y += 0.5*dt*a[i].y;
			v[i].z += 0.5*dt*a[i].z;

			r[i].x += dt*v[i].x;
			r[i].y += dt*v[i].y;
			r[i].z += dt*v[i].z;
		}
	} // fin de if
	else {
		do_part{
			v[i].x += 0.5*dt*a[i].x;
			v[i].y += 0.5*dt*a[i].y;
			v[i].z += 0.5*dt*a[i].z;
		}		
	}//fin de else
	
	
}// fin de la funcion

void ApplyBoundaryCond(){

	do_part {
		
		if (r[i].x >= 0.5*Lv.x) r[i].x -= Lv.x;
		else if (r[i].x < - 0.5*Lv.x) r[i].x += Lv.x;

		if (r[i].y >= 0.5*Lv.y) r[i].y -= Lv.y;
		else if (r[i].y < - 0.5*Lv.y) r[i].y += Lv.y;

		if (r[i].z >= 0.5*Lv.z) r[i].z -= Lv.z;
		else if (r[i].z < - 0.5*Lv.z) r[i].z += Lv.z;
	}
} 

void ComputeForces(){

	cfuerzas=1;
	
for (i=1; i<=npart; i++){
		a[i].x = 0.;
		a[i].y = 0.;
		a[i].z = 0.;
	}	

	up = 0.;
	virial = 0.;

	for (i=1; i<npart; ++i){

		for (j=i+1; j<=npart; ++j){

			dr.x = r[i].x - r[j].x;
			dr.y = r[i].y - r[j].y;
			dr.z = r[i].z - r[j].z;
			//VSub(drb,r[i],r[j]);
			/*printf("drx: %g\n",dr.x);
			printf("drxb: %g\n",drb.x);*/
			
			if (dr.x >= 0.5*Lv.x) dr.x -= Lv.x;
			else if (dr.x < - 0.5*Lv.x) dr.x += Lv.x;
			//if(i==1) printf("%g\n",dr.x);

			if (dr.y >= 0.5*Lv.y) dr.y -= Lv.y;
			else if (dr.y < - 0.5*Lv.y) dr.y += Lv.y;
			//if(i==1) printf("%g\n",dr.y);

			if (dr.z >= 0.5*Lv.z) dr.z -= Lv.z;
			else if (dr.z < - 0.5*Lv.z) dr.z += Lv.z; 
			//if(i==1) printf("%g\n",dr.z);

			rij2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
			if (rij2 < rc2 ){
				if(rij2!=0){
				ir2=1./rij2;
				ir6=ir2*ir2*ir2;
				ljr=48.*ir6*(ir6-0.5)*ir2;
				/*if(ljr==0) printf("CUIDADO!!\n");
				if(dr.x==0) printf("CUIDADO!!\n");
				if(dr.y==0) printf("CUIDADO!!\n");
				if(dr.z==0) printf("CUIDADO!!\n");*/
				
				a[i].x += ljr*dr.x;a[i].y += ljr*dr.y;a[i].z += ljr*dr.z; // fuerza sobre la part’cula i
				}
				a[j].x -= ljr*dr.x;a[j].y -= ljr*dr.y;a[j].z -= ljr*dr.z; // idem sobre la j
			
				up += 4. * ir6*(ir6-1.)+1.;
				//virial += ljr*rij2;
			} // fin de condicion if
			
		}// fin de bucle en j

	} // fin de bucle en i
	
	up /= npart;


}


void ComputeForcesCS(){
	
//	vector inva, rs, shift,rcc;
	vector inva, rs, shift;
	ivector cc, m1v, m2v, vOff[]=VECINS;
	int c, j1, j2, mi1, mi2, m1x, m1y, m1z, ivecin;
	//int iref=90;
	
	VDiv (inva, nvcelda, Lv);
	
	
	for (i=npart+1;i<=npart+VProd(nvcelda);i++) lista[i]=-1;// inicia lista[c] en -1
	 
	// este bucle crea lista de celdas
	do_part{
		
		VSAdd(rs, r[i], 0.5, Lv); // pasa coord de 0 a L, para hallar celda
		//if (i==40) {printf("r(i): %g %g %g\n", rs.x, rs.y, rs.z);printf("%g\n",inva.y);}
		//VMul(rcc, rs, inva); // numero de celda, guardado en cc
		//rcc.x=rs.x*inva.x;rcc.y=rs.y*inva.y;rcc.z=rs.z*inva.z;
		//if (i==40) printf("cc(i): %f %f %f\n", rcc.x, rcc.y, rcc.z);
		//(cc).x=(int)(rcc).x;(cc).y=(int)(rcc).y;(cc).z=(int)(rcc).z;
		VMul(cc, rs, inva); // numero de celda, guardado en cc
		//if (i==40) printf("cc(i): %d %d %d\n", cc.x, cc.y, cc.z);
		c=VLinear(cc,nvcelda)+npart+1; // numero de celda pasado a escalar
		//if (i==40) printf("c(i): %d\n", c);
		lista[i]=lista[c]; //'lista' de la 1a particula es -1, 'lista' de la siguiente es i-1
		lista[c]=i; // al final, lista[c], de la celda, guarda la ultima particula en ella
		
	} //fin bucle do_part
	
	//for (i=npart+1;i<=npart+VProd(nvcelda);i++) printf("celda: %d ultima: %d\n",i,lista[i]);
	
	do_part VZero(a[i]);
	up = 0.;
	virial = 0.;
	
	for (m1z=0; m1z < nvcelda.z; m1z++){
		for (m1y=0; m1y < nvcelda.y; m1y++){
			for (m1x=0; m1x < nvcelda.x; m1x++){
				
				VSet(m1v,m1x,m1y,m1z);
				
				mi1=VLinear(m1v,nvcelda)+npart+1;
				//if(m1z==0) if (m1y==1) if (m1x==3) printf("celda: %d, indices: %d %d %d\n",mi1,m1x,m1y,m1z);   //indice 3D a 1
				
				//if(mi1<=80) printf("mi1: %d %d %d %d CUIDADOmenos\n",VLinear(m1v,nvcelda),mi1, m1v.y, m1v.z);
				//if(mi1>VProd(nvcelda)+npart) printf("mi1: %d %d %d %d CUIDADO\n",VLinear(m1v,nvcelda),m1v.x, m1v.y, m1v.z);
				
				 for (ivecin=0;ivecin<N_VECIN;ivecin++) {//printf("%d\n",ivecin);
				
					 VAdd(m2v,m1v,vOff[ivecin]); //obtiene indice 3D de celda vecina
					 VZero(shift);
					 /*if (mi1==iref) {
						 printf("m2vantes: %d %d %d: \n",m2v.x,m2v.y,m2v.z);
						 
					 }*/
					 VCellWrapAll(); 
					/* if (mi1==iref) {
						 printf("m2v: %d %d %d: \n",m2v.x,m2v.y,m2v.z);
						 printf("shift: %g %g %g \n",shift.x,shift.y,shift.z);
					 }*/
					 mi2=VLinear(m2v,nvcelda)+npart+1; // lo pasa a formato 1D
					// if(mi1==iref) printf("mi2: %d\n",mi2);				
				 
					 do_celdas(j1,mi1){
				     do_celdas(j2,mi2){
						 
						 //if (mi1==90) printf("j2: %d\n",j2);
						 //printf("mi2: %d\n",mi2);
						 /*
						 if (j1==lista[mi1]) {// si estamos en otra celda o el par tiene indice menor (dentro de la misma celda)
							  if (mi1==iref) {if (j2==lista[mi2]) printf("j2: %d en vecino: %d ultima particula: %d\n",j2,mi2,lista[mi2]);}
						 }
						 */
						 
						 if (mi1!=mi2 || j2<j1) {
							 VSub(dr,r[j1],r[j2]);
							 VVSub(dr,shift);
							 rij2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
							 
							 if (rij2 < rc2 ){
								 if(rij2!=0){ // esta condicion probablemente no hace falta
									 ir2=1./rij2;
									 ir6=ir2*ir2*ir2;
									 ljr=48.*ir6*(ir6-0.5)*ir2;
									 
									 a[j1].x += ljr*dr.x;a[j1].y += ljr*dr.y;a[j1].z += ljr*dr.z; // fuerza sobre la part’cula i
								 }
								 a[j2].x -= ljr*dr.x;a[j2].y -= ljr*dr.y;a[j2].z -= ljr*dr.z; // idem sobre la j
								 
								 up += 4. * ir6*(ir6-1.)+1.; 
								 
							 } // fin de condicion if rij2
							 
						 } // fin if m1!=m2
						 
						 
										 
					 } //fin do_celdas(j2), mi# es fijo en cada do_celdas
					 } // fin do_celdas(j1)
				 
				 } // fin bucle ivecin
				 
				 
			} //fin bucle m1x
		}//fin bucle m1y
	} //fin bucle m1z
	
	up /= npart;

	
}


void ForcesInelastic(){
	

	
	for (i=1; i<=npart; i++){
		a[i].x = 0.; aw[i].x=0.;
		a[i].y = 0.; aw[i].y=0.;
		a[i].z = 0.; aw[i].z=0.;
	}	
	
	//up = 0.;
	
	for (i=1; i<npart; ++i){
		
		for (j=i+1; j<=npart; ++j){
			
			dr.x = r[i].x - r[j].x;
			dr.y = r[i].y - r[j].y;
			dr.z = r[i].z - r[j].z;

			// toma en cuenta las condiciones periodicas
			if (dr.x >= 0.5*Lv.x) dr.x -= Lv.x;
			else if (dr.x < - 0.5*Lv.x) dr.x += Lv.x;
			//if(i==1) printf("%g\n",dr.x);
			
			if (dr.y >= 0.5*Lv.y) dr.y -= Lv.y;
			else if (dr.y < - 0.5*Lv.y) dr.y += Lv.y;
			//if(i==1) printf("%g\n",dr.y);
			
			if (dr.z >= 0.5*Lv.z) dr.z -= Lv.z;
			else if (dr.z < - 0.5*Lv.z) dr.z += Lv.z; 
			//if(i==1) printf("%g\n",dr.z);
			// hasta aqui lo de la periodicidad
			
			rij2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
			
			if (rij2 < sigma ){
				if(rij2!=0){
					
					rij=sqrt(rij2);
					ir2=1./rij2;
					ir=1./rij;
					
					// velocidad relativa de las particulas
					vij.x=v[i].x-v[j].x;
					vij.y=v[i].y-v[j].y;
					vij.z=v[i].z-v[j].z;
					
					// vr= vijárij
					vr=vij.x*dr.x+vij.y*dr.y+vij.z*dr.z;
					// componente normal de la velocidad relativa
					vijn.x = vr*dr.x*ir2;
					vijn.y = vr*dr.y*ir2;
					vijn.z = vr*dr.z*ir2;
					// componente tangencial de la velocidad relativa
					rot.x=(w[i].y+w[j].y)*dr.z-(w[i].z+w[j].z)*dr.y;
					rot.y=(w[i].z+w[j].z)*dr.x-(w[i].x+w[j].x)*dr.z;
					rot.z=(w[i].x+w[j].x)*dr.y-(w[i].y+w[j].y)*dr.x;
					vijt.x =vij.x-vijn.x-0.5*rot.x;
					vijt.y =vij.y-vijn.y-0.5*rot.y;
					vijt.z= vij.z-vijn.z-0.5*rot.z;
					
					aNe[i].x = Y*(sigma-rij)*ir*dr.x;aNe[i].y = Y*(sigma-rij)*ir*dr.y;aNe[i].z = Y*(sigma-rij)*ir*dr.z;
					aNi[i].x = -gamman*vijn.x;aNi[i].y = -gamman*vijn.y;aNi[i].z = -gamman*vijn.z;					
					
					aT[i].x = -gammashear*vijt.x;aT[i].y = -gammashear*vijt.y;aT[i].z = -gammashear*vijt.z;
					
					atot.x= (aNe[i].x +aNi[i].x+aT[i].x);atot.y += (aNe[i].y+aNi[i].y+aT[i].y);atot.z +=(aNe[i].z+aNi[i].z+aT[i].z);
					awtot.x = dr.y*atot.z-dr.z*atot.y; awtot.y = dr.z*atot.x-dr.x*atot.z; awtot.z = dr.x*atot.y-dr.y*atot.x;
					
					a[i].x += atot.x; a[i].y += atot.y; a[i].z += atot.z;
					aw[i].x -= (awtot.x*K); aw[i].y -= (atot.y*K); aw[i].z -= (atot.z*K);
					
					
					// fuerza sobre la part’cula i
					
				}
				
				a[j].x -= atot.x;a[j].y -= atot.y;a[j].z -= atot.z;
				aw[j].x -= (awtot.x*K); aw[j].y -= (atot.y*K); aw[j].z -= (atot.z*K);
				// idem sobre la j
				
				//up += 4. * ir6*(ir6-1.)+1.;
				//virial += ljr*rij2;
			} // fin de condicion if rij
			
		}// fin de bucle en j
		
	} // fin de bucle en i
	
	up /= npart;
	
	
}




void PredictorStep()
{
	
	double cr[]={19., -10., 3.}, cv[]={27., -22., 7.}, div=24.,dtr,dtv;
//	rva *mol;
	
	dtr=Sqr (dt)/ div;
	dtv= dt / div;
	DO_MOL{
		ro[i]=r[i];
		vo[i]=v[i];
		PR(x);
		PRV(x);
		PR(y);
		PRV(y);
		PR(z);
		PRV(z);
		a2[i]=a1[i];
		a1[i]=a[i];

	}
		
}

void CorrectorStep()
{
	
	double cr[]={3., 10., -1.}, cv[]={7., 6., -1.}, div=24.,dtr,dtv;
	int i;
	
	dtr=Sqr (dt)/ div;
	dtv= dt / div;

	do_part{
		
		CR(x);
		CRV(x);
		CR(y);
		CRV(y);
		CR(z);
		CRV(z);
		
	}
	
}
