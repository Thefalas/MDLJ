/*
 *  params.h
 *  LJMD
 *
 *  Created by Juanjo Meléndez on 02/02/09.
 *  Copyright 2009. All rights reserved.
 *
 */
 
 /*	Vector operations	*/

#define VAdd(v1,v2,v3)		\
	v1.x = v2.x + v3.x,	  v1.y = v2.y + v3.y,	v1.z = v2.z + v3.z

#define VSub(v1,v2,v3)	VAdd(v1,v2,-v3)

#define VDot(v1,v2)		\
	v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

#define VSAdd(v1, v2, s3, v3)	\
	v1.x = v2.x + s3*v3.x,	v1.y = v2.y + s3*v3.y,	v1.z = v2.z + s3*v3.z

#define VSet(v, sx, sy, sz)		\
	v.x = sx,	v.y = sy,	v.z = sz

#define VSetAll(v,s)	\
	VSet(v,s,s,s)

#define VZero(v)			\
	VSetAll(v,0)

#define VVSet(v,w)			\
	v.x = w.x, v.y = w.y, v.z = w.z

#define VVSAdd(v1,s2,v2)	\
	VSAdd(v1, v1, s2, v2)

#define VLengSq(v)	\
	VDot(v,v)

#define VNorm(v)	\
	sqrt(VLengSq(v))

#define Sqr(x) ((x)*(x))
#define Cube(x) ((x)*(x)*(x))

#define VMul(v1,v2,v3)	\
(v1).x=(v2).x*(v3).x,   \
(v1).y=(v2).y*(v3).y,   \
(v1).z=(v2).z*(v3).z

#define VDiv(v1,v2,v3)	\
(v1).x=(v2).x/(v3).x,   \
(v1).y=(v2).y/(v3).y,   \
(v1).z=(v2).z/(v3).z

#define VLinear(v1, v2) ( ((v1).z*(v2).y+(v1).y)*(v2).x+(v1).x)
#define VVSub(v1,v2) VSub(v1,v1,v2)



/* funciones */

#define PCR4(r,ro,v,a,a1,a2,p) r.p=ro.p+dt*v.p+dtr*(cr[0]*a.p+cr[1]*a1.p+cr[2]*a2.p)
#define PCV4(r,ro,v,a,a1,a2,p) v.p=(r.p-ro.p)/dt+dtv*(cv[0]*a.p+cv[1]*a1.p+cv[2]*a2.p)

#define PR(p) PCR4(r[i], r[i], v[i], a[i], a1[i], a2[i], p)

#define PRV(p) PCV4(r[i], ro[i], v[i], a[i], a1[i], a2[i], p)

#define CR(p) PCR4(r[i], ro[i], vo[i], a[i], a1[i], a2[i], p)

#define CRV(p) PCV4(r[i], ro[i], v[i], a[i], a1[i], a2[i], p)

#define AllocMem(a,n,t) a=(t*)malloc((n)*sizeof(t))


/*	Procedures	*/

#define VWrap(v, t)	\
	if (v.t >= 0.5 * Lv.t) v.t -= Lv.t;	\
	else if (v.t < -0.5*Lv.t) v.t += Lv.t
  
#define VWrapAll(v)		\
	{VWrap(v, x); VWrap(v,y); VWrap(v,z);}

#define do_part for (i=1; i<=npart; i++)

#define DO_MOL for (i=1;i<=npart;i++)

#define do_celdas(j,m) for (j=lista[m]; j>0; j=lista[j])

#define VProd(v) ((v).x*(v).y*(v).z)

#define VCellWrap(t) \
	if(m2v.t>=nvcelda.t){ \
		m2v.t=0;		\
		shift.t=Lv.t;	\
}else if (m2v.t<0){		\
	m2v.t=nvcelda.t-1;		\
	shift.t=-Lv.t;		\
}

#define VCellWrapAll() {VCellWrap(x);VCellWrap(y);VCellWrap(z);}
#define VSCopy(v2,s1,v1)	\
	(v2).x=(s1)*(v1).x,	\
	(v2).y=(s1)*(v1).y,	\
	(v2).z=(s1)*(v1).z
