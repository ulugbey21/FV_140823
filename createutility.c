/*
 * Forest of oct-Trees DG implementation
 * octo-Rhodon
 *
 *  \/\/\/\/
 *   \/ / /
 *    \/ /
 *     \/
 *
 * Foundation for Research and Technology Hellas, 
 * Institute for Applied and Computational Mathematics (FORTH/IACM)
 * Heraklion Crete
 * Andreas Papoutsakis 
 * 20/05/2014
*/

#include "strdata.h"

void createutility(RUN *run){

  int i,j,v,f,ieq,neq_temp;

  run->SendS=malloc(((run->con->neq))*sizeof(double));
  run->RecvS=malloc(((run->con->neq))*sizeof(double));
  run->SendSG=malloc(run->con->neq*3*sizeof(double));
  run->RecvSG=malloc(run->con->neq*3*sizeof(double));
  if (run->con->ischeme==2){
    run->SendSF=malloc(run->con->neq*6*4*sizeof(double));
    run->RecvSF=malloc(run->con->neq*6*4*sizeof(double));
  }
    
  run->flux_adv     = malloc(run->con->neq*sizeof(double *));
  run->flux_mc     = malloc(run->con->neq*sizeof(double *));
  run->flux_viscous = malloc(run->con->neq*sizeof(double *));
  run->flux_HLLC    = malloc(run->con->neq*sizeof(double));
	for (v=0;v<run->con->neq;++v){
    run->flux_adv[v]     = malloc(3*sizeof(double));
    run->flux_viscous[v] = malloc(3*sizeof(double));
    run->flux_mc [v]    = malloc(3*sizeof(double));
  }

  run->vel_temp = malloc(3*sizeof(double));
	

	run->AUX_S=malloc((run->con->neq)*sizeof(double ));
	for(v=0;v<run->con->neq;v++){
		run->AUX_S[v]=0.0;
	}


  run->vecx0=malloc(3* sizeof(double *));
  run->vecph=malloc(3* sizeof(double *));
  for (i=0;i<3;i++){
    run->vecx0[i]=malloc(3 * sizeof(double));
    run->vecph[i]=malloc(run->con->neq * sizeof(double));
  }

  run->veccf=malloc(3* sizeof(double));
  run->xx=malloc(3* sizeof(double));

  if (run->con->iprimtv==0) { neq_temp = run->con->neq;       }
  if (run->con->iprimtv==1) { neq_temp = run->con->nprimitiv; }
  
  run->vec_reconstruct =malloc(run->con->neq * sizeof(double));
  run->wec_temp =malloc(run->con->nprimitiv * sizeof(double));
  run->vec_temp =malloc(run->con->neq * sizeof(double));
  run->delta0=malloc(neq_temp * sizeof(double));
  run->delta1=malloc(neq_temp * sizeof(double));
  run->delta2=malloc(neq_temp * sizeof(double));
  run->gradlt=malloc(neq_temp * sizeof(double));
  run->gradlt_face=malloc(6 * sizeof(double *));
  run->delta0_face=malloc(6 * sizeof(double *));
  run->delta1_face=malloc(6 * sizeof(double *));
  run->vecl=malloc(run->con->neq*sizeof(double));
  run->wecl=malloc(run->con->neq*sizeof(double));

  run->fluxs=malloc(run->con->neq*sizeof(double ));
  
  for (ieq=0;ieq<run->con->neq;ieq++) {
    run->fluxs[ieq]=0.0;
  }

  run->vdelta1=malloc(6 * sizeof(double *));
  for (f=0;f<6;f++) {
    run->vdelta1[f]=malloc(neq_temp * sizeof(double));
  }

  for (f=0;f<6;f++) {
    run->gradlt_face[f]=malloc(neq_temp * sizeof(double));
    run->delta0_face[f]=malloc(neq_temp * sizeof(double));
    run->delta1_face[f]=malloc(neq_temp * sizeof(double));
  }

  run->vel_temp    = malloc(3*sizeof(double)); // static
  run->tensor_temp = malloc(9*sizeof(double));
 
  // Auxilliary variables for hllc + hllc_br
  run->Us_temp = malloc(run->con->neq*sizeof(double));
  for(v=0;v<run->con->neq;++v){    
		run->Us_temp[v] = 0.0;        
	}	

	  
  
  run->con->avec = malloc(3*sizeof(double ));
  run->con->bvec = malloc(3*sizeof(double ));
  run->con->cvec = malloc(3*sizeof(double ));
  run->con->amat = malloc(3*sizeof(double *));
  run->con->bmat = malloc(3*sizeof(double *));
  run->con->cmat = malloc(3*sizeof(double *));
  run->con->xmat = malloc(3*sizeof(double *));
  for (i=0;i<3;i++){
    run->con->amat[i] = malloc(3*sizeof(double));
    run->con->bmat[i] = malloc(3*sizeof(double));
    run->con->cmat[i] = malloc(run->con->neq*sizeof(double));
    run->con->xmat[i] = malloc(run->con->neq*sizeof(double));
  }

	return;

}
