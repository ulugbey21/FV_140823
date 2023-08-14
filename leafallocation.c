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

void leafallocation(RUN * run,struct BRANCH *crnt){

	int m,ord,typ,iv,ifc,ing,idr;

  crnt->el=malloc(sizeof(struct LEAF));
	crnt->el->inner=0;
	crnt->el->crith=0;
	crnt->el->SN=malloc(run->con->neq *sizeof(double ));
	crnt->el->S =malloc(run->con->neq*sizeof(double ));

	for(iv=0;iv<run->con->neq;iv++){ 
	  crnt->el->S[iv]=0.0; 
	}
  for(iv=0;iv<run->con->neq;iv++){ 
	  crnt->el->SN[iv]=0.0; 
	}
      
  crnt->el->SG = malloc((run->con->neq+1)*sizeof(double *));
  for(iv=0;iv<run->con->neq+1;++iv){
    crnt->el->SG[iv] = malloc(3*sizeof(double));
    for(idr=0;idr<3;++idr){
      crnt->el->SG[iv][idr]=0.0;
    }
  }

  crnt->el->SF        = malloc(6*sizeof(double *));
  //printf(" crnt->el->SF:%d\n", crnt->el->SF);
  crnt->el->flux_face = malloc(6*sizeof(double *));
  for(ifc=0;ifc<6;++ifc){
    crnt->el->SF[ifc]        = malloc((run->con->neq)*sizeof(double ));
    crnt->el->flux_face[ifc] = malloc((run->con->neq)*sizeof(double ));
    for(iv=0;iv<run->con->neq;++iv){
      crnt->el->SF[ifc][iv]        = 0.0;
      crnt->el->flux_face[ifc][iv] = 0.0;  
    }
  }
  
	crnt->el->flux_flag = malloc(6*sizeof(int));
	for(ifc=0;ifc<6;++ifc){
    crnt->el->flux_flag[ifc] = 0;  
	}
  
  crnt->el->RHS=malloc(run->con->neq*sizeof(double ));
  for(iv=0;iv<run->con->neq;++iv){
    crnt->el->RHS[iv] = 0.0;
  }

}