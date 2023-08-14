/*
  Green-Gauss scheme with linear interpolation for face values
  ->'wrong' interpolation on different level elements 
  -> Least-square method to be implemented
*/

#include "strdata.h"

void gradcalc_V2secondorder (struct RUN* run) {
  
  int iv,f,ing;
  struct BRANCH * brch;

  brch=run->topo->locl; 
  while (brch!=NULL){

    for(iv=0;iv<run->con->neq;iv++){
      brch->el->SG[iv][0]=0.0;
      brch->el->SG[iv][1]=0.0;
      brch->el->SG[iv][2]=0.0;
    }

    for(f=0;f<brch->nlfc;f++){  // Loop for faces
      for (ing=0;ing<brch->nsfc[f];ing++) {

        facercnstr (brch,ing,f,0,0,run);
        if(brch->cl->fc[f].bc==0) {
          facercnstr (brch,ing,f,1,1,run);
        } else {
          facercnstr (brch,ing,f,0,1,run);
          bc(brch,brch->cl->fc[f].bc,f,run);
        }
	      for(iv=0;iv<run->con->neq;iv++){
          brch->el->SG[iv][0]+=0.5*(run->sol1->vec[iv]+run->sol->vec[iv])*brch->cl->nx[f]*brch->cl->Area[f][ing]/brch->cl->Vol;
	        brch->el->SG[iv][1]+=0.5*(run->sol1->vec[iv]+run->sol->vec[iv])*brch->cl->ny[f]*brch->cl->Area[f][ing]/brch->cl->Vol;
	        brch->el->SG[iv][2]+=0.5*(run->sol1->vec[iv]+run->sol->vec[iv])*brch->cl->nz[f]*brch->cl->Area[f][ing]/brch->cl->Vol;
	      }

        brch->el->SG[5][0]+=0.5*(run->sol1->p_hydro+run->sol->p_hydro)*brch->cl->nx[f]*brch->cl->Area[f][ing]/brch->cl->Vol;
	      brch->el->SG[5][1]+=0.5*(run->sol1->p_hydro+run->sol->p_hydro)*brch->cl->ny[f]*brch->cl->Area[f][ing]/brch->cl->Vol;
	      brch->el->SG[5][2]+=0.5*(run->sol1->p_hydro+run->sol->p_hydro)*brch->cl->nz[f]*brch->cl->Area[f][ing]/brch->cl->Vol;

      }
    }
    
    brch=brch->lnxt;
  }
}
