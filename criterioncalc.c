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

//---------------------------------------------------------------------------
//
//
//
//---------------------------------------------------------------------------

#include "strdata.h"

void criterioncalc (struct RUN* run) {
  int f,ing;
  
  double * criterion;
  double * criterionsum;
  struct BRANCH * brch;
  struct BRANCH * crnt;
  double grad;
  double gradx;
  double grady;
  double gradz;
  int level;
  
  brch=run->topo->locl; 
  while (brch!=NULL){
    grad=0.0;

      gradx=0.0; grady=0.0; gradz=0.0;
      for(f=0; f<brch->nlfc; f++) {  // Loop for faces
        for (ing=0;ing<brch->nsfc[f];ing++){
	   
          if(brch->cl->fc[f].bc==0){
            facevalues (brch,ing,f,0,0,run);
            facevalues (brch,ing,f,1,1,run);
	        } else {
            facevalues (brch,ing,f,0,0,run);
            facevalues (brch,ing,f,0,1,run);
            bc(brch,brch->cl->fc[f].bc,f,run);
          }
	        gradx+=0.5*(run->sol1->p_hydro+run->sol->p_hydro)*brch->cl->nx[f]*brch->cl->Area[f][ing]/brch->cl->Vol;
	        grady+=0.5*(run->sol1->p_hydro+run->sol->p_hydro)*brch->cl->ny[f]*brch->cl->Area[f][ing]/brch->cl->Vol;
	        gradz+=0.5*(run->sol1->p_hydro+run->sol->p_hydro)*brch->cl->nz[f]*brch->cl->Area[f][ing]/brch->cl->Vol;
	      }
      }

      if (run->con->ctype==1){
        grad=sqrt(gradx*gradx+grady*grady+gradz*gradz)*pow(brch->cl->Vol/brch->cl->xc,1.0/3.0);
      }
      if (run->con->ctype==2){
        grad=sqrt(gradx*gradx+grady*grady+gradz*gradz)*pow(brch->cl->Vol,1.0/3.0);
      }
      if (run->con->ctype==3){
        grad=sqrt(gradx*gradx+grady*grady+gradz*gradz)*pow(brch->cl->Vol,1.0/3.0);
      }

      //brch->el->crith=grad;
      //brch->el->crith=brch->el->SI[0];
    
      brch=brch->lnxt;
    }
}
