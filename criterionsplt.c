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

void criterionsplt (struct RUN* run) {
  
  int f,ind;
  int level;
  double rx,ry,rz;
  double rho,rho_grad_mag,u_grad_mag,v_grad_mag,p_grad_mag;
  struct BRANCH * brch;

  brch=run->topo->locl; 
  while (brch!=NULL){
    level=1;

    rx=0.0;ry=0.0;rz=0.0;
    for (ind=0;ind<brch->nlnd;ind++){
      rx+=brch->cl->nd[ind].x/((double)brch->nlnd);
      ry+=brch->cl->nd[ind].y/((double)brch->nlnd);
      rz+=brch->cl->nd[ind].z/((double)brch->nlnd);
    }	
    
    if (run->con->geo_adapth==1) {
      if(run->con->i_case==204){
        level=run->con->level;
      }

      if(run->con->i_case==2500){
        rho = brch->el->S[0];

        if (rho<50.0) {
          level=run->con->level;
        }

      }

    } else {

      if(run->con->i_case==204){
        rho_grad_mag = sqrt(pow(brch->el->SG[0][0],2.0) + pow(brch->el->SG[0][1],2.0) + pow(brch->el->SG[0][2],2.0));
        u_grad_mag   = sqrt(pow(brch->el->SG[1][0],2.0) + pow(brch->el->SG[1][1],2.0) + pow(brch->el->SG[1][2],2.0));
        v_grad_mag   = sqrt(pow(brch->el->SG[2][0],2.0) + pow(brch->el->SG[2][1],2.0) + pow(brch->el->SG[2][2],2.0));
        p_grad_mag   = sqrt(pow(brch->el->SG[5][0],2.0) + pow(brch->el->SG[5][1],2.0) + pow(brch->el->SG[5][2],2.0));

        if (rho_grad_mag>100.0) {
          level=run->con->level;
        }
      }

      if(run->con->i_case==2500){
        rho = brch->el->S[0];

        if (rho<50.0) {
          level=run->con->level;
        }

      }

    }

    brch->split=level;
    brch=brch->lnxt;
  }
}
