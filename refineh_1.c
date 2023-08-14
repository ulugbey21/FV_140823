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

void refineh (RUN * run,int imode) {
  
  struct BRANCH * brch1;
  double xnd,ynd,znd,rcrit,rvol;
  int ipass,ifc,ipt,ing,icnl;
  int nmerge;
  int ikid,kids,nkids,iel,ikeen;
  struct BRANCH * prnt;
  struct BRANCH * crnt;
  struct BRANCH * kid;

  crnt=run->topo->locl; 
  while (crnt!=NULL) {
    prnt=crnt->prnt;
    if (prnt->level>=1) {
      kids=0;
      for (ikid=0;ikid<prnt->nkids;ikid++){
        kids+=prnt->kids[ikid]->nkids;
      } 
      if (kids==0) {
        nmerge=0;
        for (ikid=0;ikid<prnt->nkids;ikid++){
          kid=prnt->kids[ikid];
          if (kid->split<kid->level&&kid->nkids==0){
            icnl=1;
            nmerge+=icnl;
          }
         }
         nkids=prnt->nkids; 
         if (nmerge==nkids) { 
           crnt=prnt->kids[prnt->nsplit-1]->lnxt;
           merge(1,run,prnt);
           run->topo->nleaves-=(prnt->nsplit-1);
         } else {
           crnt=crnt->lnxt;
          }
      }
      if (kids!=0){
        crnt=crnt->lnxt;
      }
    } else {
      crnt=crnt->lnxt;
    }
  }


  crnt=run->topo->locl; 
  while (crnt!=NULL) {
    
    if (crnt->split>crnt->level){  
      split(1,run,crnt);   
      run->topo->nleaves+=(crnt->nsplit-1);
      crnt=crnt->kids[crnt->nsplit-1]->lnxt;
    } else {
      crnt=crnt->lnxt;
    }
  }

 
}