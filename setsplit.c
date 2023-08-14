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

void setsplit(RUN *run){
struct BRANCH * crnt;

  crnt=run->topo->locl; while (crnt!=NULL){
    setsplitcalc(crnt,run);
  crnt=crnt->lnxt;}
}


void setsplitcalc(BRANCH * crnt,RUN *run){
int ifc,fcqs;
double maxprod,prod;

    maxprod=-10000.0;
    if (run->con->splitmode==4){
      for (ifc=0;ifc<crnt->nlfc;ifc++){
        prod=crnt->cl->nx[ifc]*run->con->nxquad+crnt->cl->ny[ifc]*run->con->nyquad+crnt->cl->nz[ifc]*run->con->nzquad;
        if(maxprod<prod){
          crnt->fcqsplit=ifc;
          maxprod=prod;
        }
      }
      crnt->nsplit=4;
      if (crnt->type==2&&crnt->cl->fc[crnt->fcqsplit].type==0){
        crnt->nsplit=2;
      }
    }
    if (run->con->splitmode==8){
      crnt->nsplit=8;
      crnt->fcqsplit=-1;
    }
    //printf("qsplit %d %d \n",crnt->root,crnt->fcqsplit);
}


