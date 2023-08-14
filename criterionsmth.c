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

void criterionsmth (struct RUN* run) {
  int f,ing;
  
  struct BRANCH * brch;
  struct BRANCH * prnt;
  int max_lvl,maxsplit,imerge;
  int level,ipass,ikid;

  watch(run,8,0);
  communicate_S(run,3);
  watch(run,8,1);
  for (ipass=0;ipass<(run->con->level+1);ipass++){
  brch=run->topo->locl;
  while (brch!=NULL){
    if (brch->level>1){
      prnt=brch->prnt;
      imerge=0;
      maxsplit=0;
      for(ikid=0;ikid<prnt->nkids;ikid++){
        if (prnt->kids[ikid]->nkids==0&&prnt->kids[ikid]->split<prnt->kids[ikid]->level){
          imerge++;
        }
      }
      if (imerge==prnt->nkids){
        brch->split=brch->level-1;
      }
      else{
        if (brch->level>brch->split){brch->split=brch->level;}
      }
    }
    if (brch->level<brch->split){brch->split=brch->level+1;}
  brch=brch->lnxt;}

            watch(run,8,0);
  communicate_S(run,3);
            watch(run,8,1);
  
    max_lvl=1;
    brch=run->topo->locl; 
    while (brch!=NULL){
      max_lvl=1;
      for(f=0; f<brch->nlfc; f++){  // Loop for faces
        for (ing=0;ing<brch->lfc_neigs[f];ing++){
          if (brch->cl->fc[f].bc==0){
	    max_lvl=max(brch->neigtr[f][ing]->split,max_lvl);
	  }
	}
      }
      if (max_lvl-brch->split>1){
	brch->split=max_lvl-1;
      }
    brch=brch->lnxt;}
            watch(run,8,0);
    communicate_S(run,3);
            watch(run,8,1);
  }
}
