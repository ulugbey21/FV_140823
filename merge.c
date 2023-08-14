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

void merge (int linked, RUN* run, BRANCH * brch)
{
  int ikid,ikeen,ifc,im,iv,imap,ikidnd,ikid2;
  int ix,iy,iz;
  double etakid1,etakid2,etakid3;
  double eta1,eta2,eta3;
  double of1,of2,of3;
  struct PARTICLE * drop;
  double S,V,VTOT;


  double TC1[4][4];
  double TC2[4][4];
  int IMAP[8];
  MPI_Status status;

  IMAP[0]=0;
  IMAP[1]=2;
  IMAP[2]=1;
  IMAP[3]=0;
  IMAP[4]=2;
  IMAP[5]=1;
  IMAP[6]=3;
  IMAP[7]=3;

  TC1[0][0]=-0.25;  TC2[0][0]= 0.0 ;
  TC1[1][0]= 0.25;  TC2[1][0]= 0.0 ;
  TC1[2][0]=-0.25;  TC2[2][0]= 0.5 ;
  TC1[3][0]=-0.75;  TC2[3][0]=-0.5 ;

  TC1[0][1]=-0.25;  TC2[0][1]= 0.0;
  TC1[1][1]= 0.25;  TC2[1][1]= 0.0;
  TC1[2][1]= 0.25;  TC2[2][1]= 0.5;
  TC1[3][1]= 0.75;  TC2[3][1]=-0.5;

  TC1[0][2]= 0.0 ;  TC2[0][2]= 0.0 ;
  TC1[1][2]= 1.0 ;  TC2[1][2]= 0.0 ;
  TC1[2][2]= 0.0 ;  TC2[2][2]= 0.5 ;
  TC1[3][2]= 0.0 ;  TC2[3][2]= 0.5 ;

  TC1[0][3]= 0.5 ;  TC2[0][3]= 0.0 ;
  TC1[1][3]=-0.5 ;  TC2[1][3]= 0.0 ;
  TC1[2][3]= 0.0 ;  TC2[2][3]=-0.5 ;
  TC1[3][3]= 0.0 ;  TC2[3][3]=-0.5 ;

 // printf("kid volume %d %d %e %e \n",brch->kids[0]->root,brch->kids[0]->tag,brch->kids[0]->cl->Vol,brch->kids[0]->el->S[0]);

  brch->part=brch->kids[0]->part;
  if (linked==1){
    if (run->topo->locl==brch->kids[0]){
      run->topo->locl=brch;
      brch->lnxt=brch->kids[brch->nkids-1]->lnxt;
      brch->lprv=brch->kids[0]->lprv;
    }
    else{
      brch->lnxt=brch->kids[brch->nkids-1]->lnxt;
      brch->lprv=brch->kids[0]->lprv;
    }
  }
  leafallocation(run,brch);
  brch->split=0;
  brch->el->crith=0.0;
  for (ikid=0;ikid<brch->nsplit;ikid++){
    brch->split=max(brch->split,brch->kids[ikid]->split);
    brch->el->crith=max(brch->el->crith,brch->kids[ikid]->el->crith);
  }
  normalvectorcalc(run,brch);
  tagvectorcalc(run,brch);
  volmproperties(0,brch,run);
   
  run->topo->nleavestmp++;
  run->topo->lleaves[brch->part]-=1;

  for(iv=0;iv<run->con->neq;iv++){
    brch->el->S[iv]=0.0;
    
    brch->el->SG[iv][0]=0.0;
    brch->el->SG[iv][1]=0.0;
    brch->el->SG[iv][2]=0.0;

    for (ikid=0;ikid<brch->nsplit;ikid++){
      //brch->el->S[iv]+=brch->kids[ikid]->el->S[iv]*brch->kids[ikid]->cl->Vol/brch->cl->Vol;
      brch->el->S[iv]+=brch->kids[ikid]->el->S[iv]/((double) brch->nsplit);

      brch->el->SG[iv][0]+=brch->kids[ikid]->el->SG[iv][0]*brch->kids[ikid]->cl->Vol/brch->cl->Vol;
      brch->el->SG[iv][1]+=brch->kids[ikid]->el->SG[iv][1]*brch->kids[ikid]->cl->Vol/brch->cl->Vol;
      brch->el->SG[iv][2]+=brch->kids[ikid]->el->SG[iv][2]*brch->kids[ikid]->cl->Vol/brch->cl->Vol;
    }
  }


  if (linked==1){
    if (brch->lnxt!=NULL) {brch->lnxt->lprv=brch;}
    if (brch->lprv!=NULL) {brch->lprv->lnxt=brch;}
  }
  for (ikid=0;ikid<brch->nsplit;ikid++){
    if (brch->kids[ikid]->nkids!=0){printf ("hey has kids \n");}

    leafdeallocation(run,brch->kids[ikid]);
    celldeallocation(brch->kids[ikid]);
    destroybrch(brch->kids[ikid]);
//    if(brch->kids[ikid]!=NULL){
	//    printf ("hey isnt null \n");
 //   }
 //   printf ("deallocated kid %d \n",ikid);
  } 
  free(brch->kids);
  brch->kids=NULL;
  brch->nkids=0;
}



