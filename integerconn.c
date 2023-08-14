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

int tagconn (BRANCH * brch,int ifc) {

  int tagneig;
  struct BRANCH * crnt; 
  struct BRANCH * neig; 
  int ifc2,iang,ilevel,ikid,ipnt,indnd,ifnd;
  int * adrneig;

  adrneig=malloc(10*sizeof(int));
  crnt=brch->keen[ifc];
  if (crnt!=brch->prnt||crnt->root==-1){
    tagneig=tagaddress(brch->keen[ifc]->adrs,brch->keen[ifc]->level);
    if (crnt->root==-1){brch->lfc_neigs[ifc]=0; }
  }

  if ((crnt==brch->prnt)&&(crnt->root!=-1)){
    while(crnt->level!=crnt->keen[ifc]->level&&crnt->keen[ifc]->root!=-1){
      crnt=crnt->keen[ifc];
    }
    
    if (crnt->keen[ifc]->root!=-1){
      neig=crnt->keen[ifc];
      ifc2=crnt->keenfc[ifc];
      iang=crnt->keenfcangle[ifc];
      ilevel=crnt->level;
      while((neig->nkids!=0)&&(neig->level<brch->level)){
        ilevel=neig->level;
        if (neig->nsplit==8){
	        adrneig[ilevel]=fcrefnd2elnd(ifc2,fcndrot(elnd2fcnd(ifc,brch->adrs[ilevel+1],brch->type),iang,brch->cl->fc[ifc].type),neig->type);
        }
      }
      // printf ("tagneig %d %d %d %d \n",ilevel,adrneig[0],adrneig[1],adrneig[ilevel]);
      tagneig=tagaddress(adrneig,max(ilevel-1,1));
    }
    if (crnt->keen[ifc]->root==-1){
      tagneig=-1;
    }
  }
  free(adrneig);
  return tagneig;
}


int tagconnfin (BRANCH * brch,int ifc, int ineig) {
  
  struct BRANCH * crnt;
  struct BRANCH * neig;
  int ifc2,iang,ilevel,fctype,ifnd,ifnds,ipnt,indnd,ikid;
  int tagfin;
  neig=brch->neigtr[ifc][0];
  
  if (neig->nkids!=0&&neig->root!=-1){
    ifc2=brch->neigfc[ifc];
    iang=brch->neigag[ifc];
    fctype=brch->cl->fc[ifc].type;

    if (neig->nsplit==8){
      brch->lfc_neigs[ifc]=4;
      brch->neigtr[ifc][0]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(0,iang,fctype),neig->type)];
      brch->neigtr[ifc][1]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(1,iang,fctype),neig->type)];
      brch->neigtr[ifc][2]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(2,iang,fctype),neig->type)];
      brch->neigtr[ifc][3]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(3,iang,fctype),neig->type)];
    }
  }
  return  tagfin;
}

int ielconn (BRANCH * brch,int ifc) {
  int ielneig;
  struct BRANCH * keen;
  int adrneig[10];

  // pick kin
  ielneig=0;
  keen=brch->keen[ifc];
  
  // if kin is sibling, not a ghost not a prnt
  //  printf ("ielcon A %d %d %d \n",brch->root,ifc,ielneig);
  if (keen!=brch->prnt||keen->root==-1){
    ielneig=keen->root;
    // printf ("ielcon B %d %d %d \n",brch->root,ifc,ielneig);
  }

  // if kin is parent not a ghost
  if ((keen==brch->prnt)&&(keen->root!=-1)){
    //  printf ("ielcon C %d %d %d \n",keen->root,ifc,ielneig);
    // climb to find sibling
    while(keen->level!=keen->keen[ifc]->level&&keen->keen[ifc]->root!=-1){
      keen=keen->keen[ifc];
     // printf ("ielcon C1 %d %d %d \n",keen->root,ifc,ielneig);
    }
    
    // iel is sibling root
    if (keen->keen[ifc]->root!=-1){
      ielneig=keen->keen[ifc]->root;
      //   printf ("ielcon C2 %d %d %d %d %d \n",keen->root,keen->level,brch->level,ifc,ielneig);
    }
    
    if (keen->keen[ifc]->root==-1){
      ielneig=keen->keen[ifc]->root;
      // printf ("ielcon C3 %d %d %d \n",keen->root,ifc,ielneig);
    }

  }

  return ielneig;
}