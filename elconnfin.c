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

void elconnfin (BRANCH * brch,int ifc)
{
struct BRANCH * crnt;
struct BRANCH * neig;
int ifc2,iang,ilevel,fctype,ifnd,ifnds,ineig,ipnt,indnd,ikid;

  neig=brch->neigtr[ifc][0];
 // if (brch->root==4256&&brch->tag==1){printf("in elconfin %d tag %d ifc %d ng root %d neig tag %d neig nkid %d hng %d \n",brch->root,brch->tag,ifc,neig->root,neig->tag,neig->nkids,neig->hangn);	}
  
  if (neig!=NULL&&neig->nkids!=0){
        ifc2=brch->neigfc[ifc];
        iang=brch->neigag[ifc];
        fctype=brch->cl->fc[ifc].type;

        if (neig->nsplit==8){
//	  printf("in elconfin %d ifc %d ang %d ifc2 %d typ %d \n",brch->root,ifc,iang,ifc2,fctype);	
//	  printf("elconfin var 0 %d \n",brch->root);	
//	  printf("elconfin var 1 %d \n",ifc);	
//	  printf("elconfin var 2 %d \n",iang);	
//	  printf("elconfin var 3 %d \n",ifc2);	
//	  printf("elconfin var 4 %d \n",fctype);	
//	  printf("elconfin var 5 %d \n",neig->type);	
          brch->lfc_neigs[ifc]=4;
          brch->neigtr[ifc][0]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(0,iang,fctype),neig->type)];
          brch->neigtr[ifc][1]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(1,iang,fctype),neig->type)];
          brch->neigtr[ifc][2]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(2,iang,fctype),neig->type)];
          brch->neigtr[ifc][3]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(3,iang,fctype),neig->type)];
        }
        if (brch->nsplit==4){

          brch->lfc_neigs[ifc]=2;

          ineig=0;
          for (ifnds=0;ifnds<4;ifnds++){
            if(ifnds==0){ifnd=0;}
            if(ifnds==1){ifnd=1;}
            if(ifnds==2){ifnd=3;}
            if(ifnds==3){ifnd=2;}
            if (elnd2fcnd(brch->fcqsplit,fcnd2elnd(ifc,ifnd,brch->type),brch->type)!=-1){
              ipnt=fcndrot(ifnd,iang,fctype);
              indnd=fcrefnd2elnd(ifc2,ipnt,neig->type);
              if (neig->nsplit==4){
              if (indnd==fcnd2elnd(neig->fcqsplit,0,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,0,neig->type)){ikid=0;}
              if (indnd==fcnd2elnd(neig->fcqsplit,1,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,1,neig->type)){ikid=1;}
              if (indnd==fcnd2elnd(neig->fcqsplit,2,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,2,neig->type)){ikid=2;}
              if (indnd==fcnd2elnd(neig->fcqsplit,3,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,3,neig->type)){ikid=3;}
              }


              if (neig->nsplit==2){


              if(indnd<3) {ikid=0;}
              if(indnd>=3){ikid=1;}

              }

              brch->neigtr[ifc][ineig]=neig->kids[ikid];
              brch->neignd[ifc][ineig]=ifnd;
              ineig++;
            }
          }
        }


        if (brch->nsplit==2){
          ineig=0;
          fctype=brch->cl->fc[ifc].type;
          if (fctype==0){
          brch->lfc_neigs[ifc]=2;
         
          ineig=0;
          for (ifnds=0;ifnds<2;ifnds++){
            if(ifnds==0){ifnd=0;}
            if(ifnds==1){ifnd=2;}

              ipnt=fcndrot(ifnd,iang,fctype);
              indnd=fcrefnd2elnd(ifc2,ipnt,neig->type);
              if (neig->nsplit==2){
              if(indnd<3) {ikid=0;}
              if(indnd>=3){ikid=1;}
              }


              if (neig->nsplit==4){
              if (indnd==fcnd2elnd(neig->fcqsplit,0,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,0,neig->type)){ikid=0;}
              if (indnd==fcnd2elnd(neig->fcqsplit,1,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,1,neig->type)){ikid=1;}
              if (indnd==fcnd2elnd(neig->fcqsplit,2,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,2,neig->type)){ikid=2;}
              if (indnd==fcnd2elnd(neig->fcqsplit,3,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,3,neig->type)){ikid=3;}
              }
              brch->neigtr[ifc][ineig]=neig->kids[ikid];
              ineig++;
          }
        }
          if (fctype==1){
          brch->lfc_neigs[ifc]=1;
          fctype=brch->cl->fc[ifc].type;
              if (ifc2==3){ikid=0;}
              if (ifc2==4){ikid=1;}
              brch->neigtr[ifc][0]=neig->kids[ikid];
          }
        }
}
brch->nsfc[ifc]=max(1,brch->lfc_neigs[ifc]);
}




void neigsfin (RUN * run)
{
int ifc;
struct BRANCH * crnt;

    crnt=run->topo->locl; while (crnt!=NULL){
      crnt->split=0;
      for (ifc=0;ifc<crnt->nlfc;ifc++){
        elconnfin(crnt,ifc);
	crnt->nsfc[ifc]=max(1,crnt->lfc_neigs[ifc]);
      }
    crnt=crnt->lnxt;
  }

}




