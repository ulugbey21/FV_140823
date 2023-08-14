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

void boundary(RUN *run){
int ifc,fc,i;
int ibc;
int el,iel;
int nfcbc;
FILE * fp;
struct BRANCH * crnt;
double ** BCTREES;

  BCTREES=malloc(run->topo->ntrees*sizeof(double));
  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
	   BCTREES[iel]=malloc(6*sizeof(double));
           for (ifc=0;ifc<6;ifc++){
		   BCTREES[iel][ifc]=0;
	   }
  }

  crnt=run->topo->locl; while (crnt!=NULL){
    for (ifc=0;ifc<crnt->nlfc;ifc++){crnt->cl->fc[ifc].bc=0;}
  crnt=crnt->lnxt; }
  


  for (ibc=1;ibc<5;ibc++){
    if (ibc==1) {fp=fopen(run->con->filebr1,"r");}
    if (ibc==2) {fp=fopen(run->con->filebr2,"r");}
    if (ibc==3) {fp=fopen(run->con->filebr3,"r");}
    if (ibc==4) {fp=fopen(run->con->filebr4,"r");}

    fscanf(fp,"%d \n",&nfcbc);
    for (i=0;i<nfcbc;i++){
      fscanf(fp,"%d\t%d \n",&el,&fc);
      BCTREES[el-1][fc-1]=ibc;
    }
    rewind(fp);
    fclose(fp);
  }

  crnt=run->topo->locl; while (crnt!=NULL){
    for (ifc=0;ifc<crnt->nlfc;ifc++){
      ibc=BCTREES[crnt->root][ifc];
      if (ibc==1){crnt->cl->fc[ifc].bc=1;}
      if (ibc==2){crnt->cl->fc[ifc].bc=2;}
      if (ibc==3){crnt->cl->fc[ifc].bc=3;}
      if (ibc==4){crnt->cl->fc[ifc].bc=4;}
    }
  crnt=crnt->lnxt; }

  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
	   free(BCTREES[iel]);
  }
  free(BCTREES);


}
