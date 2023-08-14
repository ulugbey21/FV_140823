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

//-------------------------------------------------------
//
// Read data from .nd (nodes) and .el (elements) and creates topology
//
//-------------------------------------------------------

// glob is the branch outside the domain, 
// neighbours to branch are boundary cells, 
// bc faces face glob, level 1 prnt is glob, 
// last list member next is NULL, first 
// list member prev is NULL

#include "strdata.h"

void spawntopo (RUN * run) {
  
  int type,typg,nlfc,nlnd;
  int iel,ielp1,istart;
  int ind;
  int node;

  struct BRANCH * crnt;
  struct BRANCH * rcnt;


  istart=0;
  crnt=run->topo->glob;
  if((run->con->verbose==1)&&(run->con->rank==0)){ printf ("MMFV: MAIN: TREES SPAWNED TO BRANCHES 0 \n");}
  for (iel=0;iel<run->topo->nleaves;iel++){ // Loop elements
    if (run->topo->drys[iel]->part==run->con->rank){
      rcnt=run->topo->drys[iel]->brch;	  
      if (istart==0){
        rcnt->lnxt=NULL;
        rcnt->lprv=run->topo->glob;
	run->topo->locl=rcnt;
	istart=1;
      }
      else{
        rcnt->lprv=crnt;
        crnt->lnxt=rcnt;
        rcnt->lnxt=NULL;
      }
      crnt=rcnt;
      crnt->prnt=run->topo->glob;  //sure?
      crnt->root=iel;
      crnt->level=1;
      crnt->adrs[0]=0;
      crnt->adrs[1]=iel;
      crnt->tag=tagaddress (crnt->adrs,crnt->level);

      for (ind=0;ind<crnt->nlnd;ind++){
        crnt->cl->nd[ind].num=run->topo->drys[iel]->vrtx[ind];
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if((run->con->verbose==1)&&(run->con->rank==0)){ printf ("MMFV: MAIN: TREES SPAWNED TO BRANCHES I \n");}
  for (iel=0;iel<run->topo->nleaves;iel++){ // Loop elements
    if (run->topo->drys[iel]->part==run->con->rank){
      crnt=run->topo->drys[iel]->brch;
      for (ind=0;ind<crnt->nlnd;ind++){
        crnt->cl->nd[ind].x=run->topo->vertex[run->topo->drys[iel]->vrtx[ind]];
        crnt->cl->nd[ind].y=run->topo->vertey[run->topo->drys[iel]->vrtx[ind]];
        crnt->cl->nd[ind].z=run->topo->vertez[run->topo->drys[iel]->vrtx[ind]];
	//printf("nodes %d %d %le %le %le \n",iel,ind,crnt->cl->nd[ind].x,crnt->cl->nd[ind].y,crnt->cl->nd[ind].z);
      }
    }
    //createlfc(crnt);
  }
	  MPI_Barrier(MPI_COMM_WORLD);
  if((run->con->verbose==1)&&(run->con->rank==0)){ printf ("MMFV: MAIN: TREES SPAWNED TO BRANCHES II \n");}
}
