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

void restart_load(RUN *run){

  int ierr,ileaf,n;
  int typ,rot,tag;
  int nleaves,ntrees,nv;
  FILE *rs;
  FILE *speciesfile;
  double buf;
  struct TREE * tree;
  struct BRANCH * crnt;

  rs=fopen("Init","r");

  ierr=fscanf(rs,"%d  %lf %d %d %d \n",&run->con->tstep,&run->con->time,&nleaves,&ntrees,&nv);
  run->con->tstep0=run->con->tstep;
  run->con->time0=run->con->time;


  for (ileaf=0;ileaf<nleaves;ileaf++){

    ierr=fscanf(rs,"%d %d ",&rot,&tag);
    
    tree=run->topo->drys[rot];

    if(tree->part==run->con->rank) {

      crnt=run->topo->drys[rot]->brch;
      typ=crnt->type;

      crnt=spawntree(run,tree,tag,rot,1);
      
      for(n=0;n<nv;n++){
        ierr=fscanf(rs,"%lf ",&crnt->el->S[n]);  
      }
    
    } else {
      for(n=0;n<nv;n++){
        ierr=fscanf(rs,"%lf ",&buf); 
      }
    }
    ierr=fscanf(rs,"\n");
    
 }
 fclose(rs);
}
