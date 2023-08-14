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

//============================================
//            Mesh Adaptation
//============================================

#include "strdata.h"

void meshadaptation (RUN * run) {

  int ipass;
  FILE *fp;
  char  fstring[50];
  
  for(ipass=0;ipass<run->con->level+1;ipass++) {

   
    gradcalc(run);
    criterionsplt(run);
    communicate_S(run,3);  
    criterionsmth(run);
    refineh(run,ipass%2);
    
    restructtopo(run);
       
    normalvector(run);
    tagvector(run);
    meshproperties(run);
    setsplit(run);
    
    communicate_S(run,0);

  }

  inner_el(run);

  if(run->con->rank==0){
    sprintf(fstring,"tot_cells.dat");
    fp=fopen(fstring,"a"); 
    fprintf(fp," %d %le %d %d \n",run->con->tstep,run->con->time,run->topo->ntrees,run->tot_leaves); 
    fclose(fp);
  }
 
}