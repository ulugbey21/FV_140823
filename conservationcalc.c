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

void conservationcalc(RUN *run){
char  fstring[50]="                    ";
int iv;
double cons[50];
double constot[50];
FILE *fp;
struct BRANCH * crnt;

  for(iv=0;iv<50;iv++) {
    cons[iv]=0.0;
  }
 
  crnt=run->topo->locl;
  while (crnt!=NULL){ 
    for(iv=0;iv<run->con->neq;iv++){
      cons[iv]+=crnt->el->S[iv]*crnt->cl->Vol;    
    }
    cons[run->con->neq]+=crnt->cl->Vol;
    crnt=crnt->lnxt;
  }

  MPI_Allreduce((&cons),(&constot),50,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if(run->con->rank==0){
    sprintf(fstring,"conservation.dat");
    fp=fopen(fstring,"a");
    fprintf(fp," %d %le   ",run->con->tstep,run->con->time);
    for(iv=0;iv<run->con->neq+1;iv++){
      fprintf(fp," %le   ",constot[iv]);
    }
    fprintf(fp," \n");
    fclose(fp);
  }

}