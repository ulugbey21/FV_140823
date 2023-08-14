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

void message(char* msg,RUN * run)
{
  int aa,ba,ca,da;
  if(run->con->verbose!=0){
    if(run->con->verbose>1){ MPI_Barrier(MPI_COMM_WORLD);}
    if(run->con->rank==0){printf (" %s ",msg);}
    if(run->con->verbose>1){ MPI_Barrier(MPI_COMM_WORLD);}
    if(run->con->verbose>2){ 
      getMemory(run,&aa, &ba, &ca, &da);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}
