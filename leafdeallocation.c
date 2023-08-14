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

void leafdeallocation(RUN * run,struct BRANCH *crnt){

  register int i;
  int m,x,y,f,j,k,isp;

  free(crnt->el->SN);
  free (crnt->el->S);
  
  
  for(m=0;m<6;m++){
    free(crnt->el->SF[m]);
    free(crnt->el->flux_face[m]);
  }
  for(m=0;m<3;m++){
    free(crnt->el->SG[m]);
  }
  
  free(crnt->el->SF);
  free(crnt->el->flux_face);
  free(crnt->el->SG);

  free(crnt->el->flux_flag);
  
  free(crnt->el->RHS);
  free(crnt->el);
	crnt->el=NULL;

}
