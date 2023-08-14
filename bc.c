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

void bc(struct BRANCH * crnt, int p, int f, struct RUN *run){
  
  switch(p){

    case 0: // no bc  
    break;
    
    case 1: // Inflow  

      bc_inlet(crnt,f,run);
      
    break;
    
    case 2: // Outflow 
   // run->sol->p_hydro=12000000;

    break;

    case 3: // Reflective (slip wall)
      bc_reflective(crnt,f,run);
    break;

    case 4: // No slip wall
      bc_wall(crnt,f,run);
    break;

    default:
    break;
  }

}//function
