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

void restructtopo (RUN * run) {
  
  struct BRANCH * crnt;
  int iel,ifc;
  int aa,ba,ca,da;
  
  //repartition(run);  // rearrange branches on processors, rearrange keen buffers
  repartition_parmetis(run);  // rearrange branches on processors, rearrange keen buffers
 
  rebuildlist(run);  // restitch local liked list
  
  rebuildprox(run);  // rearrange neig buffers
  
  rebuildconn(run);  // populate buffer
  
}
