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

int tagaddress (int * adrs,int nlvl) {
  int ilvl;
  int tag;
  tag=0;
  for (ilvl=2;ilvl<nlvl+1;ilvl++){
    tag+=((int) pow(10,ilvl-2))*(1+adrs[ilvl]);
    
  }
  return tag;
}
