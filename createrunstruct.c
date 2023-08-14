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
// Memory allocation of struck run
//
//-------------------------------------------------------

#include "strdata.h"

void  createrunstruct (struct RUN * run) {

  run->con  = malloc(sizeof(CON));
  run->topo = malloc(sizeof(FOREST));
  run->sol  = malloc(sizeof(SOL));
  run->sol1 = malloc(sizeof(SOL));
  run->par  = malloc(sizeof(PAR));

}
