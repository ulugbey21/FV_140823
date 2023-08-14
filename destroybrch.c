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

void destroybrch (BRANCH * brch)


{
int ikid,ikeen,ind,ifc,i;

    free(brch->adrs);
    brch->adrs=NULL;
    for (i=0;i<brch->nkeens;i++){
	    free(brch->neigtr[i]);
	    free(brch->neignd[i]);
    }
    free(brch->neigtr);
    brch->neigtr=NULL;
    free(brch->neignd);
    brch->neignd=NULL;
    free(brch->neigfc);
    brch->neigfc=NULL;
    free(brch->neigag);
    brch->neigag=NULL;
    free(brch->lfc_neigs);
    brch->lfc_neigs=NULL;
    free(brch->nsfc);
    brch->nsfc=NULL;
    for (i=0;i<brch->nkeens;i++){
	    free(brch->ipartbound[i]); 
    }
    free(brch->ipartbound);
    brch->ipartbound=NULL;
    free(brch->keen);
    brch->keen=NULL;
    free(brch->keenfc);
    brch->keenfc=NULL;
    free(brch->keenfcangle);
    brch->keenfcangle=NULL;
    free(brch->kids);
   
 //   if (brch->kids!=NULL) {printf("oups destroy\n");}
    free(brch);
    brch=NULL;
}
