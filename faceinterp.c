/*
 * Forest of oct-Trees FV implementation
 * octo-Rhodon
 *
 *  \/\/\/\/
 *   \/ / /
 *    \/ /
 *     \/
 *
 * City, University of London,
 * School of Mathematics, Computer Science and Engineering,
 * Department of Mechanical Engineering and Aeronautics.
 * London, UK.
 * Andreas Papoutsakis 
 * 01/09/2020
*/


#include "strdata.h"

void faceinterp(RUN *run){

	int i,j,iv,v,f,fo,ieq,Nneig,ineig;
  double dist_p1_cf,dist_p1_cc,dist_m1_cf,dist_m1_cc,dx;
  double dist_cf,dist_cc;
  struct BRANCH * crnt;

  run->T2nd=timecpu(run->T2nd,0);

  crnt=run->topo->locl;
  while (crnt!=NULL){
    
    for(f=0; f<(crnt->nlfc); ++f) {

      fo=fcopp(f,crnt->type);
       
      distance(run,crnt,f ,&dist_p1_cf,&dist_p1_cc); // center-face, center-center
      distance(run,crnt,fo,&dist_m1_cf,&dist_m1_cc); //  

      deltavalues(run,crnt,run->delta1,fo);  // i-1 , i   -> delta1
      deltavalues(run,crnt,run->delta0,f );  // i   , i+1 -> delta0
      
      reconctruct_V3(run,crnt,f,dist_p1_cf,dist_p1_cc,dist_m1_cf,dist_m1_cc); // new grad

    }
    
    crnt=crnt->lnxt;
  }
  /*
      for(i=0;i<3;i++){
      for(j=0;j<3;j++){
              printf("%d\n",G[i][j]);
          }
     }
   */

  run->T2nd=timecpu(run->T2nd,1);

  return;
	
}
