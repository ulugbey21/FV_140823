#include "strdata.h"

void reconctruct_V3(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,
                                                              double dist_m1_cf,double dist_m1_cc) {
  
  int j,iv;
  double gratio,grad,product;

  // Variable extrapolation
  
  for(iv=0;iv<run->con->neq_temp;++iv){    //variables
    run->gradlt[iv]=0.0;
		run->gradlt_face[f][iv] = 0.0;
    if (crnt->cl->fc[f].bc==0){  
      //printf("run->delta0[iv]:%f\n",run->delta0[0]);

      if (fabs(run->delta0[iv])>pow(10.0,-10.0)) {   

				product = (-run->delta1[iv]/dist_m1_cc)*(run->delta0[iv]/dist_p1_cc);

        if (product>0.0) {
					if ( fabs(-run->delta1[iv]/dist_m1_cc) > fabs(run->delta0[iv]/dist_p1_cc)) {
						grad =  run->delta0[iv]/dist_p1_cc;
					} else {
						grad = -run->delta1[iv]/dist_m1_cc;
					}
					
        } else {
          grad = 0.0;
        }

        run->gradlt[iv]         = grad;
				run->gradlt_face[f][iv] = grad;
      }
    }

		if (run->con->iprimtv==1){
			run->wec_temp[iv]   = run->sol->wec[iv] + dist_p1_cf*run->gradlt[iv];
		} else {
			//run->vec_reconstruct[iv] = run->sol->vec[iv] + dist_p1_cf*run->gradlt[iv];
      crnt->el->SF[f][iv] = run->sol->vec[iv] + dist_p1_cf*run->gradlt[iv];
    //  printf("crnt->el->SF[f][0]:%f\n",crnt->el->SF[f][iv]);
		}
        
  }

  if (run->con->iprimtv==1){
    prim2conc(run);  // run->wec_temp to run->vec_temp
    for(iv=0;iv<run->con->neq;++iv){        
      crnt->el->SF[f][iv] = run->vec_temp[iv];
    }
  }
 
  
  return;

}