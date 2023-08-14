#include "strdata.h"

void deltavalues(struct RUN *run,struct BRANCH *crnt,double * delta_temp,int f_temp) {

  int iv;
  
  facevalues_MUSCL(crnt,f_temp,0,run->sol,run); // sol

  if(crnt->cl->fc[f_temp].bc==0){
    facevalues_MUSCL(crnt,f_temp,1,run->sol1,run);  // sol 1
  } else {
    facevalues_MUSCL(crnt,f_temp,0,run->sol1,run);  // sol 1
    bc(crnt,crnt->cl->fc[f_temp].bc,f_temp,run);        // sol 1
  }
  
  if (run->con->iprimtv==1){
    cons2primtv(run->sol,run);  // vec  to wec
    cons2primtv(run->sol1,run); // vec1 to wec1
    for(iv=0;iv<run->con->nprimitiv;++iv) { delta_temp[iv]=run->sol1->wec[iv]-run->sol->wec[iv]; }
  } else {
    for(iv=0;iv<run->con->neq;++iv)       { delta_temp[iv]=run->sol1->vec[iv]-run->sol->vec[iv]; }  
  }

  return;

}