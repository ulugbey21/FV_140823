#include "strdata.h"

void RHS(struct RUN * run){
  
  struct BRANCH * crnt;

  crnt=run->topo->locl;
  while (crnt!=NULL){
    //residual(crnt,run);
    residual_faces_AMR(crnt,run);
    crnt=crnt->lnxt;
  }

}