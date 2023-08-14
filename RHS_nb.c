#include "strdata.h"

void RHS_nb(struct RUN * run,int inner_outer){
  
  struct BRANCH * crnt;

  crnt=run->topo->locl;
  while (crnt!=NULL){
    if (crnt->el->inner==inner_outer){
        //residual(crnt,run);
        residual_faces_AMR(crnt,run);
    }
    crnt=crnt->lnxt;
  }

}