#include "strdata.h"

void inner_el(RUN * run) {

	int i,iel,Nneig,ifc,f,ielneig,ineig,ikeens;
	struct BRANCH * crnt;

	

	crnt=run->topo->locl;
  while (crnt!=NULL){
    iel=crnt->root;

		crnt->el->inner=1;

    for (ifc=0; ifc<crnt->nlfc;ifc++){
			if(crnt->cl->fc[ifc].bc==0){
				ielneig=ielconn(crnt,ifc);
				if (ielneig!=-1){
					if (run->topo->drys[iel]->part!=run->topo->drys[ielneig]->part){
						crnt->el->inner=0;
					}
				}
			}
    }
    crnt=crnt->lnxt;
  }

  return;
}