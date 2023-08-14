#include "strdata.h"

void RK(int irk,RUN *run) {

	int e,f,v;
	double buf;
	struct BRANCH * crnt;
	
	crnt=run->topo->locl;  	// RUN->FOREST->TREE (points to current tree)
	while (crnt!=NULL){ 	
		
		for(f=0; f<(crnt->nlfc); f++){  // Loop faces     
			crnt->el->flux_flag[f]=0;			
		}
		
		switch(irk){

			case 0:
			for(e=0; e<run->con->neq; e++){ 	
			  crnt->el->SN[e] = crnt->el->S[e]; // SN:Solution vector 2, S:Solution vector
			}
			break;

			case 1:
			for(e=0; e<run->con->neq; e++){
				crnt->el->S[e] = crnt->el->SN[e] + run->con->dt *(1.0/crnt->cl->Vol) * crnt->el->RHS[e];	
			}
			break;

			case 2:
			for(e=0;e<run->con->neq;e++){
	      crnt->el->S[e]=(3.0/4.0)*crnt->el->SN[e]+(1.0/4.0)*crnt->el->S[e]+(1.0/4.0)*run->con->dt*crnt->el->RHS[e]/crnt->cl->Vol;
			}
			break;

			case 3:
			for(e=0;e<run->con->neq;e++){
			  crnt->el->S[e]=(1.0/3.0)*crnt->el->SN[e]+(2.0/3.0)*crnt->el->S[e]+(2.0/3.0)*run->con->dt*crnt->el->RHS[e]/crnt->cl->Vol;
			}

			break;

			case 4:
			for(e=0;e<run->con->neq;e++){
			  crnt->el->S[e]=(1.0/2.0)*crnt->el->SN[e]+(1.0/2.0)*crnt->el->S[e]+(1.0/2.0)*run->con->dt*crnt->el->RHS[e]/crnt->cl->Vol;
			}
			break;

			default:
			break;
		}

  	crnt=crnt->lnxt; 
	}

}
