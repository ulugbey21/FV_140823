#include "strdata.h"

void memallocation(RUN *run){

	struct BRANCH * crnt;
	
	solallocation(run->sol,run);
	solallocation(run->sol1,run);	

	crnt=run->topo->locl; 
	while (crnt!=NULL){
		leafallocation(run,crnt);
		crnt=crnt->lnxt;
	} //el


}

void solallocation(SOL * sol, RUN *run){
  
	int v,i;

	sol->vec       = malloc(run->con->neq*sizeof(double ));
	//printf("malloc-vel:%d\n",(run->con->neq*sizeof(double )));
	sol->wec       = malloc(run->con->nprimitiv*sizeof(double ));
	//printf("malloc-wel:%d\n",(run->con->nprimitiv*sizeof(double )));

	//sol->vec       = malloc(run->con->neq*sizeof(double ));

	sol->flux      = malloc(run->con->neq*sizeof(double* ));
	
	for (v=0;v<run->con->neq;v++){
		sol->flux[v] = malloc(3*sizeof(double));
	}

	sol->st   = malloc(3*sizeof(double *));
	for (i=0;i<3;i++){
		sol->st[i]   = malloc(3*sizeof(double));
	}

}
