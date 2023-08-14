#include "strdata.h"

void cons2primtv(struct SOL* sol,struct RUN * run) {
//void cons2primtv(struct SOL* sol_temp,struct RUN * run){
	double e_internal;
	double a;

	sol->wec[0] = sol->r;
	sol->wec[1] = sol->vec[1]/sol->r; 
	sol->wec[2] = sol->vec[2]/sol->r; 
	sol->wec[3] = sol->vec[3]/sol->r; 

 //  printf("sol->wec[1]:%f\n",a);

	
	if(run->con->model==1){
		e_internal = (sol->vec[4]/sol->r) - 0.5*(pow(sol->vec[1],2.0) + pow(sol->vec[2],2.0) + pow(sol->vec[3],2.0));
		sol->wec[4] = sol->r*(run->par->matergama[0]-1.0)*e_internal - run->par->materpinf[0];	
  }

	if(run->con->model==3 ||run->con->model==4){ // mass fraction
		sol->wec[4] = sol->vec[4]/sol->r;
	}

}
