#include "strdata.h"

void Init_S(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w,double p) {

	double e_internal;
	
	crnt->el->S[0] = rho;
	crnt->el->S[1] = rho*u;
	crnt->el->S[2] = rho*v;
	crnt->el->S[3] = rho*w;
  
	if(run->con->model==1){
		e_internal = (p + run->par->matergama[0]*run->par->materpinf[0])/(rho*(run->par->matergama[0]-1.0)); 
		crnt->el->S[4] = rho*( e_internal + 0.5*(pow(u,2.0) + pow(v,2.0) + pow(w,2.0)));
	}
	
}


void Init_S3(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w,double ymass) {

	double e_internal;
	
	crnt->el->S[0] = rho;
	crnt->el->S[1] = rho*u;
	crnt->el->S[2] = rho*v;
	crnt->el->S[3] = rho*w;

	crnt->el->S[4] = rho*ymass;
	
}
