#include "strdata.h"

void prim2conc(struct RUN * run) {

	double rho,e_internal;

	// rho
	run->vec_temp[0] = run->wec_temp[0];
	rho              = run->wec_temp[0];
	
	// ru rv rw
	run->vec_temp[1] = rho*run->wec_temp[1];
	run->vec_temp[2] = rho*run->wec_temp[2];
	run->vec_temp[3] = rho*run->wec_temp[3];
	
		// E
		e_internal = (run->wec_temp[4] + run->par->materpinf[0])/(rho*(run->par->matergama[0]-1.0));
		run->vec_temp[4] =  rho*( e_internal + 0.5*(pow(run->wec_temp[1],2.0) + pow(run->wec_temp[2],2.0) + pow(run->wec_temp[3],2.0)));
	if (run->con->model==3||run->con->model==4){
		run->vec_temp[4] = rho*run->wec_temp[4];
	}
	
}
