#include "strdata.h"

void bc_inlet(struct BRANCH * crnt,int f,struct RUN *run) {

  double rho,u,v,w,ru,rv,rw,p,e,e_internal,ymass;
  double xw, yw, zw;
  double U_B, H_C, W_N;
  crnt=run->topo->locl;  

  xw=crnt->cl->xc;
  yw=crnt->cl->yc;
  zw=crnt->cl->zc;

	if (run->con->i_case==200) { // Forward faceing step, https://amroc.sourceforge.net/examples/euler/2d/html/ffstep_n.htm
  
		rho = 1.4;
		u   = 3.0;
		v   = 0.0;
		w   = 0.0;
		ru  = rho*u;
		rv  = rho*v;
		rw  = rho*w ;
		p   = 1.0;

    e_internal = (p + run->par->matergama[0]*run->par->materpinf[0])/(rho*(run->par->matergama[0]-1.0)); 
		e = rho*( e_internal + 0.5*(pow(u,2.0) + pow(v,2.0) + pow(w,2.0)));
		
		asinvalues(run,crnt,run->sol1,rho,ru,rv,rw,e);

  }

	if (run->con->i_case==203) { // boundary layer
  
		rho = 1.0;
		u   = 68;
		v   = 0.0;
		w   = 0.0;
		ru  = rho*u;
		rv  = rho*v;
		rw  = rho*w ;
		p   = 101325.0;

    e_internal = (p + run->par->matergama[0]*run->par->materpinf[0])/(rho*(run->par->matergama[0]-1.0)); 
		e = rho*( e_internal + 0.5*(pow(u,2.0) + pow(v,2.0) + pow(w,2.0)));
		
		asinvalues(run,crnt,run->sol1,rho,ru,rv,rw,e);

  }

  	if (run->con->i_case==2900) { // injector
        U_B = 1.5625;
		H_C = 0.032;
		rho = 998.2;
		W_N = 0.001;
		u   = (9.0*U_B/4.0)*pow((1-(yw/(H_C/2))),2.0)*pow((1-(zw/(W_N/2))),2.0);
		v   = 0.0;
		w   = 0.0;
		ru  = rho*u;
		rv  = rho*v;
		rw  = rho*w ;
		ymass = 0.0;
		
		asinvalues_multiphase(run,crnt,run->sol1,rho,ru,rv,rw,e,ymass);

  }  	    
  
       if (run->con->i_case==290) { // LES channel flow edelbauer

		rho = 1000.0;
		u   = 50.00;
		v   = 0.0;
		w   = 0.0;
		ru  = rho*u;
		rv  = rho*v;
		rw  = rho*w ;
		//p   = 30000000.0;
		ymass = 1.0-run->con->ymin;
		
		asinvalues_multiphase(run,crnt,run->sol1,rho,ru,rv,rw,e,ymass);

		

  }

    	if (run->con->i_case==3400) { // injector
        rho = 722.7;
		u   = 200.0;
		v   = 0.0;
		w   = 0.0;
		p   = 1.2*pow(10.0,6.0); 
		ru  = rho*u;
		rv  = rho*v;
		rw  = rho*w ;
		ymass = 1.0-run->con->ymin;
		
		asinvalues_multiphase(run,crnt,run->sol1,rho,ru,rv,rw,e,ymass);

  }  	
  /*
  else if (run->con->i_case==2100) { // injector
		rho = 0.0170;
		u   = 6.95;
		v   = 0.0;
		w   = 0.0;
		ru  = rho*u;
		rv  = rho*v;
		rw  = rho*w ;
		ymass = 1-run->con->ymin;
		
		asinvalues_multiphase(run,crnt,run->sol1,rho,ru,rv,rw,e,ymass);

  }

*/



	return;
}
