#include "strdata.h"

void asinvalues(struct RUN *run,struct BRANCH * crnt,struct SOL* sol, double rho, double ru,double rv,double rw,double e) {

	int i,j,iv;

	sol->r = rho;
	sol->u = ru / rho;
	sol->v = rv / rho;
	sol->w = rw / rho;


  sol->vec[0] = rho;
  sol->vec[1] = ru;
  sol->vec[2] = rv;
  sol->vec[3] = rw;


  if(run->con->model==1){
    sol->e          = e;
    sol->vec[4]     = e;
    sol->e_internal = (e/rho) - 0.5*(pow(sol->u,2.0) + pow(sol->v,2.0) + pow(sol->w,2.0));
    sol->p_hydro    = sol->r*(run->par->matergama[0]-1.0)*sol->e_internal - run->par->matergama[0]*run->par->materpinf[0];
    sol->c          = sqrt(run->par->matergama[0]*(sol->p_hydro + run->par->materpinf[0])/sol->r);
  } else if(run->con->model==2) {
    barotropic_eos(crnt,0,run);
  } else if(run->con->model==3) {
    barotropic_multiphase_eos(crnt,0,run);
  }
  else if(run->con->model==4) {
    barotropic_threephase_eos(crnt,0,run);
  }
  
  sol->st[0][0] = sol->p_hydro;
  sol->st[0][1] = 0.0;
  sol->st[0][2] = 0.0;

  sol->st[1][0] = 0.0;
  sol->st[1][1] = sol->p_hydro;
  sol->st[1][2] = 0.0;

  sol->st[2][0] = 0.0;
  sol->st[2][1] = 0.0;
  sol->st[2][2] = sol->p_hydro;



  if (isnan(sol->c )==1) {
    printf("c nan | %e | %e %e | %e %e %e \n",sol->c,run->par->matergama[0],sol->p_hydro,rho,sol->e,sol->e_internal);
    printf("vel: %e %e %e \n",sol->u,sol->v,sol->w);
    int ind;
    double rx,ry,rz;
   
    rx=0.0;ry=0.0;rz=0.0;
    for (ind=0;ind<crnt->nlnd;ind++){
      rx+=crnt->cl->nd[ind].x/((double)crnt->nlnd);
      ry+=crnt->cl->nd[ind].y/((double)crnt->nlnd);
      rz+=crnt->cl->nd[ind].z/((double)crnt->nlnd);
    }

    printf("coordinates: %e %e %e \n",rx,ry,rz);
    exit(0);
  }

  return;
}
