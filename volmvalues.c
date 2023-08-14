#include "strdata.h"

void volmvalues(struct BRANCH * crnt,struct RUN *run){

  struct SOL * r_sol;
  struct LEAF * crnt_el;

  r_sol = run->sol;
	crnt_el=crnt->el;

  r_sol->r = crnt_el->S[0];
  r_sol->u = crnt_el->S[1] / r_sol->r;
  r_sol->v = crnt_el->S[2] / r_sol->r;
  r_sol->w = crnt_el->S[3] / r_sol->r;

  r_sol->vec[0] = r_sol->r;
  r_sol->vec[1] = crnt_el->S[1];
  r_sol->vec[2] = crnt_el->S[2];
  r_sol->vec[3] = crnt_el->S[3];

  if(run->con->model==1){
    r_sol->e          = crnt_el->S[4];
    r_sol->vec[4]     = crnt_el->S[4];
    r_sol->e_internal = (r_sol->e/r_sol->r) - 0.5*(pow(r_sol->u,2.0) + pow(r_sol->v,2.0) + pow(r_sol->w,2.0));
    r_sol->p_hydro    = r_sol->r*(run->par->matergama[0]-1.0)*r_sol->e_internal - run->par->matergama[0]*run->par->materpinf[0];
    r_sol->c          = sqrt(run->par->matergama[0]*(r_sol->p_hydro + run->par->materpinf[0])/r_sol->r);
  } else if(run->con->model==2) {
    barotropic_eos(crnt,0,run);
  } else if(run->con->model==3) {
    r_sol->ymass      = max(0,crnt_el->S[4]/r_sol->r);
    r_sol->vec[4]     =  max(0,crnt_el->S[4]);
    barotropic_multiphase_eos(crnt,0,run);
  }else if(run->con->model==4) {
    r_sol->ymass      = max(0,crnt_el->S[4]/r_sol->r);
    r_sol->vec[4]     =  max(0,crnt_el->S[4]);
    barotropic_threephase_eos(crnt,0,run);
  }
  
  r_sol->st[0][0] = r_sol->p_hydro; 
  r_sol->st[0][1] = 0.0;
  r_sol->st[0][2] = 0.0;
  
  r_sol->st[1][0] = 0.0;
  r_sol->st[1][1] = r_sol->p_hydro;
  r_sol->st[1][2] = 0.0;
  
  r_sol->st[2][0] = 0.0;
  r_sol->st[2][1] = 0.0;
  r_sol->st[2][2] = r_sol->p_hydro;



  if (run->con->ns==1) {  // Navier-stokes
    r_sol->ux = crnt_el->SG[1][0];
    r_sol->uy = crnt_el->SG[1][1];
    r_sol->uz = crnt_el->SG[1][2];

    r_sol->vx = crnt_el->SG[2][0];
    r_sol->vy = crnt_el->SG[2][1];
    r_sol->vz = crnt_el->SG[2][2];
    
    r_sol->wx = crnt_el->SG[3][0];
    r_sol->wy = crnt_el->SG[3][1];
    r_sol->wz = crnt_el->SG[3][2];
  }
 // printf("rho:sol1%f\n",run->sol->r);
}
