#include "strdata.h"

void facevalues(struct BRANCH * crnt,int in,int ifc,int iside, int istore,RUN *run){

  int iv,i;
  struct SOL * r_sol;
  struct SOL * r_sol1;
  struct LEAF * temp_el;
    
   if (istore==0) {       // save on sol  

    r_sol = run->sol;
    if (iside==0) {temp_el=crnt->el;}
    if (iside==1) {temp_el=crnt->neigtr[ifc][in]->el;}

  } else if (istore==1) { // save on sol1  

    r_sol  = run->sol1;
    if (iside==0) {temp_el=crnt->el;}
    if (iside==1) {temp_el=crnt->neigtr[ifc][in]->el;}

  }


  r_sol->r = temp_el->S[0];
  r_sol->u = temp_el->S[1] / r_sol->r;
  r_sol->v = temp_el->S[2] / r_sol->r;
  r_sol->w = temp_el->S[3] / r_sol->r;
  
  r_sol->vec[0] = r_sol->r;
  r_sol->vec[1] = temp_el->S[1];
  r_sol->vec[2] = temp_el->S[2];
  r_sol->vec[3] = temp_el->S[3];
  
  
  if(run->con->model==1){
    r_sol->e          = temp_el->S[4];
    r_sol->vec[4]     = temp_el->S[4];
    r_sol->e_internal = (r_sol->e/r_sol->r) - 0.5*(pow(r_sol->u,2.0) + pow(r_sol->v,2.0) + pow(r_sol->w,2.0));
    r_sol->p_hydro    = r_sol->r*(run->par->matergama[0]-1.0)*r_sol->e_internal - run->par->matergama[0]*run->par->materpinf[0];
    r_sol->c = sqrt(run->par->matergama[0]*(r_sol->p_hydro + run->par->materpinf[0])/r_sol->r);
  } else if(run->con->model==2) {
    barotropic_eos(crnt,istore,run);
  } else if(run->con->model==3) {
    r_sol->ymass = temp_el->S[4]/r_sol->r;

    r_sol->vec[4] = temp_el->S[4];
    barotropic_multiphase_eos(crnt,istore,run);
  } else if(run->con->model==4) {
    r_sol->ymass = temp_el->S[4]/r_sol->r;

    r_sol->vec[4] = temp_el->S[4];
    barotropic_threephase_eos(crnt,istore,run);
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
    r_sol->ux = temp_el->SG[1][0];
    r_sol->uy = temp_el->SG[1][1];
    r_sol->uz = temp_el->SG[1][2];

    r_sol->vx = temp_el->SG[2][0];
    r_sol->vy = temp_el->SG[2][1];
    r_sol->vz = temp_el->SG[2][2];
    
    r_sol->wx = temp_el->SG[3][0];
    r_sol->wy = temp_el->SG[3][1];
    r_sol->wz = temp_el->SG[3][2];
  }
  
} 