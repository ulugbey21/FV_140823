#include "strdata.h"

void pvrs_mc(struct BRANCH * crnt, struct RUN *run, int ifc) {

  int iv,i,j;
  double rK,uK,vK,wK,ymassK;
  double nx,ny,nz,mx,my,mz,lx,ly,lz;
  double u_star_global,v_star_global,w_star_global;

  double u_temp,v_temp,w_temp;
  double cL,cR,pL,pR,rL,rR,imp_L,imp_R;
  double p_inc,p_comp;
  double Mach,beta,a_const;

  double uL,vL,wL;
  double uR,vR,wR;
  double u_star,p_star;

  double * tensor_temp; 

  tensor_temp = malloc(9*sizeof(double));
  
  nx=crnt->cl->nx[ifc]; mx=crnt->cl->nxt1[ifc]; lx=crnt->cl->nxt2[ifc];
  ny=crnt->cl->ny[ifc]; my=crnt->cl->nyt1[ifc]; ly=crnt->cl->nyt2[ifc];
  nz=crnt->cl->nz[ifc]; mz=crnt->cl->nzt1[ifc]; lz=crnt->cl->nzt2[ifc];

  // ================================================================================
  // mc_starregion

  // =============================================================
  // Rotate sol->u,v,w from global coordinates to local
  u_temp=run->sol->u;
  v_temp=run->sol->v;
  w_temp=run->sol->w;

  uL = u_temp*nx + v_temp*ny + w_temp*nz;
  vL = u_temp*mx + v_temp*my + w_temp*mz;
  wL = u_temp*lx + v_temp*ly + w_temp*lz;

  // =============================================================
  // Rotate sol1->u,v,w from global coordinates to local
  u_temp=run->sol1->u;
  v_temp=run->sol1->v;
  w_temp=run->sol1->w;

  uR = u_temp*nx + v_temp*ny + w_temp*nz;
  vR = u_temp*mx + v_temp*my + w_temp*mz;
  wR = u_temp*lx + v_temp*ly + w_temp*lz;

  // =============================================================
  rL = run->sol->r;
  rR = run->sol1->r;
  
  cL  = run->sol->c;
  cR  = run->sol1->c;

  imp_L = (cL)*(rL);
  imp_R = (cR)*(rR);

  // =============================================================
  // Rotate stress tensors from global coordinates to local
  rotate_tensor(tensor_temp,run->sol->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                         crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                         crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                         run->con->verbose,run->con->rank);
  pL = tensor_temp[0];

  rotate_tensor(tensor_temp,run->sol1->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                          crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                          crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                          run->con->verbose,run->con->rank);
  pR = tensor_temp[0];
  // =============================================================

  u_star  =  (imp_L*uL + imp_R*uR + (pL-pR)) / (imp_L+imp_R);

  Mach    = max((abs(uL)/(cL)),(abs(uR)/(cR)));
  a_const = 10.0;
  beta    =  (1.0-exp(-a_const*Mach));
  
  p_inc    =  (imp_L*pR + imp_R*pL) / (imp_L + imp_R);
  p_comp   =  (imp_L*pR + imp_R*pL + imp_R*imp_L*(uL-uR)) / (imp_L+imp_R);
  
  p_star  =  (1.0-beta)*p_inc + beta*p_comp;
  
  if (isnan(u_star)==1){ 
    printf(" \n");
    printf(" u_star NaN: %f \n",u_star);
    printf(" imp_L,imp_R %f %f \n",imp_L,imp_R);
    printf(" uL,uR %f %f \n",uL,uR);
    printf(" pL,pR %f %f \n",pL,pR);
    printf(" rL,rR %f %f \n",rL,rR);
    printf(" \n");
    exit(0);
  }

  if (isnan(p_star)==1){ 
    printf(" \n");
    printf(" p_star NaN, %f \n",p_star);
    printf(" p_inc,p_comp %f %f \n",p_inc,p_comp); 
    printf(" Mach,beta %f %f \n",Mach,beta);
    printf(" imp_L,imp_R %f %f \n",imp_L,imp_R);
    printf(" uL,uR %f %f \n",uL,uR);
    printf(" pL,pR %f %f \n",pL,pR);
    printf(" rL,rR %f %f \n",rL,rR);
    printf(" \n");
    exit(0);

  }
  // ================================================================================

  // ================================================================================
  // fluc_mc
  if(u_star>=0.0){
    rK     = run->sol->r;
    ymassK = run->sol->ymass;
    uK = uL; 
    vK = vL; 
    wK = wL;
   } else {
    rK     = run->sol1->r;
    ymassK = run->sol1->ymass;
    uK = uL; 
    vK = vL; 
    wK = wL;
  }
  
  u_temp = u_star; //local coordinate
  v_temp = 0.0;
  w_temp = 0.0;
  
  // From local to global
  u_star_global = u_temp*nx + v_temp*mx + w_temp*lx;
  v_star_global = u_temp*ny + v_temp*my + w_temp*ly;
  w_star_global = u_temp*nz + v_temp*mz + w_temp*lz;
  
  run->flux_mc[0][0] = u_star_global *rK;
  run->flux_mc[0][1] = v_star_global *rK;
  run->flux_mc[0][2] = w_star_global *rK;

  run->flux_mc[1][0] = u_star_global *rK* uK + p_star;
  run->flux_mc[1][1] = v_star_global *rK* uK + 0.0;
  run->flux_mc[1][2] = w_star_global *rK* uK + 0.0;

  run->flux_mc[2][0] = u_star_global *rK* vK + 0.0;
  run->flux_mc[2][1] = v_star_global *rK* vK + p_star;
  run->flux_mc[2][2] = w_star_global *rK* vK + 0.0;

  run->flux_mc[3][0] = u_star_global *rK* wK + 0.0;
  run->flux_mc[3][1] = v_star_global *rK* wK + 0.0;
  run->flux_mc[3][2] = w_star_global *rK* wK + p_star;

  run->flux_mc[4][0] = u_star_global *rK* ymassK;
  run->flux_mc[4][1] = v_star_global *rK* ymassK;
  run->flux_mc[4][2] = w_star_global *rK* ymassK;

  for (i=0;i<5;i++){
    for (j=0;j<3;j++){
      if ((isnan(run->flux_mc[i][j])==1)){ 
         printf(" \n");
         printf("flux_mc NaN, %d,%d \n",i,j);
         printf(" \n");
         exit(0);
      }
    }
  }

  // ================================================================================

  free(tensor_temp);

}