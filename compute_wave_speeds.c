#include "strdata.h"

void compute_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* S_S,double* SL, double* SR,
                         double* uL, double* vL, double* wL,double* s11L,double* uR, double* vR, double* wR,double* s11R) {
                                                                              

  int iv;
  double ruL,ruR;
  double cL,cR;
  double u,v,w;
  double nx,mx,lx,ny,my,ly,nz,mz,lz;
  double p_star_comp,p_star_inc,b,d,M,p_star,CL,CR,u_star;
  double * tensor_temp; 
  
  
  tensor_temp = malloc(9*sizeof(double));

  
  u=run->sol->u;
  v=run->sol->v;
  w=run->sol->w;

  nx=crnt->cl->nx[ifc]; mx=crnt->cl->nxt1[ifc]; lx=crnt->cl->nxt2[ifc];
  ny=crnt->cl->ny[ifc]; my=crnt->cl->nyt1[ifc]; ly=crnt->cl->nyt2[ifc];
  nz=crnt->cl->nz[ifc]; mz=crnt->cl->nzt1[ifc]; lz=crnt->cl->nzt2[ifc];

  /*
  *uL = u*nx; // + v*ny + w*nz;
  *vL = 0.0; //u*mx + v*my + w*mz;
  *wL = 0.0; //u*lx + v*ly + w*lz;
  */

  *uL = u*nx + v*ny + w*nz;
  *vL = u*mx + v*my + w*mz;
  *wL = u*lx + v*ly + w*lz;

  u=run->sol1->u;
  v=run->sol1->v;
  w=run->sol1->w;

  *uR = u*nx + v*ny + w*nz;
  *vR = u*mx + v*my + w*mz;
  *wR = u*lx + v*ly + w*lz;
	

  ruL = run->sol->r* (*uL);
  ruR = run->sol1->r* (*uR);
    
  // Wave speeds, ok 
  cL = run->sol->c;
  cR = run->sol1->c;

  *SR = max( (*uL+cL) , (*uR+cR) );   // Right wave speed
  *SL = min( (*uL-cL) , (*uR-cR) );   // Left  wave speed 
    
  // Rotation of st tensors 
  rotate_tensor(tensor_temp,run->sol->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                         crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                         crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                         run->con->verbose,run->con->rank);

  *s11L = tensor_temp[0];
	//*s11L = -run->sol->p_hydro;

	//printf(" rotation 1 %e %e \n",tensor_temp[0],-run->sol->p_hydro);
	//exit(0);

  rotate_tensor(tensor_temp,run->sol1->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                          crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                          crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                          run->con->verbose,run->con->rank);
  *s11R = tensor_temp[0];
  //*s11L = -run->sol1->p_hydro;

	//printf(" rotation 1 %e %e \n",tensor_temp[0],-run->sol->p_hydro);
  

 // printf("p_star:%f\n",p_star);
 // printf("u_star:%f\n",u_star);
 //printf("*s11R:%e\n",*s11R);
  *S_S=( (run->sol->r*pow((*uL),2.0) + (*s11L)) - (run->sol1->r*pow((*uR),2.0) + (*s11R)) - ( (*SL)*ruL) + ( (*SR)*ruR)) / 
       (ruL - ruR - ( (*SL)*run->sol->r) + ( (*SR)*run->sol1->r));  

  if (isnan(*S_S)==1) {
    printf("S_S nan | %e | %e %e | %e %e | %e %e | %e %e | %e %e \n",*S_S,*SL,*SR,cL,cR,(*uL),(*uR),(*s11L),(*s11R),run->sol->r,run->sol1->r);
    exit(0);
  }

  free(tensor_temp);
}