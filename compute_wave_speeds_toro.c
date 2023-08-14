#include "strdata.h"

void compute_wave_speeds_toro(struct RUN *run, struct BRANCH * crnt,int ifc,double* S_S,double* SL, double* SR,
                         double* uL, double* vL, double* wL,double* s11L,double* uR, double* vR, double* wR,double* s11R) {
                                                                              

  int iv;
  double ruL,ruR;
  double cL,cR;
	double rho_aver,a_aver,p_vrs,p_S;
	double qL,qR,gamma;
  double u,v,w;
  double nx,mx,lx,ny,my,ly,nz,mz,lz;
  
  double * tensor_temp; 
  
  
  tensor_temp = malloc(9*sizeof(double));

  
  u=run->sol->u;
  v=run->sol->v;
  w=run->sol->w;

  nx=crnt->cl->nx[ifc]; mx=crnt->cl->nxt1[ifc]; lx=crnt->cl->nxt2[ifc];
  ny=crnt->cl->ny[ifc]; my=crnt->cl->nyt1[ifc]; ly=crnt->cl->nyt2[ifc];
  nz=crnt->cl->nz[ifc]; mz=crnt->cl->nzt1[ifc]; lz=crnt->cl->nzt2[ifc];

  *uL = u*nx; // + v*ny + w*nz;
  *vL = 0.0; //u*mx + v*my + w*mz;
  *wL = 0.0; //u*lx + v*ly + w*lz;
  
  u=run->sol1->u;
  v=run->sol1->v;
  w=run->sol1->w;

  *uR = u*nx; // + v*ny + w*nz;
  *vR = 0.0; //u*mx + v*my + w*mz;
  *wR = 0.0; //u*lx + v*ly + w*lz;

	// Rotation of st tensors 
  rotate_tensor(tensor_temp,run->sol->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                         crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                         crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                         run->con->verbose,run->con->rank);

  *s11L = tensor_temp[0];
  *s11L = run->sol->p_hydro;

  rotate_tensor(tensor_temp,run->sol1->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                          crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                          crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                          run->con->verbose,run->con->rank);
  *s11R = tensor_temp[0];
	*s11L = run->sol->p_hydro;
	//=======================================================


	rho_aver	= 0.5*(run->sol->r + run->sol1->r);
	a_aver		= 0.5*(run->sol->c + run->sol1->c);

	p_vrs = 0.5*((*s11L)+(*s11R))-0.5*(uR-uL)*rho_aver*a_aver;

	p_S = max(0.0,p_vrs);

	gamma=run->par->matergama[0];
	if (p_S<=(*s11L)){
		qL = 1.0;
	} else {
		qL = pow((1.0 + (((gamma+1)/2.0*gamma)*((p_S/(*s11L))-1.0)) ),0.5);
	}

	if (p_S<=(*s11R)){
		qR = 1.0;
	} else {
		qR = pow((1.0 + (((gamma+1)/2.0*gamma)*((p_S/(*s11R))-1.0)) ),0.5);
	}

	*SL = *uL - run->sol->c *qL;   // Left  wave speed 
	*SR = *uR + run->sol1->c*qR;   // Right wave speed
  
	ruL = run->sol->r* (*uL);
  ruR = run->sol1->r* (*uR);

	*S_S=( (*s11R) - (*s11L) + ruL*((*SL)-(*uL)) - ruR*((*SR)-(*uR))) / 
       ( run->sol->r*((*SL)-(*uL)) - run->sol1->r*((*SR)-(*uR)));  

  free(tensor_temp);
}