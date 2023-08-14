//-------------------------------------------------------
// Viscous Fluxes, page 16 Toro's book
//-------------------------------------------------------

#include "strdata.h"

void flux_viscous(struct RUN *run,struct BRANCH * crnt,struct SOL* sol,struct SOL* sol1) {
 // struct BRANCH * brch;

  crnt=run->topo->locl; 
 // struct SOL * r_sol;
  //r_sol = run->sol;
  double temp;

  double ux,uy,uz;
  double vx,vy,vz;
  double wy,wx,wz;


	double txx,txy,txz;
	double tyx,tyy,tyz;
	double tzx,tzy,tzz;

  double uf,vf,wf;
  double g_ij[3][3];
  double S_ij[3][3];

  double S_ij_d_2;
  double dist;
  int i,j,k;
  int f;
  double delta_S,C_w;
  double *dist_cf,*dist_cc;
  double g_ij_sum[3][3],g_ji_sum[3][3];
  double S_ij_d[3][3],S_ij_2;
  //double mu_turb;
  


  //
  /*
// velocites on the face
  uf = sol->u;
  vf = sol->v;
  wf = sol->w;
  // velocity gradients on the face
  ux = run->sol->ux;
  uy = run->sol->uy;
  uz = run->sol->uz;  // 0
  vx = run->sol->vx;
  vy = run->sol->vy;
  vz = run->sol->vz;  // 0
  wx = run->sol->wx;  // 0
  wy = run->sol->wy;  // 0
  wz = run->sol->wz;  // 0
  */
  // velocites on the face
  uf = 0.5*(sol->u + sol1->u);
  vf = 0.5*(sol->v + sol1->v);
  wf = 0.5*(sol->w + sol1->w);
  
  // velocity gradients on the face
  ux = 0.5*(run->sol->ux + run->sol1->ux);
  uy = 0.5*(run->sol->uy + run->sol1->uy);
  uz = 0.5*(run->sol->uz + run->sol1->uz);

  vx = 0.5*(run->sol->vx + run->sol1->vx);
  vy = 0.5*(run->sol->vy + run->sol1->vy);
  vz = 0.5*(run->sol->vz + run->sol1->vz);

  wx = 0.5*(run->sol->wx + run->sol1->wx);
  wy = 0.5*(run->sol->wy + run->sol1->wy);
  wz = 0.5*(run->sol->wz + run->sol1->wz);

  
  
if (run->con->turb==1){
  //turbulence viscosity calculations 
  S_ij[0][0] = 0.5*(ux + ux);
  S_ij[0][1]= 0.5*(uy + vx);
  S_ij[0][2]= 0.5*(uz + wx);

  S_ij[1][0]  = 0.5*(vx + uy);
  S_ij[1][1] = 0.5*(vy + vy);
  S_ij[1][2] = 0.5*(vz + wy);

  S_ij[2][0] = 0.5*(wx + uz);
  S_ij[2][1]= 0.5*(wy + vz);
  S_ij[2][2] = 0.5*(wz + wz);

   g_ij[0][0]=ux;
   g_ij[0][1]=uy;
   g_ij[0][2]=uz;
   g_ij[1][0]=vx;
   g_ij[1][1]=vy;
   g_ij[1][2]=vz;
   g_ij[2][0]=wx;
   g_ij[2][1]=wy;
   g_ij[2][2]=wz;

   for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++) {
            g_ij_sum[i][j] = 0.0;
          }
   }
   /*
      for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++) {
            g_ji_sum[i][j] = 0.0;
          }
   }
*/
  for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
          for (k = 0; k< 3; k++) {
            g_ij_sum[i][j] += g_ij[i][k] * g_ij[k][j];
          }
        }
  }
/*
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
          g_ji_sum[j][i] =  g_ij_sum[i][j];
    }
  }
  
*/
  S_ij_d[0][0]= ((g_ij_sum[0][0]+g_ij_sum[0][0])/2.0)-(g_ij_sum[0][0]+g_ij_sum[1][1]+g_ij_sum[2][2])/3.0;
  S_ij_d[0][1]= ((g_ij_sum[0][1]+g_ij_sum[1][0])/2.0)-(g_ij_sum[0][0]+g_ij_sum[1][1]+g_ij_sum[2][2])/3.0; 
  S_ij_d[0][2]= ((g_ij_sum[0][2]+g_ij_sum[2][0])/2.0)-(g_ij_sum[0][0]+g_ij_sum[1][1]+g_ij_sum[2][2])/3.0; 
  S_ij_d[1][0]= ((g_ij_sum[1][0]+g_ij_sum[0][1])/2.0)-(g_ij_sum[0][0]+g_ij_sum[1][1]+g_ij_sum[2][2])/3.0;
  S_ij_d[1][1]= ((g_ij_sum[1][1]+g_ij_sum[1][1])/2.0)-(g_ij_sum[0][0]+g_ij_sum[1][1]+g_ij_sum[2][2])/3.0;
  S_ij_d[1][2]= ((g_ij_sum[1][2]+g_ij_sum[2][1])/2.0)-(g_ij_sum[0][0]+g_ij_sum[1][1]+g_ij_sum[2][2])/3.0;
  S_ij_d[2][0]= ((g_ij_sum[2][0]+g_ij_sum[0][2])/2.0)-(g_ij_sum[0][0]+g_ij_sum[1][1]+g_ij_sum[2][2])/3.0;
  S_ij_d[2][1]= ((g_ij_sum[2][1]+g_ij_sum[1][2])/2.0)-(g_ij_sum[0][0]+g_ij_sum[1][1]+g_ij_sum[2][2])/3.0;
  S_ij_d[2][2]= ((g_ij_sum[2][2]+g_ij_sum[2][2])/2.0)-(g_ij_sum[0][0]+g_ij_sum[1][1]+g_ij_sum[2][2])/3.0;


for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        S_ij_d_2 += S_ij_d[i][j] * S_ij_d[i][j];
   //     printf("S_ij_2:%f,i:%d,j:%d\n",S_ij_d_2,i,j);
    }
}
for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        S_ij_2 += S_ij[i][j] * S_ij[i][j];
   //     printf("S_ij_2:%f,i:%d,j:%d\n",S_ij_2,i,j);
    }
}

    
  
   C_w=0.325;
  // printf("dist:%f\n",&dist);
   delta_S=C_w*pow(crnt->cl->Vol,0.3333);//*dist;
  // printf("dist:%f\n",&dist);
 if((pow(S_ij_2,2.5)+pow(S_ij_d_2,1.25))>1e-3){
  run->sol->mu_turb=((run->sol->r*pow(delta_S,2.0)*pow(S_ij_d_2,1.5))/((pow(S_ij_2,2.5)+pow(S_ij_d_2,1.25))));
//  printf("mu_turb:%f",mu_turb);
 }
 else {
 run->sol->mu_turb=0.0;
 }
}

 
 /*
 if(( crnt->cl->nd.y<4.593548831E-03)&&(crnt->cl->nd.y>4.493548831E-03)&&(crnt->cl->nd.x> 9.701266289E-01)&&(crnt->cl->nd.x< 9.711266289E-01)){ 
   printf("mu_turb:%f,density:%f,delta_S:%f,S_ij_d_2:%f,S_ij_2:%f,ux:%f,vy:%f,wz:%f\n",mu_turb,run->sol->r,delta_S,S_ij_d_2,S_ij_2,ux,vy,wz);
 }

}
 */
else if (run->con->turb==0){
  run->sol->mu_turb=0.0;
}
//if (run->sol>mu_turb>1e-4){printf("mu_turb:%f\n",run->sol>mu_turb);}

//printf("mu_tr:%f",run->sol->mu_turb);
  
/*
  temp = (2.0/3.0)*run->par->matervisc;

	txx = temp*( 2.0*ux - (vy+wz) );
	tyy = temp*( 2.0*vy - (wz+ux) );
	tzz = temp*( 2.0*wz - (ux+vy) );

	txy = run->par->matervisc * (uy+vx);
	txz = run->par->matervisc * (vz+wy);
	tyz = run->par->matervisc * (wx+uz);

	tyx = txy;
	tzy = tyz;
	tzx = txz;
*/
//Hirsch's Book page 601
 // turbulence(run,crnt,sol,sol1,mu_turb) ;
	txx = (run->par->matervisc[0]+run->sol->mu_turb)*( 2.0*ux - (2.0/3.0)* (ux+vy+wz) );
	tyy = (run->par->matervisc[0]+run->sol->mu_turb)*( 2.0*vy - (2.0/3.0)* (ux+vy+wz) );
	tzz = (run->par->matervisc[0]+run->sol->mu_turb)*( 2.0*wz - (2.0/3.0)* (ux+vy+wz) );

	txy = (run->par->matervisc[0]+run->sol->mu_turb)* (uy+vx);
	txz = (run->par->matervisc[0]+run->sol->mu_turb)*(uz+wx);
	tyz = (run->par->matervisc[0]+run->sol->mu_turb)*(vz+wy);

	tyx = txy;
	tzy = tyz;
	tzx = txz;
  // Continuity 
  run->flux_viscous[0][0] = 0.0; 
  run->flux_viscous[0][1] = 0.0;
  run->flux_viscous[0][2] = 0.0;
  
  // Momentum
  run->flux_viscous[1][0] = txx; // F
  run->flux_viscous[1][1] = tyx; // G
  run->flux_viscous[1][2] = tzx; // H
  
  run->flux_viscous[2][0] = txy; // 
  run->flux_viscous[2][1] = tyy; //
  run->flux_viscous[2][2] = tzy; // 

  run->flux_viscous[3][0] = txz; // 
  run->flux_viscous[3][1] = tyz; // 
  run->flux_viscous[3][2] = tzz; // 
  if(run->con->eos==1){
  // Total Energy
  run->flux_viscous[4][0] = uf*txx + vf*txy + wf*txz; 
  run->flux_viscous[4][1] = uf*tyx + vf*tyy + wf*tyz;  
  run->flux_viscous[4][2] = uf*tzx + vf*tzy + wf*tzz;  
  }
crnt=crnt->lnxt;

}
