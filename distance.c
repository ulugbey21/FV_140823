#include "strdata.h"

void distance(struct RUN *run,struct BRANCH *crnt,int f,double* dist_cf,double* dist_cc) {

  int ind;
  double dist_temp;
  double xc,yc,zc,nd_temp;
  double xc_neig,yc_neig,zc_neig;
  struct BRANCH * neig_temp;
  struct CELL * cl_temp;
  int i,j,k;
  double **d;
  double **dt;
  double G[3][3];



/*
  d=malloc(4 * sizeof(double*));
    for(ind=0; ind<f; ind++) {
    d[ind]=malloc(*sizeof(double*));
  }
    dt=malloc(4 * sizeof(double*));
   
     for(ind=0; ind<f; ind++) {
    dt[ind]=malloc(4*sizeof(double*));
  }

  */
  dist_temp=0.0;
  *dist_cf=0.0;
  *dist_cc=0.0;

  neig_temp=crnt->neigtr[f][0];

  if (neig_temp->level-crnt->level==0) { // Same level   

    cl_temp=crnt->cl;
    nd_temp=((double)crnt->nlnd);

    xc=0.0;yc=0.0;zc=0.0;
    for (ind=0;ind<crnt->nlnd;ind++){
      xc+=cl_temp->nd[ind].x/nd_temp;
      yc+=cl_temp->nd[ind].y/nd_temp;
      zc+=cl_temp->nd[ind].z/nd_temp;
    }

    cl_temp=neig_temp->cl;
    nd_temp=((double)neig_temp->nlnd);

    xc_neig=0.0;yc_neig=0.0;zc_neig=0.0;
    for (ind=0;ind<neig_temp->nlnd;ind++){
      xc_neig+=cl_temp->nd[ind].x/nd_temp;
      yc_neig+=cl_temp->nd[ind].y/nd_temp;
      zc_neig+=cl_temp->nd[ind].z/nd_temp;
    }
//distance matrix calculation
/*
    for (ind=0;ind<f;ind++){
      d[ind][0]=xc_neig-xc;
      d[ind][1]=yc_neig-yc;
      d[ind][2]=zc_neig-zc;
    }

    for (ind=0;ind<f;ind++){
      dt[0][ind]= d[ind][0];
      dt[1][ind]= d[ind][1];
      dt[2][ind]= d[ind][2];
      
   }
    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
        for (k=0;k<f;k++){ 
              G[i][j]+= dt[i][k]*dt[j][k];
          }
      }
   }
    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
                    printf("%d\n",G[i][j]);            
          }
     }
    */

    dist_temp=0.5*sqrt( pow((xc-xc_neig),2.0) + pow((yc-yc_neig),2.0) + pow((zc-zc_neig),2.0));
    *dist_cf=0.5*dist_temp;
    *dist_cc=dist_temp;
  }

  if (neig_temp->level-crnt->level==1) { // Neigboring element has been refined

    cl_temp=crnt->cl;
    nd_temp=((double)crnt->nlnd);

    xc=0.0;yc=0.0;zc=0.0;
    for (ind=0;ind<crnt->nlnd;ind++){
      xc+=cl_temp->nd[ind].x/nd_temp;
      yc+=cl_temp->nd[ind].y/nd_temp;
      zc+=cl_temp->nd[ind].z/nd_temp;
    }

    cl_temp=neig_temp->prnt->cl;
    nd_temp=((double)neig_temp->prnt->nlnd);

    xc_neig=0.0;yc_neig=0.0;zc_neig=0.0;
    for (ind=0;ind<neig_temp->prnt->nlnd;ind++){
      xc_neig+=cl_temp->nd[ind].x/nd_temp;
      yc_neig+=cl_temp->nd[ind].y/nd_temp;
      zc_neig+=cl_temp->nd[ind].z/nd_temp;
    }

    dist_temp=0.5*sqrt( pow((xc-xc_neig),2.0) + pow((yc-yc_neig),2.0) + pow((zc-zc_neig),2.0));
    *dist_cf=0.5*dist_temp;
    *dist_cc=0.5*dist_temp+0.25*dist_temp; 

  }

      
  if ((neig_temp->level-crnt->level==-1) && (neig_temp->level!=0)) { // Current element has been refined

    cl_temp=crnt->prnt->cl;
    nd_temp=((double)crnt->prnt->nlnd);

    xc=0.0;yc=0.0;zc=0.0;
    for (ind=0;ind<crnt->prnt->nlnd;ind++){
      xc+=cl_temp->nd[ind].x/nd_temp;
      yc+=cl_temp->nd[ind].y/nd_temp;
      zc+=cl_temp->nd[ind].z/nd_temp;
    }

    cl_temp=neig_temp->cl;
    nd_temp=((double)neig_temp->nlnd);

    xc_neig=0.0;yc_neig=0.0;zc_neig=0.0;
    for (ind=0;ind<neig_temp->nlnd;ind++){
      xc_neig+=cl_temp->nd[ind].x/nd_temp;
      yc_neig+=cl_temp->nd[ind].y/nd_temp;
      zc_neig+=cl_temp->nd[ind].z/nd_temp;
    }

    dist_temp=0.25*sqrt( pow((xc-xc_neig),2.0) + pow((yc-yc_neig),2.0) + pow((zc-zc_neig),2.0));
    *dist_cf=0.25*dist_temp;
    *dist_cc=0.25*dist_temp+0.5*dist_temp;
  } 
  
  return;

}