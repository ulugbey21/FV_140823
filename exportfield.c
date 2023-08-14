/*
 * Forest of oct-Trees DG implementation
 * octo-Rhodon
 *
 *  \/\/\/\/
 *   \/ / /
 *    \/ /
 *     \/
 *
 * Foundation for Research and Technology Hellas, 
 * Institute for Applied and Computational Mathematics (FORTH/IACM)
 * Heraklion Crete
 * Andreas Papoutsakis 
 * 20/05/2014
*/

#include "strdata.h"

void exportfield (struct RUN * run) {
  
  int k,iv,in,isp,iel,flag_rad;
  int ind;
  int ngeltot;
  int ngndtot;
  int igndtot;
  int ignd;
  int ngnd;
  int ngel;
  int type;
  int owrt;
  int imatv,jmatv;
  int TPNDS[8][8];
  int INX[8],INY[8],INZ[8];
  int INXP[8],INYP[8],INZP[8];
  int INXT[8],INYT[8],INZT[8];
  int f;
  float x,y,z,r;
  double Y[20];
  double var[100];
  int i,j;
  double xc,yc,zc,sc;
  struct BRANCH * brch;
  struct TREE * tree;
  char  fstring[200]="                                             ";
  char  fstring2[200]="                                             ";
  
  FILE * tp;
  FILE * Res;


  TPNDS[0][0]=0;
  TPNDS[1][0]=1;
  TPNDS[2][0]=3;
  TPNDS[3][0]=2;
  TPNDS[4][0]=4;
  TPNDS[5][0]=5;
  TPNDS[6][0]=7;
  TPNDS[7][0]=6;


  TPNDS[0][1]=1;
  TPNDS[1][1]=2;
  TPNDS[2][1]=2;
  TPNDS[3][1]=0;
  TPNDS[4][1]=3;
  TPNDS[5][1]=3;
  TPNDS[6][1]=3;
  TPNDS[7][1]=3;

  TPNDS[0][2]=2;
  TPNDS[1][2]=5;
  TPNDS[2][2]=3;
  TPNDS[3][2]=0;
  TPNDS[4][2]=1;
  TPNDS[5][2]=4;
  TPNDS[6][2]=4;
  TPNDS[7][2]=1;

  INX[0]=1;INY[0]=0;INZ[0]=0;
  INX[1]=1;INY[1]=1;INZ[1]=0;
  INX[2]=0;INY[2]=1;INZ[2]=0;
  INX[3]=0;INY[3]=0;INZ[3]=0;
  INX[4]=1;INY[4]=0;INZ[4]=1;
  INX[5]=1;INY[5]=1;INZ[5]=1;
  INX[6]=0;INY[6]=1;INZ[6]=1;
  INX[7]=0;INY[7]=0;INZ[7]=1;

  INXP[0]=1;INYP[0]=0;INZP[0]=0;
  INXP[1]=1;INYP[1]=1;INZP[1]=0;
  INXP[2]=0;INYP[2]=1;INZP[2]=0;
  INXP[3]=0;INYP[3]=0;INZP[3]=0;
  INXP[4]=0;INYP[4]=0;INZP[4]=1;
  INXP[5]=0;INYP[5]=1;INZP[5]=1;
  INXP[6]=0;INYP[6]=1;INZP[6]=1;
  INXP[7]=0;INYP[7]=0;INZP[7]=1;

  INXT[0]=1;INYT[0]=0;INZT[0]=0;
  INXT[1]=0;INYT[1]=1;INZT[1]=0;
  INXT[2]=0;INYT[2]=1;INZT[2]=0;
  INXT[3]=0;INYT[3]=0;INZT[3]=0;
  INXT[4]=0;INYT[4]=0;INZT[4]=1;
  INXT[5]=0;INYT[5]=0;INZT[5]=1;
  INXT[6]=0;INYT[6]=0;INZT[6]=1;
  INXT[7]=0;INYT[7]=0;INZT[7]=1;

  ngnd=0;
  ngel=0;    

  MPI_Barrier(MPI_COMM_WORLD);
  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
    tree=run->topo->drys[iel];
    if (tree->part==run->con->rank){
    brch=tree->brch;
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }

      while(brch!=NULL&&brch->root==iel){
        ngel++;
        for (ind=0; ind<brch->nlnd; ind++){
          ngnd++;
        }
      brch=brch->lnxt;}


    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Allreduce((&ngel),(&ngeltot),1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce((&ngnd),(&ngndtot),1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  ngel=ngeltot;
  ngnd=ngndtot;
  owrt=1;
  //sprintf(fstring,"field%08dprt%02d.dat", run->con->tstep,0); // Rank always 0

  if (run->con->i_case<=99) {   
    sprintf(fstring2,"Res_gnuplot.dat");   
  }

  if (run->con->rank==0) {
    //if (owrt==1) {tp=fopen(fstring, "w");}
    //if (owrt==0) {tp=fopen(fstring, "a");}

    if (run->con->i_case<=99) {
      if (owrt==1) {Res=fopen(fstring2, "w");}
      if (owrt==0) {Res=fopen(fstring2, "a");}
    }

    // Write premable
    //fprintf(tp," VARIABLES =  \"X\" \"Y\" \"Z\" \"iel\" \"part\" \"rho\" \"U\" \"V\" \"W\" \"e\" \"P\" \"C\" ");
    //fprintf(tp," \n ");
    
    //fprintf(tp, "ZONE T=\"%d\",N=%d,E=%d,F=FEPOINT,ET=BRICK,DT=(double,double,double)\n",1,ngel*8,ngel);
    //fclose(tp);
  }

 

  MPI_Barrier(MPI_COMM_WORLD);
  owrt=0;
  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (tree->part==run->con->rank){
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }

      //if (owrt==1) {tp=fopen(fstring, "w");}
      //if (owrt==0) {tp=fopen(fstring, "a");}
      
      if (run->con->i_case<=99) {
        if (owrt==1) {Res=fopen(fstring2, "w");}
        if (owrt==0) {Res=fopen(fstring2, "a");}
      }


      while(brch!=NULL&&brch->root==iel){

       // Write results
       type=brch->type;

       xc=0.0;yc=0.0;zc=0.0;
        for (ind=0;ind<brch->nlnd;ind++){
          xc+=brch->cl->nd[ind].x/((double)brch->nlnd);
          yc+=brch->cl->nd[ind].y/((double)brch->nlnd);
          zc+=brch->cl->nd[ind].z/((double)brch->nlnd);
        }


       for (in=0;in<8;in++){          // Nodes
         
         x=brch->cl->nd[TPNDS[in][type]].x; 
         y=brch->cl->nd[TPNDS[in][type]].y; 
         z=brch->cl->nd[TPNDS[in][type]].z; // Coordinates

         ind=TPNDS[in][type];                                 // Tecplot to solution vetices
         
          volmvalues(brch,run);        

          //fprintf(tp,"%e %e %e %d %d %e %e %e %e %e %e %e ",x,y,z,brch->root,brch->part,run->sol->r,run->sol->u,run->sol->v,run->sol->w,run->sol->e_internal,run->sol->p_hydro,run->sol->c); // Info

          //fprintf(tp," \n");

          if (run->con->i_case<=99)  {
            if (in==0) {
             // fprintf(Res,"%d ",1000-brch->root);
              fprintf(Res,"%d ",brch->root);
              //fprintf(Res,"%f ",x);
              if (run->con->model==1){
              fprintf(Res,"%e %e %e %e %e %e %e ",run->sol->r,run->sol->u,run->sol->v,run->sol->w,run->sol->e_internal,run->sol->p_hydro,run->sol->c); // Info
              }
              else if (run->con->model==2){
              fprintf(Res,"%e %e %e %e %e %e %e",run->sol->r,run->sol->u,run->sol->v,run->sol->w,run->sol->p_hydro,run->sol->c,run->sol->u/1450.0); // Info
              }
              else if (run->con->model==3 ||run->con->model==4){
              fprintf(Res,"%e %e %e %e %e %e %e ",run->sol->r,run->sol->u,run->sol->v,run->sol->w,run->sol->p_hydro,run->sol->c,pow(pow(run->sol->u,2.0)+pow(run->sol->v,2.0)+pow(run->sol->w,2.0),0.5)); // Info
              }
              fprintf(Res," \n");
              /*
                  for(f=0;f<brch->nlfc;f++){
  if((brch->root==499)&& (crnt->cl->nx[f]>0.6)) { 
    printf("f:%d",f);
    printf("r_sol->c:%f,brch->root:%d\n",run->sol->c,brch->root);
    

  } 
  
  if((brch->root==500)&& (crnt->cl->nx[f]<-0.6)) { 
    printf("f:%d",f);
    printf("r_sol->c:%f,brch->root:%d\n",run->sol->c,brch->root);
    exit(0);

  } 
   }
   */
              
            }
          }

       }
       brch=brch->lnxt;
      } 
      //fclose(tp);
      if (run->con->i_case<=99)  { fclose(Res); }
    } 
    MPI_Barrier(MPI_COMM_WORLD);
  }

  /*
  ignd=0;
  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (tree->part==run->con->rank){
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }

      if (owrt==1) {tp=fopen(fstring, "w");}
      if (owrt==0) {tp=fopen(fstring, "a");}
      while(brch!=NULL&&brch->root==iel){
        for (ind=0;ind<8;ind++){
          ignd++;
          fprintf(tp,"%d ",ignd);
        }
        fprintf(tp,"\n");
        brch=brch->lnxt;
      }
      fclose(tp);
    }
    MPI_Allreduce((&ignd),(&igndtot),1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    ignd=igndtot;
    MPI_Barrier(MPI_COMM_WORLD);
  }
  */
}
