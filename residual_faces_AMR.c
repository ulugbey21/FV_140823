//-------------------------------------------------------
//      Computes Right Hand Side of the equations
//-------------------------------------------------------
//-------------------------------------------------------
//      Computes Right Hand Side of the equations
//-------------------------------------------------------

#include "strdata.h"

void residual_faces_AMR(struct BRANCH * crnt, struct RUN * run) {

  int iv,v,i,j,f;
  int ing,flag_temp,f_oppo,ifc_opp,ineig,Nneig;
  int flag_con;
  double A1,A2,A3,B1,B2,B3,C1,C2,C3;
  double u_nc,v_nc,w_nc;
  double nx,ny,nz,Area_f,as,temp;
  double u_star,p_star,cL,cR,rL,rR;
  double ** Area;
  int flag_compute;
  struct BRANCH * oppo;
  
  Area=malloc(crnt->nlfc*sizeof(double *));
  for(f=0; f<(crnt->nlfc); f++) {
    Area[f]=malloc(4*sizeof(double));
  }
  

  for(iv=0;iv<run->con->neq;++iv){ 
    crnt->el->RHS[iv] = 0.0;
    run->flux_viscous[iv][0] = 0.0;
    run->flux_viscous[iv][1] = 0.0;
    run->flux_viscous[iv][2] = 0.0;
  }

 
  for(f=0; f<(crnt->nlfc); f++) {  // Loop for faces  

     
    if (crnt->cl->fc[f].bc==0) {        // Face is not boundary

      if (crnt->lfc_neigs[f]==1) {  // 1 neighbor 

        ing=0;
        Nneig=1;
        oppo=crnt->neigtr[f][ing];
        ifc_opp=crnt->neigfc[f];

        if (oppo->lfc_neigs[ifc_opp]==1) { 	// 1:1  connectivity
          flag_con=1;
          if (oppo->part==crnt->part) {  // Opposite part in current processor and its not interface cell
            
            if (crnt->el->flux_flag[f]==1) {  // Flux has been computed from the opposite face                                                  
              f_oppo=crnt->neigfc[f];
   
              for(iv=0;iv<run->con->neq;++iv){    
             //   crnt->el->flux_face[f][iv] = -oppo->el->flux_face[f_oppo][iv];
              }
              flag_compute=1;                                
              flag_temp=1;
            } else {
              f_oppo=crnt->neigfc[f];
              oppo->el->flux_flag[f_oppo]=1;    // Flag for flux calculation 
              flag_compute=1;  
            }
            
          } else {
            flag_compute=1; 
          }
        
        } else { // 2:1 connectivity
          flag_con=2;
          flag_compute=1;                             
        }
         
      }

      if (crnt->lfc_neigs[f]!=1) {  // 1:2 connectivity
        flag_con=3;
        Nneig=crnt->lfc_neigs[f];
        flag_compute=1;
      }
        
    } else { // Face is boundary
      Nneig=1;
      flag_compute=1;
    }
    
    Area_f = 0.0;
    for (ineig=0;ineig<Nneig;++ineig) {
      Area[f][ineig] = crnt->cl->Area[f][ineig];
      Area_f += crnt->cl->Area[f][ineig];
    }

    if (flag_compute==1) {
      
      for(iv=0;iv<run->con->neq;++iv){  
        crnt->el->flux_face[f][iv] = 0.0;
      } 

      ineig=0;
      if (run->con->ischeme==1){ facevalues(crnt,ineig,f,0,0,run);    } // sol
      if (run->con->ischeme==2){ facercnstr(crnt,ineig,f,0,0,run); } // sol

      for (ineig=0;ineig<Nneig;++ineig) {
  
        if (run->con->ischeme==1){
          if(crnt->cl->fc[f].bc==0){ 
            facevalues(crnt,ineig,f,1,1,run);
          }
          
          if(crnt->cl->fc[f].bc!=0){
            facevalues(crnt,ineig,f,0,1,run); //sol1
            bc(crnt,crnt->cl->fc[f].bc,f,run);
          }
        }

        if (run->con->ischeme==2){
          if(crnt->cl->fc[f].bc==0){
            facercnstr(crnt,ineig,f,1,1,run);  // sol1
          }

          if(crnt->cl->fc[f].bc!=0){
            facercnstr(crnt,ineig,f,0,1,run);  // sol1
            bc(crnt,crnt->cl->fc[f].bc,f,run);    // sol1
          }
        }

        if (run->con->ns==1) {flux_viscous(run,crnt,run->sol,run->sol1);}
        

        if ((run->con->model==1) || (run->con->model==2 || run->con->model==3)){
    
          hllc_euler(crnt,run,ineig,f); 
          for(iv=0;iv<run->con->neq;++iv){  
        
            /*
              crnt->el->flux_face[f][iv] += -(Area[f][ineig])* 
                                        ((run->flux_temp[iv][0] * crnt->cl->nx[f]  +
                                          run->flux_temp[iv][1] * crnt->cl->ny[f]  +
                                          run->flux_temp[iv][2] * crnt->cl->nz[f]) + 
                                          run->flux_HLLC[iv]); 
            */
        
            crnt->el->flux_face[f][iv] += -(Area[f][ineig])* 
                                         (((run->flux_adv[iv][0] - run->flux_viscous[iv][0]) * crnt->cl->nx[f]  +
                                           (run->flux_adv[iv][1] - run->flux_viscous[iv][1]) * crnt->cl->ny[f]  +
                                           (run->flux_adv[iv][2] - run->flux_viscous[iv][2]) * crnt->cl->nz[f]) + 
                                            run->flux_HLLC[iv]);
          }

        } else if (run->con->model==1) {
          
          pvrs_mc(crnt,run,f);
          for(iv=0;iv<run->con->neq;++iv){
      
            /*
              crnt->el->flux_face[f][iv] += -(Area[f][ineig])* 
              ((run->flux_adv[iv][0] - run->flux_viscous[iv][0]+ run->flux_mc[iv][0])* crnt->cl->nx[f]   +
              (run->flux_adv[iv][1] - run->flux_viscous[iv][1]+ run->flux_mc[iv][1]) * crnt->cl->ny[f]   +
              (run->flux_adv[iv][2] - run->flux_viscous[iv][2]+ run->flux_mc[iv][2]) * crnt->cl->nz[f]);
            */
              
            crnt->el->flux_face[f][iv] += -(Area[f][ineig])* 
                                          (((run->flux_mc[iv][0] * crnt->cl->nx[f]) +
                                            (run->flux_mc[iv][1] * crnt->cl->ny[f]) +
                                            (run->flux_mc[iv][2] * crnt->cl->nz[f])));
          }
        }
      } 
    } // if compute==1

    for(iv=0;iv<run->con->neq;++iv){
      crnt->el->RHS[iv] += crnt->el->flux_face[f][iv]; 
    }

  } // neig loop

  
    
  
  // ======================
  //      Nan - Check
  // ======================
  
  for(v=0;v<run->con->neq;v++){//variables
        
    if ((isnan(crnt->el->RHS[v])==1)){ 
      printf(" \n");
      printf("RHS NaN, eq: %d, element: %d, time-step: %d \n",v,crnt->root,run->con->tstep);
      printf(" \n");
      exit(0);
    }
  }
  
  
  for(f=0; f<(crnt->nlfc); f++) {
    free(Area[f]);
  }
  free(Area);

}







