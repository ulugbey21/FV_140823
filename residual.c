//-------------------------------------------------------
//      Computes Right Hand Side of the equations
//-------------------------------------------------------

#include "strdata.h"

void residual(struct BRANCH * crnt, struct RUN * run) {

  int iv,v,i,j,f;
  int ing,flag_temp,f_oppo,ifc_opp,ineig,Nneig;
  int flag_con;
  double u_star;
  double nx,ny,nz,Area_f,as,temp,c;
  double ** Area;
  struct BRANCH * oppo;
  
  Area=malloc(crnt->nlfc*sizeof(double *));
  for(f=0; f<(crnt->nlfc); f++) {
    Area[f]=malloc(4*sizeof(double));
  }
  
  for(iv=0;iv<run->con->neq;++iv){ 
    crnt->el->RHS[iv] = 0.0;
  }

	Nneig=1;
  for(f=0; f<(crnt->nlfc); f++) {  // Loop for faces  
    
    Area_f = 0.0;
    for (ineig=0;ineig<Nneig;++ineig) {
      Area[f][ineig] = crnt->cl->Area[f][ineig];
      Area_f += crnt->cl->Area[f][ineig];
    }
    
			ineig=0;
			if (run->con->ischeme==1){ facevalues(crnt,ineig,f,0,0,run);    } // sol
			if (run->con->ischeme==2){ facercnstr(crnt,ineig,f,0,0,run); } // sol

			if (run->con->ischeme==1){
				if(crnt->cl->fc[f].bc==0){ 
					facevalues(crnt,ineig,f,1,1,run);
				}
				
				if(crnt->cl->fc[f].bc!=0){
					facevalues(crnt,ineig,f,0,1,run);
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
      if ((run->con->model==1) || (run->con->model==2)){
			hllc_euler(crnt,run,ineig,f); 
			for(iv=0;iv<run->con->neq;++iv){  
					crnt->el->RHS[iv] += -(Area[f][ineig])* 
																			((run->flux_adv[iv][0] * crnt->cl->nx[f]  +
																				run->flux_adv[iv][1] * crnt->cl->ny[f]  +
																				run->flux_adv[iv][2] * crnt->cl->nz[f]) + 
																				run->flux_HLLC[iv]);  
       }    
			}
       /*printf("adasd");
      if (run->con->model==3){
       
        for(iv=0;iv<run->con->neq;++iv){  
        //mc_starregion(run,crnt,&u_star,&p_star);
				crnt->el->RHS[iv] += -(Area[f][ineig])* 
																			((run->flux_adv[iv][0] * crnt->cl->nx[f]   +run->flux_mc[iv][0]+
																				(run->flux_adv[iv][1]* crnt->cl->ny[f]   +run->flux_mc[iv][1]+
																				(run->flux_adv[iv][2]* crnt->cl->nz[f]   +run->flux_mc[iv][2]))));  
      
      }
    }

*/
  }
  
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
