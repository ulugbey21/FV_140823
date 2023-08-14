#include "strdata.h"

void hllc_euler (struct BRANCH * crnt,RUN *run,int in,int ifc) {

  int i,iv,flags_if;
  double S_S,SL,SR,uL,vL,wL,uR,vR,wR;
  double s11L,s11R;

  compute_wave_speeds(run,crnt,ifc,&S_S,&SL,&SR,&uL,&vL,&wL,&s11L,
                                                &uR,&vR,&wR,&s11R);
    
  flags_if=0;

  for(iv=0;iv<run->con->neq;++iv){ run->flux_HLLC[iv] = 0.0; }

  if (crnt->cl->fc[ifc].bc==3){S_S=0.0;}

  if ((0.0<=SL) && (flags_if==0)) { // Left state
    flags_if=1;
    flux_adv(run->sol,run);                                    
  }  
  
  if (( (SL<=0.0)&&(0.0<=S_S) ) && (flags_if==0)) { // Left star state  
    flags_if=2;

    star_region_euler(run,crnt,run->sol,ifc,1,S_S,SL,SR,uL,vL,wL,uR,vR,wR);
    flux_adv(run->sol,run);          

    for(iv=0;iv<run->con->neq;++iv){
      run->flux_HLLC[iv] = SL*(run->Us_temp[iv]-run->sol->vec[iv]); 
    }

  } 
  
  if (( (S_S<=0.0)&&(0.0<=SR) ) && (flags_if==0)) { // Right star state
    flags_if=3;

    star_region_euler(run,crnt,run->sol1,ifc,2,S_S,SL,SR,uL,vL,wL,uR,vR,wR);
    flux_adv(run->sol1,run);          

    for(iv=0;iv<run->con->neq;++iv){
      run->flux_HLLC[iv] = SR*(run->Us_temp[iv]-run->sol1->vec[iv]); 
    }
        
  } 
  
  if ((SR<=0.0)&& (flags_if==0)) {  // Right state
    flags_if=4;
    flux_adv(run->sol1,run);                       
  }
  
  if (flags_if==0){
    printf("stop if %e %e %e \n",S_S,SL,SR);
    exit(0);
  }
  
  /*
  // ======================================================================================
  if (crnt->cl->fc[ifc].bc==3){ // Symmetry
    facevalues (crnt,in,ifc,0,0,run);
    facevalues (crnt,in,ifc,0,1,run);
    bc(crnt,crnt->cl->fc[ifc].bc,ifc,run);
   
    flux_F1(run->sol,run);                       
    flux_F1(run->sol1,run);                        

    for(iv=0;iv<run->con->neq;iv++){       // 0->4
      run->flux_HLLC[iv] =  0.0;
      for(i=0;i<3;i++){ run->flux_adv[iv][i] = 0.0;}
    }

    for(iv=1;iv<3;iv++){       // 0->4
     for(i=0;i<3;i++){ run->flux_adv[iv][i] = 0.5*(run->sol1->flux[iv][i] + run->sol->flux[iv][i]); } 
    }

    //run->flux_adv[1][0] = run->sol->flux[1][0];
    //run->flux_adv[2][1] = run->sol->flux[2][1];

  }
  // ======================================================================================
  */

  for(iv=0;iv<run->con->neq;++iv){
    if ((isnan(run->flux_HLLC[iv])==1)){ 
      printf("HLLC nan %d \n",iv);
      exit(0);
    }
  }

}