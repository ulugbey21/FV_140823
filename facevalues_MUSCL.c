#include "strdata.h"

void facevalues_MUSCL(struct BRANCH * crnt,int ifc,int iside, struct SOL* sol_temp,struct RUN *run) {
  int iv,ing,ifc_opp;
  int ineig,Nneig;
  struct BRANCH * oppo;

  if (iside==0) {      //  run->AUX_S from current element   

    for(iv=0;iv<run->con->neq;iv++){
      sol_temp->vec[iv] = crnt->el->S[iv];   //ok
             
    }
      

  }

  if (iside==1) {     //  run->AUX_S from neighboring element

    if (crnt->lfc_neigs[ifc]==1) {        // 2:1 or 1:1 connectivity
     
      ing   = 0;   
      oppo=crnt->neigtr[ifc][ing];
      ifc_opp=crnt->neigfc[ifc];
      for(iv=0;iv<run->con->neq;iv++) {  
        sol_temp->vec[iv] = oppo->el->S[iv];
          //ok
      }
    }

    if (crnt->lfc_neigs[ifc]!=1) {          // 1:2

      for (iv=0;iv<run->con->neq;iv++) {
        sol_temp->vec[iv]=0.0;
 
      }

      Nneig=crnt->lfc_neigs[ifc];         
      
      for (ineig=0;ineig<Nneig;ineig++) {
        ing  = ineig;
        oppo = crnt->neigtr[ifc][ing];
        for(iv=0;iv<run->con->neq;iv++){  
          sol_temp->vec[iv]  += oppo->el->S[iv]/( (double) Nneig);
  
        }
      }

    }    
    
  } 

    
  return;

}