#include "strdata.h"

void pvrs(struct RUN *run, struct BRANCH * crnt,int ifc) {
   /*
  double uL,vL,wL;
  double uR,vR,wR; 
  double pL,pR,rL,rR,cL,cR;
  double r_avg,c_avg;
  double p_star,u_star,rL_star,rR_star,cL_star,cR_star,rR_shock_star;
  double *SR,*SL,*SHL,*STL,*SHR,*STR;
  double rL_fan,uL_fan,pL_fan,rR_fan,uR_fan,pR_fan,rR_star_fan;

  //Toro page 299
  uL      = run->sol->u;
  vL      = run->sol->v;
  wL      = run->sol->w;

  uR      = run->sol1->u;
  vR      = run->sol1->v;
  wR      = run->sol1->w;

  pL      = run->sol->p_hydro;
  pR      = run->sol1->p_hydro;

  rL      = run->sol->r;
  rR      = run->sol1->r;

  cL      = run->sol->c;
  cR      = run->sol1->c;
  
  r_avg   = (rL+rR)/2.0;
  c_avg   = (cL+cR)/2.0;

  p_star  = (pL+pR)/2.0+((uL-uR)/2.0)*r_avg*c_avg;
  u_star  = (uL+uR)/2.0+((pL-pR)/2.0)/(r_avg*c_avg);

  //Toro page 137

  rR_shock_star = rR* ((p_star/pR)+run->par->matergama[0]-1)/(run->par->matergama[0]+1)/(((run->par->matergama[0]-1)/(run->par->matergama[0]+1))*(p_star/pR)+1);

  rL_star = rL +(uL-u_star)*(r_avg/c_avg);
  rR_star = rR +(u_star-uR)*(r_avg/c_avg);

  cL_star = cL *pow((p_star/pL),((run->par->matergama[0]-1)/(2*run->par->matergama[0])));
  cR_star = cR *pow((p_star/pR),(run->par->matergama[0]-1)/(2*run->par->matergama[0]));

  rL_fan   = rL*pow((2.0/(run->par->matergama[0]+1.0))+((run->par->matergama[0]-1.0)/((run->par->matergama[0]+1.0)*cL))*(uL-uL),2.0/(run->par->matergama[0]-1));
  uL_fan   = (2.0/(run->par->matergama[0]+1.0))*(cL+(((run->par->matergama[0]-1.0)/2.0)*uL+uL);
  pL_fan   = pL*pow((2.0/(run->par->matergama[0]+1.0))+((run->par->matergama[0]-1.0)/((run->par->matergama[0]+1.0)*cL))*(uL-uL),2.0*run->par->matergama[0]/(run->par->matergama[0]-1));
  

  rR_fan   = rR*pow((2.0/(run->par->matergama[0]+1.0))-((run->par->matergama[0]-1.0)/((run->par->matergama[0]+1.0)*cL))*(uR-uL),2.0/(run->par->matergama[0]-1));
  uR_fan   = (2.0/(run->par->matergama[0]+1.0))*(-cR+(((run->par->matergama[0]-1.0)/2.0)*uR+uL);
  pR_fan   = pR*pow((2.0/(run->par->matergama[0]+1.0))-((run->par->matergama[0]-1.0)/((run->par->matergama[0]+1.0)*cR))*(uR-uL),2.0*run->par->matergama[0]/(run->par->matergama[0]-1));
  
  rR_star_fan = rR*pow((p_star/pR),(1.0/run->par->matergama[0]));
  *SR  = uR+cR*sqrt(((run->par->matergama[0]+1)/(2*run->par->matergama[0]))*(p_star/pR)+(matergama[0]-1)/(2*matergama[0]));   // Right wave speed
  *SL  = uL-cL*sqrt(((run->par->matergama[0]+1)/(2*run->par->matergama[0]))*(p_star/pL)+(matergama[0]-1)/(2*matergama[0]));   // Left  wave speed 
  *SHL = uL-cL;
  *STL = u_star-cL_star;
  *SHR = uR+cR;
  *STR = u_star+cR_star;

  if (uL<u_star){
    if(p_star>pL){//Left Shock
       if(uL<*SL){
        flux_adv(run->sol,run);  
         //!?
       }
       else if(*SL<uL && uL<u_star) {
        run->sol->p_hydro = p_star;
        run->sol->u       = u_star ;
        run->sol->r       = rL_star;
       }
    }
    else { //Left fan
       if (*SL<*SHL){
        flux_adv(run->sol,run);  
         //!?
       }
       else if (*SHL<uL && uL<*STL){
        run->sol->p_hydro = pL_fan;
        run->sol->u       = uL_fan ;
        run->sol->r       = rL_fan;
       }
        else if (*STL<uL && uL<u_star){
        run->sol->p_hydro = p_star;
        run->sol->u       = u_star ;
        run->sol->r       = rL_fan;
       }     
      }

    }

  if (uR>u_star){
    if(p_star>pR){//Right Shock
       if(uR>*SR){
        flux_adv(run->sol,run);  
         //!?
       }
       else if(uR<*SR && u_star<uR) {
        run->sol->p_hydro = p_star;
        run->sol->u       = u_star ;
        run->sol->r       = rR_shock_star;
       }
    }
    else { //right fan
       if (uR>*SHR){
        flux_adv(run->sol,run);  
         //!?
       }
       else if (*STR<uR && uR<*SHR){
        run->sol->p_hydro = pR_fan;
        run->sol->u       = uR_fan ;
        run->sol->r       = rR_fan;
       }
        else if (uR<*STR && u_star<uR){
        run->sol->p_hydro = p_star;
        run->sol->u       = u_star ;
        run->sol->r       = rR_star_fan;
       }     
      }

    }
    */
  }
