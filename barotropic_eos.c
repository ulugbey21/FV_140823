#include "strdata.h"
#include "math.h"
void barotropic_eos(struct BRANCH * crnt,int istore,RUN *run){

  double nexp, beta,cbar,psat,rsat;
  struct SOL * r_sol;

  nexp=7.15;
  beta=293526643.0;
  cbar=1450.0 ;
  psat=2339.0;
  rsat=998.2;
  
  if (istore==0) {
    r_sol = run->sol;     
  } else if (istore==1) {
    r_sol = run->sol1;  
  }

  if (r_sol->r>=rsat) {
    r_sol->p_hydro = psat+beta*(pow((r_sol->r/rsat),nexp)-1.0);
    r_sol->c       = sqrt(beta*nexp*(pow((r_sol->r/rsat),(nexp-1)))/rsat);
  }
  if (r_sol->r<rsat) {
    r_sol->p_hydro = psat+cbar*((1.0/rsat)-(1.0/r_sol->r));
    r_sol->c       = sqrt(cbar/((r_sol->r)*(r_sol->r)));
  }
  return;
   
}

void barotropic_multiphase_eos(struct BRANCH * crnt,int istore,struct RUN *run){

//liquid-vapour gas mixture, linear EOS 


  int f;
  int flag;

  double R_G       = 286.0;
  double T_ref     = 293.15;
  double r_sat     = 998.1618;
  double p_sat     = 2340.0;
  double c_gas     = 343.24;
  double c_liq     = 1482.5;
  double c_mixture = 1.0;

  double c_temp,r_lm,r_gas,a_gas,a_lm;
  double a_coef,b_coef,d_coef;
  double pres_1,pres_2;
  double xw,yw,zw;
  double R1;
 
  
  double check_NaN;
  
  struct SOL * r_sol;
  struct SOL * r_sol1;



  xw=crnt->cl->xc;
  yw=crnt->cl->yc;
  zw=crnt->cl->zc;

  if (istore==0){ r_sol = run->sol;  }
  if (istore==1){ r_sol = run->sol1; }   
  
  /*
  p_old = 
  if (p_old>=p_sat) {
    c_temp=c_liq;
  } else {
    c_temp=c_mixture;
  }
  */

  c_temp=c_liq;
 
  // quatradic polynomial from pressure
  a_coef = 1.0; // ok
  b_coef = pow(c_temp,2.0)*r_sat - p_sat - r_sol->r*r_sol->ymass*R_G*T_ref - r_sol->r*pow(c_temp,2.0)*(1.0-r_sol->ymass); // ok 
  d_coef = (p_sat-pow(c_temp,2.0)*r_sat)* r_sol->r*r_sol->ymass*R_G*T_ref;  // ok
  
  check_NaN=sqrt(pow(b_coef,2.0)-4.0*d_coef);
  if (isnan(check_NaN)==1) {
    c_temp=c_mixture;
      
    a_coef = 1.0; // ok
    b_coef = pow(c_temp,2.0)*r_sat - p_sat - r_sol->r*r_sol->ymass*R_G*T_ref - r_sol->r*pow(c_temp,2.0)*(1.0-r_sol->ymass); // ok 
    d_coef = (p_sat-pow(c_temp,2.0)*r_sat)* r_sol->r*r_sol->ymass*R_G*T_ref;  // ok
      
    check_NaN=sqrt(pow(b_coef,2.0)-4.0*d_coef);
    if (isnan(check_NaN)==1) {
      printf(" double nan \n");
      exit(0);
    } else {
      pres_1 = (-b_coef+sqrt(pow(b_coef,2.0)-4.0*d_coef))/2.0;
      pres_2 = (-b_coef-sqrt(pow(b_coef,2.0)-4.0*d_coef))/2.0;
  
      r_sol->p_hydro = max(pres_1,pres_2);
    }

  } else {
    pres_1 = (-b_coef+sqrt(pow(b_coef,2.0)-4.0*d_coef))/2.0;
    pres_2 = (-b_coef-sqrt(pow(b_coef,2.0)-4.0*d_coef))/2.0;
  
    r_sol->p_hydro = max(pres_1,pres_2);    
  }


  r_lm = r_sat + (1.0/pow(c_temp,2.0))*(r_sol->p_hydro-p_sat);

  // Compute vol. fractions from p_old
  r_gas=r_sol->p_hydro/(R_G*T_ref);

  a_gas = ((r_sol->ymass*r_sol->r)/r_gas); // CHECK
  a_lm  = 1.0 - a_gas;                     // CHECK


    
    
  // Wallis expression
  r_sol->c =c_temp;
  //r_sol->c=sqrt( 1.0/(r_sol->r* ( (a_lm/(r_lm*c_temp*c_temp)) + (a_gas/(r_gas*c_gas*c_gas)) ) ) );
  if ((isnan(c_temp)==1)){
    printf(" Speed of sound nan \n");
  }
}

void barotropic_threephase_eos(struct BRANCH * crnt,int istore,struct RUN *run){

  int f;
  int flag;


  double r_sat_liq = 998.16;
  double r_sat_vap = 0.0173;
  double p_sat_liq = 4664.4;
  double p_ref     = 125.0;

  double c_sat_liq = 1483.26;
  double c_sat_vap = 97.9;

  double gamma     = 1.33;
  double B         = 307100000.0;
  double N         = 1.75;

  double C_vap     = 75267;
  double a_vap,a_liq;
  double c_temp,r_lm,r_gas,a_gas,a_lm;
  double a_coef,b_coef,d_coef;
  double pres_1,pres_2;
  double p_temp;
  double xw,yw,zw;
  double R1;
 
  
  double check_NaN,check_NaN_2,check_3;
  
  struct SOL * r_sol;
  struct SOL * r_sol1;



  xw=crnt->cl->xc;
  yw=crnt->cl->yc;
  zw=crnt->cl->zc;

  if (istore==0){ r_sol = run->sol;  }
  if (istore==1){ r_sol = run->sol1; }   


  if (r_sol->r>=r_sat_liq) {
    p_temp = B * (pow((r_sol->r/r_sat_liq),N)-1.0)+p_sat_liq;
    c_temp = c_sat_liq;
  }else if (r_sol->r<r_sat_liq && r_sol->r>r_sat_vap) {
    check_NaN=(pow((c_sat_vap*r_sat_vap),2.0))-(pow((c_sat_liq*r_sat_liq),2.0));
    if (isnan(check_NaN)==1) {
      printf(" double nan \n");
      exit(0);
    }

    check_NaN_2=max(1.01,log((r_sol->r/(c_sat_liq*c_sat_liq*r_sat_liq*(r_sat_liq-r_sol->r)+c_sat_vap*c_sat_vap*r_sat_vap*(r_sol->r-r_sat_vap)))));
    check_3=((r_sol->r/(c_sat_liq*c_sat_liq*r_sat_liq*(r_sat_liq-r_sol->r)+c_sat_vap*c_sat_vap*r_sat_vap*(r_sol->r-r_sat_vap))));
    if (check_NaN_2<=1.0){
      printf (" check_3:%f,check_Nan_2:%f,r:%f,c_sat_liq:%f,r_sat_liq:%f,c_sat_vap:%f,r_sat_vap:%f\n",check_3,check_NaN_2,r_sol->r,c_sat_liq,r_sat_liq,c_sat_vap,r_sat_vap);
    }
    //(max(1.01,log(r_sol->r/(c_sat_liq*c_sat_liq*r_sat_liq*(r_sat_liq-r_sol->r)+c_sat_vap*c_sat_vap*r_sat_vap*(r_sol->r-r_sat_vap))))
    p_temp = (((pow((c_sat_vap*c_sat_liq),2.0))*r_sat_liq*r_sat_vap*(r_sat_vap-r_sat_liq))/(pow((c_sat_vap*r_sat_vap),2.0))-(pow((c_sat_liq*r_sat_liq),2.0)))* \
    check_NaN_2+p_ref;
     a_vap = (r_sol->r-r_sat_vap)/(r_sat_liq-r_sat_vap);
     a_liq=1-a_vap;
     c_temp=sqrt( 1.0/(r_sol->r* ( (a_liq/(r_sat_liq*c_sat_liq*c_sat_liq)) + (a_vap/(r_sat_vap*c_sat_vap*c_sat_vap)) ) ) );
   } else {
    p_temp = C_vap*pow(r_sol->r,gamma);
    c_temp = c_sat_vap;
  }
 
  r_sol->c =c_temp;

  r_sol->p_hydro = p_temp;
 // printf("c:%f,p:%f\n",r_sol->c,r_sol->p_hydro);
  //  r_sol->c=sqrt( 1.0/(r_sol->r* ( (a_lm/(r_sat_liq*c_sat_liq*c_sat_liq)) + (a_gas/(r_gas*c_gas*c_gas)) ) ) );
  }




