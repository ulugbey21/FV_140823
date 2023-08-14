#include "strdata.h"
void init_domain(RUN *run){

  double xw,yw,zw; 
  double rho,u,v,w,p,ymass;
  double R0,R1,dist_wall;
  double R_G,T_ref;
  struct BRANCH * crnt;
  
  crnt=run->topo->locl;
  while (crnt!=NULL){
    
    xw=crnt->cl->xc;
    yw=crnt->cl->yc;
    zw=crnt->cl->zc;
 
        
    if (run->con->i_case==1) { // Model 1, Test 1, Air, gamma=1.4, p_inf=0, Toro pg 225, table 6.2 

      if(crnt->root>300) {  // right
        
        rho = 0.125;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 0.1;
        Init_S(run,crnt,rho,u,v,w,p);

      } else {  // left

        rho = 1.0;
        u   = 0.75;
        v   = 0.0;
        w   = 0.0;
        p   = 1.0;
        Init_S(run,crnt,rho,u,v,w,p);
      }        

    }
    else if (run->con->i_case==2) { // Model 1, Test 2, Air, gamma=1.4, p_inf=0, Toro pg 225, table 6.2

      if(crnt->root<500) {  // left
        
        rho = 1.0;
        u   = -2.0;
        v   = 0.0;
        w   = 0.0;
        p   = 0.4;
        Init_S(run,crnt,rho,u,v,w,p);

      } else {  // right

        rho = 1.0;
        u   = 2.0;
        v   = 0.0;
        w   = 0.0;
        p   = 0.4;
        Init_S(run,crnt,rho,u,v,w,p);

      }        

    }
    else if (run->con->i_case==3) { // Model 1, Test 3, Air, gamma=1.4, p_inf=0, Toro pg 225, table 6.2

      if(crnt->root<500) {  // left
        
        rho = 1.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 1000.0;
        Init_S(run,crnt,rho,u,v,w,p);
        

      } else {  // right

        rho = 1.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 0.01;
        Init_S(run,crnt,rho,u,v,w,p);

      }        

    }

    else if (run->con->i_case==4) { // Model 1, Test 4, Air, gamma=1.4, p_inf=0, Toro pg 225, table 6.2

      if(crnt->root<400) {  // left
        
        rho = 5.99924;
        u   = 19.5975;
        v   = 0.0;
        w   = 0.0;
        p   = 460.894;
        Init_S(run,crnt,rho,u,v,w,p);
        

      } else {  // right

        rho = 5.99924;
        u   = -6.19633;
        v   = 0.0;
        w   = 0.0;
        p   = 46.0950;
        Init_S(run,crnt,rho,u,v,w,p);

      }        

    }

    else if (run->con->i_case==5) { // Model 1, Test 4, Air, gamma=1.4, p_inf=0, Toro pg 225, table 6.2

      if(crnt->root<800) {  // left
        
        rho = 1.0;
        u   = -19.59745;
        v   = 0.0;
        w   = 0.0;
        p   = 1000.0;
        Init_S(run,crnt,rho,u,v,w,p);
        

      } else {  // right

        rho = 1.0;
        u   = -19.59745;
        v   = 0.0;
        w   = 0.0;
        p   = 0.01;
        Init_S(run,crnt,rho,u,v,w,p);

      }        

    }
    else if (run->con->i_case==10) {  // Model 3, Case A (page 16) https://dx.doi.org/10.6084/m9.figshare.c.4388069.
      
      if(crnt->root>=5000) {   // Gas

        rho =  0.0170;
        //rho =  1.185;
        u   =  0.0; //100.0;
        v   =  0.0;
        w   =  0.0;
        //p   =  140.0;
        ymass = 1.0-run->con->ymin;
        
       // ymass = 1.0;
        Init_S(run,crnt,rho,u,v,w,p);  
       // Init_S3(run,crnt,rho,u,v,w,ymass);

        
      } else {  // Liquid
       rho = 998.206;
       // rho = 1500.202;
        u   = 0.0; //100.0;
        v   = 0.0;
        w   = 0.0;

        //ymass = run->con->ymin;
       // ymass = 0.0;
       //rho = 1500.202;
       Init_S(run,crnt,rho,u,v,w,p);  
      //  Init_S3(run,crnt,rho,u,v,w,ymass);  
      }

    }
    else if (run->con->i_case==11) {  // Model 3, Case B (page 16) https://dx.doi.org/10.6084/m9.figshare.c.4388069.
     
      if(xw<0.0) {
        rho = 1000.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 3810000.0;
        Init_S(run,crnt,rho,u,v,w,p);  
      } else {
        rho =  0.998;
        u   =  435;
        v   =  0.0;
        w   =  0.0;
        p   =  887.8;
        Init_S(run,crnt,rho,u,v,w,p);
        }
    }
    else if (run->con->i_case==12) {  // Model 3, Case C (page 16) https://dx.doi.org/10.6084/m9.figshare.c.4388069.
     
      if(xw<0.0) {
        rho = 1000.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 3810000.0;
        Init_S(run,crnt,rho,u,v,w,p);  
      } else {
        rho =  0.998;
        u   =  435;
        v   =  0.0;
        w   =  0.0;
        p   =  887.8;
        Init_S(run,crnt,rho,u,v,w,p);
        }
    }
    else if (run->con->i_case==90) {  // Model 2, shocktube case barotropic
     
      if(xw<0.0) {
        rho = 1000.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 3810000.0;
        Init_S(run,crnt,rho,u,v,w,p);  
    } else {
      rho =  0.998;
      u   =  435;
      v   =  0.0;
      w   =  0.0;
      p   =  887.8;
      Init_S(run,crnt,rho,u,v,w,p);
      }
    }
    else if (run->con->i_case==91) { // Model 2, shocktube case barotropic
      // shocktube case
      if(xw<0.5) {
        rho = 1.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 1.0;
        Init_S(run,crnt,rho,u,v,w,p);  
      } else {
        rho =  0.125;
        u   =  0;
        v   =  0.0;
        w   =  0.0;
        p   =  0.1;
        Init_S(run,crnt,rho,u,v,w,p);
      }
    }
    else if (run->con->i_case==200) { // Model 1, Forward faceing step, https://amroc.sourceforge.net/examples/euler/2d/html/ffstep_n.htm
  
      rho = 1.4;
      u   = 3.0;
      v   = 0.0;
      w   = 0.0;
      p   = 1.0;
      Init_S(run,crnt,rho,u,v,w,p);

    }
    else if (run->con->i_case==201) { // Model 1, Backward facing step, https://amroc.sourceforge.net/examples/euler/2d/html/bfstep_n.htm
  
      rho = 1.0;
      u   = 0.0;
      v   = 0.0;
      w   = 0.0;
      p   = 1.0;
      Init_S(run,crnt,rho,u,v,w,p);

    }
    else if (run->con->i_case==202) { // Model 1, Ramp, https://amroc.sourceforge.net/examples/euler/2d/html/ramp_n.htm

      if(xw+yw>0.0) { 
        rho = 1.4;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 1.0;
        Init_S(run,crnt,rho,u,v,w,p);
      } else {
        rho =  8.0;
        //u   =  8.25*cos(PI/6.0);
        //v   = -8.25*sin(PI/6.0);
        w   =  0.0;
        p   =  116.5;
        Init_S(run,crnt,rho,u,v,w,p);
      }

    }
    else if (run->con->i_case==203) { // Model 1, NS, boundary layer case for viscosity

        rho = 1.0;
        u   = 68;
        v   = 0.0;
        w   = 0.0;
        p   = 101325.0;
        Init_S(run,crnt,rho,u,v,w,p);
     
    }
    if (run->con->i_case==290) { // NS, boundary layer case for viscosity

       // run->par->materpini[2]   =  30.0*pow(10.0,5.0);
      //  run->par->materpini[1]   =  100.0*pow(10.0,5.0);
      //  run->par->materpini[0]   =  100.0*pow(10.0,5.0);
        rho =1000.0;
        u   =  0.0; //100.0;
        v   =  0.0;
        w   =  0.0;
        //p   =  140.0;
        ymass = 1.0-run->con->ymin;
        Init_S3(run,crnt,rho,u,v,w,ymass);  
     
    }

      else if (run->con->i_case==2100) {  // Model 3, Case A (page 16) https://dx.doi.org/10.6084/m9.figshare.c.4388069.
      R0  =  0.00056;
      R1  =  0.00015;
      R_G =  286.0;
      T_ref = 293.15;

      if ((pow((xw-0.00157),2.0) + pow(yw,2.0) <= pow(R0,2.0)) || xw>0.00157 ) {   // Gas
        
        p   =  pow(10.0,5.0);
        //p=1.0;
        rho =  p/(R_G*T_ref);
        
        //rho =  1.055;
        u   =  0.0; //100.0;
        v   =  0.0;
        w   =  0.0;
        //p   =  140.0;
        ymass = 1.0-run->con->ymin;
        
       // ymass = 1.0;
        if ((run->con->model==1) || (run->con->model==2)) {
            Init_S(run,crnt,rho,u,v,w,p);  
        }
        else{
             Init_S3(run,crnt,rho,u,v,w,ymass);
        }

        
      }
      else if (((pow((xw),2.0) + pow((yw),2.0)) <= pow(R1,2.0)) ) {   // Gas
        
        p   =  1.0*pow(10.0,7.0);
        rho =  p/(R_G*T_ref);
       //rho=1020.0;
        //rho =  1.055;
        u   =  0.0; //100.0;
        v   =  0.0;
        w   =  0.0;
        
        ymass = 1.0-run->con->ymin;
        
       // ymass = 1.0;
        if ((run->con->model==1) || (run->con->model==2)) {
        Init_S(run,crnt,rho,u,v,w,p);  
        }
        else{
          Init_S3(run,crnt,rho,u,v,w,ymass);
        }
      
      }
      
      
       else {  // Liquid
       rho = 998.206;
        u   = 0.0; 
        v   = 0.0;
        w   = 0.0;
        ymass = run->con->ymin;
        if ((run->con->model==1) || (run->con->model==2)) {
        Init_S(run,crnt,rho,u,v,w,p);  
        }
        else{
          Init_S3(run,crnt,rho,u,v,w,ymass);
        }
      }

    }
    
    else if (run->con->i_case==2500) { // Model 2, bubble-wall collapse
      
      R0 = 0.4/1000.0;
      dist_wall = 0.016/1000.0;
      if ( (pow((xw),2.0) + pow((yw-R0-dist_wall),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        p=2311.5;   // low pressure
        rho=49.91;    // low density
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        Init_S(run,crnt,rho,u,v,w,p);
      }
      else{
        p=1.0*pow(10.0,6); // high pressure
        rho=1002.89;               // high density
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        Init_S(run,crnt,rho,u,v,w,p);
      }
    }
    /*
    else if (run->con->i_case==2100) { // Model 3, needlesfree injection
      R0 = 0.24/1000.0;
      R1 = 0.08/1000.0;
      dist_wall = 1.64/1000.0;
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0) || xw >=dist_wall ){ // gas
        p=101325.0;
        rho=1.225;
        u  = 0.0;
        v  = 0.0;
        w  = 0.0;
        ymass = 1.0;
        Init_S3(run,crnt,rho,u,v,w,ymass);
      } else if ( (pow((xw),2.0) + pow((yw),2.0)) <= pow(R1,2.0) && (xw<=(0.08/1000))){ // 
        p=5.0*pow(10.0,7);
        rho=303.96; // 300 degrees
        u = 0.0;
        v = 0.0;
        w = 0.0;
        ymass = 1.0;
        Init_S3(run,crnt,rho,u,v,w,ymass);
      } else {
        p=101325.0;
        rho=998.17;
        u = 0.0;
        v = 0.0;
        w = 0.0;
        ymass = 0.0;
        Init_S3(run,crnt,rho,u,v,w,ymass);
      }
    }

    */
    else if (run->con->i_case==2900) { // injector

        if (xw <= 0.056) {    
          rho   =998.20;    // hight density
          u     = 0.0;
          v     = 0.0;
          w     = 0.0;
          ymass = 0.0;
          Init_S3(run,crnt,rho,u,v,w,ymass);
        }
        else{
          //p=pow(10.0,7); // high pressure
          rho   = 1.225;               // low density
          u     = 0.0;
          v     = 0.0;
          w     = 0.0;
          ymass = 1.0;
          Init_S3(run,crnt,rho,u,v,w,ymass);
        }
      }

      
    if (run->con->i_case==3100) { //edelbauer channel flow les 

        rho = 830.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 101325.0;
        ymass = 1-run->con->ymin;
        Init_S3(run,crnt,rho,u,v,w,ymass);  
      }

     else if (run->con->i_case==3400) {
        p   =   1.0*pow(10.0,5.0); // erlangen injector
        rho   = 722.64; 
        u     = 0.0;
        v     = 0.0;
        w     = 0.0;
        ymass = 1-run->con->ymin;
        Init_S3(run,crnt,rho,u,v,w,ymass);  
      }
    
    crnt=crnt->lnxt;
  }

}
