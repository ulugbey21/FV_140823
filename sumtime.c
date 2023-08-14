#include "strdata.h"
 
void sumtime(struct RUN * run) {

  /*
  run->con->stepcputime_sum += run->con->stepcputime;
  run->con->Ttot_sum        += run->con->Ttot;
  run->con->Tmeshadapt_sum  += run->con->Tmeshadapt;
  run->con->Tadv_sum        += run->con->Tadv;
  run->con->T2nd_sum        += run->con->T2nd; 
  run->con->Trhs_sum        += run->con->Trhs;
  run->con->Trelax_sum      += run->con->Trelax;
  run->con->Tcomm_sum       += run->con->Tcomm;
  run->con->Tinterfloc_sum  += run->con->Tinterfloc;
  run->con->Tepx_sum        += run->con->Tepx;
  */

  run->stepcputime_sum += run->stepcputime;
  run->Ttot_sum        += run->Ttot;
  run->Tmeshadapt_sum  += run->Tmeshadapt;
  run->Tadv_sum        += run->Tadv;
  run->T2nd_sum        += run->T2nd; 
  run->T2ndcom_sum     += run->T2ndcom;
  run->Trhs_sum        += run->Trhs;
  run->Trelax_sum      += run->Trelax;
  run->Tcomm_sum       += run->Tcomm;
  run->Tinterfloc_sum  += run->Tinterfloc;
  run->Tepx_sum        += run->Tepx;
    
}