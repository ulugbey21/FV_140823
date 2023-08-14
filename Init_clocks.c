#include "strdata.h"

void Init_clocks (struct RUN * run, int istep) {

  /*
  run->con->stepcputime=MPI_Wtime(); // Time of iteration, start clock

  run->con->Tmeshadapt = 0.0;
  run->con->T2nd       = 0.0;
  run->con->Tepx       = 0.0;
  run->con->Tcomm      = 0.0;
  run->con->Tinterfloc = 0.0;
  run->con->Tadv       = 0.0;
  run->con->Trhs       = 0.0;
  run->con->Trelax     = 0.0;
  run->con->Ttot       = 0.0;
  
  if (istep==0) {

    run->con->stepcputime_sum = 0.0;
    run->con->Ttot_sum        = 0.0;
    run->con->Tmeshadapt_sum  = 0.0;
    run->con->Tadv_sum        = 0.0;
    run->con->T2nd_sum        = 0.0;
    run->con->Trhs_sum        = 0.0;
    run->con->Trelax_sum      = 0.0;
    run->con->Tcomm_sum       = 0.0;
    run->con->Tinterfloc_sum  = 0.0;
    run->con->Tepx_sum        = 0.0;

  }
  */

  run->stepcputime=MPI_Wtime(); // Time of iteration, start clock

  run->Tmeshadapt = 0.0;
  run->T2nd       = 0.0;
  run->T2ndcom    = 0.0;
  run->Tepx       = 0.0;
  run->Tcomm      = 0.0;
  run->Tinterfloc = 0.0;
  run->Tadv       = 0.0;
  run->Trhs       = 0.0;
  run->Trelax     = 0.0;
  run->Ttot       = 0.0;

  if (istep==0) {

    run->stepcputime_sum = 0.0;
    run->Ttot_sum        = 0.0;
    run->Tmeshadapt_sum  = 0.0;
    run->Tadv_sum        = 0.0;
    run->T2nd_sum        = 0.0;
    run->T2ndcom_sum     = 0.0;
    run->Trhs_sum        = 0.0;
    run->Trelax_sum      = 0.0;
    run->Tcomm_sum       = 0.0;
    run->Tinterfloc_sum  = 0.0;
    run->Tepx_sum        = 0.0;

    run->Tdist_sum     = 0.0;
    run->Tdelta_sum    = 0.0;
    run->Trec_sum      = 0.0;
    run->T2ndrelax_sum = 0.0;
    run->T2ndface_sum  = 0.0;
  }

  run->con->iter=1;             // Flag for start of time loop

  run->con->time+=run->con->dt; // Totall time of simulation
  run->con->istep=istep;        // time steping since new launch
  run->con->tstep++;            // time steping with added restart time

}