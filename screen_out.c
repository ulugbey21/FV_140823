#include "strdata.h"
 
void screen_out(struct RUN * run,int istep) {

  double memory,memory_tot,check,time_remaining;
    
  int currRealMem; 
  int peakRealMem;
  int currVirtMem; 
  int peakVirtMem; 
  int cur_step,remaining_steps;

  if (istep%1000==0) {

    // stores each word in status file
    char buffer[1024] = "";

    // linux file contains this-process info
    FILE* file = fopen("/proc/self/status", "r");

    // read the entire file
    while (fscanf(file, " %1023s", buffer) == 1) {

        if (strcmp(buffer, "VmRSS:") == 0)  { fscanf(file, " %d", &currRealMem); }
        if (strcmp(buffer, "VmHWM:") == 0)  { fscanf(file, " %d", &peakRealMem); }
        if (strcmp(buffer, "VmSize:") == 0) { fscanf(file, " %d", &currVirtMem); }
        if (strcmp(buffer, "VmPeak:") == 0) { fscanf(file, " %d", &peakVirtMem); }
    }
    fclose(file);
   
    memory = ((double) currVirtMem)*0.001;
 
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&memory, &memory_tot, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
    
    if (run->con->istep==0) {
      time_remaining=0.0;
    } else {

      cur_step = run->con->tstep;
      remaining_steps= run->con->nstep+run->con->tstep0-cur_step;
      time_remaining=(((double)remaining_steps)*(run->cputime/cur_step))/3600.0;
    }

    if(run->con->rank==0){  printf ("\n");
      printf (">|_______________________________________________________________________________________________________|\n"                 );
      printf (">|_________________________Finite Volume on Hybrid Linked OTree Forest  CASE:%3d _______________________|\n",run->con->i_case);
      printf (">|_______________________________________________________________________________________________________|\n"                 );
      printf (">|________________   Istart|      Ifin|  Nexport|   Nadapt|        DT|       T0|     Tfin|               |\n"                 );
      printf (">|________________ %8d|  %8d| %8d| %8d|  %.2e| %.2e| %.2e|               |\n",
      run->con->tstep0,
      run->con->nstep+run->con->tstep0,
      run->con->exstep,
      run->con->adapthstep
      ,run->con->dt0
      ,run->con->time0
      ,run->con->timetot+run->con->time0);
    
      printf (">|_______________________________________________________________________________________________________|\n");
      printf (">|_______________   Ntrees|   NLeaves|   Nlevel| Nprocess| Memory usage [GB]| Expected time remaining [h]| \n");    
      printf (">|_______________ %8d|  %8d| %8d| %8d|          %f|                      %5.2f|\n",
      run->topo->ntrees,run->tot_leaves,run->con->level,run->con->size,memory_tot/1000.0,time_remaining);
      printf (">|_______________________________________________________________________________________________________|\n");
     // printf (">|Iiter  |Istep |Time     |Tcpu m| T/itr |RHS (%)|AMR (%)|WRO (%)| \n");
      printf (">|Iiter  |Istep |Time     |RT [m]|T/100[s]|AMR(s)|ADV (s)|2nd (s)|2com(s)|RHS (s)|REL (s)|COM (s)|EXP (s)| \n");
    }
    
  }

  if (run->con->istep==0) {

    if (run->con->rank==0) {

      /*
      check = (run->con->Tadv        *100.0/run->con->Ttot) + 
              (run->con->Tmeshadapt  *100.0/run->con->Ttot) +
              (run->con->Tinterfloc  *100.0/run->con->Ttot) +
              (run->con->Tepx        *100.0/run->con->Ttot);
      */
      printf (">I|%6d|%6d| %.2e|  %.2f|   %5.3f| %5.2f|  %5.2f|  %5.2f|  %5.2f|  %5.2f|  %5.2f|  %5.2f|  %5.2f|     \n",
      run->con->istep,run->con->tstep,run->con->time,
                                                            run->cputime/60.0,run->stepcputime,
                                                            run->Tmeshadapt,
                                                            run->Tadv,
                                                            run->T2nd,
                                                            run->T2ndcom,
                                                            run->Trhs,
                                                            run->Trelax,
                                                            run->Tcomm,
                                                            run->Tepx);

                                                            /*run->con->cputime/60.0,run->con->stepcputime,
                                                            run->con->Tmeshadapt,
                                                            run->con->Tadv,
                                                            run->con->T2nd,
                                                            run->con->Trhs,
                                                            run->con->Trelax,
                                                            run->con->Tcomm,
                                                            run->con->Tepx);*/
    }
  } else {
    if (run->con->rank==0) {
      /*
      check = (run->con->Tadv_sum        *100.0/run->con->Ttot_sum) + 
              (run->con->Tmeshadapt_sum  *100.0/run->con->Ttot_sum) +
              (run->con->Tinterfloc_sum  *100.0/run->con->Ttot_sum) +
              (run->con->Tepx_sum        *100.0/run->con->Ttot_sum);
      */
      //textcolor(RESET,RED,WHITE);  
      printf (">T|%6d|%6d| %.2e|  %.2f|  %5.3f| %5.2f|  %5.2f|  %5.2f|  %5.2f|  %5.2f|  %5.2f|  %5.2f|  %5.2f|    \n",
      run->con->istep,run->con->tstep,run->con->time,
                                                            run->cputime/60.0,run->stepcputime_sum,
                                                            run->Tmeshadapt_sum,
                                                            run->Tadv_sum,
                                                            run->T2nd_sum,
                                                            run->T2ndcom_sum,
                                                            run->Trhs_sum,
                                                            run->Trelax_sum,
                                                            run->Tcomm_sum,
                                                            run->Tepx_sum);

                                                            /*run->con->cputime/60.0,run->con->stepcputime_sum,
                                                            run->con->Tmeshadapt_sum,
                                                            run->con->Tadv_sum,
                                                            run->con->T2nd_sum,
                                                            run->con->Trhs_sum,
                                                            run->con->Trelax_sum,
                                                            run->con->Tcomm_sum,
                                                            run->con->Tepx_sum);*/

      // textcolor(RESET,BLACK,WHITE);
    }

    //timecpu_exp(run);

  }
  
  run->stepcputime_sum = 0.0;
  run->Tepx_sum        = 0.0;
  run->Tmeshadapt_sum  = 0.0;
  
  run->Ttot_sum        = 0.0;
  
  run->Tadv_sum        = 0.0;
  run->T2nd_sum        = 0.0;
  run->T2ndcom_sum     = 0.0;
  run->Trhs_sum        = 0.0;
  run->Trelax_sum      = 0.0;
  run->Tcomm_sum       = 0.0;
  run->Tinterfloc_sum  = 0.0;
  

  run->Tdist_sum     = 0.0;
  run->Tdelta_sum    = 0.0;
  run->Trec_sum      = 0.0;
  run->T2ndrelax_sum = 0.0;
  run->T2ndface_sum  = 0.0;
  

}