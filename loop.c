/*
 * Forest of oct-Trees DG implementation
 * octo-Rhodon
 *
 *  \/\/\/\/
 *   \/ / /
 *    \/ /
 *     \/
 *
 * Foundation for Research and Technology Hellas, 
 * Institute for Applied and Computational Mathematics (FORTH/IACM)
 * Heraklion Crete
 * Andreas Papoutsakis 
 * 20/05/2014
*/

#include "strdata.h"

void loop (RUN * run) {
  
  int istep;
  struct BRANCH * crnt;
    
  run->cputime=0.0;  
  run->adapt_flag_sensor = 0;

  watch_init(run); //LOOP
  
  //================================================================
  //                        Time Loop
  //================================================================
 
  inner_el(run);
  
  for (istep=0;istep<(run->con->nstep+1);istep++) {

    Init_clocks(run,istep);
    
    run->Ttot=timecpu(run->Ttot,0);
 
    run->Tmeshadapt=timecpu(run->Tmeshadapt,0);
    if ( (run->con->istep%run->con->adapthstep==0) && (run->con->adapth!=0) ) { 
      run->adapt_flag_sensor=1;     
      meshadaptation(run);
    }
    run->Tmeshadapt=timecpu(run->Tmeshadapt,1);
    
    advancetimeexplicit(run);
        
    run->Tepx=timecpu(run->Tepx,0);
    if (((istep+1)%run->con->exstep==0) && (istep!=0)) { 
      gradcalc(run);
      //gradcalc_V2secondorder(run);
      communicate_S(run,2); 
      exportfield_linked_V3(run);
    //  exportfield(run); 
      //export_ultra_fast(run);
    }
    if ((istep+1)%run->con->ststep==0) {
    	restart_save(run);
    }
    run->Tepx=timecpu(run->Tepx,1);
      
    run->stepcputime  = MPI_Wtime()-run->stepcputime;   // Time of iteration
    run->cputime     += run->stepcputime;               // Totall time of simulation

    run->Ttot=timecpu(run->Ttot,1);
  
    //printf(" run->con->Ttot %d %f \n",run->con->rank,Ttot);
    sumtime(run);

    if (istep%100==0) { screen_out(run,istep); }

  }

}

