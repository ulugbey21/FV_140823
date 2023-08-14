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

//-------------------------------------------------------
//
//	Counts the time between to calls with:
//		-fun=0
//		-fun=1
//
//-------------------------------------------------------

#include "strdata.h"

void watch_init(RUN * run) {
    int iw;
    for(iw=0;iw<20;iw++){
      run->con->timewatchtall[iw]=0.0;
      run->con->timewatch0[iw]=0.0;
      run->con->timewatcht[iw]=0.0;
    }
}

void watch(RUN * run, int id, int fun) {
  if (fun==0){
    run->con->timewatch0[id]=MPI_Wtime();
  }
  if (fun==1){
    run->con->timewatcht[id]+=MPI_Wtime() - run->con->timewatch0[id];

  }
}
void watch_display(RUN * run, int tag) {
  int iw,iw1,iw2;
  MPI_Allreduce(run->con->timewatcht,run->con->timewatchtall,20,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if(run->con->rank==0){
    if (tag==1){printf ("\n>W      TAG      STP  INIT      RPRT      REFN      RLST      RPRX      RCON      GEOM      COMS      COMC      EMPY");iw1=0;iw2=10;}
    if (tag==2){printf ("\n>W      TAG      STP  LOOP      ADPT      ADVC      INRP      RHSR      PRLX      COMS      COMF      COMG      IORS");iw1=10;iw2=20;}
    printf ("\n>W %8d %8d ",tag,run->con->tstep);
    for(iw=iw1;iw<iw2;iw++){
      printf (" %4.2e ",run->con->timewatchtall[iw]/(run->con->size*(run->con->istep+1)));
    }
    printf ("\n>W %8d %8d ",tag,run->con->tstep);
    for(iw=iw1;iw<iw2;iw++){
      printf (" %8.2f ",100.0*run->con->timewatchtall[iw]/run->con->timewatchtall[iw1]);
    }
    printf ("\n");
  }
}




void watch_reset(RUN * run) {
  int iw;

  for(iw=0;iw<20;iw++){
    run->con->timewatchtall[iw]=0.0;
    run->con->timewatch0[iw]=0.0;
    run->con->timewatcht[iw]=0.0;
  }




}
