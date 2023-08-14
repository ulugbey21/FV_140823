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

void normalvector(RUN *run){
  
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){
    normalvectorcalc(run,crnt);
  crnt=crnt->lnxt;}
}

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


void normalvectorcalc(RUN *run,struct BRANCH * crnt){

  int f,type,n;
  double norm;
  double xf,yf;

  type=crnt->type;
  n=2;
  for(f=0;f<crnt->nlfc;f++){//face
    norm=sqrt(pow(((crnt->cl->nd[fcnd2elnd(f,1,type)].y-crnt->cl->nd[fcnd2elnd(f,0,type)].y)*(crnt->cl->nd[fcnd2elnd(f,n,type)].z-crnt->cl->nd[fcnd2elnd(f,0,type)].z)-(crnt->cl->nd[fcnd2elnd(f,n,type)].y-crnt->cl->nd[fcnd2elnd(f,0,type)].y)*(crnt->cl->nd[fcnd2elnd(f,1,type)].z-crnt->cl->nd[fcnd2elnd(f,0,type)].z)),2.0)+pow(((crnt->cl->nd[fcnd2elnd(f,1,type)].x-crnt->cl->nd[fcnd2elnd(f,0,type)].x)*(crnt->cl->nd[fcnd2elnd(f,n,type)].z-crnt->cl->nd[fcnd2elnd(f,0,type)].z)-(crnt->cl->nd[fcnd2elnd(f,n,type)].x-crnt->cl->nd[fcnd2elnd(f,0,type)].x)*(crnt->cl->nd[fcnd2elnd(f,1,type)].z-crnt->cl->nd[fcnd2elnd(f,0,type)].z)),2.0)+pow(((crnt->cl->nd[fcnd2elnd(f,1,type)].x-crnt->cl->nd[fcnd2elnd(f,0,type)].x)*(crnt->cl->nd[fcnd2elnd(f,n,type)].y-crnt->cl->nd[fcnd2elnd(f,0,type)].y)-(crnt->cl->nd[fcnd2elnd(f,n,type)].x-crnt->cl->nd[fcnd2elnd(f,0,type)].x)*(crnt->cl->nd[fcnd2elnd(f,1,type)].y-crnt->cl->nd[fcnd2elnd(f,0,type)].y)),2.0));
    
    crnt->cl->nx[f]=((crnt->cl->nd[fcnd2elnd(f,1,type)].y-crnt->cl->nd[fcnd2elnd(f,0,type)].y)*(crnt->cl->nd[fcnd2elnd(f,n,type)].z-crnt->cl->nd[fcnd2elnd(f,0,type)].z)-(crnt->cl->nd[fcnd2elnd(f,n,type)].y-crnt->cl->nd[fcnd2elnd(f,0,type)].y)*(crnt->cl->nd[fcnd2elnd(f,1,type)].z-crnt->cl->nd[fcnd2elnd(f,0,type)].z))/norm;
    crnt->cl->ny[f]=-((crnt->cl->nd[fcnd2elnd(f,1,type)].x-crnt->cl->nd[fcnd2elnd(f,0,type)].x)*(crnt->cl->nd[fcnd2elnd(f,n,type)].z-crnt->cl->nd[fcnd2elnd(f,0,type)].z)-(crnt->cl->nd[fcnd2elnd(f,n,type)].x-crnt->cl->nd[fcnd2elnd(f,0,type)].x)*(crnt->cl->nd[fcnd2elnd(f,1,type)].z-crnt->cl->nd[fcnd2elnd(f,0,type)].z))/norm;
    crnt->cl->nz[f]=((crnt->cl->nd[fcnd2elnd(f,1,type)].x-crnt->cl->nd[fcnd2elnd(f,0,type)].x)*(crnt->cl->nd[fcnd2elnd(f,n,type)].y-crnt->cl->nd[fcnd2elnd(f,0,type)].y)-(crnt->cl->nd[fcnd2elnd(f,n,type)].x-crnt->cl->nd[fcnd2elnd(f,0,type)].x)*(crnt->cl->nd[fcnd2elnd(f,1,type)].y-crnt->cl->nd[fcnd2elnd(f,0,type)].y))/norm;

    if (crnt->cl->fc[f].bc!=0){
  
      if (run->par->bcoverwrt[0]!=0){ if (crnt->cl->nx[f]<-0.6){crnt->cl->fc[f].bc=run->par->bcoverwrt[0];} }
      if (run->par->bcoverwrt[1]!=0){ if (crnt->cl->nx[f]> 0.6){crnt->cl->fc[f].bc=run->par->bcoverwrt[1];} }

      if (run->par->bcoverwrt[2]!=0){ if (crnt->cl->ny[f]<-0.6){crnt->cl->fc[f].bc=run->par->bcoverwrt[2];} }
      if (run->par->bcoverwrt[3]!=0){ if (crnt->cl->ny[f]> 0.6){crnt->cl->fc[f].bc=run->par->bcoverwrt[3];} }
    
      if (run->par->bcoverwrt[4]!=0){ if (crnt->cl->nz[f]<-0.6){crnt->cl->fc[f].bc=run->par->bcoverwrt[4];} }
      if (run->par->bcoverwrt[5]!=0){ if (crnt->cl->nz[f]> 0.6){crnt->cl->fc[f].bc=run->par->bcoverwrt[5];} }

      if (run->con->i_case==200) {
        if ((crnt->cl->nd[fcnd2elnd(f,1,type)].x)>0.55 && (crnt->cl->nd[fcnd2elnd(f,1,type)].x<0.65)) {
          if (crnt->cl->nx[f]>0.6){crnt->cl->fc[f].bc=3;}
        }
      }
      //flat plate
      if (run->con->i_case==203) {
        if ((crnt->cl->nd[fcnd2elnd(f,0,type)].x)<0.0 ){
          if (crnt->cl->ny[f]<-0.6){crnt->cl->fc[f].bc=3;}
        }
      }
    
      if (run->con->i_case==2900) {
          if (((crnt->cl->nx[f]> 0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<(5.0*run->con->mscale))) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>(4.0*run->con->mscale)) ) {
            crnt->cl->fc[f].bc=1;
         }
        if (((crnt->cl->nx[f]<-0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<(15.0*run->con->mscale))) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>(14.0*run->con->mscale)) ) {
            crnt->cl->fc[f].bc=1;
          }
        }
//les channel flow 
      if (run->con->i_case==3000) {
          if (((crnt->cl->ny[f]<- 0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.020)) ) {
            crnt->cl->fc[f].bc=3;
         }
          if (((crnt->cl->ny[f]>0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.020)) ) {
            crnt->cl->fc[f].bc=3;
         }
          if (((crnt->cl->nz[f]<-0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.020)) ) {
            crnt->cl->fc[f].bc=3;
         }
          if (((crnt->cl->nz[f]>0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.020)) ) {
            crnt->cl->fc[f].bc=3;
         }
          if (((crnt->cl->nx[f]<- 0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.019) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<0.021) )) {
            crnt->cl->fc[f].bc=3;
         }
          if (((crnt->cl->nx[f]> 0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.011) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<0.013)) ) {
            crnt->cl->fc[f].bc=4;
         }
        }
//edelbauer les channel flow

      if (run->con->i_case==3100) {
        
        if ((crnt->cl->nz[f]<- 0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<0.0)) {
            crnt->cl->fc[f].bc=3;
        }
        if ((crnt->cl->ny[f]> 0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<0.0)){
            crnt->cl->fc[f].bc=3;
        }
        if ((crnt->cl->ny[f]> 0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.001)){
            crnt->cl->fc[f].bc=3;
        }
        if( (crnt->cl->nx[f]<-0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.0008) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<0.0012)  ){
            crnt->cl->fc[f].bc=4;
        }
        if( (crnt->cl->nx[f]>0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>-0.00008) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<0.00001)  ){
            crnt->cl->fc[f].bc=4;
        }

        }
    
      if (run->con->i_case==203) {
        if ((crnt->cl->nd[fcnd2elnd(f,1,type)].x)<0.001) {
          if (crnt->cl->ny[f]<-0.6){crnt->cl->fc[f].bc=2;}
        }
      }
      if (run->con->i_case==290) {

        if( (crnt->cl->nx[f]<-0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.0042) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<0.0047)  ){
              crnt->cl->fc[f].bc=3;
        }

/*
        if( (crnt->cl->nx[f]>0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.0142) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<0.0147)  ){
              crnt->cl->fc[f].bc=4;
        }
*/
         if( (crnt->cl->nx[f]<-0.6) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x>0.0142) && (crnt->cl->nd[fcnd2elnd(f,0,type)].x<0.0147)  ){
              crnt->cl->fc[f].bc=3;
        }

        }
 

    }
  }
}
