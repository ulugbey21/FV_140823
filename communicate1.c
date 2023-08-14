#include "strdata.h"
void communicate_C (RUN * run) {
  struct BRANCH * crnt;
  int ipt,ipr,iv,itag,iel,ilvl,nlvl;
  int tag;
  int * adr;
//  watch(run,8,0); // COMS
  MPI_Status status;
  MPI_Request request;
  adr=malloc(10*sizeof(int));
  for (ipt=0;ipt<run->con->size;ipt++){ // loop of proxies to send
    if (ipt!=run->con->rank&&run->topo->nprox[ipt]!=0){
      for (itag=0;itag<run->topo->nprox[ipt];itag++){
        iel=run->topo->proxdrys[ipt][itag];
        tag=run->topo->proxtags[ipt][itag];
        tag2adr(adr,tag,&nlvl);
        crnt=run->topo->drys[iel]->brch;
        ilvl=0;
        while (crnt->nkids!=0){
	  crnt=crnt->kids[adr[ilvl]-1];ilvl++;
	}
        MPI_Send(&(crnt->split),1,MPI_INT,ipt,itag,MPI_COMM_WORLD);
      }
    }
    for (ipr=0;ipr<run->con->size;ipr++){ // loop of proxies to send
      if (ipt==run->con->rank&&run->topo->nbuff[ipr]!=0){
        for (itag=0;itag<run->topo->nbuff[ipr];itag++){
          iel=run->topo->buffdrys[ipr][itag];
          tag=run->topo->bufftags[ipr][itag];
          tag2adr(adr,tag,&nlvl);
          crnt=run->topo->drys[iel]->brch;
          ilvl=0;
          while (crnt->nkids!=0){crnt=crnt->kids[adr[ilvl]-1];ilvl++;}
          MPI_Recv(&(crnt->split),1,MPI_INT,ipr,itag,MPI_COMM_WORLD,&status);
	}
      }
    }
//    MPI_Barrier(MPI_COMM_WORLD);
  }
  free(adr);
//  watch(run,8,1); // COMS
}

void communicate_S (RUN * run,int ifield) {

  struct BRANCH * crnt;
  int ipt,ipr,iv,itag,iel,ilvl,nlvl,ifc,idr,ing,ib;
  int ibufsize;
  int tag;
  int * adr;
  double *sendbuff;
  double *recvbuff;
  int *isendbuff;
  int *irecvbuff;


  MPI_Status status;
  MPI_Request request;
  adr=malloc(10*sizeof(int));

  for (ipt=0;ipt<run->con->size;ipt++){ // loop of proxies to send

    if (ipt==run->con->rank){

      for (ipr=0;ipr<run->con->size;ipr++){ // loop proxies to receive

        ib=0;
        if (ifield==0){ ibufsize = run->topo->nprox[ipr]*(run->con->neq);       }  // S
        if (ifield==1){ ibufsize = run->topo->nprox[ipr]*6*run->con->neq;       }  // SF
        if (ifield==2){ ibufsize = run->topo->nprox[ipr]*3*(run->con->neq+1) ;  }  // grad
        if (ifield==3){ ibufsize = run->topo->nprox[ipr];                       }  // criterion, crnt->split
        if (ifield==4){ ibufsize = run->topo->nprox[ipr];                       }  // part

        if (run->topo->nprox[ipr]!=0) {

	        if (ifield<3){ 
            sendbuff=malloc(ibufsize*sizeof(double)); 
          } else { 
            isendbuff=malloc(ibufsize*sizeof(int));
          }
        
          for (itag=0;itag<run->topo->nprox[ipr];itag++){
            iel=run->topo->proxdrys[ipr][itag];
            tag=run->topo->proxtags[ipr][itag];
            tag2adr(adr,tag,&nlvl);
            crnt=run->topo->drys[iel]->brch;
            ilvl=0;
            while (crnt->nkids!=0){
              crnt=crnt->kids[adr[ilvl]-1];ilvl++;
            }
            
            for (iv=0;iv<(run->con->neq);iv++){
              if (ifield==0) {
                sendbuff[ib]=crnt->el->S[iv];
                ib++;
              }
            }

            if (ifield==1) {
              for (ifc=0;ifc<6;ifc++){    
                for (iv=0;iv<run->con->neq;iv++){
                  sendbuff[ib]=crnt->el->SF[ifc][iv];
                  ib++;
                }
  	          }
            }
            
            if (ifield==2){
              for (iv=0;iv<(run->con->neq+1);iv++){
                for (idr=0;idr<3;idr++){
                  sendbuff[ib]=crnt->el->SG[iv][idr];
                  ib++;
                }
              }
            }

            if (ifield==3) {
              isendbuff[ib]=crnt->split;
              ib++;
            }

            if (ifield==4) {
              isendbuff[ib]=crnt->part;
              ib++;
            }
            
          }

          if (ifield<3){ 
            MPI_Send(sendbuff,ibufsize,MPI_DOUBLE,ipr,ipt,MPI_COMM_WORLD);
          } else { 
            MPI_Send(isendbuff,ibufsize,MPI_INT,ipr,ipt,MPI_COMM_WORLD);
          }
	        if (ifield<3){ 
            free(sendbuff);
          } else { 
            free(isendbuff);
          }

	      }
        //if (ifield==0){MPI_Send(run->SendS,ibufsize,MPI_DOUBLE,ipr,itag,MPI_COMM_WORLD);}
        //if (ifield==1){MPI_Send(run->SendSF,ibufsize,MPI_DOUBLE,ipr,itag,MPI_COMM_WORLD);}
        //if (ifield==2){MPI_Send(run->SendSG,ibufsize,MPI_DOUBLE,ipr,itag,MPI_COMM_WORLD);}
      }
    }

    if (ipt!=run->con->rank&&run->topo->nbuff[ipt]!=0) {

      if (ifield==0){ ibufsize = run->topo->nbuff[ipt]*(run->con->neq);      }
      if (ifield==1){ ibufsize = run->topo->nbuff[ipt]*6*run->con->neq;      }
      if (ifield==2){ ibufsize = run->topo->nbuff[ipt]*3*(run->con->neq+1);  }
      if (ifield==3){ ibufsize = run->topo->nbuff[ipt];                      }
      if (ifield==4){ ibufsize = run->topo->nbuff[ipt];                      }

	    if (ifield<3){ 
        recvbuff=malloc(ibufsize*sizeof(double));
      } else {
        irecvbuff=malloc(ibufsize*sizeof(int));
      }
      if (ifield<3){ 
        MPI_Recv(recvbuff,ibufsize,MPI_DOUBLE,ipt,ipt,MPI_COMM_WORLD,&status);
      } else { 
        MPI_Recv(irecvbuff,ibufsize,MPI_INT,ipt,ipt,MPI_COMM_WORLD,&status);
      }

      ib=0;
      for (itag=0;itag<run->topo->nbuff[ipt];itag++) {
        iel=run->topo->buffdrys[ipt][itag];
        tag=run->topo->bufftags[ipt][itag];
        tag2adr(adr,tag,&nlvl);
        crnt=run->topo->drys[iel]->brch;
        ilvl=0;
        
        while (crnt->nkids!=0){
          crnt=crnt->kids[adr[ilvl]-1];ilvl++;
        }       
        
        if (ifield==0){
          for (iv=0;iv<(run->con->neq);iv++){
            crnt->el->S[iv] = recvbuff[ib];
            ib++;
          }
        }

        
        if (ifield==1){
          for (ifc=0;ifc<6;ifc++){    
            for (iv=0;iv<run->con->neq;iv++){
              crnt->el->SF[ifc][iv] = recvbuff[ib];
              ib++;
            }
          }
        }

        if (ifield==2){
          for (iv=0;iv<(run->con->neq+1);iv++){
            for (idr=0;idr<3;idr++){
              crnt->el->SG[iv][idr] = recvbuff[ib];
              ib++;
            }
          }
        }

        if (ifield==3){ 
          crnt->split=irecvbuff[ib];
          ib++; 
        }


        if (ifield==4) {
          crnt->part=irecvbuff[ib];
          ib++;
        }

      }

      if (ifield<3){
        free(recvbuff);
      } else {
        free(irecvbuff);
      }

    }

    // MPI_Barrier(MPI_COMM_WORLD);
  }

  free(adr);
}



void tag2adr(int * n,int tag,int  * pnlvl){
  int l,t;
  int ilvl,nlvl;
  int ifc,ind,iel2;
  int match;
  
  nlvl=1;
  if (tag>0){nlvl=(((int) floor(log((double) tag)/log(10.0))+2));}
  for (ilvl=0;ilvl<nlvl;ilvl++){
    n[ilvl]=0;
  }
  t=tag;
  l=0;
  while (tag>(pow(10,l)-1)){
    n[l]=t%((int) pow(10,l+1))/((int) pow(10,l));
    t-=n[l] * ((int) pow(10,l));
    l++;
  }
  *pnlvl=nlvl;
}

