#include "strdata.h"

void restart_save_V2 (struct RUN * run) {

  int n,iel,iv;
  int owrt;
  struct BRANCH * brch;
  struct TREE * tree;
  char  fstring[100]="                      ";
  double *S_buff;
  int nbranch;
  int *tags_buff;
  MPI_Status status;

  FILE * rs;

  if(run->con->rank==0){
    sprintf(fstring,"restart%08d.dat", run->con->tstep);
    rs=fopen(fstring, "w");
    fprintf(rs,"%d %15.10e %d %d %d \n",run->con->tstep,run->con->time,run->topo->tleaves,run->topo->ntrees,run->con->neq);
  }
 
  
  owrt=0;
  for (iel=0;iel<run->topo->ntrees;iel++) { // Loop elements

    tree=run->topo->drys[iel];
    brch=tree->brch;

    if ((tree->part==run->con->rank) || (run->con->rank==0)) {
          
      if (tree->part==run->con->rank) {

        while (brch->nkids!=0){
          brch=brch->kids[0];
        }

        nbranch=0;
        while(brch!=NULL&&brch->root==iel){
          nbranch++;
          brch=brch->lnxt;
        }

        if (tree->part!=0) { MPI_Send(&nbranch,1,MPI_INT,0,0,MPI_COMM_WORLD); }
          
      }
          
      if (run->con->rank==0) {
        if (tree->part!=0) { MPI_Recv(&nbranch,1,MPI_INT,tree->part,0,MPI_COMM_WORLD,&status); }
      }
      
      S_buff = malloc((run->con->neq*nbranch)*sizeof(double));      
      tags_buff = malloc(nbranch*sizeof(int));
          
    }


    if (tree->part==run->con->rank) {
        
      brch=tree->brch;
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }

      nbranch=0;
      while(brch!=NULL&&brch->root==iel){
        nbranch++;
        tags_buff[nbranch] = brch->tag;
        for(n=0;n<run->con->neq;n++){ 
          S_buff[nbranch*run->con->neq+n] = brch->el->S[n];
        }
        
        brch=brch->lnxt;
      }
      
      if (tree->part!=0) {
        MPI_Send(&S_buff,nbranch*run->con->neq,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        MPI_Send(&tags_buff,nbranch,MPI_INT,0,0,MPI_COMM_WORLD);
      }
        
    }

    if (run->con->rank==0) {
        
      if (tree->part!=0) {
        MPI_Recv(&S_buff,nbranch*run->con->neq,MPI_DOUBLE,tree->part,0,MPI_COMM_WORLD,&status);
        MPI_Recv(&tags_buff,nbranch,MPI_INT,tree->part,0,MPI_COMM_WORLD,&status);
      }
      
      for(n=0;n<nbranch;n++){ 
        fprintf(rs,"%d %d ",iel,tags_buff[n]);
        for(n=0;n<run->con->neq;n++){ 
            fprintf(rs," %15.10e ",S_buff[nbranch*run->con->neq+n]);
        }
        fprintf(rs,"\n");
      }
            
    }

    if ((tree->part==run->con->rank) || (run->con->rank==0))  {
      free(S_buff);
      free(tags_buff);
    }

    if(run->con->rank==0){
      fclose(rs);
    }


  }

}