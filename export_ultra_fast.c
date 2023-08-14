#include "strdata.h"

void export_ultra_fast (struct RUN * run) {

  int n,iel;
  int owrt;
  struct BRANCH * brch;
  struct TREE * tree;
  char  fstring[100]="                      ";
	char  command[100]="                      ";
  FILE * rs;

	
	if (run->con->rank==0) {
		sprintf(command,"mkdir F_%08d",run->con->tstep);
		int systemRet=system(command);
		if(systemRet == -1){
			printf("System error: export_ultra_fast: mkdir F_* (not enough space?) \n");
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);


  sprintf(fstring,"F_%08d/field_%04d_%08d.dat",run->con->tstep,run->con->rank,run->con->tstep);
  rs=fopen(fstring, "w");

  fprintf(rs,"%d %15.10e %d %d %d \n",run->con->tstep,run->con->time,run->topo->tleaves,run->topo->ntrees,run->con->neq);
  fclose(rs);
  

 	owrt=0;
 	for (iel=0;iel<run->topo->ntrees;iel++) { // Loop elements

  	tree=run->topo->drys[iel];
  	brch=tree->brch;

  	if (tree->part==run->con->rank) {

    	while (brch->nkids!=0){
      	brch=brch->kids[0];
    	}

    	rs=fopen(fstring, "a");

    	while(brch!=NULL&&brch->root==iel){
      	fprintf(rs,"%d %d ",brch->root,brch->tag);
      	for(n=0;n<run->con->neq;n++){ 
        	fprintf(rs," %15.10e ",brch->el->S[n]);
      	}
      	fprintf(rs,"\n");
    		brch=brch->lnxt;
    	}
    	fclose(rs);
  	}
	
	}

}
