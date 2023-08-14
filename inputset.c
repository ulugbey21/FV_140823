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
// Called from main.c
// Read Data
//
//-------------------------------------------------------
#include "strdata.h"

void inputset(RUN *run) {

  FILE *fp;
  char *A;
  int irank;
  int i,iv,ipar;

  for (irank=0;irank<run->con->size;irank++){
    if (run->con->rank==irank){
      
      A  = malloc(1000*sizeof(char));
      fp = fopen(run->con->filecon,"r");

      if(run->con->rank==0){ printf("======================================================================\n");}

      fscanf(fp,"%s\t%d \n",A,&run->con->verbose);     if(run->con->rank==0){ printf("run->con->verbose     %d\n" ,run->con->verbose      ); }
      fscanf(fp,"%s\t%d \n",A,&run->con->i_case);      if(run->con->rank==0){ printf("run->con->i_case      %d\n" ,run->con->i_case       ); }
      fscanf(fp,"%s\t%d \n",A,&run->con->model);       if(run->con->rank==0){ printf("run->con->model       %d\n" ,run->con->model        ); }
      fscanf(fp,"%s\t%d \n",A,&run->con->ns);          if(run->con->rank==0){ printf("run->con->ns          %d\n" ,run->con->ns           ); }
   //    fscanf(fp,"%s\t%d \n",A,&run->con->turb);        if(run->con->rank==0){ printf("run->con->turb      %d\n" ,run->con->turb       ); }
      fscanf(fp,"%s\t%d \n",A,&run->con->ischeme);     if(run->con->rank==0){ printf("run->con->ischeme     %d\n" ,run->con->ischeme      ); }
      fscanf(fp,"%s\t%d \n",A,&run->con->iprimtv);     if(run->con->rank==0){ printf("run->con->iprimtv     %d\n" ,run->con->iprimtv      ); }
      fscanf(fp,"%s\t%d \n",A,&run->con->method);      if(run->con->rank==0){ printf("run->con->RK          %d\n" ,run->con->method       ); }
      fscanf(fp,"%s\t%d \n",A,&run->con->halfdt);      if(run->con->rank==0){ printf("run->con->halfdt      %d\n" ,run->con->halfdt       ); }
      fscanf(fp,"%s\t%d \n",A,&run->con->limiter);     if(run->con->rank==0){ printf("run->con->limiter     %d\n" ,run->con->limiter      ); }
      fscanf(fp,"%s\t%d \n",A,&run->con->level);       if(run->con->rank==0){ printf("run->con->level       %d\n" ,run->con->level        ); }
        if ((run->con->model==1) || (run->con->model==2)) {
        fscanf(fp,"%s\t%lf\n",A,&run->par->materpinf[0]);       if(run->con->rank==0){ printf("run->par->materpinf       %lf\n" ,run->par->materpinf[0]        ); }
        fscanf(fp,"%s\t%lf\n",A,&run->par->matergama[0]);       if(run->con->rank==0){ printf("run->par->matergama       %lf\n" ,run->par->matergama[0]        ); }
        fscanf(fp,"%s\t%lf\n",A,&run->par->matervisc[0]);       if(run->con->rank==0){ printf("run->par->matervisc       %lf\n" ,run->par->matervisc[0]        ); }
      }

      if (run->con->model==3 ||run->con->model==4) {
        
        fscanf(fp,"%s\t%d \n",A,&run->con->multimat);      if(run->con->rank==0){ printf("run->con->multimat      %d\n" ,run->con->multimat       ); }
        fscanf(fp,"%s\t%lf \n",A,&run->con->ymin);         if(run->con->rank==0){ printf("run->con->ymin          %e\n",run->con->ymin            ); }

        if(run->con->rank==0){ printf("==============================================================================\n");}
        // Reads gamma for equation of state, run->par->materpinf[] 
        if(run->con->rank==0){                                                         printf("gamma:   ");}
        fscanf(fp,"%s",A);
        for (i=0;i<run->con->multimat;i++){
          fscanf(fp,"\t%lf ",&run->par->matergama[i]);         if(run->con->rank==0){  printf(" %.2e | ",run->par->matergama[i]   ); }
        }
        fscanf(fp,"\n");                                       if(run->con->rank==0){  printf("\n");}

        // Reads p_inf for equation of state, run->par->materpinf[]
        if(run->con->rank==0){                                                         printf("p inf:   ");}
        fscanf(fp,"%s",A);
        for (i=0;i<run->con->multimat;i++){
          fscanf(fp,"\t%lf ",&run->par->materpinf[i]);         if(run->con->rank==0){  printf(" %.2e | ",run->par->materpinf[i]   ); }
        }
        fscanf(fp,"\n");                                       if(run->con->rank==0){  printf("\n");}

        // Reads initial density of the materials, run->par->matervisc[] 
        if(run->con->rank==0){                                                         printf("visc init:  ");}
        fscanf(fp,"%s",A);
        for (i=0;i<run->con->multimat;i++){
            fscanf(fp,"\t%lf ",&run->par->matervisc[i]);       if(run->con->rank==0){  printf(" %.2e | ",run->par->matervisc[i]   ); }
        }
        fscanf(fp,"\n");                                       if(run->con->rank==0){  printf("\n");}

      }
      

      fscanf(fp,"%s\t%d\n",A,&run->con->nstep);              if(run->con->rank==0){  printf("run->con->nstep          %d\n" ,run->con->nstep          ); }
      fscanf(fp,"%s\t%le\n",A,&run->con->dt);                if(run->con->rank==0){  printf("run->con->dt             %le\n",run->con->dt             ); }
      fscanf(fp,"%s\t%d\n",A,&run->con->exstep);             if(run->con->rank==0){  printf("run->con->exstep         %d\n" ,run->con->exstep         ); }
      fscanf(fp,"%s\t%d\n",A,&run->con->ststep);             if(run->con->rank==0){  printf("run->con->ststep         %d\n" ,run->con->ststep         ); }
      fscanf(fp,"%s\t%d\n",A,&run->con->geo_adapth);         if(run->con->rank==0){  printf("run->con->geo_adapth     %d\n" ,run->con->geo_adapth     ); }
      fscanf(fp,"%s\t%d\n",A,&run->con->adapth);             if(run->con->rank==0){  printf("run->con->adapth         %d\n" ,run->con->adapth         ); }
      fscanf(fp,"%s\t%d\n",A,&run->con->adapthstep);         if(run->con->rank==0){  printf("run->con->adaptstep      %d\n" ,run->con->adapthstep     ); }
      fscanf(fp,"%s\t%d\n",A,&run->con->start);              if(run->con->rank==0){  printf("run->con->start          %d\n" ,run->con->start          ); }
      fscanf(fp,"%s\t%d\n",A,&run->con->translator);         if(run->con->rank==0){  printf("run->con->translator     %d\n" ,run->con->translator     ); }

      if(run->con->rank==0){                                                         printf("bc | -x | +x | -y | +y | -z | +z | \n");printf("   ");}
      fscanf(fp,"%s",A);
      for (i=0;i<6;i++){
          fscanf(fp,"\t%d ",&run->par->bcoverwrt[i]);        if(run->con->rank==0){  printf("|  %d ",run->par->bcoverwrt[i]   );}
      }
      fscanf(fp,"\n");                                       if(run->con->rank==0){  printf("|\n");}

      fscanf(fp,"%s\t%lf\n",A,&run->con->mscale);            if(run->con->rank==0){  printf("run->con->mscale         %lf\n",run->con->mscale         ); }
      fscanf(fp,"%s\t%d\n",A,&run->con->splitmode);          if(run->con->rank==0){  printf("run->con->splitmode      %d\n" ,run->con->splitmode      ); }
      fscanf(fp,"%s\t%lf\n",A,&run->con->nxquad);            if(run->con->rank==0){  printf("run->con->nxquad         %lf\n",run->con->nxquad         ); }
      fscanf(fp,"%s\t%lf\n",A,&run->con->nyquad);            if(run->con->rank==0){  printf("run->con->nyquad         %lf\n",run->con->nyquad         ); }
      fscanf(fp,"%s\t%lf\n",A,&run->con->nzquad);            if(run->con->rank==0){  printf("run->con->nzquad         %lf\n",run->con->nzquad         ); }
    
      rewind(fp);
      fclose(fp);
      free(A);
    
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
  }


  run->con->istep=0;
  run->con->tstep=0;

  run->con->tstep0=0;
  run->con->time0=0.0;
  run->con->timetot=1000.0;

 
  if(run->con->model==1){

    run->con->nueq = 3;

    run->con->neqmass    = 1;
    run->con->neqmoment  = 3;
    run->con->neqenergy  = 1;

    run->con->neq = run->con->neqmass     +   
                    run->con->neqmoment   +    
                    run->con->neqenergy;      

    run->con->eqtypn    = malloc(run->con->nueq*sizeof(int));
    run->con->eqtypn[0] = run->con->neqmass ;      
    run->con->eqtypn[1] = run->con->neqmoment;    
    run->con->eqtypn[2] = run->con->neqenergy;      

    run->con->eqtypi    = malloc(run->con->nueq*sizeof(int));
    run->con->eqtypi[0] = 0;
    run->con->eqtypi[1] = run->con->eqtypi[0] + run->con->eqtypn[0]; 
    run->con->eqtypi[2] = run->con->eqtypi[1] + run->con->eqtypn[1]; 

    run->con->nprimitiv = run->con->eqtypn[0] + // rho
                          run->con->eqtypn[1] +  // u,v,w
                          run->con->eqtypn[2];  // p

  }

  if(run->con->model==2){

    run->con->nueq = 3;

    run->con->neqmass    = 1;
    run->con->neqmoment  = 3;

    run->con->neq = run->con->neqmass    +   
                    run->con->neqmoment;     
/*
    run->con->eqtypn    = malloc(run->con->nueq*sizeof(int));
    run->con->eqtypn[0] = run->con->neqmass ;      
    run->con->eqtypn[1] = run->con->neqmoment;     

    run->con->eqtypi    = malloc(run->con->nueq*sizeof(int));
    run->con->eqtypi[0] = 0;
    run->con->eqtypi[1] = run->con->eqtypn[0];  

    run->con->nprimitiv = run->con->eqtypn[0] + // rho
                          run->con->eqtypn[1]+
                          run->con->eqtypn[0];  // u,v,w
                          */
  run->con->eqtypn    = malloc(run->con->nueq*sizeof(int));
  run->con->eqtypn[0] = run->con->neqmass ;     // 1
  run->con->eqtypn[1] = run->con->neqmoment;      // 3
  run->con->eqtypn[2] = run->con->neqenergy;      // 1

  run->con->eqtypi    = malloc(run->con->nueq*sizeof(int));
  run->con->eqtypi[0] = 0;
  run->con->eqtypi[1] = run->con->eqtypi[0] + run->con->eqtypn[0]; // nmaterials(1)                 = 1
  run->con->eqtypi[2] = run->con->eqtypi[1] + run->con->eqtypn[1]; // nmaterials(1) + neqmoment(3)  = 4

  run->con->nprimitiv = run->con->eqtypn[0] + // rho
                        run->con->eqtypn[1] + // u,v,w
                        run->con->eqtypn[0];  // p
                      
  }

  if(run->con->model==3||run->con->model==4){

    run->con->nueq = 3;

    run->con->neqmass   = 1;
    run->con->neqmoment = 3;
    run->con->neqvf     = run->con->multimat-1;

    run->con->neq = run->con->neqmass   +   
                    run->con->neqmoment +  
                    run->con->neqvf;   

    run->con->eqtypn    = malloc(run->con->nueq*sizeof(int));
    run->con->eqtypn[0] = run->con->neqmass;      
    run->con->eqtypn[1] = run->con->neqmoment;    
    run->con->eqtypn[2] = run->con->neqvf;      

    run->con->eqtypi    = malloc(run->con->nueq*sizeof(int));
    run->con->eqtypi[0] = 0;
    run->con->eqtypi[1] = run->con->eqtypn[0]; 
    run->con->eqtypi[2] = run->con->eqtypi[1] + run->con->eqtypn[1]; 

    run->con->nprimitiv = run->con->eqtypn[0] + // rho
                          run->con->eqtypn[1] + // u,v,w
                          run->con->eqtypn[2];  // vf
  //  printf("run->con->nprimitiv:%d\n",run->con->nprimitiv);
  }
                    
  

  if (run->con->iprimtv==0) { run->con->neq_temp = run->con->neq;       }
  if (run->con->iprimtv==1) { run->con->neq_temp = run->con->nprimitiv; }

  run->con->dt0=run->con->dt;

  return;

}
