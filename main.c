/*

Models : 1: Euler
				 2: Single-phase barotropic
				 3: Multi-material barotropic
				 4: ....



*/



#include "strdata.h"

int main (int input, char **inputc) {

	int it,in,ing;

	int irank;
	int ipass;
	struct RUN  runv;
	struct RUN * run;
	struct BRANCH * crnt;
	int ifc;

	MPI_Init(&input,&inputc); 		
	PetscInitialize(&input,&inputc,0,0);

	run=&runv; 	
	createrunstruct(run);

	run->con->comm=MPI_COMM_WORLD;
	run->con->verbose=0;

	printf ("=");

	run->con->casename=inputc[1];
	run->con->runname=inputc[2];

	MPI_Comm_rank(MPI_COMM_WORLD,&run->con->rank);  
	MPI_Comm_size(MPI_COMM_WORLD,&run->con->size);  


	filenames(run); 	    
	MPI_Barrier(MPI_COMM_WORLD);

	inputset(run);
	if((run->con->verbose==1)&&(run->con->rank==0)){ printf ("MMFV: MAIN: INPUT SET \n");}
	
	if((run->con->translator==1)&&(run->con->rank==0)){ 
	  translator(run);            // writes files for mesh
	  if(run->con->verbose==1){ printf ("MMFV: MAIN: NEUTRAL FILE TRANSLATED 1 \n");}
	}
	if((run->con->translator==2)&&(run->con->rank==0)){ 
	  translator_pw(run);            // writes files for mesh
		if(run->con->verbose==1){ printf ("MMFV: MAIN: NEUTRAL FILE TRANSLATED 2 \n");}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	for (irank=0; (irank<(run->con->size)) ;irank++){// Loop all the processors
	  if (run->con->rank==irank){
	    createtopo(run); 		                 // reads mesh files
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	}
	
	
	drysconn(run->topo);
	partitiondrys(run);
	spawntopo(run);
	forestconn(run);

	MPI_Barrier(MPI_COMM_WORLD);
	run->topo->partitioned=1;

	run->con->istep=0;	
	for (irank=0;irank<run->con->size;irank++){
	  MPI_Barrier(MPI_COMM_WORLD);
    if (run->con->rank==irank){
	    boundary(run);
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	}


	createutility(run);
	watch_init(run);

	memallocation(run);

	run->con->time=0.0;
	normalvector(run);
	tagvector(run);
	meshproperties(run);

	
	setsplit(run);

  watch(run,0,0); // INIT
	run->topo->partitioned=0;

	restructtopo(run);
	run->topo->partitioned=1;

	// ToDo: Check if this is needed
	normalvector(run);
	tagvector(run);
	meshproperties(run);
	setsplit(run);

	if (run->con->start==1){
		restart_load(run);

		restructtopo(run);
       
    normalvector(run);
    tagvector(run);
    meshproperties(run);
    setsplit(run);

	}else {
		init_domain(run);
	}

	// ToDo: Check if this is needed
	if (run->con->adapth==1){
		restructtopo(run);
	  normalvector(run);
    tagvector(run);
	  meshproperties(run);
	  setsplit(run);
    communicate_S(run,0);
    criterioncalc(run);
    gradcalc(run);
    communicate_S(run,3);
	} 

  if (run->con->start==0&&run->con->geo_adapth==1){
		
    for(ipass=0;ipass<run->con->level;ipass++){
		
			if (run->con->rank==0){ printf("Geo pass %d \n",ipass); }

			init_domain(run);
			communicate_S(run,0);	

			criterionsplt(run);

			communicate_S(run,3);
			
			criterionsmth(run);
			refineh(run,ipass%2);
			restructtopo(run);
			
	    normalvector(run);
			tagvector(run);
			meshproperties(run);
	    setsplit(run);
       
   	}

		run->con->geo_adapth=0;

    init_domain(run);	
		communicate_S(run,0);
	} else {  
	  communicate_S(run,0);
  }

	//gradcalc(run);
  if (run->con->rank==0){ printf("Export 0 started \n"); }

	
	exportfield_linked_V3(run);

	//exit(0);
	//exportfield(run); 
	//exit(0);
	//export_ultra_fast(run);
  MPI_Barrier(MPI_COMM_WORLD);
 

	if (run->con->rank==0){ printf("Export 0 finished \n"); }
	loop(run);
	return 0;
}

