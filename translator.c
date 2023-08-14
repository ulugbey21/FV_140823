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
//
//
//-------------------------------------------------------

#include "strdata.h"

void translator(RUN *run){

  int nd,el,typ,numcon,c1,c2,c3,c4,c5,c6,c7,c8;
  int und1,und2,und3,bounds,topic_face,bccon;
  double xc,yc,zc;
  int i,x;
  FILE *flow_bc;

	//flows
  FILE *neu;
  FILE *nodes;
  FILE *elem; 
  FILE *bcont1;
  FILE *bcont2;
  FILE *bcont3;
  FILE *bcont4;

  char *neutral;

  neutral=(char *)malloc(sizeof(run->con->casename)+4);
  strcpy(neutral,run->con->casename);
  strcat(neutral,".neu");

  neu    = fopen(neutral,"r");
  nodes  = fopen(run->con->filend,"w");
  elem   = fopen(run->con->fileel,"w");
  bcont1 = fopen(run->con->filebr1,"w");
  bcont2 = fopen(run->con->filebr2,"w");
  bcont3 = fopen(run->con->filebr3,"w");
  bcont4 = fopen(run->con->filebr4,"w");

  char line[100];
  char text[50];
  
  for(i=0;i<100;i++){
    line[i]=0;
  }
  
  int ndnum,elnum,bcnum,dim,code1,code2;
  char CN_INFO[100]={' ',' ',' ',' ',' ',' ',' ',' ','C','O','N','T','R','O','L',' ','I','N','F','O',' ','2','.','4','.' ,'6' ,'\r','\n'};
  char ND_COOR[100]={' ',' ',' ','N','O','D','A','L',' ','C','O','O','R','D','I','N','A','T','E','S',' ','2','.','4','.','6','\r','\n'};
  char EL_CONN[100]={' ',' ',' ',' ',' ',' ','E','L','E','M','E','N','T','S','/','C','E','L','L','S',' ','2','.','4','.','6','\r','\n'};
  char BC_TEX[100]={' ','B','O','U','N','D','A','R','Y',' ','C','O','N','D','I','T','I','O','N','S',' ','2','.','4','.','6','\r','\n'};


  while((x=strncmp(CN_INFO,line,20))!=0){
    for(i=0;i<100;i++){
      line[i]=0;
    }
    fgets(line,100,neu); 
  };


  fgets(line,100,neu);
  fgets(line,100,neu);
  fgets(line,100,neu);
  fgets(line,100,neu);
  fgets(line,100,neu);

  fscanf(neu,"%d \t %d \t %d \t %d \t %d \t %d",&ndnum,&elnum,&code1,&bcnum,&dim,&code2);

  while((x=strncmp(ND_COOR,line,20))!=0){
    for(i=0;i<100;i++){
       line[i]=0;
    }
    fgets(line,100,neu); 
  };

  
  //Construction of .nd file        
  fprintf(nodes,"%d \n",ndnum);       
  for(i=0;i<ndnum;i++){
      fscanf(neu,"%d\t%le\t%le\t%le\n",&nd,&xc,&yc,&zc);
      fprintf(nodes,"%d\t%.10le\t%.10le\t%.10le\n",nd,xc,yc,zc);
    }//i
  fclose(nodes);
  
  while((x=strncmp(EL_CONN,line,20))!=0){
    for(i=0;i<100;i++){
       line[i]=0;
    }
   fgets(line,100,neu); 
  };

   fprintf(elem,"%d \n",elnum); 

  for(i=0;i<elnum;i++){
    fscanf(neu,"%d\t%d\t%d",&el,&typ,&numcon);

    if(numcon==8){  //Linear Hex
      fscanf(neu,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",&c1,&c2,&c3,&c4,&c5,&c6,&c7,&c8);
      fprintf(elem,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",el,typ,c1,c2,c3,c4,c5,c6,c7,c8);
    }

    if(numcon==6){  //Linear Prism
      fscanf(neu,"%d\t%d\t%d\t%d\t%d\t%d\n",&c1,&c2,&c3,&c4,&c5,&c6);
      fprintf(elem,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",el,typ,c1,c2,c3,c4,c5,c6);
    }

    if(numcon==5){  //Linear Pyramid
      fscanf(neu,"%d\t%d\t%d\t%d\t%d\n",&c1,&c2,&c3,&c4,&c5);
      fprintf(elem,"%d\t%d\t%d\t%d\t%d\t%d\t%d\n",el,typ,c1,c2,c3,c4,c5);
    }

    if(numcon==4){  //Linear Tetrahedral
      fscanf(neu,"%d\t%d\t%d\t%d\n",&c1,&c2,&c3,&c4);
      fprintf(elem,"%d\t%d\t%d\t%d\t%d\t%d\n",el,typ,c1,c2,c3,c4);
    }

  }//i

  fclose(elem);

   //Boundary conditions
  while((x=strncmp(BC_TEX,line,20))!=0){
    for(i=0;i<100;i++){
      line[i]=0;
    }
    fgets(line,100,neu);
  };

  for(bccon=0;bccon<bcnum;bccon++){

    fscanf(neu,"%s\t%d\t%d\t%d\t%d\n",text,&und1,&bounds,&und2,&und3);

    if(bccon==0){ fprintf(bcont1,"%d\n",bounds);}
    if(bccon==1){ fprintf(bcont2,"%d\n",bounds);}
    if(bccon==2){ fprintf(bcont3,"%d\n",bounds);}
    if(bccon==3){ fprintf(bcont4,"%d\n",bounds);}

    if(bccon==0){
      for(i=0;i<bounds;i++){                  
        fscanf(neu,"%d\t%d\t%d\n",&el,&typ,&topic_face);
        fprintf(bcont1,"%d\t%d\n",el,topic_face);
      }
    }//if 1

    if(bccon==1){
     for(i=0;i<bounds;i++){
       fscanf(neu,"%d\t%d\t%d\n",&el,&typ,&topic_face);
       fprintf(bcont2,"%d\t%d\n",el,topic_face);
     }
    }//if 2

    if(bccon==2){
      for(i=0;i<bounds;i++){
        fscanf(neu,"%d\t%d\t%d\n",&el,&typ,&topic_face);
        fprintf(bcont3,"%d\t%d\n",el,topic_face);
      }
    }//if 3

    if(bccon==3){
      for(i=0;i<bounds;i++){
        fscanf(neu,"%d\t%d\t%d\n",&el,&typ,&topic_face);
        fprintf(bcont4,"%d\t%d\n",el,topic_face);
      }
    }//if 1

    if(bccon!=bcnum-1){
        fgets(line,100,neu);
        fgets(line,100,neu);
      }
  }//bcon

  if(bcnum==3){
    fprintf(bcont4,"%d",0);
  }
  if(bcnum==2){
    fprintf(bcont4,"%d",0);
    fprintf(bcont3,"%d",0);
  }
  if(bcnum==1){
    fprintf(bcont4,"%d",0);
    fprintf(bcont3,"%d",0);
    fprintf(bcont2,"%d",0);
  }

  fclose(bcont1);
  fclose(bcont2);
  fclose(bcont3);
  fclose(bcont4);

 // printf("%d \t %d \t %d \t %d \t %d \t %d \n",ndnum,elnum,code1,bcnum,dim,code2);

   fclose(neu);  


   return;




 }
