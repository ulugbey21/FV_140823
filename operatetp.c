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
int numlfc(int typ){
  if (typ==0){return 6;}
  if (typ==1){return 4;}
  if (typ==2){return 5;}
}


int numlnd(int typ){
  if (typ==0){return 8;}
  if (typ==1){return 4;}
  if (typ==2){return 6;}
}


int numfcnd(int ifc, int typ){
  if (typ==0){return 4;}
  if (typ==1){return 3;}
  if (typ==2&&ifc<3){return 4;}
  if (typ==2&&ifc>2){return 3;}
}


int elnd2fcnd (int ifc, int ind, int typ)
{

int M[6][8];

if (typ==0){

M[0][0] = 0; M[0][1] = 1; M[0][2] =-1; M[0][3] =-1; M[0][4] = 3; M[0][5] = 2; M[0][6] =-1; M[0][7] =-1;

M[1][0] =-1; M[1][1] = 3; M[1][2] =-1; M[1][3] = 0; M[1][4] =-1; M[1][5] = 2; M[1][6] =-1; M[1][7] = 1;

M[2][0] =-1; M[2][1] =-1; M[2][2] = 0; M[2][3] = 3; M[2][4] =-1; M[2][5] =-1; M[2][6] = 1; M[2][7] = 2; 

M[3][0] = 1; M[3][1] =-1; M[3][2] = 0; M[3][3] =-1; M[3][4] = 2; M[3][5] =-1; M[3][6] = 3; M[3][7] =-1;

M[4][0] = 3; M[4][1] = 2; M[4][2] = 0; M[4][3] = 1; M[4][4] =-1; M[4][5] =-1; M[4][6] =-1; M[4][7] =-1; 

M[5][0] =-1; M[5][1] =-1; M[5][2] =-1; M[5][3] =-1; M[5][4] = 1; M[5][5] = 2; M[5][6] = 0; M[5][7] = 3;
}

if (typ==1){

M[0][0] = 1; M[0][1] = 0; M[0][2] = 2; M[0][3] =-1; M[0][4] = 3; M[0][5] =-1; M[0][6] =-1; M[0][7] =-1;

M[1][0] = 0; M[1][1] = 1; M[1][2] =-1; M[1][3] = 2; M[1][4] =-1; M[1][5] = 3; M[1][6] =-1; M[1][7] =-1;

M[2][0] =-1; M[2][1] = 0; M[2][2] = 1; M[2][3] = 2; M[2][4] =-1; M[2][5] =-1; M[2][6] = 3; M[2][7] =-1;

M[3][0] = 1; M[3][1] =-1; M[3][2] = 0; M[3][3] = 2; M[3][4] =-1; M[3][5] =-1; M[3][6] =-1; M[3][7] = 3;


}


if (typ==2){

M[0][0] = 0; M[0][1] = 1; M[0][2] =-1; M[0][3] = 3; M[0][4] = 2; M[0][5] =-1; M[0][6] =-1; M[0][7] =-1;

M[1][0] =-1; M[1][1] = 3; M[1][2] = 0; M[1][3] =-1; M[1][4] = 2; M[1][5] = 1; M[1][6] =-1; M[1][7] =-1;

M[2][0] = 0; M[2][1] =-1; M[2][2] = 3; M[2][3] = 1; M[2][4] =-1; M[2][5] = 2; M[2][6] =-1; M[2][7] =-1;

M[3][0] = 0; M[3][1] = 2; M[3][2] = 1; M[3][3] =-1; M[3][4] =-1; M[3][5] =-1; M[3][6] = 3; M[3][7] =-1; 

M[4][0] =-1; M[4][1] =-1; M[4][2] =-1; M[4][3] = 0; M[4][4] = 1; M[4][5] = 2; M[4][6] =-1; M[4][7] = 3;

M[5][0] =-1; M[5][1] =-1; M[5][2] =-1; M[5][3] =-1; M[5][4] =-1; M[5][5] =-1; M[5][6] =-1; M[5][7] =-1; 
}

//if (M[ifc][ind]==-1){printf("oups elnd2fcnd \n"); exit(-1);}

  return M[ifc][ind];
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


int fcnd2elnd (int ifc, int ifcnd,int typ)
{

int M[6][4];
if (typ==0){
M[0][0]=0;
M[0][1]=1;
M[0][2]=5;
M[0][3]=4;

M[1][0]=3;
M[1][1]=7;
M[1][2]=5;
M[1][3]=1;

M[2][0]=2;
M[2][1]=6;
M[2][2]=7;
M[2][3]=3;

M[3][0]=2;
M[3][1]=0;
M[3][2]=4;
M[3][3]=6;

M[4][0]=2;
M[4][1]=3;
M[4][2]=1;
M[4][3]=0;

M[5][0]=6;
M[5][1]=4;
M[5][2]=5;
M[5][3]=7;
}

if (typ==1){
M[0][0]=1;
M[0][1]=0;
M[0][2]=2;
M[0][3]=4;

M[1][0]=0;
M[1][1]=1;
M[1][2]=3;
M[1][3]=5;

M[2][0]=1;
M[2][1]=2;
M[2][2]=3;
M[2][3]=6;

M[3][0]=2;
M[3][1]=0;
M[3][2]=3;
M[3][3]=7;

M[4][0]=-1;
M[4][1]=-1;
M[4][2]=-1;
M[4][3]=-1;

M[5][0]=-1;
M[5][1]=-1;
M[5][2]=-1;
M[5][3]=-1;
}


if (typ==2){
M[0][0]=0;
M[0][1]=1;
M[0][2]=4;
M[0][3]=3;

M[1][0]=2;
M[1][1]=5;
M[1][2]=4;
M[1][3]=1;

M[2][0]=0;
M[2][1]=3;
M[2][2]=5;
M[2][3]=2;

M[3][0]=0;
M[3][1]=2;
M[3][2]=1;
M[3][3]=6;

M[4][0]=3;
M[4][1]=4;
M[4][2]=5;
M[4][3]=7;

M[5][0]=-1;
M[5][1]=-1;
M[5][2]=-1;
M[5][3]=-1;
}

if (M[ifc][ifcnd]==-1){printf("oups fcnd2elnd\n"); exit(-1);}

  return M[ifc][ifcnd];
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


int fcndrot (int ifcnd, int iangl,int typ)
{

int M[6][4];


if (typ==0){
M[0][0]=0;
M[1][0]=1;
M[2][0]=2;
M[3][0]=3;

M[0][1]=1;
M[1][1]=2;
M[2][1]=3;
M[3][1]=0;

M[0][2]=2;
M[1][2]=3;
M[2][2]=0;
M[3][2]=1;

M[0][3]=3;
M[1][3]=0;
M[2][3]=1;
M[3][3]=2;
}

if (typ==1){

M[0][0]=0;
M[1][0]=1;
M[2][0]=2;
M[3][0]=3;

M[0][1]=1;
M[1][1]=2;
M[2][1]=0;
M[3][1]=3;

M[0][2]=2;
M[1][2]=0;
M[2][2]=1;
M[3][2]=3;

}

  if (ifcnd>3||iangl>3||iangl<0||ifcnd<0){printf("rot %d %d %d \n",ifcnd,iangl,typ);}
  return M[ifcnd][iangl];
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


int fcopnd2elnd (int ifc, int ifcnd,int type)
{

int M[6][4];
if (type==0){
M[0][0]=2;
M[0][1]=3;
M[0][2]=7;
M[0][3]=6;

M[1][0]=2;
M[1][1]=6;
M[1][2]=4;
M[1][3]=0;

M[2][0]=0;
M[2][1]=4;
M[2][2]=5;
M[2][3]=1;

M[3][0]=3;
M[3][1]=1;
M[3][2]=5;
M[3][3]=7;

M[4][0]=6;
M[4][1]=7;
M[4][2]=5;
M[4][3]=4;

M[5][0]=2;
M[5][1]=0;
M[5][2]=1;
M[5][3]=3;
}

if (type==1){

M[0][0]=3;
M[0][1]=3;
M[0][2]=3;
M[0][3]=3;

M[1][0]=2;
M[1][1]=2;
M[1][2]=2;
M[1][3]=2;

M[2][0]=0;
M[2][1]=0;
M[2][2]=0;
M[2][3]=0;

M[3][0]=1;
M[3][1]=1;
M[3][2]=1;
M[3][3]=1;

M[4][0]=-1;
M[4][1]=-1;
M[4][2]=-1;
M[4][3]=-1;

M[5][0]=-1;
M[5][1]=-1;
M[5][2]=-1;
M[5][3]=-1;



}


if (type==2){

M[0][0]=2;
M[0][1]=1;
M[0][2]=4;
M[0][3]=5;

M[1][0]=0;
M[1][1]=3;
M[1][2]=4;
M[1][3]=1;

M[2][0]=1;
M[2][1]=4;
M[2][2]=4;
M[2][3]=1;


M[3][0]=3;
M[3][1]=5;
M[3][2]=4;
M[3][3]=7;

M[4][0]=0;
M[4][1]=1;
M[4][2]=2;
M[4][3]=6;

M[5][0]=-1;
M[5][1]=-1;
M[5][2]=-1;
M[5][3]=-1;



}


  return M[ifc][ifcnd];
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
 int  oppface (int i, int ifc,int type)
 {
 
   int fc;
   if (type==0){
     if (ifc==0){fc=2;}
     if (ifc==1){fc=3;}
     if (ifc==2){fc=0;}
     if (ifc==3){fc=1;}
     if (ifc==4){fc=5;}
     if (ifc==5){fc=4;}
   }
   
   if (type==1){
     if (ifc==0){if(i==0){fc=1;};if(i==1){fc=2;};if(i==2){fc=3;};}
     if (ifc==1){if(i==0){fc=2;};if(i==1){fc=3;};if(i==2){fc=0;};}
     if (ifc==2){if(i==0){fc=3;};if(i==1){fc=0;};if(i==2){fc=1;};}
     if (ifc==3){if(i==0){fc=0;};if(i==1){fc=1;};if(i==2){fc=2;};}
   }
   
   
   
   if (type==2){
     if (ifc==0){if(i==0){fc=1;};if(i==1){fc=2;};}
     if (ifc==1){if(i==0){fc=2;};if(i==1){fc=0;};}
     if (ifc==2){if(i==0){fc=0;};if(i==1){fc=1;};}
     if (ifc==3){if(i==0){fc=4;};}
     if (ifc==4){if(i==0){fc=3;};}
   }
   return fc;
 }
 
 int noppface (int ifc,int type)
 {
   int n;
   if (type==0){n=1;}
   if (type==1){n=3;}
   if (type==2){n=1; if(ifc<3){n=2;} }
   return n;
 }


int fcopp (int ifc,int type)
{

int M[6];
if (type==0){
M[0]=2;
M[1]=3;
M[2]=0;
M[3]=1;
M[4]=5;
M[5]=4;
}

if (type==1){
M[0]=1;
M[1]=2;
M[2]=3;
M[3]=0;
M[4]=-1;
M[5]=-1;
}



if (type==2){
M[0]=1;
M[1]=2;
M[2]=0;
M[3]=4;
M[4]=3;
M[5]=-1;
}
  return M[ifc];
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


int fcrefnd2elnd (int ifc, int ifcnd,int typ)
//reflected fc node to element nd
{

int M[6][4];

if (typ==0){

M[0][0]=0;
M[0][3]=1;
M[0][2]=5;
M[0][1]=4;

M[1][0]=3;
M[1][3]=7;
M[1][2]=5;
M[1][1]=1;

M[2][0]=2;
M[2][3]=6;
M[2][2]=7;
M[2][1]=3;

M[3][0]=2;
M[3][3]=0;
M[3][2]=4;
M[3][1]=6;

M[4][0]=2;
M[4][3]=3;
M[4][2]=1;
M[4][1]=0;

M[5][0]=6;
M[5][3]=4;
M[5][2]=5;
M[5][1]=7;
}

if (typ==1){

M[0][0]=1;
M[0][1]=2;
M[0][2]=0;
M[0][3]=4;

M[1][0]=0;
M[1][1]=3;
M[1][2]=1;
M[1][3]=5;

M[2][0]=1;
M[2][1]=3;
M[2][2]=2;
M[2][3]=6;

M[3][0]=2;
M[3][1]=3;
M[3][2]=0;
M[3][3]=7;

if (ifc<0||ifc>3||ifcnd>3||ifcnd<0){printf ("oups fcrefnd %d %d \n",ifc,ifcnd);}
//printf ("fcrefnd %d %d \n",ifc,ifcnd);
}



if (typ==2){

M[0][0]=0;
M[0][1]=3;
M[0][2]=4;
M[0][3]=1;

M[1][0]=2;
M[1][1]=1;
M[1][2]=4;
M[1][3]=5;

M[2][0]=0;
M[2][1]=2;
M[2][2]=5;
M[2][3]=3;

M[3][0]=0;
M[3][1]=1;
M[3][2]=2;
M[3][3]=6;

M[4][0]=3;
M[4][1]=5;
M[4][2]=4;
M[4][3]=7;

  if (ifc<0||ifc>4||ifcnd>3||ifcnd<0){printf ("oups fcrefnd %d %d \n",ifc,ifcnd);}



}


  return M[ifc][ifcnd];
}
