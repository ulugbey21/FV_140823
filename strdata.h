#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stddef.h>
#include<time.h>
#include<math.h>
#include<sys/types.h>
#include<sys/resource.h>
#include<float.h>
#include<mpi.h>

#include"metis.h"
#include"parmetis.h"

#include "petscsys.h"
#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscsnes.h"
#include "petscvec.h"
#include "petscts.h"


/* Define macros */
#ifndef max 
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#ifndef STOP
#define STOP 300
#endif

#ifndef EPS
#define EPS 100*DBL_EPSILON
#endif

#ifndef sign
#define sign(a,b) ((b)<0 ? -fabs(a) : fabs(a))
#endif

#define printmat3_fmt(A)	#A " =\n%g %g %g\n%g %g %g\n%g %g %g\n"

#define printmat3(A)	( printf(printmat3_fmt(A), 					\
	((double *)(A))[3*0+0], ((double *)(A))[3*1+0], ((double *)(A))[3*2+0],		\
	((double *)(A))[3*0+1], ((double *)(A))[3*1+1], ((double *)(A))[3*2+1], 	\
	((double *)(A))[3*0+2], ((double *)(A))[3*1+2], ((double *)(A))[3*2+2]) )

//PetscErrorCode (*function) (SNES,Vec,Vec,void*);

// =====================================================================
//                                GLOBAL
// =====================================================================


// =====================================================================


typedef struct ND{
    int num;
    double x;
    double y;
    double z;
}ND;

typedef struct FC{
    int num;
    int type;
    int nnds;
    int bc;
}FC;

typedef struct LEAF{

  int gnum;
  double *RHS;
  double *S;
  double *SN;
  double **SF;
  double **SG;
  double crith;


  // Face only loop
  double **flux_face;
  int *flux_flag;

  int inner;

}LEAF;

// 76 bytes x 10M x 32 tasks = 24.32 G
typedef struct TREE{
  int type;               // 4 bytes
  int * vrtx;             // 32 bytes
  int * conf;             // 24 bytes max
  int part;             // 4  bytes max
  struct BRANCH * brch;   // 8 bytes
}TREE;

typedef struct BRANCH{
  int igroup;

  int part;
  int oldpart;
  int level;
  int type;
  int nkeens;
  int nkids;
  int nlnd;
  int nlfc;
  int tot_neigs;
  int * adrs;
  int split;
  int merge;
  int tag;
  int hangn;
  struct BRANCH * prnt;   // Pointer to parent tree
  int root;     // Pointer to Ancestor tree

  struct BRANCH * next;   // Pointer to next tree (NULL if last)
  struct BRANCH * lnxt;   // 
  struct BRANCH * lprv;   // 

  struct BRANCH *** neigtr;
  int ** neignd;
  int * neigfc;
  int * neigag;
  int * lfc_neigs;
  int * nsfc;
  int ** ipartbound;

  struct BRANCH ** keen;
  int * keenfc;
  int * keenfcangle;

  int fcqsplit;
  int nsplit;
  int ctfc;

  struct BRANCH ** kids;
  struct LEAF * el;
  struct CELL * cl;

}BRANCH;

typedef struct CELL{
    double *nx;
    double *ny;
    double *nz;
    double *nxt1;
    double *nyt1;
    double *nzt1;
    double *nxt2;
    double *nyt2;
    double *nzt2;
    double Vol;
    double **Area;
    double xc;
    double yc;
    double zc;

    struct ND * nd;
    struct FC * fc;
}CELL;

typedef struct FOREST{
    int partitioned;
    int ntrees;
    int ndepth;
    int nleaves;
    int * lleaves;
    int pleaves;
    int tleaves;
    int nleavestmp;
    int ngnd;   
    int ** proxlfc1;
    int ** proxlfc2;
    int ** proxtags;
    int ** proxdrys;
    int * nprox;
    int ** bufflfc1;
    int ** bufflfc2;
    int ** bufftags;
    int ** buffdrys;
    int * nbuff;
    struct TREE ** drys; 
    struct BRANCH * locl;
    struct BRANCH * glob;
    double * vertex;   // 3 x8 bytes x10M x32 = 7.680G
    double * vertey;   // 
    double * vertez;   // 
    
}FOREST;

typedef struct PAR{        // Struck for multimaterial parameters
    double matervisc[10];
    int bcoverwrt[10];
    double materpinf[10];
    double matermush[10];
    double matergama[10];
    double materrini[10];
    double materpini[10];
    double u_init[10];
    double v_init[10];
    double w_init[10];
}PAR;

typedef struct CON{
    MPI_Comm comm;
    
    int level;
    int ns;
    int turb;
    int i_case;

    int iter;
    int nprimitiv;
    int geo_adapth;

    double ymin;
    
    double r;
    double u;
    double v;
    double w;
    double e;
    double p;
    double g;
  
    
    double ** amat;
    double ** bmat;
    double ** cmat;
    double ** xmat;
    double * avec;
    double * bvec;
    double * cvec;

    int model;
    int multimat;
    int neqvf;
    int nstep;
    double dt;
    double dt0;
    int exstep;
    int ststep;
    int adapth;
    int adapthstep;
    int start;
    int translator;
    double mscale;
    int splitmode;
    double nxquad,nyquad,nzquad;
    int eos;
    int verbose;
    
    double RE;

    int halfdt;
    int limiter;
    int istep;
    int tstep;
    int tstep0;
    double time;
    double time0;
    double timetot;
    int rank;
    int size;
    int hmaxlvl;

    double Mach;

    int neq;        // number of equations in total
    int nueq;
 
    int nfleq;
    int ndiat;
    int nreac;


    int nlimit;
    char *  casename;
    char *  runname;

    double maxM1;
    double maxM2;
    double maxM3;
    double maxM4;
    double maxM5;


    int cnamelength;

    char *filend;          //Pointer for filename which includes nodes info
    char *fileel;          //Pointer for filename which includes elements info
    char *filecon;         //Pointer for filename which includes control values
    char *filebr1;         //Pointer for filename which includes boundary conditions 1 info
    char *filebr2;         //Pointer for filename which includes boundary conditions 2 info
    char *filebr3;         //Pointer for filename which includes boundary conditions 3 info
    char *filebr4;         //Pointer for filename which includes boundary conditions 4 info

    char *filebr5;         //Pointer for filename which includes boundary conditions 4 info
    char *filesor;         //Pointer for filename which include elements with a property
    char *filev1;         //Pointer for filename which includes boundary conditions 4 info
    char *filev2;         //Pointer for filename which include elements with a property

 

    double stepcomptime;
    double stepcommtime;
    double stepadaptime;
    double timewatch0[20];
    double timewatcht[20];
    double timewatchtall[20];
    int rhsnan;

    int spadr[10];
    int sdadr[10];
    int readr[10];
    double stats[5][10][50][50];
    double statsall[5][10][50][50];
    double stycrd[50];
    double stzcrd[50];
    double stycrdall[50];
    double stzcrdall[50];
    double srcmax;
    double srcmin;


    int nmconseq;   // number of conservative equations in total
    int neqmass;
    int neqmoment;
    int neqenergy;

    int neq_temp;
    int *eqtypn;
    int *eqtypi;

    double * limitedvars;
    double umax;
    int ischeme;
    int irecons;
    int iprimtv;
    int method;
    int ctype;

    
    
    int nindm;
    int idmtp[10];
    int idmmt[10];
    double dmpar[10][10];

}CON;

typedef struct SOL{
    double r;        //density
    double u;        //velocity u
    double v;        //velocity v
    double w;        //velocity w
    double e;         //Total energy
    double e_internal;//energy
    double p_hydro;

		double u_temp;
		double v_temp;
		double w_temp;


    double **st;
    double c;

    double **flux;
    double *vec;
    double *wec;

    double ux;
    double uy;
    double uz;
    double vx;
    double vy;
    double vz;
    double wx;
    double wy;
    double wz;

    double u_vel;
    double v_vel;
    double w_vel;
    double mu_turb;
    
    double ymass;
    double vol_frac;

    double r_lm;
    double r_gas;

}SOL;


typedef struct RUN{
    

    double *SendS;
    double *RecvS;
    double *SendSG;
    double *RecvSG;
    double *SendSF;
    double *RecvSF;
    double ***val;
    double ****val4;
    int counter;

    struct SOL    * solS;
    struct SOL    * solS1;
    struct SOL    * sol;
    struct SOL    * sol1;
    struct CON    * con;
    struct INFO   * info;
    struct PAR    * par;
    struct FOREST * topo;

    double *AUX_S;         // Auxilliary variable for volmvalues function

    double * vec_reconstruct;


    double **tr2;

    double * delta0;
    double * delta1;
    double ** delta1_face;
    double ** delta0_face;
    double ** vdelta1;
    double * delta2;
    double * deltal;
    double * deltat;
    double * gradlt;
    double ** gradlt_face;
    double * vecl;
    double * wecl;
    double * fluxs;
    
  
    double ** vecx0;
    double ** vecph;
    double ** veccf;
    double * xx;

    // Auxiliary variables residual

    
    double **flux_adv;
    double **flux_viscous;
    double *flux_HLLC;
   // double *flux_mc;
    double **flux_mc;

    double * wec_temp;
    double * vec_temp;

    // Auxiliary variables for hllc + hllc_br
    int flag_ifs2;
    double r_S;
    double e_S;
    double *Us_temp;

    double * tensor_temp;
    double *vel_temp;       


    //Auxiliary for MUSCL
    double ** wec_face;
    double ** wec_face_halfdt;


    int flag_exp;
    int exp_debug;
    
    int * exp;
    int * exp_part;
    
    double * exp_RHS;
    double ** exp_flux;
    double ** exp_fl_hllc;
    int ** exp_hllc;
    int ** exp_region;
    double ** exp_Lar2;
    double ** exp_Rar2;
    double ** exp_grad;
    double ** exp_grad_a2;
    double ** exp_grad_r2;

    int deb;
    int adapt_flag_sensor;

  // Timing variables
  double cputime;
  double stepcputime;
  double Ttot;
  double Tmeshadapt;
  double Tadv;
  double T2nd;
    double Tdist;
    double Tdelta;
    double Trec;
    double T2ndrelax;
    double T2ndcom;
    double T2ndface;

    double Tdist_sum;
    double Tdelta_sum;
    double Trec_sum;
    double T2ndrelax_sum;
    double T2ndcom_sum;
    double T2ndface_sum;
    
  double Trhs;
  double Trelax;
  double Tcomm;
  double Tinterfloc;
  double Tepx;

  double stepcputime_sum;
  double Ttot_sum;
  double Tmeshadapt_sum;
  double Tadv_sum;
  double T2nd_sum;
  double Trhs_sum;
  double Trelax_sum;
  double Tcomm_sum;
  double Tinterfloc_sum;
  double Tepx_sum;
  int tot_leaves;

}RUN;


long int treetag (struct TREE * tree);


struct TREE * createtree (int type, int num, int nkids, int nlnd,int nlvs);



#define RESET		0
#define BRIGHT 		1
#define DIM		2
#define UNDERLINE 	3
#define BLINK		4
#define REVERSE		7
#define HIDDEN		8

#define BLACK 		0
#define RED		1
#define GREEN		2
#define YELLOW		3
#define BLUE		4
#define MAGENTA		5
#define CYAN		6
#define	WHITE		7

void textcolor(int attr, int fg, int bg);
void advancetimeexplicit(struct RUN * run);
void advancetimeexplicit_nb(struct RUN * run);
void boundary(RUN *run);
void calctreegm (struct RUN * run);
void statistics_calc (int id, struct TREE * tree, struct RUN * run);
void statinit (struct RUN * run);
void statsave (struct RUN * run);
void commmaster   (RUN * run);
void communicate_S (RUN * run, int ifield);
void communicate_nb_S (struct RUN * run,int ifield,int i_action);
void communicate_nb2_S (RUN * run,int ifield);
void communicate_C (RUN * run);
void createlfc (struct BRANCH * tree);
void createrunstruct (struct RUN * run);
void createtopo (RUN * run);
void restructtopo (RUN * run);
struct BRANCH * createbranch (int type, int num, int nkids, int nlnd,int nlvs);
void createutility(RUN *run);
void elconn (struct RUN * run, struct BRANCH * brch,int ifc);
int elnd2fcnd (int ifc, int ind, int typ);
void exportfield (struct RUN * run);
void exportfield_linked (struct RUN * run);
void facevalues(struct BRANCH * crnt,int in,int ifc,int iside, int istore,RUN *run);
void facevalues_MUSCL(struct BRANCH * crnt,int ifc,int iside, struct SOL* sol_temp,struct RUN *run);
void facercnstr(struct BRANCH * crnt,int in,int ifc,int iside, int istore,RUN *run);
int fcnd2elnd (int ifc, int ifcnd,int typ);
int fcndrot (int ifcnd, int iangl,int typ);
int fcopnd2elnd (int ifc, int ifcnd,int type);
int fcopp (int ifc,int type);
int noppface (int ifc,int type);
int oppface (int i,int ifc,int type);
int fcrefnd2elnd (int ifc, int ifcnd,int typ);
int fcrefnd2Prismkid (int ifc, int ifcnd);
void filenames(RUN * run);

void rotate_vector (double * vel_temp,double u,double v,double w,double nx,double ny,double nz,
                                                               double mx,double my,double mz,
                                                               double lx,double ly,double lz,int flag_write,int pros);

void rotate_tensor (double * tensor_temp,double ** st_temp,double nx,double ny,double nz,
                                                     double mx,double my,double mz,
                                                     double lx,double ly,double lz,int flag_write,int pros);

void forestconn (RUN * run);

void inputset(RUN *run);
void leafallocation(RUN * run,struct BRANCH *crnt);
void leafdeallocation(RUN * run,struct BRANCH *crnt);
void loop (RUN * run);
int main (int input, char **inputc);
void memallocation(RUN *run);
void split (int linked, RUN * run, BRANCH * tree);
void merge (int linked, RUN * run, BRANCH * tree);
void barotropic_eos(struct BRANCH * crnt,int istore,RUN *run);
void barotropic_multiphase_eos(struct BRANCH * crnt,int istore,RUN *run);
void barotropic_threephase_eos(struct BRANCH * crnt,int istore,RUN *run);
void tagvector(RUN *run);
void tagvectorcalc(RUN *run,struct BRANCH * crnt);
void textcolor(int attr, int fg, int bg);
void translator(RUN *run);
void translator_pw(RUN *run);
void brchconn (BRANCH * tree);

void init_domain(RUN *run);
void Init_S(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w,double p);
void Init_S3(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w,double ymass);

void watch(struct RUN *run, int iw, int fun);
void watch_init(struct RUN *run);
void watch_display(struct RUN *run, int tag);
void normalvector(RUN *run);
void normalvectorcalc(RUN *run,struct BRANCH * crnt);
void volmvalues(struct BRANCH * crnt,RUN *run);
void criterioncalc (struct RUN* run);
void gradcalc (struct RUN* run);
void gradcalc_V2secondorder (struct RUN* run);
void criterionsmth (struct RUN* run);
void criterionsplt (struct RUN* run);
void meshadaptation (RUN * run);
void normalvector(RUN *run);
void residuals(RUN * run);

void neigsfin (RUN * run);
void refineh (RUN * run, int imode);
void setsplit(RUN *run);
void setsplitcalc(struct BRANCH * crnt,RUN *run);
void RHS(struct RUN * run);
void RHS_nb(struct RUN * run,int inner_outer);
void RK(int v,RUN *run);
void neigs (RUN * run);
void residual(struct BRANCH * crnt,RUN * run);
void restart_load (struct RUN *run);
void restart_save (struct RUN *run);
void restart_save_V2 (struct RUN * run);

void bc(struct BRANCH * crnt, int p, int f, struct RUN *run);
void bc_inlet(struct BRANCH * crnt,int f,struct RUN *run);
void bc_reflective(struct BRANCH * crnt,int f,RUN *run);
void bc_wall(struct BRANCH * crnt,int f,struct RUN *run);

void init_domain(RUN *run);

void faceinterp_hex (struct RUN * run);
void faceinterp_linear (struct RUN * run);

//void facevalues_MUSCL(struct BRANCH * crnt,int ifc,int iside, struct SOL* sol_temp,struct RUN *run);

void faceinterp (struct RUN * run);

void distance_pp(struct RUN *run,struct BRANCH *crnt,int f,double* dx_cc,double* dy_cc,double* dz_cc);
void distance_pf(struct RUN *run,struct BRANCH *crnt,int f,double* dx_cf,double* dy_cf,double* dz_cf);

void inner_el(RUN * run);

void faceneigfc (struct RUN * run);
void facelimiting (struct RUN * run);
void facevolution (struct RUN * run);
void primitivevalues (struct SOL *, struct RUN * run);
void conservatvalues (struct SOL *, struct RUN * run);

void consrvtoprimtvJ (double ** J, struct SOL * sol, struct RUN * run);
void meshproperties(struct RUN *run);
void cellproperties(struct BRANCH * crnt, struct RUN * run);
void volmproperties(int flag, struct BRANCH * crnt, struct RUN * run);
void areaproperties(struct BRANCH * crnt, struct RUN * run);
void conservationcalc(RUN *run);
void checkconn(RUN *run);
void crammer(int n,double ** a,double ** b, double ** x,double **c);
double determinant(double **a);
int numlfc(int typ);
int numlnd(int typ);
int numfcnd(int ifc, int typ);
void drysconn (FOREST * topo);
void spawntopo (RUN * topo);
void partitiondrys (RUN * run);

void repartition (RUN * run);
void repartition_parmetis (RUN * run);
void rebuildlist (RUN * run);
void rebuildconn (RUN * run);
void rebuildprox (RUN * run);
void cellallocation (struct BRANCH * brch);
void celldeallocation (struct BRANCH *crnt);
void solallocation (SOL * sol, RUN *run);
struct BRANCH * spawntree (RUN * run,TREE * tree, int tag, int iel, int icon);
void erasetree (struct RUN * run, struct TREE * tree, int iel);

void destroybrch (BRANCH * brch);
int tagaddress (int * adrs,int nlvl);
int ielconn (BRANCH * brch,int ifc);
int tagconn (BRANCH * brch,int ifc);
int tagconnfin (BRANCH * brch,int ifc, int ineig);
void elconnfin (BRANCH * brch,int ifc);
void tag2adr(int * n, int tag, int * nlvl);
void getMemory(struct RUN * run, int* currRealMem, int* peakRealMem, int* currVirtMem, int* peakVirtMem);
void message(char* msg,RUN * run);

void Init_clocks (struct RUN * run, int istep);

void residual_faces_AMR(struct BRANCH * crnt, struct RUN * run);

void hllc_euler (struct BRANCH * crnt,struct RUN *run,int in,int ifc);
void flux_mc(struct BRANCH * crnt,RUN *run,int ifc,double u_star,double p_star, double cL, double cR, double rL, double rR);

void compute_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* S_S,double* SL, double* SR,
                         double* uL, double* vL, double* wL,double* s11L,double* uR, double* vR, double* wR,double* s11R);

void compute_wave_speeds_toro(struct RUN *run, struct BRANCH * crnt,int ifc,double* S_S,double* SL, double* SR,
                         double* uL, double* vL, double* wL,double* s11L,double* uR, double* vR, double* wR,double* s11R);

void star_region_euler(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
                            int ifc,int flag_ifs,double S_S,double SL,double SR,
                            double uL,double vL,double wL,double uR,double vR,double wR);
void mc_starregion(struct RUN *run, struct BRANCH * crnt,int ifc, double *u_star,double *cL,double *cR,double *rL,double *rR, double *p_star);

void flux_adv(struct SOL* sol,struct RUN *run);
void flux_viscous(struct RUN *run,struct BRANCH * crnt,struct SOL* sol,struct SOL* sol1);

void Init_1D(double *a,int Nx);
void Init_2D(double **a,int Nx,int Ny);
void Init_3D(double ***a,int Nx,int Ny,int Nz);

void asinvalues(struct RUN *run,struct BRANCH * crnt,struct SOL* sol, double rho, double ru,double rv,double rw,double e);

void asinvalues_multiphase(struct RUN *run,struct BRANCH * crnt,struct SOL* sol, double rho, double ru,double rv,double rw,double e, double ymass);
void rotate_vector_gtl_v2 (struct RUN *run, struct BRANCH * crnt,int *ifc, double *u,double *v,double *w);
void rotate_vector_ltg_v2 (struct RUN *run, struct BRANCH * crnt,int *ifc, double *u,double *v,double *w);


void screen_out(struct RUN * run,int istep);
double timecpu(double Ttemp, int flag);
void sumtime(struct RUN * run);


void cons2primtv(struct SOL* sol,struct RUN * run);
void prim2conc(struct RUN * run);

void residual_noamr(struct BRANCH * crnt, struct RUN * run);
void facercnstr_V2(struct BRANCH * crnt,int in,int ifc,int iside, int istore,RUN *run);

void exportfield_linked_V2 (struct RUN * run);
void exportfield_linked_V3 (struct RUN * run);
void exportfield_linked_V4 (struct RUN * run);

void distance(struct RUN *run,struct BRANCH *crnt,int f,double * dist_cf,double * dist_cc);
void deltavalues(struct RUN *run,struct BRANCH *crnt,double * delta_temp,int f_temp);
void reconctruct_V1(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,
                                                            double dist_m1_cf,double dist_m1_cc);

void reconctruct_V2(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,
                                                            double dist_m1_cf,double dist_m1_cc);        

void reconctruct_V3(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,
                                                              double dist_m1_cf,double dist_m1_cc);                                                                                                                

void halfdt(struct RUN *run,struct BRANCH *crnt,int f,double Area_temp);
void halfdt_v2(struct RUN *run,struct BRANCH *crnt,int f,double dx);


void export_ultra_fast (struct RUN * run);
