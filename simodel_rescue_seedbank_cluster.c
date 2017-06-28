/*                SI POPULATION DYNAMICS SIMULATION MODEL                */
/*           Source code: SIMO.C; Last modified on 17/06/2017           */

//* Compilation lines:

/// gcc -std=c99 -I/usr/include/gsl -Wall -o simo simodel_rescue.c vector.c -lgsl -lgslcblas -lncurses -lm
// gcc -std=c99 -I/usr/include/gsl -Wall -o simo_cluster simodel_rescue_cluster.c vector.c -lgsl -lgslcblas -ltinfo -lncurses -lm
// ./rescue -r 1 -g 500 -o 250 -n 3 -s 10 -f 0.25 -m 1.0 -v 0.0 -l 100 -x -z -k 0.2 -p 10 -t 50 -a 0.2 -b 0.1 -c 2 -d 2 -e demo -h gen -i fit -u indiv
// ./simo_cluster -r 100 -g 500 -o 250 -n 3 -s 10 -f 0.25 -m 1.0 -v 0.0 -l 100 -p 10 -t 50 -a 0.2 -b 0.1 -c 2 -d 2 -e demo_control.dat -h gen_control.dat -i fit_control.dat -u indiv_control.dat 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <curses.h>
#include <ncurses.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "vector.h"
#include <getopt.h>


/*                  PROGRAM CONSTANTS NOT IN MENUS                       */

#define IN 0          /* entry value for loops and if statements         */
#define OUT 1         /* exit value for loops and if statements          */
#define LEN 100       /* length of population grid size                  */
#define GENES 6       /* total number of loci (including the SI locus)   */ 
#define ALL 2         /* number of alleles carried at each locus         */
#define SDMAX 10      /* maximum number of seeds landing per cell        */
#define NR_END 1
#define FREE_ARG char*

/*                      GENERAL SIMULATION PARAMETERS                    */

#define RUNS 100      /* number of replicate runs per treatment          */
#define GEN 500       /* default value for number of generations         */
#define INUM 250      /* initial number of occupied cells                */
#define LINT 100      /* interval for data output (other than dynamics)  */
#define NLOC 3        /* number of alleles per neutral locus             */
#define SLOC 10       /* number of alleles per SI locus                  */
#define SITES 0.25    /* default fraction of safe sites                  */
#define LM 1.0        /* mean of lognormal distribution of ages          */
#define LV 0.0        /* variance of lognormal distribution of ages      */

/*               WITHIN POPULATION DYNAMICS PARAMETERS                   */

#define D 0.05        /* death rate for adult plants                     */
#define MINR  2       /* default age for first reproduction              */
#define PDIST 142      /* default radius for local pollen pool            */
#define SDIST 50       /* default radius for seed dispersal               */
#define EST 0.05      /* probability of seed establishment               */
#define B0 100.0       /* maximum per-capita ovule production             */
#define PLPR 0        /* probability of pollen contribution: random (0),
                         size-dependent (1), distance-dependent (2)      */
#define SITYPE 2      /* SI MODEL: gametophytic (0), codominant 
                         sporophytic (1), dominant sporophytic (2),
                         neutral (3)                                     */ 
#define SAFE_SITES_BOOL 0 /*  Implement spatial rescue */
#define DEMO_RESCUE_BOOL 0 /* Implement demographic rescue */
#define GEN_RESCUE_POOL 0 /* Implement genetic rescue */
#define SITES_INCREASE  0.1 /* Fraction of safe sites increase */  
#define POP_SIZE_INT 10  /* Population size to increase */
#define PNEW_S_ALLELE 0.2  /* Probability of new S alleles */ 
#define PNEW_N_ALLELE 0.0  /* Probability of new neutral alleles */
#define INT_THRESHOLD 50 /* Threshold of population size to intervene */
#define FIRE 0           /* Presence of fire events */
#define FIRE_SERO 1           /* Positive effect of fire on seed germination */
#define FIRE_PROB 0.5           /* Frequency of fire events */
#define FIRE_PROB_DEATH 0.075     /* Death rate by fire events */
#define FIRE_EST 0.075     /* Increased probability of seed establishment with fire  */              
#define SEEDBANK_MORT 0.1     /* Probability of seedbank seed mortality  */                

/*                         FUNCTION PROTOTYPES                           */

void free_fmatrix(float **,long int,long int,long int,long int);
void free_fvector(float *,long int,long int);
void free_imatrix(int **,long int,long int,long int,long int);
void free_ivector(int *,long int,long int);
void gene_data(int,int,int,int,float **,float **,float **,float **,float **,
      float **,float **,float **,float *,float *);
void getfile(void);
void grand_mean1(int,float **,float **,float **,int *,float **,float **,float **,
      int *, Vector *);
void grand_mean2(int,int *,float **,float **,float **,float **,float **,float **,
      float **);
void grand_mean3(int,int *,float **,float **,float **,float **,float **,float **,
      float **,float **);
void genetic_rescue(gsl_rng *, int, int, int, int, float, float, Vector *,Vector *);      
void header(FILE *,int,int,int,double,double,int,int,float,int,float,int,int,int,
      int);
void id_val(int,int,int,int *,int,int);     
void init_arrays(int,float **,float **,float **,int *,float **,float **,float **,
      float **,float **,float **,float **,float **,float **,float **,int *,
      float **,float **,float **,float **,float **,float **,float **,float **);
void init_plants(gsl_rng *,int,int,int,float,float, Vector *,Vector *);
void init_values(void);
void kill_plants(gsl_rng *,double);
void make_ovule(gsl_rng *,int,int,int *,int *);
void make_pollen_pool_r(int,int,int,int **,float *,int,int);
void make_pollen_pool_s(int,int,int,int **,float *,int,int);
void make_pollen_pool_d(int,int,int,int **,float *,int,float);
void make_seed(gsl_rng *,int,int,int *,int *,int,int,int);
void mate_system(int,int,void (**)(),void (**)());
void means1(int,float **,float **,float **,int,int,int);
void means2(int,int,float **,float **,float **,float **,float **,float **,
      float **,Vector *,Vector *);
void means3_g(int,int,int,int,int,float **,float **,float **,int *);
void means3_n(int,int,int,int,int,float **,float **,float **,int *);
void means3_scd(int,int,int,int,int,float **,float **,float **,int *);
void means3_sd(int,int,int,int,int,float **,float **,float **,int *);
void means3_sd2(int,int,int,int,int,float **,float **,float **,int *);
void means4(int,int,int,float **,float **,float **);
void means5(int,int,float **,float **,float **);
void menu(void);
void n_dads(int,int,int,int,int *,int *,float *);
void new_plants(gsl_rng *,float, float);
void nrerror(char []);     
void parms1(int *,int *,int *,int *,int *,float *,float *,float *,int *, int *, int *, int *);
void parms2(double *,int *,int *,int *,double *,float *,int *,int *, float *, int *, int *, float *,float *);
void persist(int,int,int);
void plant_data(int,int,int);
void reproduce_gr(int,int,int,int,int,float **,float **,double,double);
void reproduce_gs(int,int,int,int,int,float **,float **,double,double);
void reproduce_gd(int,int,int,int,int,float **,float **,double,double);
void reproduce_nr(int,int,int,int,int,float **,float **,double,double);
void reproduce_ns(int,int,int,int,int,float **,float **,double,double);
void reproduce_nd(int,int,int,int,int,float **,float **,double,double);
void reproduce_scdr(int,int,int,int,int,float **,float **,double,double);
void reproduce_scds(int,int,int,int,int,float **,float **,double,double);
void reproduce_scdd(int,int,int,int,int,float **,float **,double,double);
void reproduce_sdr(int,int,int,int,int,float **,float **,double,double);
void reproduce_sds(int,int,int,int,int,float **,float **,double,double);
void reproduce_sdd(int,int,int,int,int,float **,float **,double,double);
void reproduce_sdr2(int,int,int,int,int,float **,float **,double,double);
void reproduce_sds2(int,int,int,int,int,float **,float **,double,double);
void reproduce_sdd2(int,int,int,int,int,float **,float **,double,double);
void safe_sites(gsl_rng *r, float);
void SI_model(int,int,double,int,int,float,float **,float **,float **,int *,int,
      int,int,float **,float **,float **,float **,float **,float **,float **,
      float **,float **,float **,int *,float **,float **,float **,float **,
      float **,float **,float **,float **,int,void (*)(),void (*)(),double,int, int, int, float, float, float, float, int, Vector *, Vector *,Vector *, int, float, float, float, int, float);
      
int **imatrix(long int,long int,long int,long int);
int *ivector(long int,long int);
int log_ages(gsl_rng *,float,float);
int make_pollen(gsl_rng *r,int *,int **,float *,int,int *);
int nplants(int *,int *,int *,int);

float **fmatrix(long int,long int,long int,long int);
float *fvector(long int,long int);

//gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);

struct plant{
  int safe;
  int age;
  int sds;
  int gtype[GENES][ALL];
  int stype[SDMAX][GENES+1][ALL]; /*GENES+1 to keep ids of seed parents */
  int seedbankid[SDMAX][GENES+1][ALL];
  int seedbank;
  int seedbank_age[SDMAX];  
  int mom;
  int dad;
  int fsd;                        /* seeds fathered pre-dispersal */
  int msd;                        /* seeds mothered pre-dispersal */
  int fsp;                        /* seeds fathered post-dispersal */
  int msp;                        /* seeds mothered post-dispersal */
}pop[LEN][LEN];

FILE *fpw1,*fpw2,*fpw3,*fpw4;
     //char fout1[20],fout2[20],fout3[20],fout4[20];
char *fout11 = NULL,*fout22= NULL,*fout33= NULL,*fout44= NULL;
int rcnt;
float per1,per2;
const gsl_rng_type * T;
gsl_rng *r;

int main(int argc, char *argv[])
{ 

    int option = 0;
    unsigned long int randSeed;
    r = gsl_rng_alloc(gsl_rng_mt19937);
    srand(time(NULL));                    /* initialization for rand() */
    randSeed = rand();                    /* returns a non-negative integer */
    gsl_rng_set (r, randSeed);    /* seed the PRNG */  
    void (*repro)(),(*mean3)();
    int minr=MINR,pdist=PDIST,sdist=SDIST,sloc=SLOC,*nrn,*nrp,plpr=PLPR,lint=LINT;
    int rr,runs=RUNS,gen=GEN,inum=INUM,nloc=NLOC,sitype=SITYPE;  //
    float est=EST,sites=SITES,**ec1,**ec2,**ec3;
    float **mnv,**mnr,**myr,**gn1,**gn2,**gn3,**gn4,**gn5,**gn6,**gn7;
    float **ft1,**ft2,**ft3,**ft4,**ft5,**ft6,**ft7,**ft8,lm=LM,lv=LV;
    double d=D,b0=B0;
    int safe_sites_bool = SAFE_SITES_BOOL, demo_rescue_bool = DEMO_RESCUE_BOOL, gen_rescue_pool = GEN_RESCUE_POOL, pop_size_interv = POP_SIZE_INT, interv_threshold = INT_THRESHOLD, fire_flag = FIRE, fire_postive_seeds = FIRE_SERO;
    float sites_increase = SITES_INCREASE, prob_new_S_allele = PNEW_S_ALLELE, prob_new_allele = PNEW_N_ALLELE, fire_prob = FIRE_PROB, fire_prob_death = FIRE_PROB_DEATH, est_fire = FIRE_EST, seedbank_mort = SEEDBANK_MORT;
    // declare and initialize a new vector
    Vector vector_S_alleles;
    Vector vector_N_alleles;
    Vector vector_quasi_ext;
    
    

    while((option = getopt(argc, argv,"r:g:o:n:s:f:m:v:l:xyzk:p:t:a:b:c:d:j:F:A:B:C:D:e:h:i:u:")) != -1)
      {
		 switch (option) {
             case 'r' : runs = atoi(optarg);
                 break;
             case 'g' : gen = atoi(optarg);
                 break;
             case 'o' : inum = atoi(optarg); 
                 break;
             case 'n' : nloc = atoi(optarg);
                 break;
             case 's' : sloc = atoi(optarg);
                 break;
             case 'f' : sites = atof(optarg);
                 break;
             case 'm' : lm = atof(optarg); 
                 break;
             case 'v' : lv = atof(optarg);
                 break;
             case 'l' : lint = atoi(optarg);
                 break;
             case 'x' : safe_sites_bool = 1;
                 break;
             case 'y' : demo_rescue_bool = 1; 
                 break;
             case 'z' : gen_rescue_pool = 1;
                 break;
             case 'k' : sites_increase = atof(optarg);
                 break;
             case 'p' : pop_size_interv = atoi(optarg);
                 break;
             case 't' : interv_threshold = atoi(optarg); 
                 break;
             case 'a' : prob_new_S_allele = atof(optarg);             
                 break;
             case 'b' : prob_new_allele = atof(optarg); 
                 break;
             case 'c' : plpr = atoi(optarg);
                 break;
             case 'd' : sitype = atoi(optarg);
                 break;
             case 'j': fire_flag = atoi(optarg);
                 break;
             case 'F': fire_prob = atof(optarg);
                 break;
             case 'A': fire_prob_death = atof(optarg);
                 break;
             case 'B': fire_postive_seeds = atoi(optarg);
                 break;
             case 'C': est_fire = atof(optarg);
                 break;
             case 'D' : seedbank_mort = atof(optarg);
                 break;                                                                                                         
             case 'e' : fout11 = optarg;
                 break;
             case 'h' : fout22 = optarg;
                 break;
             case 'i' : fout33 = optarg;
                 break;
             case 'u' : fout44 = optarg;
                 break;                                                  
             case '?':
                 fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                 //break;                      
             //default: print_usage(); 
                 exit(EXIT_FAILURE);
             default:
                 exit(EXIT_FAILURE);
        }  
	  }
      if(fout11 == NULL){
		  printf("Initialize demography file name with -e option\n");
		  exit(2);
	  }
	  if(! (fpw1 = fopen(fout11,"w"))){
		  perror(fout11);
		  exit(2);
	  }
   // while(sflg==IN) //IN
     // {
        //system("clear");//system("cls");
	   // char choice;
       // menu();
       // scanf("%c", &choice);
        //system("clear");//system("cls");
        //refresh();
        //switch(choice)
          //{
            //case 'i':  
              //parms1(&runs,&gen,&inum,&nloc,&sloc,&sites,&lm,&lv,&lint,&safe_sites_bool,&demo_rescue_bool,&gen_rescue_pool);
              //break;
            //case 'w':  
             // parms2(&d,&minr,&pdist,&sdist,&b0,&est,&plpr,&sitype,&sites_increase,&pop_size_interv,&interv_threshold,&prob_new_S_allele,&prob_new_allele); 
             // break;
            //case 'r':
              getfile();
              rcnt = 0; 
              per1 = per2 = 0.0;
              mnv = fmatrix(0,gen,0,1); /* avg number of vegetative plants */
              mnr = fmatrix(0,gen,0,1); /* avg number of reproductive plants */
              myr = fmatrix(0,gen,0,1); /* avg plant age */
              ec1 = fmatrix(0,gen,0,1); /* avg number of empty sites/plants */
              ec2 = fmatrix(0,gen,0,1); /* avg number of potential mates/plant */
              ec3 = fmatrix(0,gen,0,1); /* avg number of compatible mates/plant */
              gn1 = fmatrix(0,gen,0,1); /* avg number of alleles/neutral locus */
              gn2 = fmatrix(0,gen,0,1); /* avg number of alleles/SI locus */
              gn3 = fmatrix(0,gen,0,1); /* avg variance in neutral allele frequencies */
              gn4 = fmatrix(0,gen,0,1); /* avg variance in SI allele frequencies */
              gn5 = fmatrix(0,gen,0,1); /* average observed heterozygosity */
              gn6 = fmatrix(0,gen,0,1); /* average expected heterozygosity */
              gn7 = fmatrix(0,gen,0,1); /* avg value of Fis */
              ft1 = fmatrix(0,gen,0,1); /* avg number of pollen donors/plant */
              ft2 = fmatrix(0,gen,0,1); /* avg variance in pollen donors/plant */
              ft3 = fmatrix(0,gen,0,1); /* avg seeds/plant pre-dispersal */
              ft4 = fmatrix(0,gen,0,1); /* avg variance: pre-dispersal seeds/father */
              ft5 = fmatrix(0,gen,0,1); /* avg variance: pre-dispersal seeds/mother */
              ft6 = fmatrix(0,gen,0,1); /* avg seeds/plant post-dispersal */
              ft7 = fmatrix(0,gen,0,1); /* avg variance: post-dispersal seeds/father */
              ft8 = fmatrix(0,gen,0,1); /* avg variance: post-dispersal seeds/mother */
              nrn = ivector(0,gen);
              nrp = ivector(0,gen);
              vector_init(&vector_quasi_ext);                  /* Quasi-extinction probability vector init */
              vector_set(&vector_quasi_ext, gen, 0); /* Quasi-extinction probability vector fill with zeroes till maximum generations per run */
              header(fpw1,runs,gen,minr,b0,d,pdist,sdist,est,inum,sites,nloc,sloc,
                     plpr,sitype);
              header(fpw2,runs,gen,minr,b0,d,pdist,sdist,est,inum,sites,nloc,sloc,
                     plpr,sitype);
              header(fpw3,runs,gen,minr,b0,d,pdist,sdist,est,inum,sites,nloc,sloc,
                     plpr,sitype);
              header(fpw4,runs,gen,minr,b0,d,pdist,sdist,est,inum,sites,nloc,sloc,
                     plpr,sitype);
              fprintf(fpw4," GEN XC YC AGE SID SI1 SI2 ID1 N11 N12 ID2 N21 N22");
              fprintf(fpw4," ID3 N31 N32 ID4 N41 N42 ID5 N51 N52\n");
              init_arrays(gen,mnv,mnr,myr,nrn,gn1,gn2,gn3,gn4,gn5,gn6,gn7,ec1,ec2,ec3,
                          nrp,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8);
              mate_system(plpr,sitype,&repro,&mean3);                                                                               
              for(rr=0;rr<runs;rr++)
                { 
				  vector_init(&vector_S_alleles);
				  vector_init(&vector_N_alleles);	
                  init_values();    
                  safe_sites(r,sites);
                  init_plants(r,inum,nloc,sloc,lm,lv,&vector_S_alleles,&vector_N_alleles);   
                  printf("Current Run: %d\n",rr+1);
                  SI_model(gen,minr,d,pdist,sdist,est,mnv,mnr,myr,nrn,sloc,nloc,lint,gn1,
                           gn2,gn3,gn4,gn5,gn6,gn7,ec1,ec2,ec3,nrp,ft1,ft2,ft3,ft4,ft5,
                           ft6,ft7,ft8,rr,repro,mean3,b0,safe_sites_bool,demo_rescue_bool,gen_rescue_pool,sites_increase,pop_size_interv,prob_new_S_allele,prob_new_allele,interv_threshold,&vector_S_alleles,&vector_N_alleles, &vector_quasi_ext, fire_flag, fire_prob, fire_prob_death, est_fire, fire_postive_seeds, seedbank_mort);
                  /// Add free_vector possibly to free memory of vector_S_alleles
                  vector_free(&vector_S_alleles);
                  vector_free(&vector_N_alleles);                                               
                }
              persist(runs,rcnt,gen);    
              grand_mean1(gen,mnv,mnr,myr,nrn,ec1,ec2,ec3,nrp,&vector_quasi_ext);
              grand_mean2(gen,nrn,gn1,gn2,gn3,gn4,gn5,gn6,gn7);
              grand_mean3(gen,nrp,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8);
              free_ivector(nrn,0,gen);
              free_ivector(nrp,0,gen);
              printf("FLAG VECTOR");
              vector_free(&vector_quasi_ext);   
              free_fmatrix(mnv,0,gen,0,1);
              free_fmatrix(mnr,0,gen,0,1);
              free_fmatrix(myr,0,gen,0,1);
              free_fmatrix(ec1,0,gen,0,1);
              free_fmatrix(ec2,0,gen,0,1);
              free_fmatrix(ec3,0,gen,0,1);      
              free_fmatrix(gn1,0,gen,0,1);
              free_fmatrix(gn2,0,gen,0,1);
              free_fmatrix(gn3,0,gen,0,1);
              free_fmatrix(gn4,0,gen,0,1);
              free_fmatrix(gn5,0,gen,0,1);
              free_fmatrix(gn6,0,gen,0,1);
              free_fmatrix(gn7,0,gen,0,1);
              free_fmatrix(ft1,0,gen,0,1);
              free_fmatrix(ft2,0,gen,0,1);
              free_fmatrix(ft3,0,gen,0,1);
              free_fmatrix(ft4,0,gen,0,1);
              free_fmatrix(ft5,0,gen,0,1);
              free_fmatrix(ft6,0,gen,0,1);
              free_fmatrix(ft7,0,gen,0,1);
              free_fmatrix(ft8,0,gen,0,1);
              fclose(fpw1);
              fclose(fpw2);
              fclose(fpw3);
              fclose(fpw4);
            //  break;
            //case 'e':
            //  sflg = OUT;
            //  break;
            //default: break;
          //}
      //}
      gsl_rng_free (r);
      return 0;
}

/* prints main menu to screen */

//void menu(void)
//  {
//    printf("MAIN MENU\n\n");
//    printf("Enter One Of The Following Choices:\n\n");
//    printf("I = initialisation and general simulation parameters\n");
//    printf("W = set parameters for within-population dynamics\n");
//    printf("R = run the simulation model\n");
//    printf("E = exit program (Return to DOS)\n");
//    return;
//  }

/* assigns pointer to appropriate functions for reproduction under a given
   mating system */

void mate_system(int plpr,int sitype,void (**repro)(),void (**mean3)())
  {
    if(sitype==0)
      {
        *mean3 = &means3_g;
        if(plpr==0)
          *repro = &reproduce_gr;
        else if(plpr==1)
          *repro = &reproduce_gs;
        else if(plpr==2)
          *repro = &reproduce_gd;
      }
    else if(sitype==1)
      {
        *mean3 = &means3_scd;
        if(plpr==0)
          *repro = &reproduce_scdr;
        else if(plpr==1)
          *repro = &reproduce_scds;
        else if(plpr==2)
          *repro = &reproduce_scdd;
      }
    else if(sitype==2)
      {
        *mean3 = &means3_sd;
        if(plpr==0)
          *repro = &reproduce_sdr;
        else if(plpr==1)
          *repro = &reproduce_sds;
        else if(plpr==2)
          *repro = &reproduce_sdd;
      }
    else if(sitype==3)
      {
        *mean3 = &means3_n;
        if(plpr==0)
          *repro = &reproduce_nr;
        else if(plpr==1)
          *repro = &reproduce_ns;
        else if(plpr==2)
          *repro = &reproduce_nd;
      }
    else if(sitype==4)
      {
        *mean3 = &means3_sd2;
        if(plpr==0)
          *repro = &reproduce_sdr2;
        else if(plpr==1)
          *repro = &reproduce_sds2;
        else if(plpr==2)
          *repro = &reproduce_sdd2;
      }  
    return;  
  }          


/* prompts user for output file names; opens files fpw1, fpw2, fpw3 */

void getfile(void)
  {
    printf("Enter output filename for demographic data\n");
    //scanf("%s",fout1);
    fpw1 = fopen(fout11,"w");
    printf("Enter output filename for genetic data\n");
    //scanf("%s",fout2);
    fpw2 = fopen(fout22,"w");
    printf("Enter output filename for fitness data\n");
    //scanf("%s",fout3);
    fpw3 = fopen(fout33,"w");
    printf("Enter output filename for individual data\n");
    //scanf("%s",fout4);
    fpw4 = fopen(fout44,"w");
    //system("clear");//system("cls");
    return;
  }

/* prints header containing input variable values to output files */

void header(FILE *fpw,int runs,int gen,int minr,double b0,double d,int pdist,
      int sdist,float est,int inum,float sites,int nloc,int sloc,int plpr,int sitype)
  {
    if(sitype==0)
      fprintf(fpw,"Mating System: Gametophytic Self Incompatibility, ");
    else if(sitype==1)
      fprintf(fpw,"Mating System: Codominant Sporophytic Self Incompatibility, ");
    else if(sitype==2)
      fprintf(fpw,"Mating System: Maternal Dominant Sporophytic Self Incompatibility, ");
    else if(sitype==3)
      fprintf(fpw,"Mating System: Neutral, ");
    else if(sitype==4)
      fprintf(fpw,"Mating System: Paternal Dominant Sporophytic Self Incompatibility, "); 
    if(plpr==0)
      fprintf(fpw,"with Random Pollination\n\n");
    else if(plpr==1)
      fprintf(fpw,"with Size-Dependent Pollination\n\n");
    else if(plpr==2)
      fprintf(fpw,"with Distance-Dependent Pollination\n\n");
    fprintf(fpw,"Runs=%d, Generations=%d, Initial PopSize=%d, Safe Sites=%.3f, First "
            "Reproduction=%d\n",runs,gen,inum,sites,minr);        
    fprintf(fpw,"Maximum Per-Capita Ovule Production=%.3lf, Seed Establishment=%.3f, "
            "Deathrate=%.3lf\n",b0,est,d);  
    fprintf(fpw,"Pollen Dispersal Range=%d, Seed Dispersal Range=%d, SI Alleles=%d, "
            "Neutral Alleles=%d\n\n",pdist,sdist,sloc,nloc);
    return;
  }

/* sets elements of output data arrays to zero */

void init_arrays(int gen,float **mnv,float **mnr,float **myr,int *nrn,float **gn1,
      float **gn2,float **gn3,float **gn4,float **gn5,float **gn6,float **gn7,
      float **ec1,float **ec2,float **ec3,int *nrp,float **ft1,float **ft2,
      float **ft3,float **ft4,float **ft5,float **ft6,float **ft7,float **ft8)
  {
    int i,j;

    for(i=0;i<=gen;i++)
      {
        nrn[i] = 0;
        nrp[i] = 0;
        for(j=0;j<ALL;j++)
          {
            mnv[i][j] = 0.0;
            mnr[i][j] = 0.0;
            myr[i][j] = 0.0;
            ec1[i][j] = 0.0;
            ec2[i][j] = 0.0;
            ec3[i][j] = 0.0;
            gn1[i][j] = 0.0;
            gn2[i][j] = 0.0;
            gn3[i][j] = 0.0;
            gn4[i][j] = 0.0;
            gn5[i][j] = 0.0;
            gn6[i][j] = 0.0;
            gn7[i][j] = 0.0;
            ft1[i][j] = 0.0;
            ft2[i][j] = 0.0;
            ft3[i][j] = 0.0;
            ft4[i][j] = 0.0;
            ft5[i][j] = 0.0;
            ft6[i][j] = 0.0;
            ft7[i][j] = 0.0;
            ft8[i][j] = 0.0;
          }
      }
    return;
  }

/* initializes elements of the population arrays for use in simod() */

void init_values(void)
  {
    int xc,yc,i,j,k;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            pop[xc][yc].safe = 0;
            pop[xc][yc].age = 0;
            pop[xc][yc].sds = 0;         
            pop[xc][yc].mom = 0;
            pop[xc][yc].dad = 0;
            pop[xc][yc].fsd = 0;
            pop[xc][yc].msd = 0;
            pop[xc][yc].fsp = 0;
            pop[xc][yc].msp = 0;
            for(i=0;i<GENES;i++)
              {
                for(j=0;j<ALL;j++)
                  pop[xc][yc].gtype[i][j] = 0;
              }
            for(i=0;i<SDMAX;i++)
              {
                for(j=0;j<=GENES;j++)
                  {
                    for(k=0;k<ALL;k++)
                      pop[xc][yc].stype[i][j][k] = 0;
                  }
              }  
          }
      }
    return;
  }

/* randomly determines a fraction of 'safe' sites (i.e. suitable for occupancy) */

void safe_sites(gsl_rng *r, float sites)
  {
    int xc,yc,cnt=0;

    sites = sites*LEN*LEN;
    if(sites==(LEN*LEN))
      {
        for(xc=0;xc<LEN;xc++)
          for(yc=0;yc<LEN;yc++)
            pop[xc][yc].safe = 1;
      }
    else if(sites<(LEN*LEN))
      {  
        while(cnt<=sites)
          {
            xc = gsl_rng_uniform_int(r,LEN-1); //g05dyc(0,LEN-1);
            yc = gsl_rng_uniform_int(r,LEN-1); //g05dyc(0,LEN-1);
            if(pop[xc][yc].safe==0)
              {
                pop[xc][yc].safe = 1;
                cnt += 1;
              }
          }
      }
    //gsl_rng_free (r);  
    return;
  }

/* Check if there is a new S allele (and neutral allele)*/

bool check_new_Sallele(int S_allele, Vector *vector){
	bool var;
	int max = vector->size; 
	for(int index=0;index < max; ++index){
	    if(S_allele == vector_get(vector, index))
	       var = true;
	    else
	       var = false; 
    }
	return var;
}


/* calculates genotypes and randomly assigns initial positions for plants; assumes
   that individuals can be heterozygous at each neutral locus */
   
   // CREATE VECTOR OF S ALLELES TO UPDATE THE LIST OF S ALLELES CREATED BY NEW MUTATIONS AT THE S LOCUS DURING THE GENETIC RESCUE

void init_plants(gsl_rng *r, int inum,int nloc,int sloc,float lm,float lv, Vector *vector_S_alleles, Vector *vector_N_alleles)
  {
	bool check;
    int i,j,xc,yc,al,cnt=0;
    while(cnt<inum)
      {
        xc = gsl_rng_uniform_int(r,LEN-1);//g05dyc(0,LEN-1);
        yc = gsl_rng_uniform_int(r,LEN-1);//g05dyc(0,LEN-1);
        if(pop[xc][yc].age==0 && pop[xc][yc].safe==1)
          {
            for(i=0;i<GENES;i++)
              {
                if(i==0)           /* SI locus - must be heterozygous */
                  {
                    al = gsl_rng_uniform_int(r,sloc);//g05dyc(1,sloc);
                    pop[xc][yc].gtype[i][0] = al;
                    if(cnt == 0){
                       vector_append(vector_S_alleles,al);
                    }else{
					   check = check_new_Sallele(al, vector_S_alleles);
                       if(check == false)
                          vector_append(vector_S_alleles,al);
					}
                       
                    while(al==pop[xc][yc].gtype[i][0])
                       al = gsl_rng_uniform_int(r,sloc);//g05dyc(1,sloc);
                    pop[xc][yc].gtype[i][1] = al;
                    check = check_new_Sallele(al, vector_S_alleles);
                    if(check == false)
                       vector_append(vector_S_alleles,al);   
                  }
                else if(i>0)       /* Neutral loci */
                  {  
                    for(j=0;j<ALL;j++)
                      {
                        al = gsl_rng_uniform_int(r,nloc);//g05dyc(1,nloc);
                        pop[xc][yc].gtype[i][j] = al;
                        check = check_new_Sallele(al, vector_N_alleles);
                        if(check == false)
                           vector_append(vector_N_alleles,al);
                      }
                  }
              }  
            pop[xc][yc].age = log_ages(r,lm,lv);  
            cnt += 1;
          }
      }
    //gsl_rng_free (r);      
    return;
  }

/* Demographic rescue --- Put N individuals in the lattice doing random mating to produce incalculates genotypes and randomly assigns initial positions for plants; assumes
   that individuals can be heterozygous at each neutral locus */

void demographic_rescue(gsl_rng *r, int inum, int minr)
  {
    int i,xc,xc_dad,yc,yc_dad,xc_mom,yc_mom,al,al2,cnt=0,count;
    int coords_mom[2] = {}, coords_dad[2] = {};
    while(cnt<inum)
      {
        xc = gsl_rng_uniform_int(r,LEN-1);//g05dyc(0,LEN-1);
        yc = gsl_rng_uniform_int(r,LEN-1);//g05dyc(0,LEN-1);
        if(pop[xc][yc].age==0 && pop[xc][yc].safe==1)
          {
			//printf(" Coord x: %d\n",xc);
			//printf(" Coord y: %d\n",yc);   
			// Find "parents":
			// Dad
			count = 0;
			while(count<1){
                xc_dad = gsl_rng_uniform_int(r,LEN-1);
                yc_dad = gsl_rng_uniform_int(r,LEN-1);
                if(pop[xc_dad][yc_dad].age>=minr){
					   coords_dad[0] = xc_dad;
					   coords_dad[1] = yc_dad;
					   count=1;
		        }
			}
			// Mom
			count = 0;
			while(count<1){
                xc_mom = gsl_rng_uniform_int(r,LEN-1);
                yc_mom = gsl_rng_uniform_int(r,LEN-1);
                if(pop[xc_mom][yc_mom].age>=minr){
					   coords_mom[0] = xc_mom;
					   coords_mom[1] = yc_mom;
					   count=1;
		        }
			} 			  
			// Assigning the genes to introduced reproductive individual from chosen individual(s) (parents) 
            for(i=0;i<GENES;i++)
              {
                if(i==0)           /* SI locus - must be heterozygous */
                  {   
                    al = pop[coords_dad[0]][coords_dad[1]].gtype[i][0]; // Assign S allele 1
                    pop[xc][yc].gtype[i][0] = al;
                    if(al != pop[coords_mom[0]][coords_mom[1]].gtype[i][1]){
                      al2 = pop[coords_mom[0]][coords_mom[1]].gtype[i][1];
                    }else{
						al2 = pop[coords_mom[0]][coords_mom[1]].gtype[i][0];
					}   
                    if(al2 == al)
                       al2 = pop[coords_dad[0]][coords_dad[1]].gtype[i][1];
                    pop[xc][yc].gtype[i][1] = al2; // Assign S allele 2
                    printf(" S allele 1: %d\n",pop[xc][yc].gtype[i][0]);
                    printf(" S allele 2: %d\n",pop[xc][yc].gtype[i][1]);                    
                    //if(pop[xc][yc].gtype[i][0] > 100 || pop[xc][yc].gtype[i][1] > 100){
                       //printf(" S allele 1: %d\n",pop[xc][yc].gtype[i][0]);
                       //printf(" S allele 2: %d\n",pop[xc][yc].gtype[i][1]);
                       //printf(" S allele 1 dad: %d\n",pop[coords_dad[0]][coords_dad[1]].gtype[i][0]);
			           //printf(" S allele 2 dad: %d\n",pop[coords_dad[0]][coords_dad[1]].gtype[i][1]); 
			           //printf(" S allele 1 mom: %d\n",pop[coords_mom[0]][coords_mom[1]].gtype[i][0]);
			           //printf(" S allele 2 mom: %d\n",pop[coords_mom[0]][coords_mom[1]].gtype[i][1]);                       
				   //}
                  }
                else if(i>0)       /* Neutral loci */
                  {  
                      //Assign neutral allele, nloc: number of neutral alleles
                      pop[xc][yc].gtype[i][0] = pop[coords_dad[0]][coords_dad[1]].gtype[i][0];
                      pop[xc][yc].gtype[i][1] = pop[coords_mom[0]][coords_mom[1]].gtype[i][1];
                  }
              }
            coords_dad[0] = 0;
            coords_dad[1] = 0;
            coords_mom[0] = 0;
            coords_mom[1] = 0;    
            pop[xc][yc].age = 1;//log_ages(r,lm,lv);  
            cnt += 1;
          }
      }
        
    return;
  }

/* Genetic rescue --- Put N individuals with 20% different genetic diversity in the lattice doing random mating to produce incalculates genotypes and randomly assigns initial positions for plants; assumes
   that individuals can be heterozygous at each neutral locus */

void genetic_rescue(gsl_rng *r, int inum, int minr, int sloc, int nloc, float prob_new_S_allele, float prob_new_allele, Vector *Stype, Vector *vector_N_alleles)
  {
    int i,xc,xc_dad,yc,yc_dad,xc_mom,yc_mom,al,al2,cnt=0,count,update_sloc=sloc;
    int coords_mom[2] = {}, coords_dad[2] = {};
    bool check;
    while(cnt<inum)
      {
        xc = gsl_rng_uniform_int(r,LEN-1);//g05dyc(0,LEN-1);
        yc = gsl_rng_uniform_int(r,LEN-1);//g05dyc(0,LEN-1);
        if(pop[xc][yc].age==0 && pop[xc][yc].safe==1)
          {
			//printf(" Coord x: %d\n",xc);
			//printf(" Coord y: %d\n",yc);   
			// Find "parents":
			// Dad
			count = 0;
			while(count<1){
                xc_dad = gsl_rng_uniform_int(r,LEN-1);
                yc_dad = gsl_rng_uniform_int(r,LEN-1);
                if(pop[xc_dad][yc_dad].age>=minr){
					   coords_dad[0] = xc_dad;
					   coords_dad[1] = yc_dad;
					   count=1;
		        }
			}
			// Mom
			count = 0;
			while(count<1){
                xc_mom = gsl_rng_uniform_int(r,LEN-1);
                yc_mom = gsl_rng_uniform_int(r,LEN-1);
                if(pop[xc_mom][yc_mom].age>=minr){
					   coords_mom[0] = xc_mom;
					   coords_mom[1] = yc_mom;
					   count=1;
		        }
			} 			  
			// Assigning the genes to introduced reproductive individual from chosen individual(s) (parents) 
            for(i=0;i<GENES;i++)
              {
                if(i==0)           /* SI locus - must be heterozygous */
                  { 
					if(prob_new_S_allele >  gsl_rng_uniform(r)){    
                      al = gsl_rng_uniform_int(r,sloc) + 10; // Assign S allele 1
                      update_sloc ++;
				    }else{
					  al = pop[coords_dad[0]][coords_dad[1]].gtype[i][0];
				    } 
                    pop[xc][yc].gtype[i][0] = al;
                    if(prob_new_S_allele >  gsl_rng_uniform(r)){
						al2 = gsl_rng_uniform_int(r,sloc) + 10; // Assign S allele 2
					}else{
						al2 = pop[coords_mom[0]][coords_mom[1]].gtype[i][1];
					}
					
                    if(al == al2){
						al2 = pop[coords_mom[0]][coords_mom[1]].gtype[i][0];
					}else{
						update_sloc ++;
					}   
                    if(al == al2){
						if(prob_new_S_allele >  gsl_rng_uniform(r)){
						    al2 = gsl_rng_uniform_int(r,sloc) + 10; //////// --------------------> rethink this, may be gsl_rng_uniform_int(r,HUGE_NUMBER) or different distribution
                        }else{ 
                            al2 = pop[coords_dad[0]][coords_dad[1]].gtype[i][1];
						}
					}
                    pop[xc][yc].gtype[i][1] = al2; // Assign S allele 2
                    //printf(" S allele 1: %d\n",pop[xc][yc].gtype[i][0]);
                    //printf(" S allele 2: %d\n",pop[xc][yc].gtype[i][1]);
                    check = check_new_Sallele(pop[xc][yc].gtype[i][0], Stype);
                    if(check == false)
                       vector_append(Stype,pop[xc][yc].gtype[i][0]);
                    check = check_new_Sallele(pop[xc][yc].gtype[i][1], Stype);
                    if(check == false)
                       vector_append(Stype,pop[xc][yc].gtype[i][1]);                                          
                    //if(pop[xc][yc].gtype[i][0] > 100 || pop[xc][yc].gtype[i][1] > 100){
                       //printf(" S allele 1: %d\n",pop[xc][yc].gtype[i][0]);
                       //printf(" S allele 2: %d\n",pop[xc][yc].gtype[i][1]);
                       //printf(" S allele 1 dad: %d\n",pop[coords_dad[0]][coords_dad[1]].gtype[i][0]);
			           //printf(" S allele 2 dad: %d\n",pop[coords_dad[0]][coords_dad[1]].gtype[i][1]); 
			           //printf(" S allele 1 mom: %d\n",pop[coords_mom[0]][coords_mom[1]].gtype[i][0]);
			           //printf(" S allele 2 mom: %d\n",pop[coords_mom[0]][coords_mom[1]].gtype[i][1]);                       
				   //}
                  }
                else if(i>0)       /* Neutral loci */
                  {  
                      //Assign neutral allele, nloc: number of neutral alleles
                      if(prob_new_allele >  gsl_rng_uniform(r)){
						  pop[xc][yc].gtype[i][0] = gsl_rng_uniform_int(r,nloc) + 5;
					  }else{
                          pop[xc][yc].gtype[i][0] = pop[coords_dad[0]][coords_dad[1]].gtype[i][0];
					  }
					  if(prob_new_allele >  gsl_rng_uniform(r)){
						  pop[xc][yc].gtype[i][1] = gsl_rng_uniform_int(r,nloc) + 5;
					  }else{
                          pop[xc][yc].gtype[i][1] = pop[coords_mom[0]][coords_mom[1]].gtype[i][1];
					  }
					                      check = check_new_Sallele(pop[xc][yc].gtype[i][0], Stype);
					  // Check if there are new neutral alleles, if yes add to the vector of neutral alleles                    
                      if(check == false)
                         vector_append(vector_N_alleles,pop[xc][yc].gtype[i][0]);
                      check = check_new_Sallele(pop[xc][yc].gtype[i][1], vector_N_alleles);
                      if(check == false)
                         vector_append(vector_N_alleles,pop[xc][yc].gtype[i][1]); 
                  }
              }
            coords_dad[0] = 0;
            coords_dad[1] = 0;
            coords_mom[0] = 0;
            coords_mom[1] = 0;    
            pop[xc][yc].age = 1;//log_ages(r,lm,lv);  
            cnt += 1;
          }
      }
        
    return;
  }


/* assigns ages from a lognormal distribution, with mean (lm) and variance (lv) */

int log_ages(gsl_rng *r,float lm, float lv)
  {
    int flg=OUT;
    float nm,nv,sd,tmp;

    nm = log(lm) - 0.5*log(lv/(lm*lm) + 1.0);
    nv = log(lv/(lm*lm) + 1.0);
    sd = sqrt(nv);
    while(flg==OUT)
      {
        tmp= exp((gsl_ran_gaussian(r,sd) + nm));//exp(g05ddc(nm,sd));
        ((tmp-floor(tmp))<0.5)?(tmp=floor(tmp)):(tmp=ceil(tmp));
        if(tmp>=1)
          flg = IN;
      }
    //gsl_rng_free (r);    
    return (int)tmp;
  }

/* self-incompatibility population dynamical model */

void SI_model(int gen,int minr,double d,int pdist,int sdist,float est,float **mnv,
      float **mnr,float **myr,int *nrn,int sloc,int nloc,int lint,float **gn1,
      float **gn2,float **gn3,float **gn4,float **gn5,float **gn6,float **gn7,
      float **ec1,float **ec2,float **ec3,int *nrp,float **ft1,float **ft2,
      float **ft3,float **ft4,float **ft5,float **ft6,float **ft7,float **ft8,int rr,
      void (*repro)(),void (*mean3)(),double b0, int safe_sites_bool, int demo_rescue_bool, int gen_rescue_pool,
      float sites_increase, float pop_size_interv, float prob_new_S_allele, float prob_new_allele, int interv_threshold, 
      Vector *Stype, Vector *vector_N_alleles, Vector *vector_quasi_ext,
      int fire_flag, float fire_prob, float fire_prob_death, float est_fire, int fire_postive_seeds, float seedbank_mort)
  {
    int g,i,plants,rep,veg,age,cnt=0, number_interv=0, rcnt_fake, update_quasi_prob = 0;
    //int pop_size_interv = 100, new_sloc = sloc + 10, new_nloc = nloc + 5;
    //float sites_increase = 0.2;
    //float prob_new_S_allele = 0.5, prob_new_allele = 0.5;
    for(g=0;g<=gen;g++)
      {
		//printf(" FLAG 4 \n");  
        if((plants=nplants(&veg,&rep,&age,minr))>1)
          {
            means1(g,mnv,mnr,myr,veg,rep,age);
            printf(" FLAG 5 \n"); 
            means2(g,plants,gn1,gn2,gn3,gn4,gn5,gn6,gn7,Stype,vector_N_alleles);//means2(g,plants,new_sloc,new_nloc,gn1,gn2,gn3,gn4,gn5,gn6,gn7);
            printf(" FLAG 6 \n"); 
            (*mean3)(g,rep,minr,pdist,sdist,ec1,ec2,ec3,nrp);
            printf(" FLAG 7 \n"); 
            (*repro)(g,rep,minr,pdist,sdist,ft1,ft2,b0,d);
            printf(" FLAG 8 \n");    
            means4(g,rep,minr,ft3,ft4,ft5);
            if(cnt==0 && rr==0)
              plant_data(g,sloc,nloc);
            (cnt==lint-1)?(cnt=0):(cnt+=1);
            printf(" Population size: %d\n",plants);
            printf(" --------------------------------------------------------------------Generation: %d\n",g);
            if(plants < interv_threshold && number_interv == 0){
              if(safe_sites_bool == 1)
                   safe_sites(r, sites_increase);
              if(demo_rescue_bool == 1)     
                   demographic_rescue(r, pop_size_interv, minr);
              if(gen_rescue_pool == 1)
                   genetic_rescue(r, pop_size_interv, minr, sloc, nloc, prob_new_S_allele, prob_new_allele, Stype, vector_N_alleles);
              plants=nplants(&veg,&rep,&age,minr);
              printf(" Population size after rescue: %d\n",plants);
              printf(" Genetic rescue: %d\n",gen_rescue_pool);
              number_interv=1;
		    }
		    // Filling vector to calculate quasi-extinction probability risk curve 
		    if(plants < interv_threshold){
				update_quasi_prob = vector_get(vector_quasi_ext,g) + 1;
				vector_set(vector_quasi_ext,g, update_quasi_prob);
			}
		    //printf(" Genetic rescue: %d\n",gen_rescue_pool);
		    //printf(" Demographic rescue: %d\n",demo_rescue_bool);
            //printf(" Spatial rescue: %d\n",safe_sites_bool);
            /// FIRE
            if(fire_flag == 1){ /// Flags if there will be fire events in the simulations or not at all
				if(fire_prob >  gsl_rng_uniform(r)){ /// prob_fire_freq: frequency of fire events along one simulation 
					kill_plants(r,fire_prob_death);   /// fire_death_prob: increased (compared to "d") probability of death by a fire event
					printf("             ----->   FIRE!\n");
					if(fire_postive_seeds == 1)   new_plants(r,est_fire, seedbank_mort); /// Add positive effect of fire in seed germination probability when species is serotinous
					else  new_plants(r,est,seedbank_mort); /// No serotinous
				}else{
					kill_plants(r,d);
				}
			}else{ /// NO FIRE
				kill_plants(r,d);
				new_plants(r,est,seedbank_mort);
			}
            
            printf(" FLAG 2 \n"); 
            means5(g,rep,ft6,ft7,ft8);
            printf(" FLAG 3 \n");
          } 
        else break;
      }
    for(i=0;i<g;i++)
      nrn[i] += 1;
    rcnt_fake = rcnt;                /* number of runs with plants present */
    (g<gen+1)?(rcnt=rcnt_fake):(rcnt+=1);
    per1 += (g-1);
    per2 += (g-1)*(g-1);
    return;
  }

/* calculates number of vegetatives, number of reproductives and total ages; returns
   the total population size */

int nplants(int *veg,int *rep,int *age,int minr)
  {
    int xc,yc;

    *veg = 0;
    *rep = 0;
    *age = 0;
    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>0)
              {
                *age += pop[xc][yc].age;  
                if(pop[xc][yc].age<minr)
                  *veg += 1;
                else if(pop[xc][yc].age>=minr)
                  *rep += 1;
              }
          }
      }
    return (*veg + *rep);
  }  

/* sums and squared sums of demographic variables for average dynamics */

void means1(int g,float **mnv,float **mnr,float **myr,int veg,int rep,int age)
  {
    mnv[g][0] += (float)veg;
    mnv[g][1] += (float)veg*veg;
    mnr[g][0] += (float)rep;
    mnr[g][1] += (float)rep*rep;
    myr[g][0] += (float)age/(veg+rep);
    myr[g][1] += (float)(age/(veg+rep))*(age/(veg+rep));    
    return;
  }

/* calculates number of alleles and number of heterozygous individuals at each locus
   across the entire population; calls function gene_data() */

void means2(int g,int plants,float **gn1,float **gn2,float **gn3,
      float **gn4,float **gn5,float **gn6,float **gn7,Vector  *Stype, Vector *vector_N_alleles)
  {
    int xc,yc,i,j;
    float *si_gene,*heteros,**n_genes;
    int sloc = Stype->size;
    int nloc = vector_N_alleles->size;
    n_genes = fmatrix(1,GENES-1,1,nloc);
    si_gene = fvector(1,sloc);
    heteros = fvector(1,GENES-1);
    for(i=1;i<GENES;i++)
      {
        heteros[i] = 0.0;
        for(j=1;j<=nloc;j++)
          n_genes[i][j] = 0.0;
      }
    for(i=1;i<=sloc;i++)
      si_gene[i] = 0.0;  
    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>0)
              {
                for(i=1;i<GENES;i++)
                  {
                    if(pop[xc][yc].gtype[i][0]!=pop[xc][yc].gtype[i][1])
                      heteros[i] += 1.0;
                    n_genes[i][pop[xc][yc].gtype[i][0]] += 1.0;
                    n_genes[i][pop[xc][yc].gtype[i][1]] += 1.0;
                  }
                si_gene[pop[xc][yc].gtype[0][0]] += 1.0;
                si_gene[pop[xc][yc].gtype[0][1]] += 1.0;
              }
          }
      }
    gene_data(g,plants,nloc,sloc,gn1,gn2,gn3,gn4,gn5,gn6,gn7,n_genes,si_gene,heteros);
    printf(" Flag means2 \n"); 
    free_fmatrix(n_genes,1,GENES-1,1,nloc);
    printf(" Flag means3 \n"); 
    free_fvector(si_gene,1,sloc);
    printf(" Flag means4 \n");
    free_fvector(heteros,1,GENES-1);
    printf(" Flag means5 \n");  
    return;
  }        

/* sums and squared sums of additional demographic variables for average dynamics
   assuming gametophytic self-incompatibility */

void means3_g(int g,int rep,int minr,int pdist,int sdist,float **ec1,float **ec2,
      float **ec3,int *nrp)
  {
    int i,j,k,l,xc,yc,dads=0,open=0,tmp;
    float mate=0.0;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)
              {
                for(i=xc-pdist;i<=xc+pdist;i++)
                  {
                    for(j=yc-pdist;j<=yc+pdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age>=minr)
                              {
                                tmp = 0;
                                dads += 1;
                                for(k=0;k<ALL;k++)
                                  for(l=0;l<ALL;l++)
                                    if(pop[xc][yc].gtype[0][k]!=pop[i][j].gtype[0][l])
                                      tmp += 1;
                                if(tmp==4)
                                  mate += 1.0;
                                else if(tmp==3)
                                  mate += 0.5;
                              }
                          }
                      }
                  }
                for(i=xc-sdist;i<=xc+sdist;i++)
                  {
                    for(j=yc-sdist;j<=yc+sdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age==0)
                              open += 1;
                          }
                      }
                  }  
              } 
          }
      }
    if(rep>0)
      {
        nrp[g] += 1;            /* number of runs with reproductives present */
        ec1[g][0] += open/rep;  /* empty sites w/in seed dispersal range */
        ec1[g][1] += (open/rep)*(open/rep);
        ec2[g][0] += dads/rep;  /* reproductive plants w/in pollen dispersal range */ 
        ec2[g][1] += (dads/rep)*(dads/rep);
        ec3[g][0] += mate/rep;  /* compatible males w/in pollen dispersal range */
        ec3[g][1] += (mate/rep)*(mate/rep);
      } 
    return;
  }

/* sums and squared sums of additional demographic variables for average dynamics
   assuming neutral self-incompatibility */

void means3_n(int g,int rep,int minr,int pdist,int sdist,float **ec1,float **ec2,
      float **ec3,int *nrp)
  {
    int i,j,xc,yc,dads=0,open=0;
    float mate=0.0;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)
              {
                for(i=xc-pdist;i<=xc+pdist;i++)
                  {
                    for(j=yc-pdist;j<=yc+pdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age>=minr)
                              dads += 1;
                          }
                      }
                  }
                for(i=xc-sdist;i<=xc+sdist;i++)
                  {
                    for(j=yc-sdist;j<=yc+sdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age==0)
                              open += 1;
                          }
                      }
                  }  
              } 
          }
      }
    mate = (float)dads;  
    if(rep>0)
      {
        nrp[g] += 1;            /* number of runs with reproductives present */
        ec1[g][0] += open/rep;  /* empty sites w/in seed dispersal range */
        ec1[g][1] += (open/rep)*(open/rep);
        ec2[g][0] += dads/rep;  /* reproductive plants w/in pollen dispersal range */ 
        ec2[g][1] += (dads/rep)*(dads/rep);
        ec3[g][0] += mate/rep;  /* compatible males w/in pollen dispersal range */
        ec3[g][1] += (mate/rep)*(mate/rep);
      } 
    return;
  }

/* sums and squared sums of additional demographic variables for average dynamics
   assuming codominant sporophytic self-incompatibility */

void means3_scd(int g,int rep,int minr,int pdist,int sdist,float **ec1,float **ec2,
      float **ec3,int *nrp)
  {
    int i,j,k,l,xc,yc,dads=0,open=0,flg;
    float mate=0.0;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)
              {
                for(i=xc-pdist;i<=xc+pdist;i++)
                  {
                    for(j=yc-pdist;j<=yc+pdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age>=minr)
                              {
                                dads += 1;
                                flg = IN;
                                for(k=0;k<ALL;k++)
                                  for(l=0;l<ALL;l++)
                                    if(pop[xc][yc].gtype[0][k]==pop[i][j].gtype[0][l])
                                      flg = OUT;
                                if(flg==IN)
                                  mate += 1.0;
                              }
                          }
                      }
                  }
                for(i=xc-sdist;i<=xc+sdist;i++)
                  {
                    for(j=yc-sdist;j<=yc+sdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age==0)
                              open += 1;
                          }
                      }
                  }  
              } 
          }
      }
    if(rep>0)
      {
        nrp[g] += 1;            /* number of runs with reproductives present */
        ec1[g][0] += open/rep;  /* empty sites w/in seed dispersal range */
        ec1[g][1] += (open/rep)*(open/rep);
        ec2[g][0] += dads/rep;  /* reproductive plants w/in pollen dispersal range */ 
        ec2[g][1] += (dads/rep)*(dads/rep);
        ec3[g][0] += mate/rep;  /* compatible males w/in pollen dispersal range */
        ec3[g][1] += (mate/rep)*(mate/rep);
      } 
    return;
  }

/* sums and squared sums of additional demographic variables for average dynamics
   assuming maternal dominant sporophytic self-incompatibility */

void means3_sd(int g,int rep,int minr,int pdist,int sdist,float **ec1,float **ec2,
      float **ec3,int *nrp)
  {
    int i,j,k,xc,yc,dads=0,open=0,flg,al;
    float mate=0.0;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)
              {
                for(i=xc-pdist;i<=xc+pdist;i++)
                  {
                    for(j=yc-pdist;j<=yc+pdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age>=minr)
                              {
                                dads += 1;
                                flg = IN;
                                if((al=pop[xc][yc].gtype[0][0])<pop[xc][yc].gtype[0][1])
                                  al = pop[xc][yc].gtype[0][1];
                                for(k=0;k<ALL;k++)
                                  if(pop[i][j].gtype[0][k]==al)
                                    flg = OUT;
                                if(flg==IN)
                                  mate += 1.0;
                              }
                          }
                      }
                  }
                for(i=xc-sdist;i<=xc+sdist;i++)
                  {
                    for(j=yc-sdist;j<=yc+sdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age==0)
                              open += 1;
                          }
                      }
                  }  
              } 
          }
      }
    if(rep>0)
      {
        nrp[g] += 1;            /* number of runs with reproductives present */
        ec1[g][0] += open/rep;  /* empty sites w/in seed dispersal range */
        ec1[g][1] += (open/rep)*(open/rep);
        ec2[g][0] += dads/rep;  /* reproductive plants w/in pollen dispersal range */ 
        ec2[g][1] += (dads/rep)*(dads/rep);
        ec3[g][0] += mate/rep;  /* compatible males w/in pollen dispersal range */
        ec3[g][1] += (mate/rep)*(mate/rep);
      } 
    return;
  }

/* sums and squared sums of additional demographic variables for average dynamics
   assuming paternal dominant sporophytic self-incompatibility */

void means3_sd2(int g,int rep,int minr,int pdist,int sdist,float **ec1,float **ec2,
      float **ec3,int *nrp)
  {
    int i,j,k,xc,yc,dads=0,open=0,flg,al;
    float mate=0.0;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)
              {
                for(i=xc-pdist;i<=xc+pdist;i++)
                  {
                    for(j=yc-pdist;j<=yc+pdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age>=minr)
                              {
                                dads += 1;
                                flg = IN;                        
                                if((al=pop[i][j].gtype[0][0])<pop[i][j].gtype[0][1])
                                  al = pop[i][j].gtype[0][1];
                                for(k=0;k<ALL;k++)
                                  if(pop[xc][yc].gtype[0][k]==al)
                                    flg = OUT;
                                if(flg==IN)
                                  mate += 1.0;
                              }
                          }
                      }
                  }
                for(i=xc-sdist;i<=xc+sdist;i++)
                  {
                    for(j=yc-sdist;j<=yc+sdist;j++)
                      {
                        if(i>=0 && i<LEN && j>=0 && j<LEN)
                          {
                            if(pop[i][j].age==0)
                              open += 1;
                          }
                      }
                  }  
              } 
          }
      }
    if(rep>0)
      {
        nrp[g] += 1;            /* number of runs with reproductives present */
        ec1[g][0] += open/rep;  /* empty sites w/in seed dispersal range */
        ec1[g][1] += (open/rep)*(open/rep);
        ec2[g][0] += dads/rep;  /* reproductive plants w/in pollen dispersal range */ 
        ec2[g][1] += (dads/rep)*(dads/rep);
        ec3[g][0] += mate/rep;  /* compatible males w/in pollen dispersal range */
        ec3[g][1] += (mate/rep)*(mate/rep);
      } 
    return;
  }

/* sums and squared sums of fitness variables (pre-dispersal) for average dynamics */

void means4(int g,int rep,int minr,float **ft3,float **ft4,float **ft5)
  {
    int xc,yc;
    float fsm1=0.0,fsm2=0.0,msm1=0.0,msm2=0.0,vfs,vms;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)
              {
                fsm1 += pop[xc][yc].fsd;
                fsm2 += (pop[xc][yc].fsd)*(pop[xc][yc].fsd);
                msm1 += pop[xc][yc].msd;
                msm2 += (pop[xc][yc].msd)*(pop[xc][yc].msd);
              }
            pop[xc][yc].fsd = 0;
            pop[xc][yc].msd = 0;  
          }
      }
    if(rep>0)
      {  
        ft3[g][0] += fsm1/rep;     /* mean seeds/parent before dispersal */
        ft3[g][1] += (fsm1/rep)*(fsm1/rep);
        (rep==1)?(vfs=0.0):(vfs=(fsm2-(fsm1*fsm1)/rep)/(rep-1));
        ft4[g][0] += vfs;           /* variance in seeds/father before dispersal */
        ft4[g][1] += vfs*vfs;
        (rep==1)?(vms=0.0):(vms=(msm2-(msm1*msm1)/rep)/(rep-1));
        ft5[g][0] += vms;           /* variance in seeds/mother before dispersal */
        ft5[g][1] += vms*vms;
      }
    return;
  }      

/* sums and squared sums of fitness variables (post-dispersal) for average dynamics */

void means5(int g,int rep,float **ft6,float **ft7,float **ft8)
  {
    int xc,yc,xn,yn;
    float fsp1=0.0,fsp2=0.0,msp1=0.0,msp2=0.0,vfs,vms;
    double ival,fval,len=LEN;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age==1)
              {
                fval = modf(pop[xc][yc].dad/len,&ival);
                xn = ival;
                yn = fval*LEN;
                pop[xn][yn].fsp += 1;
                fval = modf(pop[xc][yc].mom/len,&ival);
                xn = ival;
                yn = fval*LEN;
                pop[xn][yn].msp += 1;
              }
          }
      }
    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            fsp1 += pop[xc][yc].fsp;
            fsp2 += (pop[xc][yc].fsp)*(pop[xc][yc].fsp);
            msp1 += pop[xc][yc].msp;
            msp2 += (pop[xc][yc].msp)*(pop[xc][yc].msp);
            pop[xc][yc].fsp = 0;
            pop[xc][yc].msp = 0;  
          }
      }
    if(rep>0)
      {  
        ft6[g][0] += fsp1/rep;     /* mean seeds/parent after dispersal */
        ft6[g][1] += (fsp1/rep)*(fsp1/rep);
        (rep==1)?(vfs=0.0):(vfs=(fsp2-(fsp1*fsp1)/rep)/(rep-1));
        ft7[g][0] += vfs;           /* variance in seeds/father after dispersal */
        ft7[g][1] += vfs*vfs;
        (rep==1)?(vms=0.0):(vms=(msp2-(msp1*msp1)/rep)/(rep-1));
        ft8[g][0] += vms;           /* variance in seeds/mother after dispersal */
        ft8[g][1] += vms*vms;
      }  
    return;
  }

/* sums and squared sums of genetic variables for average dynamics */

void gene_data(int g,int plants,int nloc,int sloc,float **gn1,float **gn2,
      float **gn3,float **gn4,float **gn5,float **gn6,float **gn7,float **n_genes,
      float *si_gene,float *heteros)
  {
    int i,j;
    float NallL,TNall=0.0,SIall=0.0,Nvar=0.0,SIvar=0.0,Nmfrq,SIfrq=0.0,Ho=0.0;
    float He=0.0;
    
    for(i=1;i<GENES;i++)
      {
        NallL = 0.0;
        Nmfrq = 0.0;
        Ho += heteros[i]/plants;
        for(j=1;j<=nloc;j++)
          {
            if(n_genes[i][j]>0.0)
              {
                NallL += 1.0;
                Nmfrq += (n_genes[i][j]/(2.0*plants))*(n_genes[i][j]/(2.0*plants));
              }
          }
        He += 1.0 - Nmfrq;  
        if(NallL>1.0)  
          Nvar += (Nmfrq-1.0/NallL)/(NallL-1.0);
        TNall += NallL;
      }
    for(i=1;i<=sloc;i++)
      {
        if(si_gene[i]>0.0)
          {
            SIall += 1.0;
            SIfrq += (si_gene[i]/(2.0*plants))*(si_gene[i]/(2.0*plants));
          }
      }
    if(SIall>1.0)
      SIvar = (SIfrq-1.0/SIall)/(SIall-1.0); 
    gn1[g][0] += TNall/(GENES-1);   /* mean number of alleles per neutral locus */
    gn1[g][1] += (TNall/(GENES-1))*(TNall/(GENES-1));
    gn2[g][0] += SIall;             /* mean number of alleles per SI locus */
    gn2[g][1] += SIall*SIall;
    gn3[g][0] += Nvar/(GENES-1);    /* variance in neutral allele frequencies */
    gn3[g][1] += (Nvar/(GENES-1))*(Nvar/(GENES-1));
    gn4[g][0] += SIvar;             /* variance in SI allele frequencies */
    gn4[g][1] += SIvar*SIvar;
    gn5[g][0] += Ho/(GENES-1);      /* observed heterozygosity (neutral loci) */
    gn5[g][1] += (Ho/(GENES-1))*(Ho/(GENES-1));
    gn6[g][0] += He/(GENES-1);      /* expected heterozygosity (neutral loci) */
    gn6[g][1] += (He/(GENES-1))*(He/(GENES-1));
    if(He>0.0)
      {
        gn7[g][0] += 1-(Ho/He);     /* Fis (neutral loci) */
        gn7[g][1] += (1-(Ho/He))*(1-(Ho/He));
      }
    return;
  }   
          
/* assumes gametophytic self-incompatibility and random pollen donors: for each
   reproductive adult, calls functions that determines the local pollen pool, create
   ovules and pollen grains, makes seeds and disperses them to empty sites; also
   calculates information used for fitness data output */

void reproduce_gr(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_r(xc,yc,pdist,xyval,probs,minr,dads);   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;       
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    for(j=0;j<ALL;j++)
                      if(ptype[0]==pop[xc][yc].gtype[0][j])
                        flg=OUT;
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes neutral self-incompatibility and random pollen donors: for each
   reproductive adult, calls functions that determines the local pollen pool, create
   ovules and pollen grains, makes seeds and disperses them to empty sites; also
   calculates information used for fitness data output */

void reproduce_nr(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_r(xc,yc,pdist,xyval,probs,minr,dads);   
                for(i=0;i<ovn;i++)
                  {  
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    donor[pd] = IN;
                    pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                    pop[xc][yc].msd += 1;
                    make_ovule(r,xc,yc,otype,&mid);
                    make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes co-dominant sporophytic self-incompatibility and random pollen donors: for
   each reproductive adult, calls functions that determines the local pollen pool,
   create ovules and pollen grains, makes seeds and disperses them to empty sites;
   also calculates information used for fitness data output */

void reproduce_scdr(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,k,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_r(xc,yc,pdist,xyval,probs,minr,dads);   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;       
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    for(j=0;j<ALL;j++)
                      {
                        for(k=0;k<ALL;k++)
                          if(pop[xyval[pd][0]][xyval[pd][1]].gtype[0][j]==\
                             pop[xc][yc].gtype[0][k])
                            flg=OUT;
                      }
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes maternally dominant sporophytic self-incompatibility and random pollen
   donors: for each reproductive adult, calls functions that determines the local
   pollen pool, create ovules and pollen grains, makes seeds and disperses them to
   empty sites; also calculates information used for fitness data output */

void reproduce_sdr(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd,md;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_r(xc,yc,pdist,xyval,probs,minr,dads);
                //printf("FLAG 1\n");   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;
                    //printf("FLAG POLLEN\n");
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    //printf("FLAG 2\n");       
                    if((md=pop[xc][yc].gtype[0][0])<pop[xc][yc].gtype[0][1])
                      md = pop[xc][yc].gtype[0][1];
                      //printf("MD: %d\n", md);
                    for(j=0;j<ALL;j++)
                      if(pop[xyval[pd][0]][xyval[pd][1]].gtype[0][j]==md)
                        flg=OUT; 
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        //printf("FLAG 3\n"); 
                        make_ovule(r,xc,yc,otype,&mid);
                        //printf("FLAG 4\n"); 
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                        //printf("FLAG 5\n"); 
                      }
                  }
                //printf("DADS %d\n", dads);  
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                //printf("FLAG 6\n");    
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes paternally dominant sporophytic self-incompatibility and random pollen
   donors: for each reproductive adult, calls functions that determines the local
   pollen pool, create ovules and pollen grains, makes seeds and disperses them to
   empty sites; also calculates information used for fitness data output */

void reproduce_sdr2(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd,md;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_r(xc,yc,pdist,xyval,probs,minr,dads);   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);       
                    md = pop[xyval[pd][0]][xyval[pd][1]].gtype[0][0];    
                    if(md<pop[xyval[pd][0]][xyval[pd][1]].gtype[0][1])    
                      md = pop[xyval[pd][0]][xyval[pd][1]].gtype[0][1];
                    for(j=0;j<ALL;j++)
                      if(pop[xc][yc].gtype[0][j]==md)
                        flg=OUT;    
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes gametophytic self-incompatibility and size-dependent pollen donors: for each
   reproductive adult, calls functions that determines the local pollen pool,
   create ovules and pollen grains, makes seeds and disperses them to empty sites;
   also calculates information used for fitness data output */

void reproduce_gs(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_s(xc,yc,pdist,xyval,probs,minr,tage);      
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;       
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    for(j=0;j<ALL;j++)
                      if(ptype[0]==pop[xc][yc].gtype[0][j])
                        flg=OUT;
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes neutral self-incompatibility and size-dependent pollen donors: for each
   reproductive adult, calls functions that determines the local pollen pool,
   create ovules and pollen grains, makes seeds and disperses them to empty sites;
   also calculates information used for fitness data output */

void reproduce_ns(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_s(xc,yc,pdist,xyval,probs,minr,tage);      
                for(i=0;i<ovn;i++)
                  {    
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    donor[pd] = IN;
                    pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                    pop[xc][yc].msd += 1;
                    make_ovule(r,xc,yc,otype,&mid);
                    make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes co-dominant sporophytic self-incompatibility and size-dependent pollen
   donors: for each reproductive adult, calls functions that determines the local pollen
   pool, create ovules and pollen grains, makes seeds and disperses them to empty sites;
   also calculates information used for fitness data output */

void reproduce_scds(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,k,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_s(xc,yc,pdist,xyval,probs,minr,tage);   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;       
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    for(j=0;j<ALL;j++)
                      {
                        for(k=0;k<ALL;k++)
                          if(pop[xyval[pd][0]][xyval[pd][1]].gtype[0][j]==\
                             pop[xc][yc].gtype[0][k])
                            flg=OUT;
                      }
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes maternally dominant sporophytic self-incompatibility and size-dependent
   pollen donors: for each reproductive adult, calls functions that determines the
   local pollen pool, create ovules and pollen grains, makes seeds and disperses them
   to empty sites; also calculates information used for fitness data output */

void reproduce_sds(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd,md;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_s(xc,yc,pdist,xyval,probs,minr,tage);   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);       
                    if((md=pop[xc][yc].gtype[0][0])<pop[xc][yc].gtype[0][1])
                      md = pop[xc][yc].gtype[0][1];
                    for(j=0;j<ALL;j++)
                      if(pop[xyval[pd][0]][xyval[pd][1]].gtype[0][j]==md)
                        flg=OUT;
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes paternally dominant sporophytic self-incompatibility and size-dependent
   pollen donors: for each reproductive adult, calls functions that determines the
   local pollen pool, create ovules and pollen grains, makes seeds and disperses them
   to empty sites; also calculates information used for fitness data output */

void reproduce_sds2(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd,md;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_s(xc,yc,pdist,xyval,probs,minr,tage);   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);       
                    md = pop[xyval[pd][0]][xyval[pd][1]].gtype[0][0];    
                    if(md<pop[xyval[pd][0]][xyval[pd][1]].gtype[0][1])    
                      md = pop[xyval[pd][0]][xyval[pd][1]].gtype[0][1];
                    for(j=0;j<ALL;j++)
                      if(pop[xc][yc].gtype[0][j]==md)
                        flg=OUT;  
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes gametophytic self-incompatibility and distance-dependent pollen donors: for
   each reproductive adult, calls functions that determines the local pollen pool,
   create ovules and pollen grains, makes seeds and disperses them to empty sites; also
   calculates information used for fitness data output */

void reproduce_gd(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_d(xc,yc,pdist,xyval,probs,minr,tdst);      
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;       
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    for(j=0;j<ALL;j++)
                      if(ptype[0]==pop[xc][yc].gtype[0][j])
                        flg=OUT;
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes neutral self-incompatibility and distance-dependent pollen donors: for
   each reproductive adult, calls functions that determines the local pollen pool,
   create ovules and pollen grains, makes seeds and disperses them to empty sites; also
   calculates information used for fitness data output */

void reproduce_nd(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_d(xc,yc,pdist,xyval,probs,minr,tdst);      
                for(i=0;i<ovn;i++)
                  {     
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    donor[pd] = IN;
                    pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                    pop[xc][yc].msd += 1;
                    make_ovule(r,xc,yc,otype,&mid);
                    make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes co-dominant sporophytic self-incompatibility and distance-dependent pollen
   donors: for each reproductive adult, calls functions that determines the local pollen
   pool, create ovules and pollen grains, makes seeds and disperses them to empty sites;
   also calculates information used for fitness data output */

void reproduce_scdd(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,k,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_d(xc,yc,pdist,xyval,probs,minr,tdst);   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;       
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);
                    for(j=0;j<ALL;j++)
                      {
                        for(k=0;k<ALL;k++)
                          if(pop[xyval[pd][0]][xyval[pd][1]].gtype[0][j]==\
                             pop[xc][yc].gtype[0][k])
                            flg=OUT;
                      }
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes maternally dominant sporophytic self-incompatibility and distance-dependent
   pollen donors: for each reproductive adult, calls functions that determines the local
   pollen pool, create ovules and pollen grains, makes seeds and disperses them to empty
   sites; also calculates information used for fitness data output */

void reproduce_sdd(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd,md;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_d(xc,yc,pdist,xyval,probs,minr,tdst);
                //printf("FLAG_sdd\n");   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);       
                    if((md=pop[xc][yc].gtype[0][0])<pop[xc][yc].gtype[0][1])
                      md = pop[xc][yc].gtype[0][1];
                    for(j=0;j<ALL;j++)
                      if(pop[xyval[pd][0]][xyval[pd][1]].gtype[0][j]==md)
                        flg=OUT;
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* assumes paternally dominant sporophytic self-incompatibility and distance-dependent
   pollen donors: for each reproductive adult, calls functions that determines the local
   pollen pool, create ovules and pollen grains, makes seeds and disperses them to empty
   sites; also calculates information used for fitness data output */

void reproduce_sdd2(int g,int rep,int minr,int pdist,int sdist,float **ft1,
     float **ft2,double b0,double d)
  {
    int i,j,flg,xc,yc,otype[GENES],ptype[GENES],dads,tage,**xyval,did,mid;
    int *donor,pd=0,lpd,md;
    float ovn,*probs,tpd1=0.0,tpd2=0.0,vpd,tdst;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>=minr)         
              {
                lpd = 0;
                (pop[xc][yc].age<=(1.0/d))?(ovn=b0*d*pop[xc][yc].age):(ovn=b0);
                ((ovn-floor(ovn))<0.5)?(ovn=floor(ovn)):(ovn=ceil(ovn));
                n_dads(xc,yc,pdist,minr,&dads,&tage,&tdst);
                xyval = imatrix(0,dads-1,0,1); 
                probs = fvector(0,dads-1);
                donor = ivector(0,dads-1);
                for(i=0;i<dads;i++)
                  donor[i] = OUT;
                make_pollen_pool_d(xc,yc,pdist,xyval,probs,minr,tdst);   
                for(i=0;i<ovn;i++)
                  {
                    flg = IN;
                    pd = make_pollen(r,ptype,xyval,probs,dads,&did);       
                    md = pop[xyval[pd][0]][xyval[pd][1]].gtype[0][0];    
                    if(md<pop[xyval[pd][0]][xyval[pd][1]].gtype[0][1])    
                      md = pop[xyval[pd][0]][xyval[pd][1]].gtype[0][1];
                    for(j=0;j<ALL;j++)
                      if(pop[xc][yc].gtype[0][j]==md)
                        flg=OUT; 
                    if(flg==IN)
                      {
                        donor[pd] = IN;
                        pop[xyval[pd][0]][xyval[pd][1]].fsd += 1;
                        pop[xc][yc].msd += 1;
                        make_ovule(r,xc,yc,otype,&mid);
                        make_seed(r,xc,yc,ptype,otype,sdist,did,mid);
                      }
                  }
                for(i=0;i<dads;i++) /* count for local pollen donors */
                  if(donor[i]==IN)
                    lpd += 1;
                tpd1 += lpd;   
                tpd2 += lpd*lpd;       
                free_imatrix(xyval,0,dads-1,0,1);
                free_fvector(probs,0,dads-1);
                free_ivector(donor,0,dads-1);
              }
          }
      }
    if(rep>0)
      {
        ft1[g][0] += tpd1/rep;     /* average number of pollen donors */
        ft1[g][1] += (tpd1/rep)*(tpd1/rep);
        (rep==1)?(vpd=0.0):(vpd=(tpd2-(tpd1*tpd1)/rep)/(rep-1));  
        ft2[g][0] += vpd;           /* variance in number of pollen donors */
        ft2[g][1] += vpd*vpd; 
      } 
    return;
  }

/* calculates # potential pollen donors within a local area determined by pdist; sums
   ages for calculating probabilities of size-dependent pollen dispersal; sums distances
   for calculating probabilities of distance-dependent pollen dispersal */

void n_dads(int xc,int yc,int pdist,int minr,int *dads,int *tage,float *tdst)
  {
    int i,j;

    *dads = 0;
    *tage = 0;
    *tdst = 0.0;
    for(i=xc-pdist;i<=xc+pdist;i++)
      {
        for(j=yc-pdist;j<=yc+pdist;j++)
          {
            if(i>=0 && i<LEN && j>=0 && j<LEN)
              {
                if(pop[i][j].age>=minr)
                  {
                    *dads += 1;
                    *tage += pop[i][j].age;
                    *tdst += sqrt((xc-i)*(xc-i)+(yc-j)*(yc-j));
                  }
              }
          }
      }
    return;
  }

/* makes an ovule by randomly choosing one of the alleles at each locus of a
   maternal parent (stored in otype[]) */
  
void make_ovule(gsl_rng *r, int xc,int yc,int *otype,int *mid)
  {
    int i,al;

    for(i=0;i<GENES;i++)
      {
        al = gsl_rng_uniform_int(r,2);//g05dyc(0,1);
        otype[i] = pop[xc][yc].gtype[i][al];
      }
    *mid = LEN*xc + yc;
    //gsl_rng_free (r);    
    return;
  }

/* constructs the pool of pollen donors and calculates their relative probabilities
   of reproductive success assuming no dependence on size/age or distance */

void make_pollen_pool_r(int xc,int yc,int pdist,int **xyval,float *probs,int minr,
                      int dads)
  {
    int i,j,cnt=0;
    
    for(i=xc-pdist;i<=xc+pdist;i++)
      {
        for(j=yc-pdist;j<=yc+pdist;j++)
          {        
            if(i>=0 && i<LEN && j>=0 && j<LEN) 
              {
                if(pop[i][j].age>=minr)
                  {
                    probs[cnt] = (float)1.0/dads;
                    if(cnt>0)
                      probs[cnt] += probs[cnt-1];
                    xyval[cnt][0] = i;
                    xyval[cnt][1] = j;
                    cnt += 1;
                  }
              }
          }
      }
    return;
  }

/* constructs the pool of pollen donors and calculates their relative probabilities
   of reproductive success assuming size/age-dependence */

void make_pollen_pool_s(int xc,int yc,int pdist,int **xyval,float *probs,int minr,
                      int tage)
  {
    int i,j,cnt=0;
    
    for(i=xc-pdist;i<=xc+pdist;i++)
      {
        for(j=yc-pdist;j<=yc+pdist;j++)
          {        
            if(i>=0 && i<LEN && j>=0 && j<LEN) 
              {
                if(pop[i][j].age>=minr)
                  {
                    probs[cnt] = (float)pop[i][j].age/tage;
                    if(cnt>0)
                      probs[cnt] += probs[cnt-1];
                    xyval[cnt][0] = i;
                    xyval[cnt][1] = j;
                    cnt += 1;
                  }
              }
          }
      }
    return;
  }

/* constructs the pool of pollen donors and calculates their relative probabilities
   of reproductive success assuming distance-dependence */

void make_pollen_pool_d(int xc,int yc,int pdist,int **xyval,float *probs,int minr,
                      float tdst)
  {
    int i,j,cnt=0,dst;
    
    for(i=xc-pdist;i<=xc+pdist;i++)
      {
        for(j=yc-pdist;j<=yc+pdist;j++)
          {        
            if(i>=0 && i<LEN && j>=0 && j<LEN) 
              {
                if(pop[i][j].age>=minr)
                  {
                    dst = sqrt((xc-i)*(xc-i)+(yc-j)*(yc-j));
                    probs[cnt] = 1.0 - dst/tdst;
                    if(cnt>0)
                      probs[cnt] += probs[cnt-1];
                    xyval[cnt][0] = i;
                    xyval[cnt][1] = j;
                    cnt += 1;
                  }
              }
          }
      }
    return;
  }

/* randomly chooses a pollen donor and makes a pollen grain using one of the alleles
   at each locus of a paternal parent (stored in ptype[]) */

int make_pollen(gsl_rng *r, int *ptype,int **xyval,float *probs,int dads, int *did)
  {
    int i,j,al,flg=0;
    double rnum = 1.0;//gsl_rng_uniform(r) - 0.001;
    
    if(rnum == 1.0){
	   //printf(" For Rnum1= %.3f\n",rnum);
       while(rnum == 1.0){
          rnum = gsl_rng_uniform(r) - 0.001;//g05cac();
          //printf(" For Rnum1= %.3f\n",rnum);
	   }
    }
	
    for(i=0;i<=dads && !flg;i++){
      //printf("FLG = %d\n",flg);
      flg = rnum <= probs[i];
      //printf(" For dads= %d\n",flg);
      //printf(" For Probs= %.3f\n",probs[i]);
      //printf(" For Rnum2= %.3f\n",rnum);
    }
    for(j=0;j<GENES;j++)
      {
        al = gsl_rng_uniform_int(r,2);//g05dyc(0,1);
        //printf("xyval x: %d", xyval[i-1][0]);
        //printf(" y: %d\n", xyval[i-1][1]);
        ptype[j] = pop[xyval[i-1][0]][xyval[i-1][1]].gtype[j][al];
        //printf("HIT ! gen=%d\n",j);
      }
    *did = LEN*xyval[i-1][0] + xyval[i-1][1];
    //gsl_rng_free (r);      
    return i-1;
  }

/* makes a seed and disperses it into an empty site if possible */

void make_seed(gsl_rng *r, int xc,int yc,int *ptype,int *otype,int sdist,int did,int mid)
  {
    int i,xn,yn;

    xn = gsl_rng_uniform_int(r,sdist+1) + xc;//g05dyc(xc-sdist,xc+sdist);
    yn = gsl_rng_uniform_int(r,sdist+1) + yc;//g05dyc(yc-sdist,yc+sdist);
    if(xn>=0 && xn<LEN && yn>=0 && yn<LEN)
      {
        if(pop[xn][yn].safe==1 && pop[xn][yn].age==0 && pop[xn][yn].sds<SDMAX) 
          {
            for(i=0;i<GENES;i++)
              {
                pop[xn][yn].stype[pop[xn][yn].sds][i][0] = otype[i];
                pop[xn][yn].stype[pop[xn][yn].sds][i][1] = ptype[i];
              }
            pop[xn][yn].stype[pop[xn][yn].sds][GENES][0] = mid;
            pop[xn][yn].stype[pop[xn][yn].sds][GENES][1] = did;
            pop[xn][yn].sds += 1;
          }
      }
    //gsl_rng_free (r);        
    return;
  }

/* determines mortality for adult plants assuming constant death rate */

void kill_plants(gsl_rng *r, double d)
  {
    int i,j,xc,yc;
    //double rnum;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].age>0)
              {
                if((gsl_rng_uniform(r))<=d) //g05cac()
                  {
                    pop[xc][yc].age = 0;
                    pop[xc][yc].mom = 0;
                    pop[xc][yc].dad = 0;          
                    for(i=0;i<GENES;i++)
                      {
                        for(j=0;j<ALL;j++)
                          pop[xc][yc].gtype[i][j] = 0;
                      }
                  }
                else if(gsl_rng_uniform(r)>d) //rnum>d
                  pop[xc][yc].age += 1;  
              }
          }
      }
    //gsl_rng_free (r);     
    return;
  }  

/* adds new seedlings at time t+1 to the population (after adult mortality) */
/* Function to add probability of germination from seed bank */

void new_plants(gsl_rng *r,float est, float  mort_seed_seedbank_prob)
  {
    int i,j,k,xc,yc,sn;
    double rnum, seedbank_prob = 0.5;

    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].sds>0 || pop[xc][yc].seedbank>0) // sds is number of seeds per individual if(pop[xc][yc].sds>0 || pop[xc][yc].seedbank[]>0)
              {
                if((rnum=gsl_rng_uniform(r))<=est) //g05cac()
                  {
					  if(pop[xc][yc].seedbank>0){
					       //printf("SDS %d \n", pop[xc][yc].sds);
					       //printf("Seedbank %d \n", pop[xc][yc].seedbank);
					       if(seedbank_prob < gsl_rng_uniform(r) && pop[xc][yc].sds>0){   
                               sn = gsl_rng_uniform_int(r,pop[xc][yc].sds);//g05dyc(0,pop[xc][yc].sds-1); Pick one seed randomly to survive
                               for(i=0;i<GENES;i++)
                               {
                                  for(j=0;j<ALL;j++)
                                     pop[xc][yc].gtype[i][j] = pop[xc][yc].stype[sn][i][j];
                               }
                               pop[xc][yc].mom = pop[xc][yc].stype[sn][GENES][0];
                               pop[xc][yc].dad = pop[xc][yc].stype[sn][GENES][1]; 
                               pop[xc][yc].age = 1;
					       }else{
						       sn = gsl_rng_uniform_int(r,pop[xc][yc].seedbank);// Pick one seed randomly to survive
						       for(i=0;i<GENES;i++)
                              {
                                 for(j=0;j<ALL;j++)
                                    pop[xc][yc].gtype[i][j] = pop[xc][yc].seedbankid[sn][i][j];
                              }
                              pop[xc][yc].mom = pop[xc][yc].seedbankid[sn][GENES][0];
                              pop[xc][yc].dad = pop[xc][yc].seedbankid[sn][GENES][1]; 
                              pop[xc][yc].age = 1;
					       }
					  }else{
						      sn = gsl_rng_uniform_int(r,pop[xc][yc].sds);//g05dyc(0,pop[xc][yc].sds-1); Pick one seed randomly to survive
                              for(i=0;i<GENES;i++)
                               {
                                  for(j=0;j<ALL;j++)
                                     pop[xc][yc].gtype[i][j] = pop[xc][yc].stype[sn][i][j];
                               }
                               pop[xc][yc].mom = pop[xc][yc].stype[sn][GENES][0];
                               pop[xc][yc].dad = pop[xc][yc].stype[sn][GENES][1]; 
                               pop[xc][yc].age = 1;
					  }
                  }else{ // If no seed germinates then they become part of the seed bank, next gen they need to be considered as part of the pool of seeds to germinate
					  pop[xc][yc].seedbank += pop[xc][yc].sds;
					  for(i=0;i<SDMAX;i++){
						  for(j=0;j<=GENES;j++)
                           {
                              for(k=0;k<ALL;k++)
                                  pop[xc][yc].seedbankid[i][j][k] = pop[xc][yc].stype[i][j][k]; // Given genetic identity to seed bank seeds
                           }
                           pop[xc][yc].seedbank_age[i] += 1;  
					  }
				  }
                pop[xc][yc].sds = 0;  // Reset everything (seeds disappear) to zero in that cell
                for(i=0;i<SDMAX;i++)
                  {
                    for(j=0;j<=GENES;j++)
                      {
                        for(k=0;k<ALL;k++)
                          pop[xc][yc].stype[i][j][k] = 0;
                      }
                  }  
              }
              // Seedbank seeds that are older than 10 years old die
              for(i=0;i<SDMAX;i++){
				  if(pop[xc][yc].seedbank_age[i] > 10){
					  if(pop[xc][yc].seedbank > 0){
					      pop[xc][yc].seedbank -= 1;
					      //printf("Seedbank reset");
					  }
					  pop[xc][yc].seedbank_age[i] = 0;
					  for(j=0;j<=GENES;j++){
                              for(k=0;k<ALL;k++)
                                 pop[xc][yc].seedbankid[i][j][k] = 0;
					  }
					  
				  }else{ //Seedbank seeds mortality filter
					  if(mort_seed_seedbank_prob > gsl_rng_uniform(r)){
					     if(pop[xc][yc].seedbank > 0){
					          pop[xc][yc].seedbank -= 1;
					      //printf("Seedbank reset");
					     }
					     pop[xc][yc].seedbank_age[i] = 0;
					     for(j=0;j<=GENES;j++){
                              for(k=0;k<ALL;k++)
                                 pop[xc][yc].seedbankid[i][j][k] = 0;
					      }
					  }
				  }
			  }
          }
      }
    //gsl_rng_free (r);      
    return;
  }

/* allocates a float vector */

float *fvector(long nl,long nh)
  {
    float *v;

    v = (float *) malloc((size_t)((nh-nl+1+NR_END)*sizeof(float)));
    if(!v)
      nrerror("allocation failure in fvector()");
    return v-nl+NR_END;
  }

/* frees a float vector allocated by fvector() */

void free_fvector(float *v,long int nl,long int nh)
  {
    free((FREE_ARG) (v+nl-NR_END));
  }

/* allocates an int vector */

int *ivector(long nl,long nh)
  {
    int *v;

    v = (int *) malloc((size_t)((nh-nl+1+NR_END)*sizeof(int)));
    if(!v)
      nrerror("allocation failure in ivector()");
    return v-nl+NR_END;
  }

/* frees an int vector allocated by ivector() */

void free_ivector(int *v,long int nl,long int nh)
  {
    free((FREE_ARG) (v+nl-NR_END));
  }

/* allocates an int matrix with subscript range m[nrl..nrh][ncl..nch] */

int **imatrix(long int nrl, long int nrh, long int ncl, long int nch)
  {
    long int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
    int **m;

    m = (int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
    if(!m)
      nrerror("allocation failure 1 in imatrix()");
    m += NR_END;
    m -= nrl;
    m[nrl] = (int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
    if(!m[nrl])
      nrerror("allocation failure 2 in imatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++)
      m[i] = m[i-1] + ncol;
    return m;
  }

/* frees an int matrix allocated by imatrix() */

void free_imatrix(int **m, long int nrl, long int nrh, long int ncl,
                  long int nch)
  {
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
  }

/* allocates a float matrix with subscript range m[nrl..nrh][ncl..nch] */

float **fmatrix(long int nrl, long int nrh, long int ncl, long int nch)
  {
    long int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
    float **m;

    m = (float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
    if(!m)
      nrerror("allocation failure 1 in fmatrix()");
    m += NR_END;
    m -= nrl;
    m[nrl] = (float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
    if(!m[nrl])
      nrerror("allocation failure 2 in fmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++)
      m[i] = m[i-1] + ncol;
    return m;
  }

/* frees a float matrix allocated by fmatrix() */

void free_fmatrix(float **m, long int nrl, long int nrh, long int ncl,
                  long int nch)
  {
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
  }

/* Numerical Recipes standard error handler */

void nrerror(char error_text[])
  {
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
  }

/* calculates means and standard errors for demographic data; prints to file fpw1 */

void grand_mean1(int gen,float **mnv,float **mnr,float **myr,int *nrn,float **ec1,
      float **ec2,float **ec3,int *nrp, Vector *vector_quasi_ext)
  {
    int g;
    float se;

    fprintf(fpw1,"GEN\tNRUN\tTVEG\tERR1\tTREP\tERR2\tAGE\tERR3\t");
    fprintf(fpw1,"EMPTY\tERR4\tLREP\tERR5\tMATE\tERR6\tQUASI\n");
    for(g=0;g<=gen;g++)
      {
        fprintf(fpw1,"%d\t%d\t",g,nrn[g]);
        if(nrn[g]==0)
          fprintf(fpw1,".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n");    
        else if(nrn[g]==1)
          {
            fprintf(fpw1,"%.1f\t.\t%.1f\t.\t%.1f\t.\t",mnv[g][0],mnr[g][0],myr[g][0]);
            if(nrp[g]==0)
              fprintf(fpw1,".\t.\t.\t.\t.\t.\n");
            else if(nrp[g]==1)
              fprintf(fpw1,"%.1f\t.\t%.1f\t.\t%.1f\t%3d\t\n",ec1[g][0],ec2[g][0],ec3[g][0],vector_get(vector_quasi_ext,g)); 
          }
        else
          {
            fprintf(fpw1,"%.1f\t",(mnv[g][0]/nrn[g]));
            mnv[g][0] *= mnv[g][0];
            se = (sqrt(fabs((mnv[g][1]-(mnv[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw1,"%.3f\t",se);
            fprintf(fpw1,"%.1f\t",(mnr[g][0]/nrn[g]));
            mnr[g][0] *= mnr[g][0];
            se = (sqrt(fabs((mnr[g][1]-(mnr[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw1,"%.3f\t",se);
            fprintf(fpw1,"%.1f\t",(myr[g][0]/nrn[g]));
            myr[g][0] *= myr[g][0];
            se = (sqrt(fabs((myr[g][1]-(myr[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw1,"%.3f\t",se);    
            if(nrp[g]==0)
              fprintf(fpw1,".\t.\t.\t.\t.\t.\n");
            else if(nrp[g]==1)
              fprintf(fpw1,"%.1f\t.\t%.1f\t.\t%.1f\t%3d\t\n",ec1[g][0],ec2[g][0],ec3[g][0],vector_get(vector_quasi_ext,g)); 
            else
              {    
                fprintf(fpw1,"%.1f\t",(ec1[g][0]/nrp[g]));    
                ec1[g][0] *= ec1[g][0];
                se=(sqrt(fabs((ec1[g][1]-(ec1[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
                fprintf(fpw1,"%.3f\t",se);  
                fprintf(fpw1,"%.1f\t",(ec2[g][0]/nrp[g]));
                ec2[g][0] *= ec2[g][0];
                se=(sqrt(fabs((ec2[g][1]-(ec2[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
                fprintf(fpw1,"%.3f\t",se);
                fprintf(fpw1,"%.1f\t",(ec3[g][0]/nrp[g]));
                ec3[g][0] *= ec3[g][0];
                se=(sqrt(fabs((ec3[g][1]-(ec3[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
                fprintf(fpw1,"%.3f\t",se);
                fprintf(fpw1,"%3d\n",vector_get(vector_quasi_ext,g));
              }  
          }
      }
    return;
  }

/* calculates means and standard errors for genetic data; prints to file fpw2 */

void grand_mean2(int gen,int *nrn,float **gn1,float **gn2,float **gn3,float **gn4,
      float **gn5,float **gn6,float **gn7)
  {
    int g;
    float se;

    fprintf(fpw2,"GEN\tNRUN\tNGENE\tERR1\tSGENE\tERR2\tNVAR\tERR3\tSVAR\tERR4\t");
    fprintf(fpw2,"HO\tERR5\tHE\tERR6\tFIS\tERR7\n");
    for(g=0;g<=gen;g++)
      {
        fprintf(fpw2,"%d\t%d\t",g,nrn[g]);
        if(nrn[g]==0)
          fprintf(fpw2,".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n");    
        else if(nrn[g]==1)
          {
            fprintf(fpw2,"%.3f\t.\t%.3f\t.\t%.3f\t.\t%.3f\t.\t%.3f\t.\t",gn1[g][0],
                    gn2[g][0],gn3[g][0],gn4[g][0],gn5[g][0]);
            fprintf(fpw2,"%.3f\t.\t%.3f\t.\n",gn6[g][0],gn7[g][0]);
          }               
        else
          {
            fprintf(fpw2,"%.3f\t",(gn1[g][0]/nrn[g]));
            gn1[g][0] *= gn1[g][0];
            se = (sqrt(fabs((gn1[g][1]-(gn1[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw2,"%.3f\t",se);
            fprintf(fpw2,"%.3f\t",(gn2[g][0]/nrn[g]));
            gn2[g][0] *= gn2[g][0];
            se = (sqrt(fabs((gn2[g][1]-(gn2[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw2,"%.3f\t",se);
            fprintf(fpw2,"%.3f\t",(gn3[g][0]/nrn[g]));
            gn3[g][0] *= gn3[g][0];
            se = (sqrt(fabs((gn3[g][1]-(gn3[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw2,"%.3f\t",se);
            fprintf(fpw2,"%.3f\t",(gn4[g][0]/nrn[g]));
            gn4[g][0] *= gn4[g][0];
            se = (sqrt(fabs((gn4[g][1]-(gn4[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw2,"%.3f\t",se);                
            fprintf(fpw2,"%.3f\t",(gn5[g][0]/nrn[g]));
            gn5[g][0] *= gn5[g][0];
            se = (sqrt(fabs((gn5[g][1]-(gn5[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw2,"%.3f\t",se);                
            fprintf(fpw2,"%.3f\t",(gn6[g][0]/nrn[g]));
            gn6[g][0] *= gn6[g][0];
            se = (sqrt(fabs((gn6[g][1]-(gn6[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw2,"%.3f\t",se);
            fprintf(fpw2,"%.3f\t",(gn7[g][0]/nrn[g]));
            gn7[g][0] *= gn7[g][0];
            se = (sqrt(fabs((gn7[g][1]-(gn7[g][0]/nrn[g]))/(nrn[g]-1))))/sqrt(nrn[g]);
            fprintf(fpw2,"%.3f\n",se);      
          }
      }
    return;
  }

/* calculates means and standard errors for fitness data; prints to file fpw3 */

void grand_mean3(int gen,int *nrp,float **ft1,float **ft2,float **ft3,float **ft4,
      float **ft5,float **ft6,float **ft7,float **ft8)
  {
    int g;
    float se;

    fprintf(fpw3,"GEN\tNPOLL\tERR1\tVPOLL\tERR2\tSDS1\tERR3\tVFSD1\tERR4\t");
    fprintf(fpw3,"VMSD1\tERR5\tSDS2\tERR6\tVFSD2\tERR7\tVMSD2\tERR8\n");
    for(g=0;g<=gen;g++)
      {
        fprintf(fpw3,"%d\t",g);
        if(nrp[g]==0)
          fprintf(fpw3,".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n");    
        else if(nrp[g]==1)
          {
            fprintf(fpw3,"%.3f\t.\t%.3f\t.\t%.3f\t.\t",ft1[g][0],ft2[g][0],ft3[g][0]);
            fprintf(fpw3,"%.3f\t.\t%.3f\t.\t%.3f\t.\t%.3f\t.\t%.3f\t.\n",ft4[g][0],
                    ft5[g][0],ft6[g][0],ft7[g][0],ft8[g][0]);
          }               
        else
          {
            fprintf(fpw3,"%.3f\t",(ft1[g][0]/nrp[g]));
            ft1[g][0] *= ft1[g][0];
            se = (sqrt(fabs((ft1[g][1]-(ft1[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
            fprintf(fpw3,"%.3f\t",se);
            fprintf(fpw3,"%.3f\t",(ft2[g][0]/nrp[g]));
            ft2[g][0] *= ft2[g][0];
            se = (sqrt(fabs((ft2[g][1]-(ft2[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
            fprintf(fpw3,"%.3f\t",se);
            fprintf(fpw3,"%.3f\t",(ft3[g][0]/nrp[g]));
            ft3[g][0] *= ft3[g][0];
            se = (sqrt(fabs((ft3[g][1]-(ft3[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
            fprintf(fpw3,"%.3f\t",se);
            fprintf(fpw3,"%.3f\t",(ft4[g][0]/nrp[g]));
            ft4[g][0] *= ft4[g][0];
            se = (sqrt(fabs((ft4[g][1]-(ft4[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
            fprintf(fpw3,"%.3f\t",se);
            fprintf(fpw3,"%.3f\t",(ft5[g][0]/nrp[g]));
            ft5[g][0] *= ft5[g][0];
            se = (sqrt(fabs((ft5[g][1]-(ft5[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
            fprintf(fpw3,"%.3f\t",se);
            fprintf(fpw3,"%.3f\t",(ft6[g][0]/nrp[g]));
            ft6[g][0] *= ft6[g][0];
            se = (sqrt(fabs((ft6[g][1]-(ft6[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
            fprintf(fpw3,"%.3f\t",se);
            fprintf(fpw3,"%.3f\t",(ft7[g][0]/nrp[g]));
            ft7[g][0] *= ft7[g][0];
            se = (sqrt(fabs((ft7[g][1]-(ft7[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
            fprintf(fpw3,"%.3f\t",se);
            fprintf(fpw3,"%.3f\t",(ft8[g][0]/nrp[g]));
            ft8[g][0] *= ft8[g][0];
            se = (sqrt(fabs((ft8[g][1]-(ft8[g][0]/nrp[g]))/(nrp[g]-1))))/sqrt(nrp[g]);
            fprintf(fpw3,"%.3f\n",se);                        
          }
      }
    return;
  }

/* calculates mean and standard error for persistence; prints to file fpw1 */

void persist(int runs,int rcnt,int gen)
  {
    float se;

    fprintf(fpw1,"Number of Runs Persisting to %d Generations: %d\n",gen,rcnt);    
    fprintf(fpw1,"Average Persistence: %.1f\t",(per1/runs));
    if(runs>1)
      {
        per1 *= per1;
        se = (sqrt(fabs((per2-(per1/runs))/(runs-1))))/sqrt(runs);
        fprintf(fpw1,"Standard Error: %.3f\n\n",se);   
      }
    else
      fprintf(fpw1,"Standard error: N/A\n\n");  
    return;
  }

/* prints genetic and fitness data for individual plants to file fpw4 */

void plant_data(int g,int sloc,int nloc)
  {
    int i,id,xc,yc;
 
    for(xc=0;xc<LEN;xc++)
      {
        for(yc=0;yc<LEN;yc++)
          {
            if(pop[xc][yc].safe==1)
              { 
                fprintf(fpw4,"%4d %2d %2d %3d",g,xc,yc,pop[xc][yc].age);
                if(pop[xc][yc].age>0)
                  {
                    for(i=0;i<GENES;i++)
                      {
                        id_val(xc,yc,i,&id,sloc,nloc);
                        fprintf(fpw4," %3d %3d %3d",id,pop[xc][yc].gtype[i][0],
                                pop[xc][yc].gtype[i][1]);
                        printf(" %3d %3d %3d\n",id,pop[xc][yc].gtype[i][0],
                                pop[xc][yc].gtype[i][1]); 
                      }
                  }
                else if(pop[xc][yc].age==0)           
                  for(i=0;i<GENES;i++)
                    fprintf(fpw4," %3d %3d %3d",0,0,0);
                fprintf(fpw4,"\n");      
              }
          }  
      }
    return;
  }

/* calculates genotype id values for each locus */

void id_val(int xc,int yc,int i,int *id,int sloc,int nloc)
  {
    int j,sum=0,loc;

    (i==0)?(loc=sloc):(loc=nloc);
    if(pop[xc][yc].gtype[i][0]<=pop[xc][yc].gtype[i][1])
      {
        for(j=0;j<pop[xc][yc].gtype[i][0];j++)
          sum += j;
        *id = pop[xc][yc].gtype[i][1]+(loc*(pop[xc][yc].gtype[i][0]-1))-sum;
      }
    else if(pop[xc][yc].gtype[i][0]>pop[xc][yc].gtype[i][1])
      {
        for(j=0;j<pop[xc][yc].gtype[i][1];j++)
          sum += j;
        *id = pop[xc][yc].gtype[i][0]+(loc*(pop[xc][yc].gtype[i][1]-1))-sum;
      }
    return;
  }
