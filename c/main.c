//gcc -fopenmp main.c -lgsl -lgslcblas -lm -o CVW.out
// valgrind --leak-check=yes --error-limit=no --track-origins=yes --log-file=CVWvalgrind.log ./CVW_dbg.out &

#ifdef _MPI_USE
#include <mpi.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "utils.h"
#include <math.h>
#include <time.h>
#include <nlopt.h>


#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_interp.h>


// No solve to debug the optimzation routine
int nosolve = 0;


// declare some parameters that will be global scope
int const JJ = 4;     // number of occupations
int const NS = 2;     // number of specific skill types
int const NZ = 9;     // 9 number of occ match-quality types
int const NE = 8;     // 8 number of firm epsilon match-quality types
int const NP = 5;     // 5 number of occupation-specific productivies
int const NA = 2;     // 2 number of aggregate productivities

int NN ,NUN;

int verbose = 3;
int print_lev = 3;
int opt_print = 0;

int maxiter = 1000; //1000
double vftol = 1e-4;
double rhotightening = 1.;
double caltol = 1e-3;
double smth_flows = 0.25;

// THESE are set in main, possibly by command line
int cal_now, cf_now,jac_now;


int const TT      = 12*13;    // periods per simulation path
int const burnin  = 12*3;     // number of periods to throw away each time
int        TTT ;
int const Npaths  = 120;//120      // number of simulation paths to draw
int const Nsim    = 500; // 500

int const Npwave  = 4;
int const Anan    = 1;
int const Nqtls   = 5; // number of quantiles that will use to compare distributions
double     qtlgrid[] = {0.1,0.25,0.5,0.75,0.9};
double     edgeqtls[] = {0.025, 0.05, 0.5, 0.95,0.975};
int        Nparams = 16;
int        Ntargets= 8;

int        Ncluster = 3;
int        Npar_cluster[4] ={6,0,0,0}; // first the flows parameters, then the dist params then the cyclcical parameters
int        Ntgt_cluster[4] ={6,0,0,0};

int nflows=0; // number of net flows across occupations, will be useful for allocating parameters, etc
int ntgtavgflow = 9; //overall flows, like J2J rate, finding rate, separation rate etc.


int eps_2emg = 1; //should we use a double-exponentially modified gaussian, or just normal
int const nstarts = 1;  // how many starts per node

double beta	= 0.997;		// discount factor
double b 	= 0.; 		// unemployment benefit
double wage_lev = 1.;        // will be a shifter so the average wage is >0
double occ_wlevs[4] = {0.,-0.2657589,-0.4975667,-0.2189076}; // wage levels for each occupation
double occ_size_dat[] = {0.2636133, 0.3117827, 0.1493095, 0.2752945};

double urt_avg = .055;     // average separation rate

// specifics for the calibration routine.
FILE * calhist;
char calhi_f[] = "calhistX.csv";
FILE * parhist;
char parhi_f[] = "paramhistX.csv";
FILE * par_opt;
char parop_f[] = "param_optX.csv";

char exper_f[] = "noXX"; // the label for all the output that depends on a counter-factual.
double * solver_state;
int solX =0;

char * parnames_clu0[] = {"alpha0","lambdaUS0","lambdaUM0","lambdaES","lambdaEM","delta","zloss","alphaE1","alphaU1",
                                      "alphaf0","alphaf1","alphaf2","alphaf3"};
                          //"alpha01","alpha02","alpha03","alpha12","alpha13","alpha23"};
char * parnames_clu1[] = {"update_z","scale_z","shape_z","var_pe","autop","var_a","autoa","gdfather",
						  "stwupdt","var_eps","ltl_eps","rtl_eps"};
char * parnames_clu2[] = {"Delta_A","lamEM_A","lamES_A","lamUS_A","lamUM_A","zloss_A"
                          ,"z_Acoef","z_Amag","eps_Acoef","eps_Amag","occ1","occ2","occ3"};
char * parnames_cluall[] = {"alpha0","lambdaUS0","lambdaUM0","lambdaES","lambdaEM","delta","zloss","alphaE1","alphaU1",
                                        "alphaf0","alphaf1","alphaf2","alphaf3",
                            //"alpha01","alpha02","alpha03","alpha12","alpha13","alpha23",
                            "update_z","scale_z","shape_z","var_pe","autop","var_a","autoa","gdfather","stwupdt","var_eps","ltl_eps","rtl_eps",
                            "Delta_A","lamEM_A","lamES_A","lamUS_A","lamUM_A","zloss_A",
                            "z_Acoef","z_Amag","eps_Acoef","eps_Amag","occ1","occ2","occ3"};
char ** parnames[] = {parnames_clu0,parnames_clu1,parnames_clu2,parnames_cluall};

char * tgtnames_clu0[]   = {"J2J_err","fnd_err","sep_err","swEE_err","swU_err","swSt_err","dur_ratio","nrmflowE","nrmflowU",
                            "flow0","flow1","flow2","flow3"};
                            //"flow01","flow02","flow03","flow12","flow13","flow23"};

char * tgtnames_clu1;
char * tgtnames_clu2;
char * tgtnames_cluall;

char ** tgtnames;

struct cal_params{
	int cluster_hr, rank;
	double gdfthr;      // probability of a J2J without a choice
	double lambdaEM0, lambdaES0, lambdaUS0, lambdaUM0; //probability of ability to make a move
	double alpha0; 	    // scale of alpha function
	double alphaU1;		// concavity of alpha function
	double alphaE1;		// concavity of alpha function
	double **alpha_nf; //scale on net flow of switching from l to d
	double kappa;		// cost of switching
	double autoa;		// persistence of aggregate shock
	double autop;		// persistence of occ-specific shock
	double var_ae;		// innovations to aggregate shock
	double var_pe;		// innovations to occ-specific shock
	double var_ze;		// innovations to match-quality shock
	double update_z;    // re-draw match quality shock
	double scale_z;     // scale parameter of match-quality
	double shape_z;     // shape parameter of match-quality
	double scale_eps,shape_eps ; //scale and shape of firm-match quality
    double eps_Acoef,z_Acoef; // coefficient for distribution convolution to make more skewed in recession
    double eps_Amag,z_Amag;   // coefficient for magnitude of the mixed in convolution to make more skewed in recession

	double var_eps;     // variance of epsilon
	double lshape_eps;  // left-tail skewness of epsilon (exponential parameter)
	double rshape_eps;  // right-tail skewness of epsilon
	double wage_curve;  // curviness of utility function over wages
    double delta_Acoef;
	double zloss_Acoef;   // elasticity wrt zloss prob
    double delta_avg;     // average separation rate
	double lambdaUS_Acoef,lambdaUM_Acoef,lambdaES_Acoef,lambdaEM_Acoef; //How the move chance changes over the cycle
    double zloss;
	double stwupdate; // update rate for wages of stayers
    double * xparopt, *xsolopt;

    double * param_lb; double * param_ub;  // the upper and lower bounds for all parameters
    double ** cluster_lb; double** cluster_ub;

	gsl_vector * AloadP; //loading on A for each P
	gsl_vector * Plev;
	gsl_vector * Alev;
	gsl_vector * xSlev;
	gsl_vector * zlev;
	gsl_vector * epslev;

	gsl_matrix **Ptrans; // markov transition matrix for P, for each J
	gsl_matrix * Atrans;
	gsl_matrix * xStrans;
//	gsl_matrix * ztrans[NA]; // ztrans, one for each A
	gsl_vector **zprob; // distribution from which to draw when z first drawn
	gsl_vector **epsprob; // iid, so just a vector of probabilities
	gsl_vector * jprob;
};

struct valfuns{
	gsl_matrix * WE; // NN X JJ - value function for employed workers
	gsl_matrix * RE; // NN X JJ - value of search for employed workers
	gsl_matrix * WU; // NUN X JJ - value function for unemployed workers
	gsl_matrix * RU; // NUN X JJ - value of search for unemployed workers
	gsl_matrix * WEdist; // distance between two WE's
};

struct polfuns{
	gsl_matrix * mE;    // NN X JJ - to switch, employed
	gsl_matrix * mU;    // NUN X JJ - to switch, unemployed
	gsl_matrix ** sE;    // NN X JJ X JJ  - where to search, employed
	gsl_matrix ** sU;    // NUN X JJ X JJ - where to search, unemployed

};

struct shocks{
	gsl_matrix ** xGsel;
	gsl_matrix ** xSsel;
	gsl_matrix ** zsel;
	gsl_matrix ** epssel;
	gsl_matrix ** jsel;
	gsl_matrix ** dsel;
	gsl_matrix ** msel;
	gsl_matrix ** lambdaUSsel;
    gsl_matrix ** lambdaUMsel;
	gsl_matrix ** lambdaEMsel;
	gsl_matrix ** lambdaESsel;
	gsl_matrix ** Psel;
	gsl_vector ** Asel;


};

struct hists{
	gsl_matrix_int ** uhist;
	gsl_matrix     ** whist;
	gsl_matrix_int**  xShist;
    gsl_matrix_int**  epshist;
    gsl_matrix_int**  zhist;

    gsl_matrix_int ** jhist;
	gsl_vector_int ** Ahist;
	gsl_matrix_int ** Phist;
	gsl_matrix_int ** J2Jhist;
	gsl_matrix     *** sijhist;

};

struct stats{
	double swProb_U, swProb_st, swProb_EE;
	double doubleswE, doubleswU;
	double J2Jprob;
	double unrate;
	double findrate;
	double seprate;
	double udur_nosw;
	double udur_sw;
	double corrEE_wgocc, corrEUE_wgocc;
	double occ_ten;

	double seprt_ratio;
	double fndrt_ratio;
	double EE_ratio;
	double swProb_EE_ratio; // cyclical EE - ratio of expansion over recession
	double swProb_U_ratio;

	double var_wg; //variance of wages
	double var_wg_rec[2]; //variance of wages in expansions and recessions
	double * MVsw_qtls_ratio;
	double * MVns_qtls_ratio;

	double *UEsw_qtls; // will be Nqtls long
	double *UEns_qtls;
	double *EEsw_qtls;
	double *EEns_qtls;
	double *EUsw_qtls;
	double *EUns_qtls;
	double *stsw_qtls;
	double *stns_qtls;

	double *all_qtls, *edge_qtls;
	double ** all_qtls_rec, **edge_qtls_rec;

	double *occ_netflow; //  will have JJ ! / (JJ-2)!2! points.
	double *occ_netflowU; //  will have JJ ! / (JJ-2)!2! points.
	double *occ_netflowE; //  will have JJ ! / (JJ-2)!2! points.

	double *occ_margflow; //will have JJ points

	double nrmflows,nrmflowsE,nrmflowsU;
	double varGflowsE,varGflowsU;

};

void allocate_pars( struct cal_params * par );
void free_pars( struct cal_params * par);
void allocate_mats(struct valfuns * vf, struct polfuns * pf, struct hists * ht, struct shocks * sk);
void alloc_valfuns(struct valfuns *vf );
void alloc_hists( struct hists *ht );
void alloc_shocks(struct shocks * sk);
void memcpy_shocks(struct shocks * sk_dest , struct shocks * sk_orig);

void alloc_qtls( struct stats *st );
void free_qtls( struct stats *st );

void free_mats(struct valfuns * vf, struct polfuns * pf, struct hists *ht, struct shocks * sk);
void free_valfuns(struct valfuns *vf);
void free_hists( struct hists *ht);
void free_shocks(struct shocks * sk);

void init_pf( struct polfuns *pf ,struct cal_params * par);
void init_vf( struct valfuns *vf ,struct cal_params * par);
void memcpy_pf(struct polfuns *pf_dest, struct polfuns * pf_orig );
void memcpy_vf(struct valfuns *vf_dest, struct valfuns * vf_orig);
void convprobs(gsl_vector* probs, gsl_vector * lev,double coef, double magn);

int draw_shocks(struct shocks * sk);
// solve, simulate and compute summary statistics
int sol_dyn( struct cal_params * par, struct valfuns * vf, struct polfuns * pf, struct shocks * sk);
int sim( struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk );
int sum_stats(   struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk, struct stats *st  );

double param_dist( double * x, struct cal_params *par, int Npar, double * err_vec , int Nerr);
void compute_derivs(double ** cal_Jacobian, struct cal_params * par0, int printJac );

void fix_wagedists(struct stats *dat, struct stats *mod, struct hists * ht );

void shock_cf(  struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk ,struct stats *st);
void ion_permute( int nshocks, int ** idx_ion );

void set_dat( struct stats * );
void set_params( double * x, int n, struct cal_params * par,int ci);
void set_xpspace( double *x, struct cal_params *par );
void print_params(double *x, int n, struct cal_params * par);
void print_targets(char *name, struct stats * st, struct stats * dat);
void read_params(char* name, struct cal_params *par);
int w_qtls( double* vec, int stride, int len, double * qtlgrid, double * qtls_out );


// wrapper for nlopt algorithms
double f_wrapper_nlopt(unsigned n, const double * x, double * grad, void * par);

// interface for dfbols
void dfovec_iface_(double * f, double * x, int * n, int* mv);
// global parameters structure:
struct cal_params * glb_par;


// some helper functions to reduce typing
int static inline ggi_get( gsl_matrix_int * mat, int ri ,int ci){
	return gsl_matrix_int_get( mat, (size_t) ri, (size_t) ci );
}
double static inline gg_get( gsl_matrix * mat, int ri ,int ci){
	return gsl_matrix_get( mat, (size_t) ri, (size_t) ci );
}
void static inline ggi_set( gsl_matrix_int * mat, int ri ,int ci, int val){
	gsl_matrix_int_set( mat, (size_t) ri, (size_t) ci , val);
}
void static inline gg_set( gsl_matrix * mat, int ri ,int ci, double val){
	gsl_matrix_set( mat, (size_t) ri, (size_t) ci , val);
}

int main(int argc,char *argv[] ) {

	int i,ii,ji,j, ci,ai;
	int success, rank,nnodes;

	struct cal_params par;
	struct valfuns vf;
	struct polfuns pf;
	struct shocks sk;
	struct hists ht;

	struct stats st;
    
	// actually take this from argc
	if(argc >1)
		cal_now = atoi(argv[1]);
	else
		cal_now = 0;

	if(argc>2)
		cf_now = atoi(argv[2]);
	else
		cf_now = 0;
	if(argc>3)
	    jac_now = atoi(argv[3]);
	else
	    jac_now =0;

	sprintf(exper_f,"%s","");
	nflows = fact_int(JJ)/(fact_int(2)*fact_int(JJ-2));

	Npar_cluster[0] = ntgtavgflow + JJ; // was +nflows
	Npar_cluster[1] = 10;
	if(eps_2emg==1) Npar_cluster[1] +=2;
	Npar_cluster[2] = 10+JJ-1; //cyclical parameters
	Nparams = 0;
	for(i=0;i<Ncluster;i++)Nparams += Npar_cluster[i];

	Ntgt_cluster[0] = ntgtavgflow + JJ; // was + nflows
	Ntgt_cluster[1] = Nqtls*8; // 2 occwg and wage correlations. Plus, wage change distributions
	Ntgt_cluster[2] = 5+Nqtls*2; // cyclicality stuff
	Ntargets =0;
	for(i=0;i<Ncluster;i++) Ntargets += Ntgt_cluster[i];

	Npar_cluster[Ncluster] = 0;
	for(ci=0;ci<Ncluster;ci++) Npar_cluster[Ncluster] +=Npar_cluster[ci];
	Ntgt_cluster[Ncluster] = 0;
	for(ci=0;ci<Ncluster;ci++) Ntgt_cluster[Ncluster] +=Ntgt_cluster[ci];


	if(verbose>1){
		printf("Program version June 23d, 2020\n");

		printf("Solving for 3 clusters, sized %d,%d,%d\n", Npar_cluster[0],Npar_cluster[1],Npar_cluster[2]);
		printf("Targeting 3 clusters, sized %d,%d,%d\n", Ntgt_cluster[0],Ntgt_cluster[1],Ntgt_cluster[2]);
#ifdef _MPI_USE
        printf("Using MPI \n");
#endif
#ifdef _DFBOLS_USE
        printf("Using DFBOLS \n");
#endif
	}
	alloc_qtls(&st);

	allocate_mats(&vf,&pf,&ht,&sk);
	allocate_pars(&par);

	init_pf(&pf,&par);
	init_vf(&vf,&par);

    for(i=0;i<JJ;i++){
		gsl_matrix_set_all(par.Ptrans[i], 0.05/( (double)NP-1. ));
		for(ii=0;ii<NP;ii++)gg_set(par.Ptrans[i],ii,ii, 0.95);
		gg_set(par.Ptrans[i],0,NP-1, 0.);gg_set(par.Ptrans[i],NP-1,0, 0.);

		int ri,ci;
		for(ri=0;ri<NP;ri++){
			double rsum =0.;
			for(ci=0;ci<NP;ci++)rsum += gg_get(par.Ptrans[i],ri,ci);
			gsl_vector_view Pr = gsl_matrix_row( par.Ptrans[i],ri);
			gsl_vector_scale(& Pr.vector,1./rsum);
		}
	}


	// calibrate specific human capital growth for 1 year of employer tenure and 2 years of occupational tenure
	gsl_matrix_set_zero(par.xStrans);
	for(ii=0;ii<NS;ii++) gg_set(par.xStrans,ii,ii,1. - 1./12.);//two years to get tenured
	for(ii=0;ii<NS-1;ii++) gg_set(par.xStrans,ii,ii+1, 1./12.);//two years to get tenured
	gg_set(par.xStrans,NS-1,NS-1,1.);


	// just to initialize
	// will over-write this
	gsl_matrix_set_all(par.Atrans, 0.025/( (double)NA-1. ));
	for(ii=0;ii<NA;ii++) gg_set(par.Atrans,ii,ii,0.975);

	par.scale_z = 1.;
	par.shape_z = 1.;
	par.update_z = 0.05;
	par.zloss  = 0.025;
	par.zloss_Acoef =0.;


	par.autoa = 0.95;par.var_ae = 0.01*0.01;
	if(  NA>2)
		rouwenhorst(par.autoa,pow(par.var_ae,0.5),par.Atrans,par.Alev);
	else{
		double probRec = .21; double probExp=(1.-probRec);
		double qAtrans = (probExp/probRec+par.autoa)/(probExp/probRec+1);
		double pAtrans = par.autoa + 1. - qAtrans;
		double zAlev   = pow( par.var_ae ,0.5);
		gg_set(par.Atrans,0,0,   pAtrans);gg_set(par.Atrans,0,1,1.-pAtrans);
		gg_set(par.Atrans,1,0,1.-qAtrans);gg_set(par.Atrans,1,1,   qAtrans);
		gsl_vector_set(par.Alev,0,-zAlev);gsl_vector_set(par.Alev,1,zAlev);
	}
	gsl_vector * ergAprob = gsl_vector_calloc(NA);
	ergod_dist(par.Atrans,ergAprob);
	gsl_vector_free(ergAprob);
	par.autop = 0.95; par.var_pe = 0.02*0.2;
	rouwenhorst(par.autop,pow(par.var_pe,0.5),par.Ptrans[0],par.Plev);
	for(i=1;i<JJ;i++){
		gsl_matrix_memcpy(par.Ptrans[i],par.Ptrans[0]);
	}
	double zlev1[NZ]; double zprob1[NZ];
	success = disc_Weibull( zprob1,zlev1,NZ, 0., par.scale_z, par.shape_z );
	for(ai=0;ai<NA;ai++){
		for(ii=0;ii<NZ;ii++){
			gsl_vector_set( par.zlev,ii,zlev1[ii] );
			gsl_vector_set( par.zprob[ai],ii,zprob1[ii] );
		}
        gsl_vector_set_all(par.epsprob[ai], 1./(double) NE);
    }
	gsl_vector_set_all(par.AloadP,1.0);


	ai = 0;int pi = 0;int gi = 0;int si = 0;int zi = 0;	int thi = 0; ji=0;

	par.alphaE1 = 0.5;
	par.alphaU1 = 0.5;
    par.alpha0  = 0.05*pow((double)(JJ-1),-par.alphaE1);
    par.lambdaUS0  = 0.2 / 0.5;
    par.lambdaUM0  = 0.2 / 0.5;
    par.lambdaES0 = 0.01;
	par.lambdaEM0 = 0.8;
	par.kappa     = 0.00 ;
	par.gdfthr    = 0.5 ;
	par.wage_curve= 0.0 ;
	par.delta_avg = 0.01;
	par.delta_Acoef = 0.;
	par.lambdaUS_Acoef = 0.0;
    par.lambdaUM_Acoef = 0.0;
    par.lambdaES_Acoef = 0.0;
	par.lambdaEM_Acoef = 0.0;
	par.eps_Acoef = -0.25;
    par.eps_Amag = 0.0;
    par.z_Acoef = -0.25;
    par.z_Amag = 0.0;

	par.var_eps = 0.1;
	par.lshape_eps = 1.0;par.rshape_eps = 1.0;

	par.scale_eps=1.; par.shape_eps = 1.;

    double wage_lev0 = exp(par.AloadP->data[ji] * par.Alev->data[ai] +
                               par.Plev->data[pi] +
                               par.epslev->data[thi] +
                               par.zlev->data[zi] +
                               par.xSlev->data[si]);
    wage_lev = 1.;


    // parameter space:
	// alpha0 , lambdaUS, lambdaUM, lambdaES, lambdaEM, delta_avg, zloss_prob
	ii =0;
	par.param_lb[ii] = 0.001; par.param_ub[ii] = 0.10;ii++;

	par.param_lb[ii] = 0.001; par.param_ub[ii] = 0.80;ii++;
    par.param_lb[ii] = 0.001; par.param_ub[ii] = 0.80;ii++;
	par.param_lb[ii] = 0.001; par.param_ub[ii] = 0.05;ii++;
	par.param_lb[ii] = 0.001; par.param_ub[ii] = 0.99;ii++;
	par.param_lb[ii] = 0.001; par.param_ub[ii] = 0.02;ii++;
	par.param_lb[ii] = 0.001; par.param_ub[ii] = 0.02;ii++;

	//, alphaE_1, alphaU_1
	par.param_lb[ii] = 0.010; par.param_ub[ii] = 0.90;ii++; //potentially these are not <1
	par.param_lb[ii] = 0.010; par.param_ub[ii] = 0.90;ii++; //potentially these are not <1

	// alpha_nf matrix
	for(i=0;i<JJ;i++){ //was for(i=0;i<nflows;i++)
		par.param_lb[i+ii] = -0.50; par.param_ub[i+ii] = 0.50;}

	if(Ncluster>1){
		ii = Npar_cluster[0];

		// update_z, scale_z, shape_z,
		par.param_lb[ii] = 0.001; par.param_ub[ii] = 0.10;  ii++;
		par.param_lb[ii] = 0.500; par.param_ub[ii] = 4.000; ii++;
		par.param_lb[ii] = 0.500; par.param_ub[ii] = 7.000; ii++;
		//var_pe, autop, var_ae,autoa, gdfather, stwupdate
		par.param_lb[ii]= 0.001;  par.param_ub[ii]= 0.010;  ii++;
		par.param_lb[ii]= 0.500;  par.param_ub[ii]= 0.999;  ii++;
		par.param_lb[ii]= 0.001;  par.param_ub[ii]= 0.010;  ii++;
		par.param_lb[ii]= 0.750;  par.param_ub[ii]= 0.999;  ii++;
		par.param_lb[ii]= 0.01;   par.param_ub[ii]= 0.99;   ii++;
		par.param_lb[ii]= 0.01;   par.param_ub[ii]= 0.50;   ii++;
		//var_eps, lshape_eps, rshape_eps -- larger values of epsilon mean shorter tails:
		par.param_lb[ii]= 0.01;  par.param_ub[ii]= 0.10; ii++; //std of 0.5 as upper limit
		if(eps_2emg == 1){
			par.param_lb[ii]= 0.100;  par.param_ub[ii]= 4.00; ii++;
			par.param_lb[ii]= 0.100;  par.param_ub[ii]= 4.00; ii++;
		}

	}
	if(Ncluster>2){
		ii = Npar_cluster[1]+Npar_cluster[0];
		//delta_Acoef, lambdaEM, lambdaES, lambdaUS, lambdaUM, zloss (can up-to double)
		par.param_lb[ii]=-1.999;   par.param_ub[ii] = 0.999; ii++;
		par.param_lb[ii]=-0.999;   par.param_ub[ii] = 2.999; ii++;
		par.param_lb[ii]=-0.999;   par.param_ub[ii] = 1.999; ii++;
		par.param_lb[ii]=-0.999;   par.param_ub[ii] = 0.999; ii++;
        par.param_lb[ii]=-1.999;   par.param_ub[ii] = 1.999; ii++;
		par.param_lb[ii]=-1.999;   par.param_ub[ii] = 0.999; ii++;
        // z_Acoef, z_Amag, eps_Acoef, eps_Amag
        par.param_lb[ii]=-0.999;   par.param_ub[ii] = 0.499; ii++;
        par.param_lb[ii]= 0.   ;   par.param_ub[ii] = 0.499; ii++;
        par.param_lb[ii]=-0.999;   par.param_ub[ii] = 0.499; ii++;
        par.param_lb[ii]= 0.   ;   par.param_ub[ii] = 0.499; ii++;

        for(ji=0;ji<JJ;ji++){
			if(ji>0){
				par.param_lb[ii]=0.01;
				par.param_ub[ii] = 1.99;
				ii++;
			}

		}
	}


	par.cluster_lb = malloc(sizeof(double*)*(Ncluster+1));
	par.cluster_ub = malloc(sizeof(double*)*(Ncluster+1));
	par.xparopt = malloc(sizeof(double)*Nparams);
	par.xsolopt = malloc(sizeof(double)*Nparams);

	double * err = malloc(sizeof(double)*Ntargets);
	double *x0 = malloc(sizeof(double)*Nparams);
	double *x0_clu;double *err_clu;



	// branch for MPI here
    rank = 0;
    nnodes = 1;
#ifdef _MPI_USE
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status mpistatus;
#endif

	printf("\n On rank %d with %d nodes\n",rank,nnodes);

	par.rank = rank;
	double * x0starts,*caldist,*calx;
	double * x0starts_j = malloc(sizeof(double)*nstarts*Nparams);
	double * caldist_j = malloc(sizeof(double)*nstarts);
	double * calx_j    = malloc(sizeof(double)*nstarts*Nparams);

	if(rank==0){
		// rank 0 should scatter this out to everyone:
		x0starts = malloc(sizeof(double)*nnodes*nstarts*Nparams);
		caldist = malloc(sizeof(double)*nnodes*nstarts);
		calx    = malloc(sizeof(double)*nnodes*nstarts*Nparams);

		if(nnodes*nstarts>1){
			//pseudo-random starting point draws
            gsl_qrng *qrng = gsl_qrng_alloc(gsl_qrng_sobol,Nparams);
			for(ii=0;ii<10;ii++)
                gsl_qrng_get(qrng,x0); // burn it in with 10 draws

			for(i=0;i<nstarts*nnodes;i++){
				gsl_qrng_get(qrng,x0);
				for(ii=0;ii<Nparams;ii++) x0starts[i*Nparams+ii] = x0[ii];
			}
			gsl_qrng_free(qrng);
		}
		else{
			for(i=0;i<Nparams;i++) x0starts[i] = 0.5;
		}
	}

	if(cal_now==1){
	    printf("\n In the calibration \n");
		sprintf(calhi_f,"calhist%d.csv",rank);
		calhist = fopen(calhi_f,"w+");
		fprintf(calhist,"Value,");
		for(i=0;i<Ntgt_cluster[0];i++){
		    fprintf(calhist,"%s,",tgtnames_clu0[i]);
		}
		if(Ntargets>Ntgt_cluster[0]){

			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," stns%f, ",qtlgrid[ii]); //for(ii=0;ii<Nqtls;ii++) fprintf(calhist," stns%f, ",qtlgrid[ii]);
			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," stsw%f, ",qtlgrid[ii]);
			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," EEns%f, ",qtlgrid[ii]);
			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," EEsw%f, ",qtlgrid[ii]);
			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," EUns%f, ",qtlgrid[ii]);
			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," EUsw%f, ",qtlgrid[ii]);
			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," UEns%f, ",qtlgrid[ii]);
			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," UEsw%f, ",qtlgrid[ii]);
		}
		if(Ncluster>2){
			fprintf(calhist,"swPrUratio,swPrEEratio,frtratio,seprtratio,EEratio,");
			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," MVsw%f, ",edgeqtls[ii] );
			for(ii=0;ii<Nqtls;ii++) fprintf(calhist," MVns%f, ",edgeqtls[ii] );
		}
		for(i=0;i<Nparams;i++){
			fprintf(calhist,"%s,",parnames_cluall[i]);
		}

		fprintf(calhist,"\n");
		fclose(calhist);



		sprintf(parhi_f,"paramhist%d.csv",rank);
		parhist = fopen(parhi_f,"w+");
		for(i=0;i<Nparams;i++){
			fprintf(parhist,"%s,",parnames_cluall[i]);
		}
		fprintf(parhist,"\n");
		fclose(parhist);


        #ifdef _MPI_USE
		int nsend = Nparams*nstarts;
		MPI_Scatter( x0starts , nsend,MPI_DOUBLE,x0starts_j,nsend,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(rank==0){
            printf("\n From rank 0, sending x0starts, length %d\n", nsend);
        }else{
            printf("\n Received starts at rank %d rank. x0[0]: %f \n. ", rank, x0starts_j[0]);
        }
        #endif

        #ifndef _MPI_USE
        for(i=0;i<nstarts*Nparams;i++){
            x0starts_j[i] = x0starts[i];
        }
        // ad hoc adjustment:
        x0starts_j[ntgtavgflow] = .1; x0starts_j[ntgtavgflow+1] = .9;
        #endif

        double mindist = 1e6;

		i=0;
		while(i<nstarts){

			for(ii=0;ii<Nparams;ii++) x0[ii] = x0starts_j[i*Nparams+ii] ;
			for(ii=0;ii<Nparams;ii++)
				par.xparopt[ii] = x0[ii] * (par.param_ub[ii]-par.param_lb[ii])
			        +par.param_lb[ii];
			set_params(par.xparopt,Nparams,&par,Ncluster); //sets all of the parameters to the initial guess

			double dist;
			//if(verbose>0)
			//	printf("Beginning to evaluate a DFBOLS start point \n");

			glb_par = &par;
            // iterate on cluster loops
            int clu_iter = 0;
            int max_cluiter = 1;
            for(clu_iter=0;clu_iter<max_cluiter;clu_iter++) {

                // loop over clusters
                int ci0 =0<clu_iter && clu_iter<max_cluiter-1 ? 0        : Ncluster;
                int ciN =0<clu_iter && clu_iter<max_cluiter-1 ? Ncluster : Ncluster+1;
                for (ci = ci0; ci < ciN; ci++) {
                    printf("\n Doing first set of parameters, to check, x[0]=%f alpha0=%f: \n", x0[0], par.xparopt[0]);
                    set_params(par.xparopt, Nparams, &par, Ncluster); //sets all of the parameters to the initial guess
                    par.cluster_hr = ci;
                    if (cal_now == 0) {
                        ci = Ncluster;
                        par.cluster_hr = ci;
                    }
                    solver_state = malloc(sizeof(double) * (Npar_cluster[ci] + 1));
                    solver_state[0] = 1e4;
                    // allocations for DFBOLS
                    x0_clu = malloc(Npar_cluster[ci] * sizeof(double));
                    int npt = 2 * Npar_cluster[ci] + 1;
                    double rhobeg = pow(.5,(double)clu_iter/2.+1.);//  0.33*pow((double)(nnodes*nstarts),-1./(double)Npar_cluster[ci]);
                    rhobeg = rhobeg/ (double)nnodes;
                    double rhoend = clu_iter <max_cluiter-1 ? rhobeg/100 : rhobeg / 1000;
                    //if(ci==Ncluster) rhobeg /= 3;

                    //iterations choice in DFBOLS
                    int maxfun = 10 * (2 * Npar_cluster[ci] + 1); //approximate it 5 times
                    //int maxfun = ci<Ncluster ? 4*Npar_cluster[ci]+3  : 50*(2*Npar_cluster[ci]+1);
                    //int maxfun = 2;

                    double *wspace = calloc(
                            (npt + 5) * (npt + Npar_cluster[ci]) + 3 * Npar_cluster[ci] * (Npar_cluster[ci] + 5) / 2,
                            sizeof(double));

                    double *dfbols_lb, *dfbols_ub;
                    dfbols_lb = calloc(Npar_cluster[ci], sizeof(double));
                    dfbols_ub = calloc(Npar_cluster[ci], sizeof(double));
                    par.cluster_lb[ci] = malloc(Npar_cluster[ci] * sizeof(double));
                    par.cluster_ub[ci] = malloc(Npar_cluster[ci] * sizeof(double));
                    for (ii = 0; ii < Npar_cluster[ci]; ii++) {
                        dfbols_lb[ii] = 0.;
                        dfbols_ub[ii] = 1.;
                    }
                    int par_hr = 0;
                    // set up the bounds as a subset of the whole set of bounds
                    if (ci < Ncluster) {
                        for (ii = 0; ii < ci; ii++) { par_hr += Npar_cluster[ii]; }
                    }
                    for (ii = 0; ii < Npar_cluster[ci]; ii++) {
                        x0_clu[ii] = x0[ii + par_hr];
                        par.cluster_lb[ci][ii] = par.param_lb[par_hr + ii];
                        par.cluster_ub[ci][ii] = par.param_ub[par_hr + ii];
                    }
                    if (verbose > 1 && rank==0) {
                        printf("Bounds are: \n    ");
                        for (ii = 0; ii < Npar_cluster[ci]; ii++) { printf("%8s,", parnames[ci][ii]); }
                        printf("\n");
                        printf("lb: ");
                        for (ii = 0; ii < Npar_cluster[ci]; ii++) { printf("%8.6f,", par.cluster_lb[ci][ii]); }
                        printf("\n");
                        printf("ub: ");
                        for (ii = 0; ii < Npar_cluster[ci]; ii++) { printf("%8.6f,", par.cluster_ub[ci][ii]); }
                        printf("\n");
                    }

                    int dfbols_printlev = print_lev > 3 ? 3 : print_lev;
                    dfbols_printlev = print_lev < 0 ? 0 : print_lev;
                    if (cal_now == 1) {
                        #ifdef _DFBOLS_USE
                        bobyqa_h_(&(Npar_cluster[ci]),&npt,x0_clu,dfbols_lb,dfbols_ub,&rhobeg,&rhoend,&dfbols_printlev ,&maxfun,wspace,&(Ntgt_cluster[ci]));
                        #endif
                        #ifndef _DFBOLS_USE
                        err_clu = malloc(sizeof(double) * Ntgt_cluster[ci]);
                        dfovec_iface_(err_clu, x0_clu, &(Npar_cluster[ci]), &(Ntgt_cluster[ci]));
                        free(err_clu);
                        #endif

                        for (ii = 0; ii < Npar_cluster[ci]; ii++) par.xparopt[ii + par_hr] = solver_state[ii + 1];
                        for (ii = 0; ii < Npar_cluster[ci]; ii++) {
                            x0[ii + par_hr] = (par.xparopt[ii + par_hr] - par.cluster_lb[ci][ii]) /
                                              (par.cluster_ub[ci][ii] - par.cluster_lb[ci][ii]);
                        }
                    } else {
                        ci = Ncluster;
                        par.cluster_hr = ci;
                        for (ii = 0; ii < Nparams; ii++)
                            par.xparopt[ii] = x0[ii] * (par.param_ub[ii] - par.param_lb[ii])
                                              + par.param_lb[ii];

                    }
                    int print_lev_old = print_lev;
                    print_lev = 2;
                    sprintf(exper_f, "");
                    opt_print = 1;
                    par.cluster_hr = Ncluster;
                    dist = param_dist(par.xparopt, &par, Nparams, err, Ntargets);
                    print_lev = print_lev_old;

                    if (verbose > 1) {
                        printf("error is %f at vector:  (", dist);
                        for (ii = 0; ii < Ntgt_cluster[ci] - 1; ii++)
                            printf("%f,", err[ii]);
                        printf("%f)\n", err[Ntgt_cluster[ci] - 1]);
                        printf("evaluated at (");
                        for (ii = 0; ii < Npar_cluster[ci] - 1; ii++)
                            printf("%f,", par.xparopt[par_hr + ii]);
                        printf("%f)\n", par.xparopt[par_hr + Npar_cluster[ci] - 1]);
                    }
                    opt_print = 0;
                    free(wspace);
                    free(dfbols_lb);
                    free(dfbols_ub);
                    free(par.cluster_ub[ci]);
                    free(par.cluster_lb[ci]);
                    free(x0_clu);
                }//ci=0:Ncluster

            } // loop over all of the clusters several times
			free(solver_state);
			if(verbose>0)
				printf("Evaluated a DFBOLS start point \n");

			caldist_j[i] = dist;
			for(ii=0;ii<Nparams;ii++) calx_j[i*Nparams + ii] = x0[ii];

	        #ifdef _MPI_USE
			int ngather = nstarts*Nparams;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gather(caldist_j,nstarts,MPI_DOUBLE,caldist,nstarts,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Gather(calx_j,ngather ,MPI_DOUBLE,calx,ngather ,MPI_DOUBLE,0,MPI_COMM_WORLD);
            #endif

			#ifndef _MPI_USE
			for(ii=0;ii<nstarts;ii++) caldist[ii] = caldist_j[ii];
			for(ii=0;ii<nstarts*Nparams;ii++) calx[ii] = calx_j[ii];
	        #endif
			if(rank==0){
				if(verbose>0)
					printf("Looking for best point among starts \n");
				for(j=0;j<nnodes;j++){
					dist = caldist[j+i*nnodes];

					if (dist<mindist){
						mindist = dist;
						for(ii=0;ii<Nparams;ii++) par.xparopt[ii] = calx[j*Nparams+ ii]*(par.param_ub[ii]- par.param_lb[ii]) +par.param_lb[ii];
					}
				}
			}
			//check if less than tolerance and then broadcast a stop of the nstarts loop
	        #ifdef _MPI_USE
			MPI_Bcast(&mindist,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	        #endif
			if(dist<caltol)
				i=nstarts;
			else
				i++;
		}
		char parhi_old[15];
		snprintf(parhi_old,15,"%s",parhi_f); //strcpy(parhi_old,parhi_f);
		strcpy(parhi_f , "param_opt.csv");
        remove(parhi_f);
		print_params(calx, Nparams, &par);
        snprintf(parhi_f,15,"%s",parhi_old); //strcpy(parhi_f,parhi_old);
	}
	if(jac_now==1) {
        int xi;
        read_params("param_opt.csv", &par);
        int print_lev_old = print_lev;
        print_lev = 2;
        opt_print = 1;
        sprintf(exper_f,"");
        par.cluster_hr= Ncluster;
        //baseline, print where we are
        double dist = param_dist(par.xparopt, & par ,Nparams,err,Ntargets);
        print_lev = 0;
        opt_print= 0;
        printf("calculating derivatives\n");
        double **cal_Jacobian = malloc(sizeof(double *) * Nparams);
        for (xi = 0; xi < Nparams; xi++)
            cal_Jacobian[xi] = malloc(sizeof(double) * Ntargets);
        compute_derivs(cal_Jacobian, &par, 1);
        xi = 0;
        print_lev = print_lev_old;
        free(cal_Jacobian);
    }

    if( cf_now == 1 ){
		read_params("param_opt.csv", &par);
		int print_lev_old = print_lev;
		print_lev = 2;
		par.cluster_hr= Ncluster;
		set_params(par.xparopt,Nparams,&par,Ncluster );
		// solve once at the normal set of points
		init_pf(&pf,&par);
		init_vf(&vf,&par);
		//set the markov transition matrix for P
		rouwenhorst(par.autop,pow(par.var_pe,0.5),par.Ptrans[0],par.Plev);
		for(i=1;i<JJ;i++){
			gsl_matrix_memcpy(par.Ptrans[i],par.Ptrans[0]);
		}
		// set the markov transition matrix for A
		if(NA>2){
			rouwenhorst(par.autoa,pow(par.var_ae,0.5),par.Atrans,par.Alev);
		} else{
			double probRec = .21; double probExp=(1.-probRec);
			double qAtrans = (probExp/probRec+par.autoa)/(probExp/probRec+1);
			double pAtrans = par.autoa + 1. - qAtrans;
			double zAlev   = pow( par.var_ae ,0.5);
			gg_set(par.Atrans,0,0,   pAtrans);gg_set(par.Atrans,0,1,1.-pAtrans);
			gg_set(par.Atrans,1,0,1.-qAtrans);gg_set(par.Atrans,1,1,   qAtrans);
			gsl_vector_set(par.Alev,0,-zAlev);gsl_vector_set(par.Alev,1,zAlev);
		}

		// setup the z matrix
		double zlev1[NZ]; double zprob1[NZ];
		success = disc_Weibull( zprob1,zlev1,NZ, 0., par.scale_z, par.shape_z );
		for(ai=0;ai<NA;ai++){
			for(ii=0;ii<NZ;ii++){
				gsl_vector_set( par.zlev,ii,zlev1[ii] );
				gsl_vector_set( par.zprob[ai],ii,zprob1[ii] );
			}
		}
        convprobs(par.zprob[0],par.zlev,par.z_Acoef,par.z_Amag);

		//success  = disc_Weibull(par.epsprob->data, par.epslev->data, NE,0.,par.scale_eps,par.shape_eps);
		for(ai=0;ai<NA;ai++)
            success = disc_2emg(par.epsprob[ai]->data, par.epslev->data, (int) par.epsprob[ai]->size,
                                0., par.var_eps, par.lshape_eps, par.rshape_eps);
        convprobs(par.epsprob[0],par.epslev,par.eps_Acoef,par.eps_Amag);

		if(success > 0){
			printf(" Did not compute the distribution properly");
		}
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		for(i=0;i<NS;i++) gsl_vector_set(par.xSlev,i, .010984* (double)i/(double) (NS-1)); // really just 0 or 1--- matches 1% increase due to employer tenure in KM2009
		memcpy( par.jprob->data, occ_size_dat,sizeof(occ_size_dat));


		// ensure don't voluntarily quit, except for when z is the lowest value
		double w_hr;
		for(ji=0;ji<JJ;ji++){
			for (ii = 0; ii < NN; ii++) {
				int ai = ii / (NP * NS * NZ * NE);
				int pi = (ii - ai * NP * NS * NZ * NE) / (NS * NZ * NE);
				int si = (ii - ai * NP * NS * NZ * NE - pi * NS * NZ * NE) / (NZ * NE);
				int zi =
						(ii - ai * NP * NS * NZ * NE - pi * NS * NZ * NE - si * NZ * NE) / NE;
				int ti = ii - ai * NP * NS * NZ * NE - pi * NS * NZ * NE - si * NZ * NE -
				         zi * NE;
				w_hr = exp(par.AloadP->data[ji] * par.Alev->data[ai] +
						           par.Plev->data[pi] +
						           par.epslev->data[ti] +
						           par.zlev->data[zi] +
						           par.xSlev->data[si] +
						           occ_wlevs[ji]) + wage_lev;
				if( w_hr< b && w_hr>0 && zi>0) b = w_hr;
			}
		}

		success = draw_shocks(&sk);
		if(verbose>2 && success != 0) printf("Problem drawing shocks\n");
		success = sol_dyn(&par, &vf, &pf, &sk);
		if (verbose > 2 && success != 0) printf("Problem solving model\n");

		success = sim(&par, &vf, &pf, &ht, &sk);
		if (verbose > 2 && success != 0) printf("Problem simulating model\n");

		success = sum_stats(&par, &vf, &pf, &ht, &sk, &st);

		shock_cf( &par, &vf, &pf , &ht, &sk , &st);
		print_lev = print_lev_old;

	}

	free(x0starts_j);
	free(caldist_j);free(calx_j);

	free(x0starts);
	free(caldist);free(calx);

	free(err); free(x0);
	free(par.xparopt);

	free_mats(&vf,&pf,&ht,&sk);
	free_pars(&par);
	free_qtls(&st);

    return success;
}

int sol_dyn( struct cal_params * par, struct valfuns * vf, struct polfuns * pf, struct shocks * sk ){
	int ii,ji,viter;

	int success=0;

	struct valfuns vf0;
	alloc_valfuns(&vf0);

	double**wagevec;
	wagevec = malloc(sizeof(double*)*NN);
	for(ii=0;ii<NN;ii++)wagevec[ii] = malloc( sizeof(double)*JJ );

	for(ji=0;ji<JJ;ji++){
		for (ii = 0; ii < NN; ii++) {
			int ai = ii / (NP * NS * NZ * NE);
			int pi = (ii - ai * NP * NS * NZ * NE) / (NS * NZ * NE);
			int si = (ii - ai * NP * NS * NZ * NE - pi * NS * NZ * NE ) / (NZ * NE);
			int zi =
					(ii - ai * NP * NS * NZ * NE - pi * NS * NZ * NE - si * NZ * NE) / NE;
			int ti = ii - ai * NP * NS * NZ * NE - pi * NS * NZ * NE - si * NZ * NE -
			         zi * NE;

			wagevec[ii][ji] = exp(par->AloadP->data[ji] * par->Alev->data[ai] +
			                  par->Plev->data[pi]   +
			                  par->epslev->data[ti] +
			                  par->zlev->data[zi]   +
			                  par->xSlev->data[si]  +
			                  occ_wlevs[ji] ) + wage_lev;
		}
	}
    gsl_vector * ergPprobs = gsl_vector_calloc(NP);
    ergod_dist(par->Ptrans[0],ergPprobs);


	for(ji=0;ji<JJ;ji++) {
		for (ii = 0; ii < NN; ii++) {
			if( par->wage_curve<0.0001 && par->wage_curve>-0.0001 )
				gg_set(vf0.WE, ii, ji,
				       wagevec[ii][ji] / (1. - beta));
			else
				gg_set(vf0.WE, ii, ji,
			               pow(wagevec[ii][ji], 1.-par->wage_curve) / (1. - par->wage_curve) / (1. - beta));
		}
	}

	for(ii=0;ii<NUN;ii++){
		if( par->wage_curve< 0.0001 && par->wage_curve> -0.0001 )
			gg_set(vf0.WU,ii,0,
			       (1.-par->lambdaUS0)*b
			       /(1.-beta) +
			       par->lambdaUS0*(wagevec[ii*NE+NE/2][0])
			       /(1.-beta));
		else
			gg_set(vf0.WU,ii,0,
		               (1.-par->lambdaUS0)*pow(b, 1.-par->wage_curve) / (1. - par->wage_curve)
		               /(1.-beta) +
		               par->lambdaUS0*pow(wagevec[ii*NE+NE/2][0], 1.-par->wage_curve)/(1.-par->wage_curve)
		               /(1.-beta));
		for(ji=1;ji<JJ;ji++)  gg_set(vf0.WU,ii,ji, gg_get(vf0.WU,ii,0));
	}

	double maxdist;
	for(viter = 0;viter<maxiter;viter++){
        //figuring out value of employment
		for(ji=0;ji<JJ;ji++){
			#pragma omp parallel for private(ii) firstprivate(ji)
			for(ii=0;ii<NN;ii++){
				int ai = ii/(NP*NS*NZ*NE);
				int pi = (ii - ai*NP*NS*NZ*NE)/(NS*NZ*NE);
				int si = (ii - ai*NP*NS*NZ*NE - pi*NS*NZ*NE )/(NZ*NE) ;
				int zi = (ii - ai*NP*NS*NZ*NE - pi*NS*NZ*NE  - si*NZ*NE)/NE ;
				int ti = ii -  ai*NP*NS*NZ*NE - pi*NS*NZ*NE  - si*NZ*NE - zi*NE;

				int iU = ai*NP*NS*NZ + pi*NS*NZ + si*NZ +zi;

                double lambdaEMhr = par->lambdaEM0 * exp( par->lambdaEM_Acoef* (double)ai);
                double lambdaEShr = par->lambdaES0 * exp( par->lambdaES_Acoef*(double)ai);

				double delta_hr = par->delta_avg * exp( par->delta_Acoef * (double)ai );
                double zloss_hr = par->zloss     * exp( par->zloss_Acoef * (double)ai );

                if( zloss_hr + delta_hr >1 ){ //shouldn't get to this situation
                    double uprob =zloss_hr + delta_hr;
                    zloss_hr /= uprob;
                    delta_hr /= uprob;
                }
				//compute expectations over A, Pt
				double EAPWE = 0.;
				int aai, ppi;
				int ssi;
				for(ssi=0;ssi<NS;ssi++){
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++)  EAPWE += gsl_max(gg_get( vf0.WE,  aai*NP*NS*NZ*NE + ppi*NS*NZ*NE+ ssi*NZ*NE +zi*NE+ ti ,ji),
					                                          gg_get( vf0.WU,  aai*NP*NS*NZ    + ppi*NS*NZ   + ssi*NZ    +zi        ,ji))*
							gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi)*gg_get(par->xStrans,si,ssi);
				}
				}
				EAPWE = gsl_max(EAPWE, gg_get(vf0.WU,iU,ji));

				double EtTWE = 0.;
				int tti,zzi;
				double ttiprod = 0.;
				ssi =0; //if drawing new epsilon, then also getting ssi=0
				for(tti=ti;tti<NE;tti++) ttiprod+=gsl_vector_get(par->epsprob[ai],tti); //make sure tpobs sum to 1
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++){
						for(tti=ti;tti<NE;tti++)
							EtTWE +=  gsl_max(gg_get(vf0.WE, aai*NP*NS*NZ*NE + ppi*NS*NZ*NE+ ssi*NZ*NE +zi*NE+ tti ,ji) ,
							                  gg_get(vf0.WU, aai*NP*NS*NZ    + ppi*NS*NZ   + ssi*NZ    +zi         ,ji))*
									gsl_vector_get(par->epsprob[ai],tti)/ttiprod *gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
					}
				}
				EtTWE = gsl_max(EtTWE,gg_get(vf0.WU,iU,ji));

				double EtWE = 0.;
				ssi =0; //if drawing new epsilon, then also getting ssi=0
				for(aai=0;aai<NA;aai++) {
					for (ppi = 0; ppi < NP; ppi++) {
						for (tti = 0; tti < NE; tti++)
							EtWE += gsl_max(gg_get(vf0.WE, aai*NP*NS*NZ*NE + ppi*NS*NZ*NE + ssi*NZ*NE + zi*NE + tti, ji) ,
									        gg_get(vf0.WU, aai*NP*NS*NZ    + ppi*NS*NZ    + ssi*NZ    + zi         , ji))*
							        par->epsprob[ai]->data[tti]*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
					}
				}
				EtWE = gsl_max(EtWE,gg_get(vf0.WU,iU,ji));
                int jji;
                double EzWU = 0.; //value if switch occ and lose job
                for(jji=0;jji<JJ;jji++) {
                    for (aai = 0; aai < NA; aai++) {
                        for (ppi = 0; ppi < NP; ppi++) {
                            for (zzi = 0; zzi < NZ; zzi++)
                                EzWU += gg_get(vf0.WU,
                                               aai * NP * NS * NZ + ppi * NS * NZ + zzi, jji) * occ_size_dat[jji] *
                                        gsl_vector_get(par->zprob[ai], zzi) * gg_get(par->Atrans, ai, aai) * gsl_vector_get(ergPprobs, ppi);
                        }
                    }
                }



				double REhr = -par->kappa;
				double *EzWE =malloc(JJ*sizeof(double)); //value if switch occ
                double *EztWE=malloc(JJ*sizeof(double)); //value if switch occ and job

                double dirjE[JJ];
				ssi =0; //if drawing new occupation, then also getting ssi=0

				for(jji=0;jji<JJ;jji++){
					EzWE[jji]  = 0.;
					EztWE[jji] = 0.;

					for(aai=0;aai<NA;aai++){
						for(ppi=0;ppi<NP;ppi++) {
							for (zzi = 0; zzi < NZ; zzi++)
								EzWE[jji] += gsl_max(
								        gg_get(vf0.WE,aai*NP*NS*NZ*NE + ppi*NS*NZ*NE +zzi*NE + ti, jji) ,EzWU)*
								        gsl_vector_get(par->zprob[ai], zzi)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[jji],pi,ppi);
						}
					}
					//may get zloss shock
					EzWE[jji] = zloss_hr*EzWU + (1. - zloss_hr)*EzWE[jji];

					for(aai=0;aai<NA;aai++){
						for(ppi=0;ppi<NP;ppi++) {
							for(tti=0;tti<NE;tti++){
								for(zzi=0;zzi<NZ;zzi++)
									EztWE[jji] += gsl_max(
									        gg_get(vf0.WE,aai*NP*NS*NZ*NE + ppi*NS*NZ*NE+ zzi*NE + tti,jji) ,EzWU)*
											gsl_vector_get(par->zprob[ai], zzi)*(par->epsprob[ai]->data[tti])*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
								// should this be EztWE =max(EzWE,EztWE)
							}
						}
					}
                 //   EztWE[jji] = zloss_hr*EzWU + (1. - zloss_hr)*EztWE[jji];
					if(jji!=ji){
						// constructing RE
                        dirjE[jji] =
								(lambdaEMhr*EztWE[jji] +
								(1.- lambdaEMhr)*(par->stwupdate*EzWE[jji] +(1.-par->stwupdate)*EAPWE  ) );
                        // REhr is needed for whether we want to switch or not. That should have optimal search directions in it
						REhr += dirjE[jji]*gsl_min(par->alpha0*
                                 exp(par->alpha_nf[ji][jji]*par->alphaE1)*
                                 pow( gg_get( pf->sE[jji],ii,ji),1.-par->alphaE1),1.);
					}else{ dirjE[jji] = 0.; }
				}
				double totalphaS = 0;
				for(jji=0;jji<JJ;jji++){
					if(jji!=ji) totalphaS += gsl_min(par->alpha0*
					        exp(par->alpha_nf[ji][jji]*par->alphaE1)*
							pow(gg_get(pf->sE[jji],ii,ji),1.-par->alphaE1),1.);
				}
				REhr +=   gsl_max(1. - totalphaS,0. )* EAPWE;
				if(gsl_finite(REhr)==0){
					//printf("Uhoh. Bad REhr");
				}
				gg_set( vf->RE,ii,ji,REhr);
				double EWnom=(lambdaEShr * (par->gdfthr * EtWE + (1. - par->gdfthr) * EtTWE) +
                              (1. - lambdaEShr) *// potentially update z:
                              ((1.-par->update_z)*EAPWE + par->update_z*EzWE[ji]));

				double mhr;
                mhr = exp((REhr-EWnom)/rhotightening)/
                             ( exp((REhr-EWnom)/rhotightening)+ 1.);
				if( isinf(mhr) | isnan(mhr) ){
		    		mhr = REhr >= EAPWE ? 1. : 0. ;
				}

				gg_set( pf->mE,ii,ji, mhr );

				//set search direction for next iteration:
				double sEjiDenom = 0.;
				double sEnorm = 0.;//-par->kappa+ lambdaEMhr*EztWE[ji] + (1.-lambdaEMhr)*EzWE[ji]- EAPWE;
				for(jji=0;jji<JJ;jji++) sEnorm = jji != ji ? sEnorm + (dirjE[jji]-par->kappa)/( (double)JJ -1) : sEnorm;
				sEnorm = sEnorm>0 && gsl_finite(sEnorm)==1 ? sEnorm : 1.;
				for(jji=0;jji<JJ;jji++){
					double dirreturn = dirjE[jji]-par->kappa; //-par->kappa+ lambdaEMhr*EztWE[jji] + (1.-lambdaEMhr)*EzWE[jji]- EAPWE;
					dirreturn = dirreturn < 0. ? 0. : dirreturn;
					if(jji!=ji) sEjiDenom += exp(par->alpha_nf[ji][jji])*
					        pow(dirreturn/sEnorm , 1./par->alphaE1);
				}
				if(sEjiDenom >0){
					for(jji=0;jji<JJ;jji++){
						if(jji!=ji){
							double dirreturn =  dirjE[jji] -par->kappa ; //-par->kappa+ lambdaEMhr*EztWE[jji] + (1.-lambdaEMhr)*EzWE[jji]- EAPWE;
							dirreturn = dirreturn< 0 ? 0. : dirreturn;
							gg_set(pf->sE[jji] ,ii,ji,
                                   exp(par->alpha_nf[ji][jji])*
                                   pow(dirreturn /sEnorm, 1./par->alphaE1)
                                   /sEjiDenom);
						}else{ // this is kind of redundant because it should have initialized to 0
							gg_set(pf->sE[jji] ,ii,ji, 0.);
						}
					}
				}else{
					for(jji=0;jji<JJ;jji++){
						if(jji!=ji){ gg_set(pf->sE[jji],ii,ji,1./(double)(JJ-1) );
						}else
							gg_set(pf->sE[jji],ii,ji,0. );
					}
				}
				// make sure it adds to 1
				sEnorm = 0.;
				for(jji=0;jji<JJ;jji++)
					sEnorm += gg_get(pf->sE[jji],ii,ji);
				for(jji=0;jji<JJ;jji++)
					gg_set(pf->sE[jji],ii,ji, gg_get(pf->sE[jji],ii,ji)/sEnorm);


				// update the value function
				double WEhr = 0.;
				if(par->wage_curve<0.0001 && par->wage_curve>-0.0001) {
                    WEhr = (wagevec[ii][ji]) +
                           beta * delta_hr * gg_get(vf0.WU, iU, ji) + beta * zloss_hr * EzWU +
                           beta * (1. - delta_hr - zloss_hr) * (
                                   gg_get(pf->mE, ii, ji) * gg_get(vf->RE, ii, ji) +
                                   (1. - gg_get(pf->mE, ii, ji)) *EWnom );
                }else {
                    WEhr = pow(wagevec[ii][ji], 1. - par->wage_curve) / (1. - par->wage_curve) +
                           beta * delta_hr * gg_get(vf0.WU, iU, ji) + beta * zloss_hr * EzWU +
                           beta * (1. - delta_hr - zloss_hr) * (
                                   gg_get(pf->mE, ii, ji) * gg_get(vf->RE, ii, ji) +
                                   (1. - gg_get(pf->mE, ii, ji)) * EWnom);
				}
				gg_set( vf->WE, ii,ji,WEhr);

				gg_set( vf-> WEdist, ii,ji, gg_get(vf->WE,ii,ji) - gg_get(vf0.WE,ii,ji) );

				free(EztWE);free(EzWE);
			} // OMP loop over state ii
		}
		gsl_matrix_memcpy(vf0.WE,vf->WE);

		// now do the value of unemployment
		for(ji=0;ji<JJ;ji++){
			#pragma omp parallel for private(ii) firstprivate(ji)
			for(ii=0;ii<NUN;ii++) {
				int ai =  ii/(NP*NS*NZ);
				int pi = (ii - ai*NP*NS*NZ)/(NS*NZ);
				int si = (ii - ai*NP*NS*NZ - pi*NS*NZ)/(NZ) ;
				int zi =  ii - ai*NP*NS*NZ - pi*NS*NZ -  si*NZ ;

				int jji,zzi,tti,aai,ppi;

				double lambdaUShr = par->lambdaUS0*exp( par->lambdaUS_Acoef* (double)ai);
                double lambdaUMhr = par->lambdaUM0*exp( par->lambdaUM_Acoef* (double)ai);
                double zloss_hr   = par->zloss    *exp( par->zloss_Acoef   * (double)ai); //relevant because might get a really bad z draw

                double EzDelWU = 0.; //value if in the forced-switch position
                for(jji=0;jji<JJ;jji++) {
                    for (aai = 0; aai < NA; aai++) {
                        for (ppi = 0; ppi < NP; ppi++) {
                            for (zzi = 0; zzi < NZ; zzi++)
                                EzDelWU += gg_get(vf0.WU,
                                               aai * NP * NS * NZ + ppi * NS * NZ + zzi, jji) * occ_size_dat[jji] *
                                        gsl_vector_get(par->zprob[ai], zzi) * gg_get(par->Atrans, ai, aai) * gsl_vector_get(ergPprobs, ppi);
                        }
                    }
                }


                // staying in the same occupation through U
				double EAPWU = 0.;
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++) {
						EAPWU += gg_get(vf0.WU, aai*NP*NS*NZ+ppi*NS*NZ + si*NZ + zi,ji)*
								gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
					}
				}
				double EzWU[JJ];
				double RUhr = -par->kappa;
				double dirjU[JJ];
				// if one changes occupation through U
				for(jji=0;jji<JJ;jji++){
					EzWU[jji]=0;
					double EtWE_jji=0;

					for(aai=0;aai<NA;aai++){
						for(ppi=0;ppi<NP;ppi++) {
							for(zzi=0;zzi<NZ;zzi++){
								for(tti=0;tti<NE;tti++){
									// new draws on A,P, epsilon. have si=0
									int iihr = aai*NP*NS*NZ*NE+ppi*NS*NZ*NE+0*NZ*NE+zzi*NE+tti;
									int iUhr = aai*NP*NS*NZ   +ppi*NS*NZ   +0*NZ   +zzi ;
									EtWE_jji += gsl_max(gg_get(vf0.WE, iihr, jji), gg_get(vf0.WU,iUhr,jji) )*
											gsl_vector_get(par->zprob[ai],zzi)*gsl_vector_get(par->epsprob[ai],tti)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[jji], pi, ppi);}
							}
						}
					}
					EtWE_jji = zloss_hr*EzDelWU + (1.-zloss_hr)*EtWE_jji;
					for(aai=0;aai<NA;aai++){
						for(ppi=0;ppi<NP;ppi++) {
							for( zzi=0;zzi<NZ;zzi++ )
								EzWU[jji] += gg_get(vf0.WU, aai*NP*NS*NZ+ppi*NS*NZ + zzi,jji)*
										gsl_vector_get(par->zprob[ai],zzi)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[jji],pi,ppi);
						}
					}
                    EzWU[jji] =zloss_hr*EzDelWU + (1.-zloss_hr)*EzWU[jji];
					if(jji!=ji){
						dirjU[jji] =
								(EzWU[jji]*(1.-lambdaUMhr)+lambdaUMhr*EtWE_jji);
						RUhr += dirjU[jji]*gsl_min(
                                par-> alpha0*exp(par->alpha_nf[ji][jji]*par->alphaU1)
                                *pow(gg_get(pf->sU[jji],ii,ji),1.-par->alphaU1 ),1.);
					}else{dirjU[jji] =0;} //just to not have uninitialized stuff
				}
				double totalalphaS = 0;
				for(jji=0;jji<JJ;jji++)
					if(jji!=ji) totalalphaS += gsl_min(par-> alpha0*
					        exp(par->alpha_nf[ji][jji] *par->alphaU1)
							*pow(gg_get(pf->sU[jji],ii,ji),1.-par->alphaU1),1.);
				RUhr += gsl_max(1.-totalalphaS,0. )*EAPWU;
				gg_set(vf->RU, ii,ji, RUhr);

				double EtWE=0;
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++) {
						for(tti=0;tti<NE;tti++){
							int iihr = aai*NP*NS*NZ*NE+ppi*NS*NZ*NE+ si*NZ*NE+zi*NE+tti;
							EtWE += gsl_max(gg_get(vf0.WE, iihr, ji), vf0.WU->data[ii*vf0.WU->tda+ji] )*
							        gsl_vector_get(par->epsprob[ai],tti)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji], pi, ppi);}
					}
				}
				double EWnom = (1.-lambdaUShr)*EAPWU+ lambdaUShr*gsl_max(EtWE,EAPWU)  ;
				double mhr = exp((RUhr-EWnom)/rhotightening)/
				               (exp((RUhr-EWnom)/rhotightening)+1.);
				if( isinf(mhr) | isnan(mhr) ){
					mhr = RUhr > (1.-lambdaUShr)*EAPWU+ lambdaUShr*gsl_max(EtWE,EAPWU) ?
							1. : 0. ;
				}
				gg_set(pf->mU,ii,ji, mhr );
				//search dir for next iterations
				double sUdenom = 0.;
				double sUnorm = 0.; //-par->kappa + EzWU[ji]- EAPWU;
				for(jji=0;jji<JJ;jji++) sUnorm = jji != ji ? sUnorm + (dirjU[jji]-par->kappa)/( (double)JJ -1) : sUnorm;
				sUnorm = sUnorm >0 && gsl_finite(sUnorm)==1 ? sUnorm : 1.;
				for(jji=0;jji<JJ;jji++){
					double dirreturn = dirjU[jji]-par->kappa; //-par->kappa + EzWU[jji]- EAPWU;
					dirreturn =dirreturn<0 ? 0. : dirreturn;
					if(jji!=ji) sUdenom += exp(par->alpha_nf[ji][jji])*
					        pow(dirreturn/sUnorm ,1./par->alphaU1);
				}
				if(sUdenom > 0 ) {
					for (jji = 0; jji < JJ; jji++) {
						if (jji != ji) {
							double dirreturn = dirjU[jji]-par->kappa; //-par->kappa + EzWU[jji] - EAPWU;
							dirreturn = dirreturn < 0 ? 0. : dirreturn;
							gg_set(pf->sU[jji], ii, ji,
							        exp(par->alpha_nf[ji][jji])*
							        pow(dirreturn/sUnorm, 1. / par->alphaU1)
							        / sUdenom);
						} else {//this is kind of redundant, because it should have initialized to 0
							gg_set(pf->sU[jji], ii, ji, 0.);
						}
					}
				}else{
					for(jji=0;jji<JJ;jji++){
						if(jji!=ji){ gg_set(pf->sU[jji],ii,ji,1./(double)(JJ-1) );
						}else
							gg_set(pf->sU[jji],ii,ji,0.);
					}
				}
				sUnorm = 0.;
				for(jji=0;jji<JJ;jji++)
					sUnorm += gg_get(pf->sU[jji],ii,ji);
				for(jji=0;jji<JJ;jji++)
					gg_set(pf->sU[jji],ii,ji,gg_get(pf->sU[jji],ii,ji)/sUnorm);
				double WUhr;
				if( par->wage_curve<0.00001 & par->wage_curve>-0.0001 )
					WUhr = b +
					       beta*( pf->mU->data[ii*pf->mU->tda + ji]*vf->RU->data[ii*vf->RU->tda+ji] +
					              (1.-pf->mU->data[ii*pf->mU->tda + ji])*( (1.-lambdaUShr)*EAPWU + lambdaUShr*EtWE ) ) ;

				else
					WUhr=pow(b,1.-par->wage_curve)/(1.-par->wage_curve) +
											beta*( pf->mU->data[ii*pf->mU->tda + ji]*vf->RU->data[ii*vf->RU->tda+ji] +
											(1.-pf->mU->data[ii*pf->mU->tda + ji])*( (1.-lambdaUShr)*EAPWU + lambdaUShr*EtWE ) ) ;

				WUhr = (1.-par->update_z)*WUhr + par->update_z * EzWU[ji];

				gg_set( vf->WU, ii, ji, WUhr);
			}
		}

		gsl_matrix_memcpy(vf0.WU,vf->WU);
		maxdist = gsl_max( gsl_matrix_max(vf->WEdist),-gsl_matrix_min(vf->WEdist));
 		if(viter % 200 == 0 && verbose >1 )  printf("Max distance is %f on iteration %d \n", maxdist,viter);
		double tolhr = pow(viter+1,-1)*1e-2+vftol;
		if( maxdist < tolhr ){
			success = 0;
			printf("Solved after %d iterations with error %f \n", viter,maxdist);
			break;
		}else{

			success = (int) maxdist;
		}

	}
	printf("Exited w/o convergence after %d iterations with error %f \n", viter,maxdist);

	for(ii=0;ii<NN;ii++)free(wagevec[ii]);
	free(wagevec);
	free_valfuns(&vf0);
    gsl_vector_free(ergPprobs);

    return success;

}

void ion_permute( int nshocks, int ** idx_ion ){
    int ion1, ion2,idx=0; // will loop over scenarios
    int i;

    int nscenario = 0;
    for(i=1;i<nshocks;i++){
        nscenario += fact_int(nshocks)/fact_int(i)/fact_int(nshocks - i);
    }

    gsl_permutation * shockpermutations = gsl_permutation_alloc (nshocks);
    gsl_permutation_init (shockpermutations);
    int npermut = nshocks;
    int non,permut_idx;
    for(i=1;i<nshocks;i++) npermut *= i;

    int plist[npermut][nshocks];
    permut_idx =0;
    do{
        for(i=0;i<nshocks;i++) plist[permut_idx][i] = (int) gsl_permutation_get(shockpermutations,i);
        permut_idx ++;
    }while (gsl_permutation_next(shockpermutations) == GSL_SUCCESS);
    gsl_permutation_free(shockpermutations);

    idx =0;int ionL=0;
    for(non=1;non<nshocks;non++) {
        int **plist_n = malloc(sizeof(int*)*npermut );
        for (ion2 = 0; ion2 < npermut; ion2++)plist_n[ion2] = malloc(sizeof(int) * non);
        for (ion1 = 0; ion1 < npermut; ion1++) {
            for (ion2 = 0; ion2 < non; ion2++) plist_n[ion1][ion2] = plist[ion1][ion2];
            gsl_sort_int(plist_n[ion1], 1, non);
        }
        ion1 = 0;
        for (ion2 = 0; ion2 < nshocks; ion2++)idx_ion[idx][ion2]=0; //init to 0
        for (ion2 = 0; ion2 < non; ion2++)
            idx_ion[idx][plist_n[ion1][ion2]] = 1;  // set to 1 some of them
        for (ion1 = 1; ion1 < npermut; ion1++) {
            int pnew = 1;
            // check if this plist_n is different from any we've already done in plist_n
            for(ionL=0;ionL<ion1;ionL++){
                int psame =0;
                for (ion2 = 0; ion2 < non; ion2++)
                    psame = (plist_n[ion1][ion2] == plist_n[ionL][ion2]) ? psame+1 : psame;
                pnew = (psame==non) ? 0 : pnew;
            }
            if (pnew == 1) {
                idx++;
                for (ion2 = 0; ion2 < nshocks; ion2++)idx_ion[idx][ion2]=0;
                for (ion2 = 0; ion2 < non; ion2++)
                    idx_ion[idx][plist_n[ion1][ion2]] = 1;
            }
        }
        idx++; // to setup for the next non
        for(ion2=0;ion2<npermut;ion2++)free(plist_n[ion2]);
        free(plist_n);
    }
    if(idx> nscenario ){
        printf("TOO MANY SCENARIOS!!! ");
    }

}

void shock_cf(struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk,struct stats * st ) {

	// simulate without z shocks, epsilon shocks, A shocks and P shocks

	struct shocks sk_noflow;
	struct hists ht_noflow;
	struct valfuns vf_noflow;
	struct polfuns pf_noflow;

	struct stats st_noflow; alloc_qtls(&st_noflow);
	int ll,i,it, wi;
	// the number of combinations of scenarios:
	// noz ye yp ya, noz noe yp ya, noz noe nop ya, noz noe nop na;; yz ne yp ya, yz ne np ya, yz ne np na;; yz ye np ya, yz ye np na;; yz ye yp na
	int nshocks = 4;
	int nscenario = 0;
	for(i=1;i<nshocks;i++){
		nscenario += fact_int(nshocks)/fact_int(i)/fact_int(nshocks - i);
	}
	struct shocks  sk_noX[nscenario];
	struct hists   ht_noX[nscenario];
	struct valfuns vf_noX[nscenario];
	struct polfuns pf_noX[nscenario];
	struct stats   st_noX[nscenario];
	for(i=0;i<nscenario;i++)
		alloc_qtls( &st_noX[i] );

	int pl_0 = print_lev;
	int vbs_0 = verbose;
    int noff,ri,qi;
    FILE * scenariokey = fopen( "scenariokey.csv","w+" );
    fprintf(scenariokey,"Num, Z off, Eps off, P off, A off\n");
    int ion1, ion2,idx=0; // will loop over scenarios
    int **idx_ion;

	if(print_lev>1 && par->rank==0){
        printvec("zlev.csv", par->zlev,0);
        int ai;for(ai=0;ai<NA;ai++) printvec("zprobs.csv",par->zprob[ai],0);
        printvec("Alev.csv", par->Alev,0);
		printmat("Atrans.csv",par->Atrans,0);
		printmat("Ptrans.csv",par->Ptrans[0],0);
		printvec("Plev.csv", par->Plev,0);
        for(ai=0;ai<NA;ai++) printvec("epsprob.csv", par->epsprob[ai],0);
		printvec("epslev.csv", par->epslev,0);
	}
	/* make the edge quantiles for the baseline case:
	int idx_all=0;
	int idx_cyc[2];
	double *w_all = malloc(sizeof(double)*Npaths*Nsim*TT/Npwave);
	double **w_cyc = malloc(sizeof(double*)*2);
	for(wi=0;wi<2;wi++){
		w_cyc[wi]= malloc(sizeof(double)*Npaths*Nsim*TT/Npwave);
		idx_cyc[wi]=0;
	}

	for(ll=0;ll<Npaths;ll++){
		for(i=0;i<Nsim;i++){
			for(wi=0;wi<(TT/Npwave);wi++){
				it = wi*Npwave;
				// no need to also look at each month per wave b/c nothing changes
				int ri = gsl_vector_int_get(ht->Ahist[ll],it);
				// identify if we're in a recession or not
				double wlast=0.,wnext=0.;
				if( wi>3 && wi< TT/Npwave-3 && Anan==1  ) {
					int si = 0;
					for (si = 1; si < Npwave * 3 + 1; si++)
						wlast += gg_get(ht->whist[ll], i, it - si);
					for (si = 0; si < Npwave * 3; si++)
						wnext += gg_get(ht->whist[ll], i, it + si);
					if( wnext>0. && wlast>0. && ggi_get(ht->uhist[ll],i,it) ==0){
						double wdif = wnext - wlast;
						w_all[idx_all] = wdif;
						w_cyc[ri][idx_cyc[ri]] = wdif;
						idx_all++;
						idx_cyc[ri]++;
					}
				}
			}
		}
	}
	w_qtls(w_all,1,idx_all,edgeqtls,st->edge_qtls);
	for(wi=0;wi<2;wi++)
		w_qtls(w_cyc[wi],1,idx_cyc[wi],edgeqtls,st->edge_qtls_rec[wi]);
	free(w_all);
	for(wi=0;wi<2;wi++)free(w_cyc[wi]);
	free(w_cyc);
*/

	gsl_vector * zlev_0 = gsl_vector_calloc(NZ);
	gsl_vector_memcpy(zlev_0,par->zlev);
	double zloss0 = par->zloss;

	double delta0 = par->delta_avg;
	double gdfthr0 = par->gdfthr;

	gsl_vector * eps_0 = gsl_vector_calloc(NE);
	gsl_vector_memcpy(eps_0,par->epslev);

	gsl_vector * Plev_0 = gsl_vector_calloc(NP);
	gsl_vector_memcpy(Plev_0,par->Plev);
	gsl_vector * ergPprobs = gsl_vector_calloc(NP);
	ergod_dist(par->Ptrans[0],ergPprobs);

	gsl_vector * Alev_0 = gsl_vector_calloc(NA);
	gsl_vector_memcpy(Alev_0,par->Alev);
    gsl_vector * ergAprobs = gsl_vector_calloc(NA);
    ergod_dist(par->Atrans,ergAprobs);

/*

	idx_ion = malloc(sizeof(int*)*nscenario);
	for(idx=0;idx<nscenario;idx++)
		idx_ion[idx] = malloc(sizeof(int)*nshocks);

    ion_permute( nshocks, idx_ion );
	// ++++++++++++++++++++++++++++++++++++++++++++++++++
    //cycle through scenarios turning on and off:
    for(idx=0;idx<nscenario;idx++){

        sprintf(exper_f, "_%d",idx);
	    fprintf( scenariokey, "%d, ", idx );

        allocate_mats(&vf_noX[idx],&pf_noX[idx],&ht_noX[idx],&sk_noX[idx]);

        memcpy_shocks(&sk_noX[idx],sk);
        memcpy_pf(&pf_noX[idx],pf);
        memcpy_vf(&vf_noX[idx],vf);

        if(idx_ion[idx][0] ==1 ){
            // then we're turning off z

			double zmean = 0.;//gsl_stats_mean(par->zlev->data,par->zlev->stride,par->zlev->size);
            int zi;
            for(zi=1;zi<NZ;zi++) zmean+=gsl_vector_get(par->zlev,zi)/((double)NZ - 1. ) ;
            for(zi=0;zi<NZ;zi++) gsl_vector_set(par->zlev,zi,zmean + 0.0001*((double)zi-(double)NZ/2.)/(double)NZ);
	        printf("Setting z shocks to their mean %f\n",zmean);

	        par->zloss=0.;

            fprintf(scenariokey," 1,");
        }else
            fprintf(scenariokey," 0,");

        if( idx_ion[idx][1]==1 ){
            // turning off epsilon
            double epsmean = gsl_stats_mean(par->epslev->data,par->epslev->stride,par->epslev->size);
            int ei;
            for(ei=0;ei<NE;ei++)
                gsl_vector_set(par->epslev,ei,epsmean+ 0.0001*((double)ei-(double)NE/2. )/(double)NE);

            par->gdfthr =0.;
            par->delta_avg=0.;

	        printf("Setting epsilon shocks to their mean %f\n",epsmean);

	        fprintf(scenariokey," 1,");
        }else
            fprintf(scenariokey," 0,");

        if( idx_ion[idx][2]==1){
            //turning off P

            double Pmean = 0.; int ip,ji;
            for(ji=0;ji<JJ;ji++){
                for(ip=0;ip<NP;ip++)
                    Pmean += gsl_vector_get(Plev_0,ip)*gsl_vector_get(ergPprobs,ip)*occ_size_dat[ji];
            }
            for(ip=0;ip<NP;ip++)
                gsl_vector_set(par->Plev,ip, 0. + 0.00001* ( (double)ip-(double)NP/2. )/(double)NP );
			printf("setting P shocks to their mean %f \n", Pmean);
            fprintf(scenariokey," 1,");
        }else
            fprintf(scenariokey," 0,");
        if( idx_ion[idx][3]==1 ){
            //turning off A

            double Amean = 0.; int ia;
            for(ia=0;ia<NA;ia++)
                Amean += gsl_vector_get(Alev_0,ia)*gsl_vector_get(ergAprobs,ia);
            for(ia=0;ia<NA;ia++)
                gsl_vector_set(par->Alev, ia,  0. + 0.00001* ( (double)ia-(double)NA/2. )/(double)NA );
	        printf("setting A shocks to their mean %f \n", Amean);

            fprintf(scenariokey," 1 \n");
        }else
            fprintf(scenariokey," 0 \n");


        sol_dyn( par, &vf_noX[idx], &pf_noX[idx], &sk_noX[idx] );
        print_lev = 2; verbose =vbs_0;
        sim( par, &vf_noX[idx], &pf_noX[idx], &ht_noX[idx], &sk_noX[idx] );
        sum_stats( par, &vf_noX[idx],&pf_noX[idx],&ht_noX[idx],&sk_noX[idx], &st_noX[idx]);
        print_lev = pl_0;

	    par->zloss = zloss0;
	    par->gdfthr =gdfthr0;
	    par->delta_avg= delta0;

	    gsl_vector_memcpy(par->Plev, Plev_0);
        gsl_vector_memcpy(par->Alev,Alev_0);
        gsl_vector_memcpy(par->zlev,zlev_0);
        gsl_vector_memcpy(par->epslev, eps_0);
    }
	fclose(scenariokey);

	double varcontrib[nshocks*(nscenario+1)];
	for(i=0;i<(nshocks*(nscenario+1));i++) varcontrib[i]=0;
	double varcontrib_rec[2][nshocks*(nscenario+1)];
	for(i=0;i<(nshocks*(nscenario+1));i++) varcontrib_rec[0][i]=0;
	for(i=0;i<(nshocks*(nscenario+1));i++) varcontrib_rec[1][i]=0;
	double qtlcontrib[Nqtls][nshocks*(nscenario+1)];
	double qtlcontrib_rec[Nqtls*2][nshocks*(nscenario+1)];



	double qtl_avg_rec[Nqtls*2][nshocks];
	for(qi=0;qi<2*Nqtls;qi++){
		for(ion1=0;ion1<nshocks;ion1++)
			qtl_avg_rec[qi][ion1] =0.;
	}
	double qtl_avg[Nqtls][nshocks];
	for(qi=0;qi<(Nqtls);qi++){
		for(ion1=0;ion1<nshocks;ion1++)
		qtl_avg[qi][ion1] = 0.;}
	double nexpers[nshocks];
	for(i=0;i<nshocks;i++) nexpers[i]=0;
	//do the decomposition, find the buddy of each:
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	for( idx=0;idx<nscenario;idx++ ){
	// loop through each scenario, and get the variance contribution by looking at noff-1 in which that one is turned off.
		noff =0;
		for(ion1=0;ion1<nshocks;ion1++) noff+= idx_ion[idx][ion1]; // the number that are on in this scenario
		for(ion1=0;ion1<nshocks;ion1++){ // turn this one off
			if( idx_ion[idx][ion1] ==1 ){ //only turn it on if it's off
				nexpers[ion1] += 1.;
				double var_on=0.;
				double var_off = st_noX[idx].var_wg;
				double var_on_rec[2];
				for(ri=0;ri<2;ri++) var_on_rec[ri] = 0.;
				double var_off_rec[2];
				for(ri=0;ri<2;ri++) var_off_rec[ri] = st_noX[idx].var_wg_rec[ri];
				double qtl_on[Nqtls];
				double qtl_off[Nqtls];
				for(qi=0;qi<Nqtls;qi++) qtl_off[qi] = st_noX[idx].all_qtls[qi];
				double qtl_on_rec[2][Nqtls];
				for(ri=0;ri<2;ri++) var_on_rec[ri]= 0.;
				double qtl_off_rec[2][Nqtls];
				for(ri=0;ri<2;ri++){ for(qi=0;qi<Nqtls;qi++)qtl_off_rec[ri][qi] = st_noX[idx].all_qtls_rec[ri][qi];}
				// find the buddy where it's on:
				if(noff==1){
					var_on =st->var_wg;
					for(ri=0;ri<2;ri++){
						var_on_rec[ri] = st->var_wg_rec[ri];
						for(qi=0;qi<Nqtls;qi++)qtl_on_rec[ri][qi] = st->all_qtls_rec[ri][qi];
					}
					for(qi=0;qi<Nqtls;qi++) qtl_on[qi] = st->all_qtls[qi] ;
				}else{
					int idx_ion_1off[nshocks];
					for(ion2=0;ion2<nshocks;ion2++)	idx_ion_1off[ion2] = idx_ion[idx][ion2];
					idx_ion_1off[ion1] = 0;
					int matches=0, idx_match=0;
					//do{
					for(idx_match=0;idx_match<nscenario;idx_match++){
						matches=0;
						for(ion2=0;ion2<nshocks;ion2++)
							matches= idx_ion[idx_match][ion2] == idx_ion_1off[ion2]? matches+1 : matches;
						if(matches == nshocks){
							var_on = st_noX[idx_match].var_wg;
							for(ri=0;ri<2;ri++){
								var_on_rec[ri] = st_noX[idx_match].var_wg_rec[ri];
								for(qi=0;qi<Nqtls;qi++)qtl_on_rec[ri][qi] = st_noX[idx_match].all_qtls_rec[ri][qi];
							}
							for(qi=0;qi<Nqtls;qi++) qtl_on[qi] = st_noX[idx_match].all_qtls[qi];

							printf("match is %d for scenario %d,index %d \n", idx_match, idx, ion1);
						}
						if(idx_match >=nscenario){
							printf("couldn't find a match for the variance contribution at %d, noff=%d \n",ion1,noff);
							break;
						}
					}
				}
				varcontrib[idx*nshocks+ion1] = var_on-var_off;
				for(ri=0;ri<2;ri++)varcontrib_rec[ri][idx*nshocks+ion1] = var_on_rec[ri]-var_off_rec[ri];
				for(qi=0;qi<Nqtls;qi++) qtlcontrib[qi][idx*nshocks+ion1] = qtl_on[qi] - qtl_off[qi];
				for(ri=0;ri<2;ri++){for(qi=0;qi<Nqtls;qi++) qtlcontrib_rec[qi+ri*Nqtls][idx*nshocks+ion1] = qtl_on_rec[ri][qi] - qtl_off_rec[ri][qi];}
				for(qi=0;qi<Nqtls;qi++) qtl_avg[qi][ion1] += qtlcontrib[qi][idx*nshocks+ion1];
			}
		}
	}
	//also do the ones where we just turn one on:
	for( idx=0;idx<nscenario;idx++ ){
		noff =0;
		for(ion1=0;ion1<nshocks;ion1++) noff+= idx_ion[idx][ion1]; // the number that are on in this scenario
		if(noff==nshocks-1){
			// will compare to the baseline with all on
			ion1=0;
			for(ion2=0;ion2<nshocks;ion2++) ion1 = idx_ion[idx][ion2]==0 ? ion2 : ion1;
			nexpers[ion1]+=1;
			double var_off = 0.;
			double var_on  = st_noX[idx].var_wg;
			double qtl_on[Nqtls],qlt_off[Nqtls];
			for(qi=0;qi<Nqtls;qi++)  qtl_on[qi]= st_noX[idx].all_qtls[qi];

			double var_on_rec[2] = {st_noX[idx].var_wg_rec[0],st_noX[idx].var_wg_rec[1]};
			for(ri=0;ri<2;ri++) varcontrib_rec[ri][(nscenario)*nshocks+ion1] = var_on_rec[ri];

			qtlcontrib[qi][nscenario*nshocks+ion1] = qtl_on[qi]- 0.;
			for(qi=0;qi<Nqtls;qi++) qtl_avg[qi][ion1] += qtlcontrib[qi][nscenario*nshocks+ion1];
		}
	}

	printf("Had number of scenarios: z %f2, eps %f2, P %f2, A %f2 \n", nexpers[0],nexpers[1],nexpers[2],nexpers[3]);


	gsl_matrix_view  vcontrib_mat = gsl_matrix_view_array(varcontrib, nscenario+1,nshocks );
	printmat("var_contrib.csv",&(vcontrib_mat.matrix),0);
	gsl_matrix_view vcontrib_e_mat = gsl_matrix_view_array(varcontrib_rec[0], nscenario+1,nshocks );
	printmat("var_contrib_exp.csv",&(vcontrib_e_mat.matrix),0);
	gsl_matrix_view vcontrib_r_mat = gsl_matrix_view_array(varcontrib_rec[1], nscenario+1,nshocks );
	printmat("var_contrib_rec.csv",&(vcontrib_r_mat.matrix),0);

	gsl_matrix*  qcontrib_mat = gsl_matrix_calloc(Nqtls,nscenario);;
	for(qi=0;qi<Nqtls;qi++){
		for(ion1=0;ion1<nshocks;ion1++) gg_set(qcontrib_mat,qi,ion1, qtlcontrib[qi][ion1] ); }
	printmat("qtl_contrib.csv",qcontrib_mat,0);

    for(i=0;i<nscenario;i++){
		free_qtls(&st_noX[i]);
        free_mats(&vf_noX[i],&pf_noX[i],&ht_noX[i],&sk_noX[i]);
		free(idx_ion[i]);
	}
    for(i=0;i<nscenario;i++)
        free(idx_ion[i]);
	free(idx_ion);

    ////////////////////////////////////////////////////////
    /// turn off flows
*/
    printf("+++++++++++++++++++\n Going into experiments turning off flows");
	sprintf(exper_f,"noff") ;
	allocate_mats( &vf_noflow,&pf_noflow,&ht_noflow, &sk_noflow);

	memcpy_pf(&pf_noflow,pf);
	memcpy_vf(&vf_noflow, vf);
	memcpy_shocks(&sk_noflow,sk);
	//change some parameter values:
	double Aexp = 1.; //gsl_vector_get(par->Alev,1);

	double delA0 = par->delta_Acoef;
	double delm0 = par->delta_avg;
	par->delta_Acoef =0.;
	par->delta_avg = delm0*exp(delA0);

	double zlossA0 = par->zloss_Acoef;
	double zlossm0 = par->zloss;
	par->zloss_Acoef =0.;
	par->zloss = zlossm0*exp(zlossA0);

	double lamEMA0 = par->lambdaEM_Acoef;
	double lamEMm0 = par->lambdaEM0;
	par->lambdaEM0 = lamEMm0*exp(lamEMA0);
	par->lambdaEM_Acoef = 0.;

	double lamESA0 = par->lambdaES_Acoef;
	double lamESm0 = par->lambdaES0;
	par->lambdaES_Acoef = 0.;
	par->lambdaES0 = lamESm0*exp(lamESA0);

	double lamUSA0 = par->lambdaUS_Acoef;
	double lamUSm0 = par->lambdaUS0;
	par->lambdaUS0 = lamUSm0*exp(lamUSA0);
	par->lambdaUS_Acoef = 0.;

    double lamUMA0 = par->lambdaUM_Acoef;
    double lamUMm0 = par->lambdaUM0;
    par->lambdaUM0 = lamUMm0*exp(lamUMA0);
    par->lambdaUM_Acoef = 0.;

    par->zloss =0.;
	par->gdfthr =0.;

	print_lev=0;verbose =0;
	sol_dyn(par, &vf_noflow, &pf_noflow, &sk_noflow);
	//print_lev = 2;verbose = vbs_0;
	sim( par, &vf_noflow, &pf_noflow, &ht_noflow, &sk_noflow);
	sum_stats( par, &vf_noflow, &pf_noflow,&ht_noflow,&sk_noflow, &st_noflow);


	par->delta_Acoef =delA0;
	par->zloss_Acoef = zlossA0;
	par->lambdaEM_Acoef = lamEMA0;
	par->lambdaES_Acoef =lamESA0;
	par->lambdaUS_Acoef = lamUSA0;
    par->lambdaUM_Acoef = lamUMA0;

	par->delta_avg =delm0;
	par->zloss     = zlossm0;
	par->lambdaEM0 = lamEMm0;
	par->lambdaES0 =lamESm0;
	par->lambdaUS0 = lamUSm0;
    par->lambdaUM0 = lamUMm0;

	par->zloss =zloss0;
	par->gdfthr =gdfthr0;

	printf("The additional variance due to flows is: %f \n", st->var_wg - st_noflow.var_wg);
	printf("The additional variance due to flows is: %f in recession and %f in expansion \n",
	       st->var_wg_rec[1] - st_noflow.var_wg_rec[1],st->var_wg_rec[0] - st_noflow.var_wg_rec[0]);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// now do experiments with flows off:
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int print_old = print_lev;
    print_lev = 0;
    int ncycflows = 6; //each cyclical sensitivity
	nscenario = 0;
	for(i=1;i<ncycflows;i++){
		nscenario += fact_int(ncycflows)/fact_int(i)/fact_int(ncycflows - i);
	}
	struct shocks  sk_nofX[nscenario];
	struct hists   ht_nofX[nscenario];
	struct valfuns vf_nofX[nscenario];
	struct polfuns pf_nofX[nscenario];
	struct stats   st_nofX[nscenario];
	for(i=0;i<nscenario;i++)
		alloc_qtls( &st_nofX[i] );
	double edge_qtls_all[nscenario][Nqtls];
	double edge_qtls_cyc[2][nscenario][Nqtls];


	scenariokey = fopen( "fscenariokey.csv","w+" );
	fprintf(scenariokey,"Num, delta off, zloss off, lamUM off, lamUS off, lamES off, lamEM off \n");
	idx_ion = malloc(sizeof(int*)*nscenario);
	for(idx=0;idx<nscenario;idx++)
		idx_ion[idx] = malloc(sizeof(int)*ncycflows);

	ion_permute( ncycflows, idx_ion );

	for(idx=0;idx<nscenario;idx++){
		sprintf(exper_f,"fl%02d",idx);
		fprintf( scenariokey, "%02d, ", idx );

		allocate_mats(&vf_nofX[idx],&pf_nofX[idx],&ht_nofX[idx],&sk_nofX[idx]);

		memcpy_shocks(&sk_nofX[idx],sk);
		memcpy_pf(&pf_nofX[idx],pf);
		memcpy_vf(&vf_nofX[idx],vf);


		if(idx_ion[idx][0]==1){
			//then we're turning off delta_Acoef
            par->delta_Acoef =0.;
            par->delta_avg = delm0*exp(delA0);
            fprintf( scenariokey," 1,");
		}else
            fprintf( scenariokey," 0,");
		if(idx_ion[idx][1]==1){
			//then we're turning off zloss_Acoef
            par->zloss_Acoef =0.;
            par->zloss = zlossm0*exp(zlossA0);
            fprintf( scenariokey," 1,");
		} else
		    fprintf(scenariokey, " 0,");
        if(idx_ion[idx][2]==1){
            //then we're turning off lambdaUS_Acoef
            par->lambdaUS0 = lamUSm0*exp(lamUSA0);
            par->lambdaUS_Acoef = 0.;
            fprintf( scenariokey," 1,");
        } else
            fprintf(scenariokey, " 0,");
        if(idx_ion[idx][3]==1){
            //then we're turning off lambdaUM_Acoef
            par->lambdaUM0 = lamUMm0*exp(lamUMA0);
            par->lambdaUM_Acoef = 0.;
            fprintf( scenariokey," 1,");
        } else
            fprintf(scenariokey, " 0,");
        if(idx_ion[idx][4]==1){
            //then we're turning off lambdaES_Acoef
            par->lambdaES_Acoef =0.;
            par->lambdaES0 = lamESm0*exp(lamESA0);
            fprintf( scenariokey," 1,");
        } else
            fprintf(scenariokey, " 0,");
        if(idx_ion[idx][5]==1){
			//then we're turning off lambdaEM_Acoef
            par->lambdaEM0 = lamEMm0*exp(lamEMA0);
            par->lambdaEM_Acoef = 0.;
            fprintf( scenariokey," 1,");
		} else
            fprintf(scenariokey, " 0,");


        fprintf(scenariokey, " \n");
        sol_dyn( par, &vf_nofX[idx], &pf_nofX[idx], &sk_nofX[idx] );
        verbose =vbs_0;
        sim( par, &vf_nofX[idx], &pf_nofX[idx], &ht_nofX[idx], &sk_nofX[idx] );
        int cf_old = cf_now; cf_now=1; // pretty sure should already have cf_now=1
        sum_stats( par, &vf_nofX[idx],&pf_nofX[idx],&ht_nofX[idx],&sk_nofX[idx], &st_nofX[idx]);
        print_lev = pl_0; cf_now = cf_old;

/*
        int si;
        int idx_all=0;
        int idx_cyc[2];
        double *w_all = malloc(sizeof(double)*Npaths*Nsim*TT/Npwave);
        double **w_cyc = malloc(sizeof(double*)*2);
        for(ri=0;ri<2;ri++){
        	w_cyc[ri]= malloc(sizeof(double)*Npaths*Nsim*TT/Npwave);
        	idx_cyc[ri]=0;
        }

		for(ll=0;ll<Npaths;ll++){
            for(i=0;i<Nsim;i++){
                for(wi=0;wi<(TT/Npwave);wi++){
                    it = wi*Npwave;
                    // no need to also look at each month per wave b/c nothing changes
                    ri = gsl_vector_int_get(ht_nofX[idx].Ahist[ll],it);
                    // identify if we're in a recession or not
                    double wlast=0.,wnext=0.;

                    if( wi>3 && wi< TT/Npwave-3 && Anan==1  ) {
                        int si = 0;
                        for (si = 1; si < Npwave * 3 + 1; si++)
                            wlast += gg_get(ht_nofX[idx].whist[ll], i, it - si);
                        for (si = 0; si < Npwave * 3; si++)
                            wnext += gg_get(ht_nofX[idx].whist[ll], i, it + si);
                        if( wnext>0. && wlast>0. && ggi_get(ht_nofX[idx].uhist[ll],i,it) ==0){
                            double wdif = wnext - wlast;
                        	w_all[idx_all] = wdif;
	                        w_cyc[ri][idx_cyc[ri]] = wdif;

                        	idx_all++;
                            idx_cyc[ri]++;

                        }
                    }
                }
            }
        }
		w_qtls(w_all,1,idx_all,edgeqtls,edge_qtls_all[idx]);
		for(ri=0;ri<2;ri++)
			w_qtls(w_cyc[ri],1,idx_cyc[ri],edgeqtls,edge_qtls_cyc[ri][idx]);


        free(w_all);
        for(ri=0;ri<2;ri++)free(w_cyc[ri]);
        free(w_cyc);
*/
        for(qi=0;qi<Nqtls;qi++) edge_qtls_all[idx][qi] = st_nofX[idx].edge_qtls[qi];
        for(ri=0;ri<2;ri++) {
            for (qi = 0; qi < Nqtls; qi++) edge_qtls_cyc[ri][idx][qi] = st_nofX[idx].edge_qtls_rec[ri][qi];
        }
        free_mats( &vf_nofX[idx],&pf_nofX[idx],&ht_nofX[idx],&sk_nofX[idx] );

		// put it all back
        par->delta_Acoef =delA0;
        par->zloss_Acoef = zlossA0;
        par->lambdaUS_Acoef = lamUSA0;
        par->lambdaUM_Acoef = lamUMA0;
        par->lambdaES_Acoef =lamESA0;
        par->lambdaEM_Acoef = lamEMA0;

        par->delta_avg = delm0;
        par->zloss     = zlossm0;
        par->lambdaEM0 = lamEMm0;
        par->lambdaES0 =lamESm0;
        par->lambdaUS0 = lamUSm0;
        par->lambdaUM0 = lamUMm0;

        par->zloss =zloss0;
        par->gdfthr =gdfthr0;
	}

	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Do the decomposition, finding the buddy and averaging
	double flow_vcontrib[ncycflows*nscenario]; double flow_vcontrib_rec[2][ncycflows*nscenario];
	double flow_qtlcontrib_rec[2][ncycflows*nscenario][Nqtls];
	double flow_qtlcontrib[ncycflows*nscenario][Nqtls];
	double flow_qtl_avg[ncycflows][Nqtls];
	for(ion1=0;ion1<ncycflows;ion1++){
		for(qi=0;qi<Nqtls;qi++) flow_qtl_avg[ion1][qi]=0.;
	}


    double cycflow_qtl_rec[Nqtls*2][ncycflows];
    for(qi=0;qi<2*Nqtls;qi++){
        for(ion1=0;ion1<ncycflows;ion1++)
            cycflow_qtl_rec[qi][ion1] =0.;
    }
    double cycflow_qtl[Nqtls][ncycflows];
    for(qi=0;qi<(Nqtls);qi++){
        for(ion1=0;ion1<ncycflows;ion1++)
            cycflow_qtl[qi][ion1] = 0.;}
    double cycflow_expers[ncycflows];
    for(i=0;i<ncycflows;i++) cycflow_expers[i]=0;
    //do the decomposition, find the buddy of each:
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    FILE* fqtlcontrib = fopen("flows_qtlscontrib.csv","w+");
    FILE* fqtlcontrib_rec = fopen("flows_qtlscontrib_rec.csv","w+");
    fclose(fqtlcontrib);fclose(fqtlcontrib_rec);
    for( idx=0;idx<nscenario;idx++ ){

        FILE* fqtlcontrib = fopen("flows_qtlscontrib.csv","a+");
        FILE* fqtlcontrib_rec = fopen("flows_qtlscontrib_rec.csv","a+");

        // loop through each scenario, and get the variance contribution by looking at noff-1 in which that one is turned off.
        noff =0;
        for(ion1=0;ion1<ncycflows;ion1++) noff+= idx_ion[idx][ion1]; // the number that are on in this scenario
        if(noff==nshocks-1) continue;
        for(ion1=0;ion1<ncycflows;ion1++){ // flip this one from off to on
            if( idx_ion[idx][ion1] ==1 ){ //only turn it on if it's off
                fprintf(fqtlcontrib,"%d,",idx);
                fprintf(fqtlcontrib_rec,"%d,",idx);

                cycflow_expers[ion1] += 1.;
                double var_on=0.;
                double var_off = st_nofX[idx].var_wg;
                double var_on_rec[2];
                for(ri=0;ri<2;ri++) var_on_rec[ri] = 0.;
                double var_off_rec[2];
                for(ri=0;ri<2;ri++) var_off_rec[ri] = st_nofX[idx].var_wg_rec[ri];
                double qtl_on[Nqtls];
                double qtl_off[Nqtls];
                for(qi=0;qi<Nqtls;qi++) qtl_off[qi] = edge_qtls_all[idx][qi];
                double qtl_on_rec[2][Nqtls];
                for(ri=0;ri<2;ri++) var_on_rec[ri]= 0.;
                double qtl_off_rec[2][Nqtls];
                for(ri=0;ri<2;ri++){
                	for(qi=0;qi<Nqtls;qi++)qtl_off_rec[ri][qi] = edge_qtls_cyc[ri][idx][qi];
                }
                // find the buddy where it's on:
                if(noff==1){
                    var_on =st->var_wg;
                    for(ri=0;ri<2;ri++){
                        var_on_rec[ri] = st->var_wg_rec[ri];
                        for(qi=0;qi<Nqtls;qi++)qtl_on_rec[ri][qi] = st->edge_qtls_rec[ri][qi];
                    }
                    for(qi=0;qi<Nqtls;qi++) qtl_on[qi] = st->edge_qtls[qi] ;
                }else{
                    int idx_ion_1off[ncycflows];
                    for(ion2=0;ion2<ncycflows;ion2++)	idx_ion_1off[ion2] = idx_ion[idx][ion2];
                    idx_ion_1off[ion1] = 0;
                    int matches=0, idx_match=0;
                    //do{
                    for(idx_match=0;idx_match<nscenario;idx_match++){
                        matches=0;
                        for(ion2=0;ion2<ncycflows;ion2++)
                            matches= idx_ion[idx_match][ion2] == idx_ion_1off[ion2]? matches+1 : matches;
                        if(matches == ncycflows){
                            var_on = st_nofX[idx_match].var_wg;
                            for(ri=0;ri<2;ri++){
                                var_on_rec[ri] = st_nofX[idx_match].var_wg_rec[ri];
                                for(qi=0;qi<Nqtls;qi++)qtl_on_rec[ri][qi] = edge_qtls_cyc[ri][idx_match][qi];
                            }
                            for(qi=0;qi<Nqtls;qi++) qtl_on[qi] = edge_qtls_all[idx_match][qi];

                            printf("match is %d for scenario %d,index %d \n", idx_match, idx, ion1);
                        }
                    }
                }
                flow_vcontrib[idx*ncycflows+ion1] = var_on-var_off;
                for(ri=0;ri<2;ri++)flow_vcontrib_rec[ri][idx*ncycflows+ion1] = var_on_rec[ri]-var_off_rec[ri];
                for(qi=0;qi<Nqtls;qi++) flow_qtlcontrib[idx*ncycflows+ion1][qi] = qtl_on[qi] - qtl_off[qi];
                for(ri=0;ri<2;ri++){
                    for(qi=0;qi<Nqtls;qi++) flow_qtlcontrib_rec[ri][idx*ncycflows+ion1][qi] = (qtl_on[qi] - qtl_off[qi]);
                }
                for(qi=0;qi<Nqtls;qi++) flow_qtl_avg[ion1][qi] += flow_qtlcontrib[idx*ncycflows+ion1][qi];

                for(qi=0;qi<Nqtls;qi++) fprintf(fqtlcontrib,"%f, ",qtl_on[qi] - qtl_off[qi]);fprintf(fqtlcontrib,"\n");

                printf("Scenario %d,",idx); for(qi=0;qi<Nqtls;qi++) printf("%f, ",qtl_on[qi] - qtl_off[qi]); printf("\n");

                for(ri=0;ri<2;ri++) {
                    for (qi = 0; qi < Nqtls; qi++)
                        fprintf(fqtlcontrib_rec, "%f, ", qtl_on_rec[ri][qi] - qtl_off_rec[ri][qi]);
                }
                fprintf(fqtlcontrib_rec,"\n");

            }
        }
        fclose(fqtlcontrib);fclose(fqtlcontrib_rec);
    }
    /*also do the ones where we just turn one on:
    for( idx=0;idx<nscenario;idx++ ){
        noff =0;
        for(ion1=0;ion1<ncycflows;ion1++) noff+= idx_ion[idx][ion1]; // the number that are on in this scenario
        if(noff==nshocks-1){
            fprintf(fqtlcontrib,"%d,",idx);
            fprintf(fqtlcontrib_rec,"%d,",idx);
            // will compare to the baseline with all on
            ion1=0;
            for(ion2=0;ion2<nshocks;ion2++) ion1 = idx_ion[idx][ion2]==0 ? ion2 : ion1;
            cycflow_expers[ion1]+=1;
            double var_off = 0.;
            double var_on  = st_nofX[idx].var_wg;
            double qtl_on[Nqtls],qtl_off[Nqtls]; double qtl_on_rec[2][Nqtls],qtl_off_rec[2][Nqtls];
            for(qi=0;qi<Nqtls;qi++) qtl_off[qi] = edge_qtls_all[idx][qi];
            for(qi=0;qi<Nqtls;qi++)  qtl_on[qi] = st_nofX[idx].edge_qtls[qi];

            double var_on_rec[2] = {st_nofX[idx].var_wg_rec[0],st_nofX[idx].var_wg_rec[1]};
            for(ri=0;ri<2;ri++) {
                flow_vcontrib_rec[ri][(nscenario) * ncycflows + ion1] = var_on_rec[ri];
                for(qi=0;qi<Nqtls;qi++)qtl_on_rec[ri][qi]  = st->edge_qtls_rec[ri][qi];
                for(qi=0;qi<Nqtls;qi++)qtl_off_rec[ri][qi] =     edge_qtls_cyc[ri][idx][qi];
            }

            flow_qtlcontrib[qi][nscenario*ncycflows+ion1] = qtl_on[qi]- qtl_off[qi];
            for(ri=0;ri<2;ri++) flow_qtlcontrib_rec[ri][qi][nscenario*ncycflows+ion1] = qtl_on_rec[ri][qi] - qtl_off_rec[ri][qi];
	        for(qi=0;qi<Nqtls;qi++) flow_qtl_avg[ion1][qi] += flow_qtlcontrib[nscenario*ncycflows+ion1][qi];
	        for(qi=0;qi<Nqtls;qi++) fprintf(fqtlcontrib,"%f,",flow_qtlcontrib[qi][nscenario*ncycflows+ion1]);
            for(ri=0;ri<2;ri++) {
                for (qi = 0; qi < Nqtls; qi++) fprintf(fqtlcontrib_rec, "%f,",qtl_on_rec[ri][qi] - qtl_off_rec[ri][qi]);
            }
            fprintf(fqtlcontrib,"\n");
            fprintf(fqtlcontrib_rec,"\n");

        }
    }*/
    for(ion1=0;ion1<ncycflows;ion1++) {
        for(qi=0;qi<Nqtls;qi++)flow_qtl_avg[ion1][qi] = flow_qtl_avg[ion1][qi]/cycflow_expers[ion1];
    }


    for(idx=0;idx<nscenario;idx++)
        free_qtls(&st_nofX[idx]);


	print_lev = pl_0;

    printarray("mod_wgqtls_rec.csv",st->all_qtls,Nqtls,0);
    for(ri=0;ri<2;ri++)
        printarray("mod_wgqtls_rec.csv",st->all_qtls_rec[ri],Nqtls,1);

    printarray("mod_nof_wgqtls_rec.csv",st_noflow.all_qtls,Nqtls,0);
    for(ri=0;ri<2;ri++)
        printarray("mod_nof_wgqtls_rec.csv",st_noflow.all_qtls_rec[ri],Nqtls,1);

    printarray("mod_edgqtls_rec.csv",st->edge_qtls,Nqtls,0);
    for(ri=0;ri<2;ri++)
        printarray("mod_edgqtls_rec.csv",st->edge_qtls_rec[ri],Nqtls,1);

    printarray("mod_nof_edgqtls_rec.csv",st_noflow.edge_qtls,Nqtls,0);
    for(ri=0;ri<2;ri++)
        printarray("mod_nof_edgqtls_rec.csv",st_noflow.edge_qtls_rec[ri],Nqtls,1);

    //for each flow:
    printarray("flow_qtl_avg.csv", flow_qtl_avg[0],Nqtls,0);
    for(idx=1;idx<ncycflows;idx++)
        printarray("flow_qtl_avg.csv", flow_qtl_avg[idx],Nqtls,1);


	free_mats( &vf_noflow,&pf_noflow, &ht_noflow, &sk_noflow);

	gsl_vector_free(Alev_0);gsl_vector_free(Plev_0);gsl_vector_free(zlev_0);gsl_vector_free(eps_0);
	gsl_vector_free(ergPprobs);gsl_vector_free(ergAprobs);

}


int sim( struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk ){

	int i,ti,ll,ji,ai;
	int distiter;

	double cmzprob[NA][NZ], cmzprob0[NA][NZ];
	double cmepsprob[NA][NE],cmepsprob0[NA][NE];
	double cmjprob[JJ];
	double cmAtrans[NA][NA];
	double cmPtrans[JJ][NP][NP];
	double cmxStrans[NS][NS];

	double cmxSprob[NS];

	gsl_matrix ** swprob_hist = malloc(sizeof(gsl_matrix *)*Npaths);
	gsl_matrix ** WE_hist = malloc(sizeof(gsl_matrix *)*Npaths);
	gsl_matrix ** WU_hist = malloc(sizeof(gsl_matrix *)*Npaths);
	gsl_matrix ** epssel_hist = malloc(sizeof(gsl_matrix *)*Npaths);

	double quitwage[Npaths] , noquitwage[Npaths], Nquit[Npaths], Nnoquit[Npaths], quitVFdiff[Npaths];
	double quitWE[Npaths],quitWU[Npaths],noquitWE[Npaths],noquitWU[Npaths];
	long quitepsval[Npaths][NE];
	long quitPval[Npaths][NP];
	long quitzval[Npaths][NZ];long noquitzval[Npaths][NZ];


	for(ll=0;ll<Npaths;ll++){
		swprob_hist[ll] = gsl_matrix_calloc(Nsim,TT);
		WE_hist[ll] = gsl_matrix_calloc(Nsim,TT);
		WU_hist[ll] = gsl_matrix_calloc(Nsim,TT);
		epssel_hist[ll] = gsl_matrix_calloc(Nsim,TTT);
		gsl_matrix_int_set_all(ht->jhist[ll],-1);
	}

	cmjprob[0] = par->jprob->data[0];
	for(ji=0;ji<JJ-1;ji++) cmjprob[ji+1] = par->jprob->data[ji+1] + cmjprob[ji];

	for(ai=0;ai<NA;ai++){
        cmepsprob[ai][0] = par->epsprob[ai]->data[0];
        for(i=0;i<NE-1;i++) cmepsprob[ai][i+1] = par->epsprob[ai]->data[i+1] + cmepsprob[ai][i];
        for(i=0;i<NE;i++) cmepsprob0[ai][i]=cmepsprob[ai][i];
        cmzprob[ai][0] = par->zprob[ai]->data[0];
		for(i=0;i<NZ-1;i++) cmzprob[ai][i+1] = par->zprob[ai]->data[i+1] + cmzprob[ai][i];
		for(i=0;i<NZ;i++) cmzprob0[ai][i] = cmzprob[ai][i];
	}

	gsl_matrix * xStrans2 = gsl_matrix_calloc(NS,NS);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,par->xStrans,par->xStrans,0.,xStrans2);
	gsl_matrix * xSprob = gsl_matrix_calloc(NS,NS);
	for(i=0;i<NS*NS*NS*NS*NS*NS*NS*NS;i++) {
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,par->xStrans,xStrans2,0.,xSprob);
		gsl_matrix_memcpy(xStrans2,xSprob);
	}
	cmxSprob[0]=gg_get(xSprob,0,0);
	for(i=0;i<NS-1;i++) cmxSprob[i+1] = gg_get(xSprob,0,i+1) + cmxSprob[i];
	gsl_matrix_free(xSprob);gsl_matrix_free(xStrans2);

	for(ji=0;ji<JJ;ji++){
		for(i=0;i<NP;i++){
			for(ti=0;ti<NP;ti++) cmPtrans[ji][i][ti] = ti>0 ? cmPtrans[ji][i][ti-1] + gg_get(par->Ptrans[ji],i,ti) : gg_get(par->Ptrans[ji],i,ti);
		}
	}
	for(i=0;i<NA;i++){
		for(ti=0;ti<NA;ti++) cmAtrans[i][ti] = ti>0 ? gg_get(par->Atrans,i,ti) + cmAtrans[i][ti-1]: gg_get(par->Atrans,i,ti) ;
	}

	for(i=0;i<NS;i++){
		for(ti=0;ti<NS;ti++) cmxStrans[i][ti] = ti>0 ? gg_get(par->xStrans,i,ti) + cmxStrans[i][ti-1] : gg_get(par->xStrans,i,ti);
	}
	double topepsprob;
	//for(distiter = 0;distiter<maxiter;distiter++){
	for(distiter = 0;distiter<4;distiter++){
		#pragma omp parallel for private(i,ti,ll,ji) firstprivate(cmzprob,cmepsprob,cmzprob0,cmepsprob0,cmjprob,cmAtrans,cmPtrans,cmxStrans)
		for(ll=0;ll<Npaths;ll++){
			int * xSt,*xStm1; //xS.z.eps
            int * zt,*ztm1;int * epst,*epstm1;
			int At,Atm1;
			int Pt[JJ],Ptm1[JJ];
			int * jt,*jtm1;
			int * ut,*utm1;

			quitwage[ll] = 0. ; noquitwage[ll]=0.; Nquit[ll] =0.; Nnoquit[ll]=0.;  quitVFdiff[ll]=0.;
			quitWE[ll]=0.;quitWU[ll]=0.;noquitWE[ll]=0.;noquitWU[ll]=0.;
			for(i=0;i<NE;i++)quitepsval[ll][i]=0;for(i=0;i<NZ;i++)quitzval[ll][i]=0;for(i=0;i<NZ;i++)noquitzval[ll][i]=0;
			for(i=0;i<NP;i++)quitPval[ll][i]=0;

			xSt = malloc(sizeof(int)*Nsim);
            zt = malloc(sizeof(int)*Nsim);
            epst = malloc(sizeof(int)*Nsim);
            jt = malloc(sizeof(int)*Nsim);
			xStm1 = malloc(sizeof(int)*Nsim);
            epstm1 = malloc(sizeof(int)*Nsim);
            ztm1 = malloc(sizeof(int)*Nsim);

			jtm1 = malloc(sizeof(int)*Nsim);
			ut = malloc(sizeof(int)*Nsim);
			utm1 = malloc(sizeof(int)*Nsim);

			int ai,gi,si,zi,thi,xi,ii,iU,jji;

			// initial productivity:
			At = NA>2 ?  NA/2 : 1;
			for(ji=0;ji<JJ;ji++) Pt[ji] = NP/2;

			// initial allocation for all these guys
			for(i=0;i<Nsim;i++) {
				jt[i] = 0; //gsl_interp_bsearch(cmjprob,gg_get(sk->jsel[ll],i,0),0,JJ-1); //
				for (ji = 0; ji < JJ; ji++) if (gg_get(sk->jsel[ll], i, 0) > cmjprob[ji]) ++ jt[i] ;
				xSt[i] = 0;
				for (si = 0; si < NS; si++) if (gg_get(sk->xSsel[ll], i, 0) > cmxSprob[si]) ++ xSt[i] ;
				zt[i] =0;
				for(zi=0;zi<NZ;zi++) if( gg_get(sk->zsel[ll],i,0) > cmzprob0[At][zi] ) ++ zt[i] ;
				epst[i] =0;
				for(thi=0;thi<NE;thi++) if( gg_get(sk->epssel[ll],i,0) > cmepsprob0[At][thi] ) ++ epst[i];

				if( gg_get( sk->lambdaUSsel[ll],i,0) <  urt_avg ){
					ut[i] =1;
				}else{
					ut[i] =0;
				}
                // save this guy's state
                jtm1[i] = jt[i];
                utm1[i] = ut[i];
                xStm1[i] = xSt[i];
                ztm1[i] = zt[i];
                epstm1[i] = epst[i];
			}



            double Athist[TTT];double Pjhist[JJ][TTT];
			double lambdaEMhist[TTT],lambdaEShist[TTT],lambdaUShist[TTT],lambdaUMhist[TTT];
			double deltahist[TTT],zlosshist[TTT];
			for(ti=0;ti<TTT;ti++) {
                // increment shocks
                Atm1 = At;
                At = 0;
                for (ai = 0; ai < NA; ai++) if (gsl_vector_get(sk->Asel[ll], ti) > cmAtrans[Atm1][ai]) ++At;
                Athist[ti] = At;
                if (ti >= burnin) gsl_vector_int_set(ht->Ahist[ll], ti - burnin, At);
                for (ji = 0; ji < JJ; ji++) {
                    Ptm1[ji] = Pt[ji];
                    Pt[ji] = 0;
                    for (ai = 0; ai < NP; ai++) if (gg_get(sk->Psel[ll], ti, ji) > cmPtrans[ji][Ptm1[ji]][ai]) ++Pt[ji];
                    Pjhist[ji][ti] = Pt[ji];
                }

                lambdaEMhist[ti] = (par->lambdaEM0 * exp(par->lambdaEM_Acoef * (double) At));//par->Alev->data[At]));
                lambdaEShist[ti] = (par->lambdaES0 * exp(par->lambdaES_Acoef * (double) At));//par->Alev->data[At]));
                lambdaUShist[ti] = (par->lambdaUS0 * exp(par->lambdaUS_Acoef * (double) At));//par->Alev->data[At]));
                lambdaUMhist[ti] = (par->lambdaUM0 * exp(par->lambdaUM_Acoef * (double) At));//par->Alev->data[At]));
                deltahist[ti]    =  par->delta_avg * exp(par->delta_Acoef * (double) At);//par->Alev->data[At]);
                zlosshist[ti]    =  par->zloss     * exp(par->zloss_Acoef * (double) At); //par->Alev->data[At]
            }

            for(i=0;i<Nsim;i++){
                for(ti=0;ti<TTT;ti++) {

                    double lambdaEM_hr = lambdaEMhist[ti] ;
                    double lambdaES_hr = lambdaEShist[ti];
                    double lambdaUM_hr = lambdaUMhist[ti];
                    double lambdaUS_hr = lambdaUShist[ti];
                    double delta_hr = deltahist[ti];
                    double zloss_hr = zlosshist[ti];

                    for(ji=0;ji<JJ;ji++) Pt[ji] = Pjhist[ji][ti];
                    At = Athist[ti];

                    // increment individual specific shocks:
                    xSt[i] = 0;
                    for (si = 0; si < NS; si++) if (gg_get(sk->xSsel[ll], i, ti) > cmxStrans[xStm1[i]][si]) ++xSt[i];
                    // no longer incrementing z:
                    if (gg_get(sk->xGsel[ll], i, ti) > (1. -
                                                        par->update_z)) { //shouldn't technically use xGsel, but it's basically uncorrelated from all the other stuff I care about with this decision.
                        zt[i] = 0;
                        for (zi = 0; zi < NZ; zi++) if (gg_get(sk->zsel[ll], i, ti) > cmzprob[At][zi]) ++zt[i];
                    } else {
                        zt[i] = ztm1[i];
                    }
                    // no need to redraw epsilon unless later there's a job switch
                    epst[i] = epstm1[i];
                    ut[i] = utm1[i]; // this is not strictly necessary except in the first period?

                    ii = At * NP * NS * NZ * NE + Pt[jt[i]] * NS * NZ * NE + xSt[i] * NZ * NE + zt[i] * NE + epst[i];
                    iU = At * NP * NS * NZ + Pt[jt[i]] * NS * NZ + xSt[i] * NZ + zt[i];

                    //evaluate decision rules for the worker:
                    if (utm1[i] == 0) {
                        //ji = jt[i];
                        // employed workers' choices
                        // separate if zi =0 or unemployment shock or want to separate
                        if (delta_hr > gg_get(sk->dsel[ll], i, ti)
                            || zloss_hr > gg_get(sk->xGsel[ll], i, ti)
                            || gg_get(vf->WU, iU, jt[i]) > gg_get(vf->WE, ii, jt[i])) {
                            // mark unemployed
                            ut[i] = 1;

                            // record some stuff about whether voluntary or not
                            if (gg_get(vf->WU, iU, jt[i]) > gg_get(vf->WE, ii, jt[i]) &&
                                delta_hr <= gg_get(sk->dsel[ll], i, ti) && gg_get(sk->xGsel[ll], i, ti) >= zloss_hr) {
                                Nquit[ll] += 1.;
                                quitVFdiff[ll] += gg_get(vf->WU, iU, jt[i]) - gg_get(vf->WE, ii, jt[i]);
                                quitepsval[ll][epst[i]]++;
                                quitPval[ll][Pt[jt[i]]]++;
                                quitzval[ll][zt[i]]++;
                                quitWE[ll] += gg_get(vf->WE, ii, jt[i]);
                                quitWU[ll] += gg_get(vf->WU, iU, jt[i]);
                            } else {
                                noquitzval[ll][zt[i]]++;
                                Nnoquit[ll] += 1;
                                noquitWE[ll] += gg_get(vf->WE, ii, jt[i]);
                                noquitWU[ll] += gg_get(vf->WU, iU, jt[i]);
                            }


                        } else { //did not separate
                            ut[i] = 0;
                            // stay or go?
                            if (ti >= burnin) gg_set(swprob_hist[ll], i, ti - burnin, gg_get(pf->mE, ii, jt[i]));
                            if (gg_get(sk->msel[ll], i, ti) <= gg_get(pf->mE, ii, jt[i])) {
                                //RE
                                double cumsij[JJ];
                                double alphasij[JJ];
                                double sEhr;
                                for (jji = 0; jji < JJ; jji++) {
                                    sEhr = gg_get(pf->sE[jji], ii, jt[i]);
                                    alphasij[jji] = gsl_min(par->alpha0 *
                                                            exp(par->alpha_nf[jt[i]][jji] * par->alphaE1)
                                                            * pow(sEhr, 1. - par->alphaE1), 1.);
                                }
                                cumsij[0] = alphasij[0];
                                for (jji = 1; jji < JJ; jji++) {
                                    sEhr = gg_get(pf->sE[jji], ii, jt[i]);
                                    cumsij[jji] = cumsij[jji - 1] + alphasij[jji];
                                }//record sij
                                if (ti >= burnin) {
                                    for (jji = 0; jji < JJ; jji++)
                                        gg_set(ht->sijhist[ll][jji], i, ti - burnin, alphasij[jji]);
                                }
                                if (gg_get(sk->jsel[ll], i, ti) > cumsij[JJ - 1]) {
                                    jt[i] = jtm1[i];
                                } else {
                                    // successfully switched: if( gg_get(sk->jsel[ll],i,ti) <  alpha(sij[JJ]))
                                    jt[i] = 0;
                                    for (jji = 0; jji < JJ; jji++) {
                                        if (gg_get(sk->jsel[ll], i, ti) > cumsij[jji]) {
                                            jt[i]++;
                                        }
                                    }
                                }

                                if (jt[i] != jtm1[i]) { // switchers
                                    // draw a new z:
                                    int zipp = 0;
                                    for (zi = 0; zi < NZ; zi++)
                                        if (gg_get(sk->zsel[ll], i, ti) > cmzprob[At][zi]) ++zipp;
                                    // also switch employers?
                                    if (gg_get(sk->lambdaEMsel[ll], i, ti) < lambdaEM_hr) {
                                        int xti1 = 0; // lose specific skill
                                        // draw a new epsilon
                                        gg_set(epssel_hist[ll], i, ti, gg_get(sk->epssel[ll], i, ti));
                                        int epspp = 0;
                                        for (thi = 0; thi < NE; thi++) {
                                            if (gg_get(sk->epssel[ll], i, ti) > cmepsprob[At][thi])++epspp;
                                        }
                                        epspp = epspp > NE - 1 ? NE - 1 : epspp;
                                        int ip = At * NP * NS * NZ * NE + Pt[jt[i]] * NS * NZ * NE + xti1 * NZ * NE +
                                                 zipp * NE + epspp; // new epislon
                                        int ipp = At * NP * NS * NZ * NE + Pt[jt[i]] * NS * NZ * NE + xti1 * NZ * NE +
                                                  zipp * NE + epst[i]; //oldepsilon
                                        // either you want to or get godfather shock (using lambdaESsel b/c wasn't relevant here, so that's an uncorrelated draw)
                                        if (gg_get(vf->WE, ip, jt[i]) >= gg_get(vf->WE, ipp, jt[i]) ||
                                            gg_get(sk->lambdaESsel[ll], i, ti) < par->gdfthr) {
                                            if (ti >= burnin) ggi_set(ht->J2Jhist[ll], i, ti - burnin, 1);
                                            xSt[i] = xti1;
                                            zt[i] = zipp;
                                            epst[i] = epspp;
                                        } else {//don't switch employer
                                            if (ti >= burnin) ggi_set(ht->J2Jhist[ll], i, ti - burnin, 0);
                                            xSt[i] = xti1;
                                            zt[i] = zipp;
                                        }
                                    } else if (gg_get(sk->lambdaUMsel[ll], i, ti) < par->stwupdate) {
                                        // occ switch and no job switch
                                        zt[i] = zipp;
                                        xSt[i] = 0;
                                    } else {
                                        zt[i] = ztm1[i]; // no occupation switch
                                        jt[i] = jtm1[i];
                                    }
                                }// else nothing happens b/c/ did not successfully switch occupations
                            } else {
                                //WEs
                                if (gg_get(sk->lambdaESsel[ll], i, ti) < lambdaES_hr) {
                                    gg_set(epssel_hist[ll], i, ti, gg_get(sk->epssel[ll], i, ti));
                                    if (gg_get(sk->lambdaEMsel[ll], i, ti) <
                                        par->gdfthr) { //godfather (gamma) shock? lambdaEMsel hasn't been used for this guy
                                        int epspp = 0;
                                        for (thi = 0; thi < NE; thi++) {
                                            if (gg_get(sk->epssel[ll], i, ti) > cmepsprob[At][thi])++epspp;
                                        }
                                        epspp = epspp > NE - 1 ? NE - 1 : epspp;
                                        // want to quit?
                                        iU = At * NP * NS * NZ + Pt[jt[i]] * NS * NZ + xSt[i] * NZ + zt[i];
                                        int ipp = At * NP * NS * NZ * NE + Pt[jt[i]] * NS * NZ * NE + xSt[i] * NZ * NE +
                                                  zt[i] * NE + epst[i];
                                        if(gg_get(vf->WE, ipp, jt[i]) >= gg_get(vf->WU, iU, jt[i])) {
                                            epst[i] = epspp;
                                            if (ti >= burnin) ggi_set(ht->J2Jhist[ll], i, ti - burnin, 1);
                                        }else{
                                            ut[i] =1;
                                        }

                                    } else { //climbing the ladder!
                                        if ((gg_get(sk->epssel[ll], i, ti) > cmepsprob[At][epst[i]])) {
                                            if (ti >= burnin) ggi_set(ht->J2Jhist[ll], i, ti - burnin, 1);

                                            epst[i] = 0;
                                            for (thi = 0; thi < NE; thi++) {
                                                if (gg_get(sk->epssel[ll], i, ti) > cmepsprob[At][thi])++epst[i];
                                            }
                                            epst[i] = epst[i] > NE - 1 ? NE - 1 : epst[i];
                                        }
                                    }

                                }
                            }
                            // impose wage stickiness for stayers. Doing ti> burnin to start at period 1, not 0
                            if (ti > burnin) {
                                if (ggi_get(ht->J2Jhist[ll], i, ti - burnin) == 0 && jt[i] == jtm1[i]) {
                                    if (gg_get(sk->lambdaUSsel[ll], i, ti) <
                                        1. - par->stwupdate) //lambdaUsel not ever used before here
                                        gg_set(ht->whist[ll], i, ti - burnin,
                                               gg_get(ht->whist[ll], i, ti - burnin - 1));
                                }
                            }
                        }

                    }
                    if (ut[i] == 1) {  // ut ==1 unemployed. Can search even if just lost job
                        ii = At * NP * NS * NZ + Pt[jt[i]] * NS * NZ + xSt[i] * NZ + zt[i];

                        if (ti >= burnin) { gg_set(ht->whist[ll], i, ti - burnin, 0.); }
                        // stay or go?
                        if (ti >= burnin) gg_set(swprob_hist[ll], i, ti - burnin, gg_get(pf->mU, ii, jt[i]));
                        if (gg_get(pf->mU, ii, jt[i]) > gg_get(sk->msel[ll], i, ti) ||
                            // switch voluntarily or involuntarily
                            zloss_hr > gg_get(sk->xGsel[ll], i, ti)) {
                            //RU
                            double cumsij[JJ]; //cumul match prob for occupational placement
                            double alphasijhr[JJ];
                            double sUhr;
                            for (jji = 0; jji < JJ; jji++) {
                                sUhr = gg_get(pf->sU[jji], ii, jt[i]);
                                alphasijhr[jji] = gsl_min(par->alpha0 *
                                                          exp(par->alpha_nf[jt[i]][jji] * par->alphaU1) *
                                                          pow(sUhr, 1. - par->alphaU1), 1.);
                            }
                            cumsij[0] = alphasijhr[0];
                            for (jji = 1; jji < JJ; jji++)
                                cumsij[jji] = cumsij[jji - 1] + alphasijhr[jji];
                            // if forced to move, make the probabilities sum to 1:
                            if (zloss_hr > gg_get(sk->xGsel[ll], i, ti)) {
                                for (jji = 0; jji < JJ; jji++)
                                    alphasijhr[jji] = alphasijhr[jji] / cumsij[JJ - 1];
                            }
                            //record sij
                            if (ti >= burnin) {
                                for (jji = 0; jji < JJ; jji++)
                                    gg_set(ht->sijhist[ll][jji], i, ti - burnin, alphasijhr[jji]);
                            }
                            if (gg_get(sk->jsel[ll], i, ti) > cumsij[JJ - 1]) {
                                jt[i] = jtm1[i];
                            } else {
                                // successfully switched: if( gg_get(sk->jsel[ll],i,ti) <  sij[JJ])
                                jt[i] = 0;
                                for (jji = 0; jji < JJ; jji++)
                                    if (gg_get(sk->jsel[ll], i, ti) > cumsij[jji]) jt[i]++;
                            }
                        }
                        if (jt[i] != jtm1[i]) { // switchers
                            xSt[i] = 0; // lose specific skill
                            // draw a new z:
                            zt[i] = 0;
                            for (zi = 0; zi < NZ; zi++) {
                                if (gg_get(sk->zsel[ll], i, ti) > cmzprob[At][zi]) ++zt[i];
                            }
                            if (gg_get(sk->lambdaUMsel[ll], i, ti) < lambdaUM_hr) {
                                // draw a new epsilon
                                epst[i] = 0;
                                gg_set(epssel_hist[ll], i, ti, gg_get(sk->epssel[ll], i, ti));
                                for (thi = 0; thi < NE; thi++)
                                    if (gg_get(sk->epssel[ll], i, ti) > cmepsprob[At][thi])
                                        ++epst[i];
                                int iM = At * NP * NS * NZ * NE + Pt[jt[i]] * NS * NZ * NE + xSt[i] * NZ * NE +
                                         zt[i] * NE + epst[i];
                                int iU = At * NP * NS * NZ + Pt[jt[i]] * NS * NZ + xSt[i] * NZ + zt[i];
                                if (gg_get(vf->WE, iM, jt[i]) < gg_get(vf->WU, iU, jt[i])) {
                                    ut[i] = 1;
                                    epst[i] = epstm1[i];
                                } else {
                                    ut[i] = 0;
                                }
                            }
                        } else {
                            if (gg_get(sk->lambdaUSsel[ll], i, ti) < lambdaUS_hr) { // found a job??
                                epst[i] = 0;
                                gg_set(epssel_hist[ll], i, ti, gg_get(sk->epssel[ll], i, ti));
                                for (thi = 0; thi < NE; thi++)
                                    if (gg_get(sk->epssel[ll], i, ti) > cmepsprob[At][thi])
                                        ++epst[i];
                                int iM = At * NP * NS * NZ * NE + Pt[jt[i]] * NS * NZ * NE + xSt[i] * NZ * NE +
                                         zt[i] * NE + epst[i];
                                int iU = At * NP * NS * NZ + Pt[jt[i]] * NS * NZ + xSt[i] * NZ + zt[i];
                                if (gg_get(vf->WE, iM, jt[i]) < gg_get(vf->WU, iU, jt[i])) {
                                    ut[i] = 1;
                                    epst[i] = epstm1[i];
                                } else {
                                    ut[i] = 0;
                                }
                            }
                        }
                    }
                    //put the shocks in the history matrix
                    if (ti >= burnin) {
                        ggi_set(ht->xShist[ll], i, ti - burnin, xStm1[i]);
                        ggi_set(ht->zhist[ll], i, ti - burnin, ztm1[i]);
                        ggi_set(ht->epshist[ll], i, ti - burnin, epstm1[i]);
                        ggi_set(ht->uhist[ll], i, ti - burnin, ut[i]);
                        ggi_set(ht->jhist[ll], i, ti - burnin, jt[i]);
                        ggi_set(ht->Phist[ll], i, ti - burnin, Pt[jt[i]]);

                        gg_set(WE_hist[ll], i, ti - burnin, gg_get(vf->WE, ii, jt[i]));
                        gg_set(WU_hist[ll], i, ti - burnin, gg_get(vf->WU, iU, jt[i]));
                    }

                    // save this guy's state
                    jtm1[i] = jt[i];
                    utm1[i] = ut[i];
                    xStm1[i] = xSt[i];
                    ztm1[i] = zt[i];
                    epstm1[i] = epst[i];
                    // save wages
                    double wagehr = exp(par->AloadP->data[jt[i]] * par->Alev->data[At] +
                                        par->Plev->data[Pt[jt[i]]] +
                                        par->epslev->data[epst[i]] +
                                        par->zlev->data[zt[i]] +
                                        par->xSlev->data[xSt[i]] +
                                        occ_wlevs[jt[i]]) + 1;
                    if (ti >= burnin)
                        gg_set(ht->whist[ll], i, ti - burnin, wagehr);

                }// t=1:TTT
            } // i=1:NSim


			quitVFdiff[ll] /= Nquit[ll];


            free(xSt);free(xStm1);free(zt);free(ztm1);free(epst);free(epstm1);
            free(jt);free(jtm1);free(ut);free(utm1);
		} // end omp loop over ll
		double nemp =0., nemp_At[NA];
		topepsprob = cmepsprob0[1][NE-1] - cmepsprob0[1][NE-2];

		double topzprob[NA] ;
		for(ai=0;ai<NA;ai++)
			topzprob[ai] = (cmzprob0[ai][NZ-1] - cmzprob0[ai][NZ-2]);

		for(ai=0;ai<NA;ai++){
			nemp_At[ai]=0.;
			for(i=0;i<NZ;i++)cmzprob0[ai][i]=0.;
            for(i=0;i<NE;i++)cmepsprob0[ai][i]=0.;
		}
		for(ll=0;ll<Npaths;ll++){
			for(ti=0;ti<TT;ti++){
				int Ahr = gsl_vector_int_get( ht->Ahist[ll],ti );
				for(i=0;i<Nsim;i++){
					if( gsl_matrix_int_get(ht->uhist[ll],i,ti)==0 ){
						int epshr = gsl_matrix_int_get( ht->epshist[ll],i,ti);
						cmepsprob0[Ahr][ epshr ] += 1.;
						int zhr =  gsl_matrix_int_get( ht->zhist[ll],i,ti);
						cmzprob0[Ahr][zhr] += 1.;
						nemp += 1.;
						nemp_At[Ahr] +=1.;
					}
				}
			}
		}
		for(ai=0;ai<NA;ai++){
            for(i=0;i<NE;i++)cmepsprob0[ai][i] *= 1./nemp_At[ai];
            for(i=0;i<NZ;i++)cmzprob0[ai][i] *= 1./nemp_At[ai];
		}

		for(ai=0;ai<NA;ai++){
            for(i=1;i<NE;i++) cmepsprob0[ai][i] += cmepsprob0[ai][i-1];
            for(i=1;i<NZ;i++) cmzprob0[ai][i] += cmzprob0[ai][i-1];
		}
		topepsprob -= (cmepsprob0[1][NE-1] - cmepsprob0[1][NE-2] );
		topepsprob  = topepsprob < 0. ? - topepsprob : topepsprob;
		for(ai=0;ai<NA;ai++){
			topzprob[ai]   -= (cmzprob0[ai][NZ-1] - cmzprob0[ai][NZ-2]);
			topzprob[ai]    = topzprob[ai] < 0. ? - topzprob[ai] : topzprob[ai];
		}
		if( (gsl_min(topzprob[1],topepsprob) <1e-5 && gsl_max(topzprob[1],topepsprob) < 1e-4)  ) // && distiter >0
			break;
	}


	if(verbose>2){
		double Nquit_pr=0,Nnoquit_pr=0,quitVFdiff_pr=0;
		double noquitWE_pr=0.,noquitWU_pr=0.,quitWE_pr=0.,quitWU_pr=0.;
		long quitepsval_pr[NE],quitzval_pr[NZ],noquitzval_pr[NZ],quitPval_pr[NP];
		for(i=0;i<NE;i++)quitepsval_pr[i] =0;
		for(i=0;i<NZ;i++)quitzval_pr[i] =0;
		for(i=0;i<NZ;i++)noquitzval_pr[i] =0;
		for(i=0;i<NP;i++)quitPval_pr[i] =0;
		for(ll=0;ll<Npaths;ll++){
			quitVFdiff_pr += quitVFdiff[ll];
			Nquit_pr +=Nquit[ll];
			Nnoquit_pr += Nnoquit[ll];
			noquitWE_pr+= noquitWE[ll];
			noquitWU_pr+= noquitWU[ll];
			quitWE_pr+= quitWE[ll];
			quitWU_pr+= quitWU[ll];
			for(i=0;i<NE;i++)quitepsval_pr[i] += quitepsval[ll][i];
			for(i=0;i<NZ;i++)quitzval_pr[i] += quitzval[ll][i];
			for(i=0;i<NZ;i++)noquitzval_pr[i] += noquitzval[ll][i];
			for(i=0;i<NP;i++)quitPval_pr[i] += quitPval[ll][i];
		}
		printf(" The fraction of quits out of separations was %f. The VF difference was %f \n", Nquit_pr/(Nquit_pr+Nnoquit_pr), quitVFdiff_pr);
		printf(" Location is (");for(i=0;i<NE;i++) printf(" %ld,", quitepsval_pr[i]); printf(") in epsilon \n");
		printf(" Location is (");for(i=0;i<NZ;i++) printf(" %ld,", quitzval_pr[i]); printf(") in z \n");
		printf(" No quits location is (");for(i=0;i<NZ;i++) printf(" %ld,", noquitzval_pr[i]); printf(") in z \n");
		printf(" Location is (");for(i=0;i<NP;i++) printf(" %ld,", quitPval_pr[i]); printf(") in P \n");
		printf("WE in quit is %f and no quit %f \n", quitWE_pr/Nquit_pr, noquitWE_pr/Nnoquit_pr );
		printf("WU in quit is %f and no quit %f \n", quitWU_pr/Nquit_pr, noquitWU_pr/Nnoquit_pr );
		printf("Number of distribution iterations is %d, distance %f \n",distiter,
		       gsl_max(topepsprob,-topepsprob) );
		// put in the WE and the WU for the levels of z Is WE coming down with z, or WU going up too fast?
	}



    if(print_lev >=2){
		int append =0;
		char matname[20];
		for(ll=0;ll<Npaths;ll++){
			sprintf(matname, "whist%s.csv",exper_f);
			printmat(matname,ht->whist[ll],append );
			sprintf(matname, "uhist%s.csv",exper_f);
			printmat_int(matname,ht->uhist[ll],append );
			sprintf(matname, "Ahist%s.csv",exper_f);
			printvec_int(matname,ht->Ahist[ll],append );
			sprintf(matname, "Phist%s.csv",exper_f);
			printmat_int(matname,ht->Phist[ll],append );
			sprintf(matname, "jhist%s.csv",exper_f);
			printmat_int(matname,ht->jhist[ll],append );
			sprintf(matname, "xShist%s.csv",exper_f);
			printmat_int(matname,ht->xShist[ll],append );
			sprintf(matname, "zhist%s.csv",exper_f);
			printmat_int(matname,ht->zhist[ll],append );
			sprintf(matname, "epshist%s.csv",exper_f);
			printmat_int(matname,ht->epshist[ll],append );
			sprintf(matname, "swpr_hist%s.csv",exper_f);
			printmat(matname, swprob_hist[ll],append );
			sprintf(matname, "WUhist%s.csv",exper_f);
			printmat(matname, WU_hist[ll],append );
			sprintf(matname, "WEhist%s.csv",exper_f);
			printmat(matname, WE_hist[ll],append );
			sprintf(matname, "epsselhist%s.csv",exper_f);
			printmat(matname, epssel_hist[ll],append );
			sprintf(matname, "J2Jhist%s.csv",exper_f);
			printmat_int(matname, ht->J2Jhist[ll],append );

			append = 1;
		}
	}

	for(ll=0;ll<Npaths;ll++){
		gsl_matrix_free(swprob_hist[ll]);gsl_matrix_free(epssel_hist[ll]);
		gsl_matrix_free(WE_hist[ll]);gsl_matrix_free(WU_hist[ll]);
	}
	free(swprob_hist);
	free(epssel_hist);
	free(WE_hist);
	free(WU_hist);

}

int sum_stats(   struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk, struct stats *st  ){

	int ll, i,ti,wi,si,ji,jji;
	int ri;

	int Nemp=0, Nunemp = 0, Nnosep=0, Nfnd =0, NswU = 0, NswE=0, NJ2J=0, Nspell=0, NswSt=0,Ndur_sw=0,Ndur_nosw=0,Nsep=0,NErisksep=0;
	int NdoubleswU =0 , NdoubleswE=0;
	int Nemp_rec[2]={0,0}, Nunemp_rec[2] = {0,0}, Nnosep_rec[2]={0,0}, Nfnd_rec[2]={0,0}, NswU_rec[2]= {0,0}, NswE_rec[2]={0,0}, NJ2J_rec[2]={0,0}, Nspell_rec[2]={0,0}, NswSt_rec[2]={0,0},Ndur_sw_rec[2]={0,0},Ndur_nosw_rec[2]={0,0};
	int Nsep_rec[2]={0,0},NErisksep_rec[2]={0,0};
	double Nocc_ten=0.;
	double * occwg = calloc(JJ*Npaths,sizeof(double));
	int    * occsz = calloc(JJ*Npaths,sizeof(int)   );
	double **occ_gflow= malloc(JJ*sizeof(double*)); //the realized gross flow
	double **occ_sijflow= malloc(JJ*sizeof(double*));  // the probabilistically defined flows
    double **occ_gflowE= malloc(JJ*sizeof(double*)); //the realized gross flow for EE
    double **occ_sijflowE= malloc(JJ*sizeof(double*));  // the probabilistically defined flows for EE
    double **occ_gflowU= malloc(JJ*sizeof(double*)); //the realized gross flow for EU
    double **occ_sijflowU= malloc(JJ*sizeof(double*));  // the probabilistically defined flows for EU

    for(ji=0;ji<JJ;ji++){
		occ_gflow[ji] = calloc(JJ, sizeof(double));
		occ_sijflow[ji] = calloc(JJ, sizeof(double));
		for(i=0;i<JJ;i++){
			occ_gflow[ji][i] = 0.;
			occ_sijflow[ji][i] = 0.;
		}
        occ_gflowE[ji] = calloc(JJ, sizeof(double));
        occ_sijflowE[ji] = calloc(JJ, sizeof(double));
        for(i=0;i<JJ;i++){
            occ_gflowE[ji][i] = 0.;
            occ_sijflowE[ji][i] = 0.;
        }
        occ_gflowU[ji] = calloc(JJ, sizeof(double));
        occ_sijflowU[ji] = calloc(JJ, sizeof(double));
        for(i=0;i<JJ;i++){
            occ_gflowU[ji][i] = 0.;
            occ_sijflowU[ji][i] = 0.;
        }
	}
	int Nfndmo =0, Nfndmo_denom=0, Nsepmo =0,Nsepmo_denom=0,NJ2Jmo=0;

	double invlaid_wval = -1e6;
	gsl_vector * Atrans_ergod = gsl_vector_calloc(NA);
	int recIndic[NA];double cumAtrans_ergod[NA];
	ergod_dist(par->Atrans,Atrans_ergod);
	cumAtrans_ergod[0] = Atrans_ergod->data[0];
	for(i=1;i<NA;i++) cumAtrans_ergod[i] =Atrans_ergod->data[i] + cumAtrans_ergod[i-1];
	for(i=0;i<NA;i++) recIndic[i] = 1-i; //1 when a recession, so this is the opposite of A's index

	if(verbose>1){
		int recFrac =0;
		for( ll=0;ll<Npaths;ll++){
			for(ti=0;ti<TT;ti++)
				recFrac += recIndic[gsl_vector_int_get(ht->Ahist[ll],ti)];
		}
		printf(" The fraction of periods in recession is %f, with %d periods \n", (double)recFrac /(double)(TT*Npaths) , recFrac);
	}

	for(ri=0;ri<3;ri++){ //run this 3 times for (0) expansion, (1) recession and (2) pooled
		Nemp     = 0;
		Nunemp   = 0;
		Nnosep   = 0;
		Nfnd     = 0;
		NswU     = 0;
		NswE     = 0;
		NJ2J     = 0;
		Nspell   = 0;
		NswSt    = 0;
		Ndur_sw  = 0;
		Ndur_nosw= 0;
		Nsep     = 0;
		NErisksep= 0;
		Nocc_ten = 0.;
		NdoubleswU =0 ; NdoubleswE=0;
		#pragma omp parallel for private(ll,ti,i,wi,si,ji,jji) \
		reduction( +: Nemp ) reduction( +:Nunemp) reduction( +:Nnosep) \
		reduction( +: Nfnd) reduction( +: NswU) reduction( +: NswE) reduction( +: Nspell) \
		reduction( +: NswSt) reduction( +: Ndur_sw) reduction( +: Ndur_nosw) reduction(+: Nsep) \
		reduction(+:NErisksep) reduction(+:Nocc_ten) reduction(+:NdoubleswU) reduction(+:NdoubleswE)
		for(ll=0;ll<Npaths;ll++){

			int count_wi;

			int Nemp_ll = 0, Nunemp_ll = 0,Nnosep_ll=0, Nfnd_ll =0, NswU_ll = 0, NswE_ll=0, NJ2J_ll=0,Nspell_ll=0, NswSt_ll =0, Ndur_sw_ll=0, Ndur_nosw_ll=0,Nsep_ll = 0,NErisksep_ll=0;
			int NdoubleswU_ll =0 , NdoubleswE_ll=0;
			int sw_spell=-1,tsep=-1;
			int Nocc_ll = 0;
			for(i=0;i<Nsim;i++){
				sw_spell = -1;tsep=-1;
				int Noccs_i=0;
				int lastswU_wi =0,lastswE_wi =0;
				for(wi=0;wi<(TT/Npwave);wi++){//first loop over waves (wi), then loop over reference month (si)
					count_wi=0;
	                int Nemp_wi = 0,Nspell_wi=0,Nfnd_wi=0,NJ2J_wi=0,NswE_wi=0,NswSt_wi=0,Nnosep_wi=0,Nunemp_wi=0,NswU_wi=0,Ndur_sw_wi=0,Ndur_nosw_wi=0;
	                int Nsep_wi = 0,NErisksep_wi=0;
					int NdoubleswU_wi =0 , NdoubleswE_wi=0;
				    for(si=0;si<Npwave;si++){
	                    ti = si + wi * Npwave;
						// should I count this observation? for ri=2, always true. For ri=0 only true when expansion and for ri=1 only true when recession
					    if(ri>= recIndic[ gsl_vector_int_get(ht->Ahist[ll],ti)] && count_wi==0)
					        count_wi = 1;
	                    if(ti<TT-1 && ti>0){
	                        if( ggi_get(ht->uhist[ll],i,ti) ==0 ){
	                            Nemp_wi=1;
		                        NErisksep_wi=1;
	                            occwg[ ll*JJ + ggi_get(ht->jhist[ll],i,ti)] += gg_get(ht->whist[ll],i ,ti);
		                        occsz[ ll*JJ + ggi_get(ht->jhist[ll],i,ti)] ++;
	                            if(ri==2) Nsepmo_denom ++;
		                        if( ggi_get(ht->uhist[ll],i,ti+1) ==1 ){
			                        if(ri==2) Nsepmo ++;
	                            	sw_spell = ggi_get(ht->jhist[ll],i,ti-1);
	                            	if(sw_spell >0 && sw_spell<JJ){
			                            tsep = ti;
			                            Nnosep_wi =0;
			                            Nsep_wi = 1;
	                            	}else{
			                            sw_spell = -1;
	                            	}
	                            }else{
	                                Nnosep_wi =1;
		                            if( ggi_get(ht->J2Jhist[ll],i,ti) ==1 ){
		                                NJ2J_wi = 1;
			                            if(ri==2) NJ2Jmo ++;
		                                if( ggi_get(ht->jhist[ll],i,ti-1) !=ggi_get(ht->jhist[ll],i,ti) ){
		                                	occ_gflow[ggi_get(ht->jhist[ll],i,ti-1)][ggi_get(ht->jhist[ll],i,ti)]++;
                                            occ_gflowE[ggi_get(ht->jhist[ll],i,ti-1)][ggi_get(ht->jhist[ll],i,ti)]++;
		                                	for(jji=0;jji<JJ;jji++){
				                                occ_sijflow[ggi_get(ht->jhist[ll],i,ti-1)][jji] += gg_get(ht->sijhist[ll][jji],i,ti);
                                                occ_sijflowE[ggi_get(ht->jhist[ll],i,ti-1)][jji] += gg_get(ht->sijhist[ll][jji],i,ti);
		                                	}
		                                	NswE_wi =1;
		                                	if( lastswE_wi ==1 ) NdoubleswE_wi = 1;
			                                if( lastswU_wi ==1 ) NdoubleswU_wi = 1;

			                                Noccs_i ++;
		                                }
		                            }else{  // not EE
		                                if( ggi_get(ht->jhist[ll],i,ti-1) !=ggi_get(ht->jhist[ll],i,ti) ) NswSt_wi = 1;
		                            }
	                            }
	                        }else{ //uhist = 1
	                            if(sw_spell>-1 && tsep>-1){
	                            	Nunemp_wi =1;
		                            if(ri==2) Nfndmo_denom ++;
	                            }
	                            if( ggi_get(ht->uhist[ll],i,ti+1) ==0 && sw_spell>-1 && tsep>-1 ){
	                                Nfnd_wi =1;
	                                Nspell_wi = 1;
		                            if(ri==2) Nfndmo ++;
		                            if( ggi_get(ht->jhist[ll],i,ti+1) != sw_spell ){
		                            	occ_gflow[sw_spell][ggi_get(ht->jhist[ll],i,ti+1)] ++;
                                        occ_gflowU[sw_spell][ggi_get(ht->jhist[ll],i,ti+1)] ++;
			                            for(jji=0;jji<JJ;jji++){
				                            occ_sijflow[sw_spell][jji] += gg_get(ht->sijhist[ll][jji],i,ti);
                                            occ_sijflowU[sw_spell][jji] += gg_get(ht->sijhist[ll][jji],i,ti);
			                            }
	                                    NswU_wi = 1; //only count switches at the end of the spell.
			                            if( lastswE_wi ==1 ) NdoubleswE_wi = 1;
			                            if( lastswU_wi ==1 ) NdoubleswU_wi = 1;
	                                    Ndur_sw_wi = ti - tsep+1;
	                                    Noccs_i += 1;
	                                }else{
		                                Ndur_nosw_wi = ti - tsep+1;
	                                }
	                            }
	                        }
	                    }
				    }
					if(NswSt_wi ==1 && (NswE_wi==1 || Nunemp_wi ==1) ) NswSt_wi=0; // stayer is exclusive

				    //sum over waves. Keeping s.t. count any transition only once during wave
				    if(count_wi==1){
		                Nemp_ll  += Nemp_wi;
		                Nnosep_ll+= Nnosep_wi;
		                NJ2J_ll  += NJ2J_wi;
		                NswE_ll  += NswE_wi;
		                NswSt_ll += NswSt_wi;
		                Nunemp_ll+= Nunemp_wi;
		                Nfnd_ll  += Nfnd_wi;
		                Nspell_ll+= Nspell_wi;
		                NswU_ll  += NswU_wi;
		                Ndur_sw_ll  += Ndur_sw_wi;
						Ndur_nosw_ll  += Ndur_nosw_wi;
						Nsep_ll += Nsep_wi;
						NErisksep_ll += NErisksep_wi;
						NdoubleswE_ll += NdoubleswE_wi;
					    NdoubleswU_ll += NdoubleswE_wi;
				    }
				    lastswE_wi = NswE_wi >0 ? 1 : 0;
					lastswU_wi = NswU_wi >0 ? 1 : 0;
				}
				Nocc_ll += Noccs_i; //counting the number of occupations for each individual.
			}

			Nemp += Nemp_ll;
			Nunemp += Nunemp_ll;
			Nnosep += Nnosep_ll;
			Nfnd += Nfnd_ll;
			NswU += NswU_ll;
			NswE += NswE_ll;
			NJ2J += NJ2J_ll;
			Nspell += Nspell_ll;
			NswSt  += NswSt_ll;
			Ndur_sw  += Ndur_sw_ll;
			Ndur_nosw  += Ndur_nosw_ll;
			Nsep += Nsep_ll;
			NErisksep += NErisksep_ll;
			Nocc_ten += Nocc_ll;
			NdoubleswE += NdoubleswE_ll;
			NdoubleswU += NdoubleswE_ll;
		}
		if(ri<2){
			Nemp_rec[ri] = Nemp;
			Nunemp_rec[ri] = Nunemp;
			Nnosep_rec[ri] = Nnosep;
			Nfnd_rec[ri] = Nfnd;
			NswU_rec[ri] = NswU;
			NswE_rec[ri] = NswE;
			NJ2J_rec[ri] = NJ2J;
			Nspell_rec[ri] = Nspell;
			NswSt_rec[ri]  = NswSt;
			Ndur_sw_rec[ri]  = Ndur_sw;
			Ndur_nosw_rec[ri]  = Ndur_nosw;
			Nsep_rec[ri] += Nsep;
			NErisksep_rec[ri] += NErisksep;
		}
	}
    st->J2Jprob  = Nnosep>0? (double) NJ2J / (double) Nnosep : 0.;
    st->findrate = Nunemp>0 ? (double) Nfnd / (double) Nunemp : 1.;
    st->seprate  = Nemp >0 ? (double) Nsep / (double) NErisksep : 1.; // can use Nfnd because only count separations that find again
    st->swProb_EE = NJ2J>0 ? (double) NswE / (double) NJ2J : 0.;
    st->swProb_U = Nspell>0 ? (double) NswU / (double) Nspell : 0.;
	st->swProb_st = Nnosep - NJ2J > 0 ? (double) NswSt / (double) (Nnosep- NJ2J): 0.;
    st->unrate =  (double) Nunemp / (double) ( Nunemp + Nemp );
    st->udur_nosw = Nspell - NswU>0 ? (double) Ndur_nosw/ (double)(Nspell - NswU): 0.;
	st->udur_sw = NswU>0 ? (double) Ndur_sw/ (double)NswU : 0.;
	st->occ_ten = (double)TT / ((double)Nocc_ten/ (double)Npaths / (double)Nsim);

	st->doubleswE = NswE>0 ? (double)NdoubleswE / (double)NswE : 0.;
	st->doubleswU = NswU>0 ? (double)NdoubleswU / (double)NswU : 0.;

	if(verbose>3) {
        printf("Conditional prob double switch, U: %f \n", st->doubleswU);
        printf("Conditional prob double switch, E: %f \n", st->doubleswE);
    }
    printf("swProb_st : %f, swProb_EE: %f,  swProb_U: %f \n", st->swProb_st , st->swProb_EE, st->swProb_U );
    printf("Average number of occupations per person: %f and occupational tenure %f \n",
            (double)Nocc_ten/ (double)Npaths / (double)Nsim, st->occ_ten);

	double occ_gtot = 0.; double occ_stot = 0.; double varGflow=0.;
	double occ_Ltot[JJ]; for(jji=0;jji<JJ;jji++) occ_Ltot[jji]=0.;
    double occ_marg[JJ]; for(jji=0;jji<JJ;jji++) occ_marg[jji]=0.;
    double occ_smarg[JJ]; for(jji=0;jji<JJ;jji++) occ_smarg[jji]=0.;
	for(jji=0;jji<JJ;jji++){
		for(ji=0;ji<JJ;ji++){
			occ_gtot+= occ_gflow[jji][ji];
			occ_Ltot[jji] += occ_gflow[jji][ji];
			occ_stot+= occ_sijflow[jji][ji];
		}
	}
	for(jji=0;jji<JJ;jji++){
	    for(ji=0;ji<JJ;ji++){
            occ_marg[ji]        += occ_gflow[jji][ji];
            occ_smarg[ji]       += occ_sijflow[jji][ji];
            occ_gflow[jji][ji]   =  occ_gtot>0 ? occ_gflow[jji][ji]/occ_gtot  : 0.;
			occ_sijflow[jji][ji] =  occ_stot>0 ? occ_sijflow[jji][ji]/occ_stot: 0.;
		}
        occ_Ltot[jji] = occ_gtot>0 ? occ_Ltot[jji]/occ_gtot : 0.;
    }
    for(jji=0;jji<JJ;jji++) occ_marg[jji]  = occ_marg[jji]/occ_gtot;
    for(jji=0;jji<JJ;jji++) occ_smarg[jji] = occ_smarg[jji]/occ_stot;
    for(jji=0;jji<JJ;jji++) st->occ_margflow[jji] = (1.-smth_flows)*occ_marg[jji] + smth_flows*occ_smarg[jji];

        for(jji=0;jji<JJ;jji++){
        si =0;
        double gflow_hr[JJ - 1];
        for (ji = 0; ji < JJ; ji++) {
            if (jji != ji){ gflow_hr[si] = occ_gflow[jji][ji]/occ_Ltot[jji]; si++; }
        }
        double vGflow_j =gsl_stats_variance(gflow_hr,1,JJ-1) * occ_Ltot[jji];
        varGflow = gsl_finite(vGflow_j) ?varGflow+ vGflow_j : varGflow;
    }
    si=0;
	for(jji=0;jji<JJ;jji++){

		for(ji=jji+1;ji<JJ;ji++){
				st->occ_netflow[si] = (1.-smth_flows)* (occ_gflow[jji][ji] - occ_gflow[ji][jji])
						+ smth_flows*(occ_sijflow[jji][ji] - occ_sijflow[ji][jji]);
				si++;
		}

	}
	// now compute flows for EE transitions
	occ_gtot=0.;occ_stot=0.;
    for(jji=0;jji<JJ;jji++) occ_Ltot[jji]=0.;
    for(jji=0;jji<JJ;jji++){
        for(ji=0;ji<JJ;ji++){
            occ_gtot+= occ_gflowE[jji][ji];
            occ_Ltot[jji] += occ_gflowE[jji][ji];
            occ_stot+= occ_sijflowE[jji][ji];
        }
    }
    for(jji=0;jji<JJ;jji++){
        for(ji=0;ji<JJ;ji++){
            occ_gflowE[jji][ji] = occ_gtot>0 ? occ_gflowE[jji][ji]/occ_gtot: 0.;
            occ_sijflowE[jji][ji] = occ_stot>0 ? occ_sijflowE[jji][ji]/occ_stot: 0.;
        }
        occ_Ltot[jji] = occ_gtot>0 ? occ_Ltot[jji]/occ_gtot : 0.;
    }
    if(occ_gtot <= 0.){
        printf("+++++++++++++++++++++++++++ \n occ_gtotE ==0 \n");
        for(jji=0;jji<JJ;jji++){
            for(ji=0;ji<JJ;ji++){ printf("%f,", occ_gflowE[jji][ji]);}
            printf("\n");
        }
    }
    if(occ_stot <= 0.){
        printf("+++++++++++++++++++++++++++ \n occ_stotE ==0 \n");
        for(jji=0;jji<JJ;jji++){
            for(ji=0;ji<JJ;ji++){ printf("%f,", occ_sijflowE[jji][ji]);}
            printf("\n");
        }
    }
    st->varGflowsE = 0.;
    for(jji=0;jji<JJ;jji++) {
        si =0;
        double gflow_hr[JJ - 1];
        for (ji = 0; ji < JJ; ji++) {
            if (jji != ji){ gflow_hr[si] = occ_gflowE[jji][ji]/occ_Ltot[jji]; si++; }
        }
        double vGflowE_j = gsl_stats_variance(gflow_hr,1,JJ-1) * occ_Ltot[jji];
        st->varGflowsE = gsl_finite(vGflowE_j)?vGflowE_j+st->varGflowsE :st->varGflowsE;
    }
    si =0;
    for(jji=0;jji<JJ;jji++){
        for(ji=jji+1;ji<JJ;ji++){
            st->occ_netflowE[si] = (1.-smth_flows)* (occ_gflowE[jji][ji] - occ_gflowE[ji][jji])
                                  + smth_flows*(occ_sijflowE[jji][ji] - occ_sijflowE[ji][jji]);
            si++;
        }
    }
// now compute flows for EU transitions
    occ_gtot=0.;occ_stot=0.;
    for(jji=0;jji<JJ;jji++) occ_Ltot[jji]=0.;
    for(jji=0;jji<JJ;jji++){
        for(ji=0;ji<JJ;ji++){
            occ_gtot+= occ_gflowU[jji][ji];
            occ_Ltot[jji] += occ_gflowU[jji][ji];
            occ_stot+= occ_sijflowU[jji][ji];
        }
    }
    for(jji=0;jji<JJ;jji++){
        for(ji=0;ji<JJ;ji++){
            occ_gflowU[jji][ji] =  occ_gtot>0 ? occ_gflowU[jji][ji]/occ_gtot : 0.;
            occ_sijflowU[jji][ji] =  occ_stot>0 ?  occ_sijflowU[jji][ji]/occ_stot: 0.;
        }
        occ_Ltot[jji] = occ_gtot>0 ? occ_Ltot[jji]/occ_gtot : 0.;
    }
    if(occ_gtot <= 0.){
        printf("+++++++++++++++++++++++++++ \n occ_gtotU ==0 \n");
        for(jji=0;jji<JJ;jji++){
            for(ji=0;ji<JJ;ji++){ printf("%f,", occ_gflowU[jji][ji]);}
            printf("\n");
        }
    }
    if(occ_stot <= 0.){
        printf("+++++++++++++++++++++++++++ \n occ_stotU ==0 \n");
        for(jji=0;jji<JJ;jji++){
            for(ji=0;ji<JJ;ji++){ printf("%f,", occ_sijflowU[jji][ji]);}
            printf("\n");
        }
    }
    st->varGflowsU = 0.;
    for(jji=0;jji<JJ;jji++) {
        si =0;
        double gflow_hr[JJ - 1];
        for (ji = 0; ji < JJ; ji++) {
            if (jji != ji){
                gflow_hr[si] =occ_Ltot[jji]>0? occ_gflowU[jji][ji]/occ_Ltot[jji] : 0.; si++; }
        }
        double vGflowU_j = gsl_stats_variance(gflow_hr,1,JJ-1) * occ_Ltot[jji];
        st->varGflowsU = gsl_finite(vGflowU_j)?st->varGflowsU+ vGflowU_j: st->varGflowsU;
    }
    si =0;
    for(jji=0;jji<JJ;jji++){
        for(ji=jji+1;ji<JJ;ji++){
            st->occ_netflowU[si] = (1.-smth_flows)* (occ_gflowU[jji][ji] - occ_gflowU[ji][jji])
                                   + smth_flows*(occ_sijflowU[jji][ji] - occ_sijflowU[ji][jji]);
            si++;
        }
    }

	st->nrmflows = 0;
	for(i=0;i<nflows;i++) st->nrmflows = st->occ_netflow[i] >0 ? st->occ_netflow[i]+st->nrmflows : - (st->occ_netflow[i])+st->nrmflows;
	for(i=0;i<nflows;i++) st->occ_netflow[i] = st->occ_netflow[i]/st->nrmflows;
    st->nrmflowsE = 0;
    for(i=0;i<nflows;i++) st->nrmflowsE = st->occ_netflowE[i] >0 ? st->occ_netflowE[i]+st->nrmflowsE : - (st->occ_netflowE[i])+st->nrmflowsE;
    st->nrmflowsU = 0;
    for(i=0;i<nflows;i++) st->nrmflowsU = st->occ_netflowU[i] >0 ? st->occ_netflowU[i]+st->nrmflowsU : - (st->occ_netflowU[i])+st->nrmflowsU;

    double varnetflowsU = gsl_stats_variance(st->occ_netflowU,1,nflows);
    double varnetflowsE = gsl_stats_variance(st->occ_netflowE,1,nflows);
    if(verbose>1) {
        printf(" \n +++++++++++++++ \n model varnetflowsU: %f  model varnetflowsE: %f \n",
               varnetflowsU, varnetflowsE);
        printf(" \n +++++++++++++++ \n model absdevnetflowsU: %f  model absdevnetflowsE: %f \n",
               st->nrmflowsU, st->nrmflowsE);
        printf(" \n +++++++++++++++ \n model varGflowsU: %f  model varGflowsE: %f \n",
               st->varGflowsU, st->varGflowsE);
    }
	double fndrt_rec0 = Nunemp_rec[0] >0 ? (double) Nfnd_rec[0] / (double) Nunemp_rec[0] : 1.;
	double fndrt_rec1 = Nunemp_rec[1] >0 ? (double) Nfnd_rec[1] / (double) Nunemp_rec[1] : 1.;
	st->fndrt_ratio = fndrt_rec1 >0 ? fndrt_rec0/fndrt_rec1 : 0.;

	double seprt_rec0 = Nunemp_rec[0] >0 ? (double) Nsep_rec[0] / (double) NErisksep_rec[0] : 1.;
	double seprt_rec1 = Nunemp_rec[1] >0 ? (double) Nsep_rec[1] / (double) NErisksep_rec[1] : 1.;
	st->seprt_ratio = seprt_rec1 > 0 ? seprt_rec0/seprt_rec1 : 0.;

	double swProb_EE_rec0 = NJ2J_rec[0] >0 ? (double) NswE_rec[0] / (double) NJ2J_rec[0] : 0.;
	double swProb_EE_rec1 = NJ2J_rec[1] >0 ? (double) NswE_rec[1] / (double) NJ2J_rec[1] : 0.;
	st->swProb_EE_ratio = swProb_EE_rec1>0 ? swProb_EE_rec0/swProb_EE_rec1 : 0.;

	double swProb_U_rec0 = Nspell_rec[0] >0 ? (double) NswU_rec[0] / (double) Nspell_rec[0] : 0.;
	double swProb_U_rec1 = Nspell_rec[1] >0 ? (double) NswU_rec[1] / (double) Nspell_rec[1] : 0.;
	st->swProb_U_ratio = swProb_U_rec1>0 ? swProb_U_rec0/swProb_U_rec1: 0;

	double EErt_rec0 = Nnosep_rec[0] >0 ? (double) NJ2J_rec[0] / (double) Nnosep_rec[0] : 1.;
	double EErt_rec1 = Nnosep_rec[1] >0 ? (double) NJ2J_rec[1] / (double) Nnosep_rec[1] : 1.;
	st->EE_ratio = EErt_rec1 > 0 ? EErt_rec0/EErt_rec1 : 0.;

	if(verbose > 1){
		printf( "Wave UE rate, %f, Mo UE rate %f \n" , st->findrate , (double)Nfndmo/(double)Nfndmo_denom);
		printf( "Wave EU rate, %f, Mo EU rate %f \n" , st->seprate  , (double) Nsepmo/Nsepmo_denom);
		printf( "Wave EE rate, %f, Mo EE rate %f \n" , st->J2Jprob  , (double) NJ2Jmo/Nsepmo_denom);
	}

	for(ll=0;ll<Npaths;ll++){
		for(ji=0;ji<JJ;ji++)
			occwg[ll*JJ + ji ] =occsz[ll*JJ + ji ]>0 ?  occwg[ll*JJ + ji ] / (double)occsz[ll*JJ + ji ] : 0.;
	}
	if(verbose>1){// occsz is:
	    double occfr[JJ];double tot=0.;
        for(ji=0;ji<JJ;ji++) occfr[ji] = 0.;
	    for(ll=0;ll<Npaths;ll++){
	        for(ji=0;ji<JJ;ji++) occfr[ji] +=occsz[ ll*JJ+ji ];
	    }
        for(ji=0;ji<JJ;ji++) tot += occfr[ji] ;
        printf("Occupation sizes: \n");
        for(ji=0;ji<JJ;ji++) printf( "%f, " ,occfr[ji]/tot);
        printf(" \n");
	}

    double * w_stns = Nemp>0   ? malloc(sizeof(double) * Nemp  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave );
	double * w_stsw = NswSt>0  ? malloc(sizeof(double) * NswSt ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_EEsw = NswE>0   ? malloc(sizeof(double) * NswE  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );

	double * w_EEns = NJ2J>0   ? malloc(sizeof(double) * NJ2J  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_EUsw = NswU>0   ? malloc(sizeof(double) * NswU  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_EUns = Nspell>0 ? malloc(sizeof(double) * Nspell): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_UEsw = NswU>0   ? malloc(sizeof(double) * NswU  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_UEns = Nspell>0 ? malloc(sizeof(double) * Nspell): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );

	double * w_all = malloc(sizeof(double)*Npaths*Nsim*TT/Npwave);
	double **w_all_rec = malloc(sizeof(double*)*2);
	for(ri=0;ri<2;ri++) w_all_rec[ri] = malloc(sizeof(double)*Npaths*Nsim*TT/Npwave);

	double * wwv_UEsw    = NswU>0   ? malloc(sizeof(double) * NswU  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * wwv_UEoccsw = NswU>0   ? malloc(sizeof(double) * NswU  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );

	double * wwv_EEsw    = NswE>0   ? malloc(sizeof(double) * NswE  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * wwv_EEoccsw = NswE>0   ? malloc(sizeof(double) * NswE  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );


	double ** w_MVswrec = malloc(sizeof(double*)*2);
	double ** w_MVnsrec = malloc(sizeof(double*)*2);
	w_MVswrec[0] = malloc((NJ2J + Nspell*2)*sizeof(double) );
	w_MVswrec[1] = malloc((NJ2J + Nspell*2)*sizeof(double) );

	w_MVnsrec[0] = malloc((NJ2J + Nspell*2)*sizeof(double) );
	w_MVnsrec[1] = malloc((NJ2J + Nspell*2)*sizeof(double) );

	int idx_EUns =0, idx_UEns =0, idx_EEns=0, idx_EUsw=0,idx_UEsw=0,idx_EEsw=0;
	int idx_UEwvsw=0,idx_EEwvsw=0;
	int idx_stns =0, idx_stsw =0;
	int idx_all=0;
	int valid_w =0;int valid_EUE=0;

	int idx_MVswrec[2]; idx_MVswrec[0]=0; idx_MVswrec[1]=0;
	int idx_MVnsrec[2]; idx_MVnsrec[0]=0; idx_MVnsrec[1]=0;
	int idx_all_rec[2]; idx_all_rec[0]=0; idx_all_rec[1]=0;

	// in its current architecture, this can't be parallelized
//	#pragma omp parallel for private(ll,ti,i,wi,si,ji,jji) need to make idx' w/in each branch
	for(ll=0;ll<Npaths;ll++){
		int sw_spell=-1;
		double wlast =invlaid_wval, wnext=invlaid_wval, wlast_wave=invlaid_wval, wnext_wave=invlaid_wval,wlast_waveEU=invlaid_wval, sw_occwg=0, sw_occwgEU=0;
		double w_EU=0.;
		for(i=0;i<Nsim;i++){
			for(wi=0;wi<TT/Npwave;wi++){
				double w_EEsw_wi =0., w_EEns_wi=0., w_stsw_wi=0.,w_stns_wi=0.,w_UEsw_wi=0.,w_UEns_wi=0.,w_EUsw_wi=0.,w_EUns_wi=0.;
				double  wwv_EEsw_wi=0.,wwv_EEoccsw_wi =0., wwv_UEsw_wi=0., wwv_UEoccsw_wi=0. ;
				int    I_EEsw_wi =0,  I_EEns_wi=0 , I_stsw_wi=0 ,I_stns_wi=0 ,I_UEsw_wi=0 ,I_UEns_wi=0 ,I_EUsw_wi=0 ,I_EUns_wi=0 ;
				int    I_UEwvsw_wi=0,I_EEwvsw_wi=0;
				int    sep=0;
				int    ri_wv = 0;
				wlast=0;wnext=0;
				ti = wi*Npwave;

				ri_wv = recIndic[ gsl_vector_int_get(ht->Ahist[ll],ti)]==0? 0 : 1; // ri_wv = 0 means expansion, ri_wv=1 means recession
				if( wi>3 && wi< TT/Npwave-3 && Anan==1  ){
					int ri=0;
					for(ri=1;ri<Npwave*3+1;ri++)
						wlast += gg_get(ht->whist[ll], i, ti - ri) ;
					for(ri=0;ri<Npwave*3;ri++)
						wnext += gg_get( ht->whist[ll] ,i,ti+ri) ;
				}else if( ti>3 && ti< TT-Npwave || Anan==0 ){
					wlast = gg_get( ht->whist[ll] ,i,ti-1) + gg_get( ht->whist[ll] ,i,ti-2)+ gg_get( ht->whist[ll] ,i,ti-3);
					wnext = gg_get( ht->whist[ll] ,i,ti+Npwave) + gg_get( ht->whist[ll] ,i,ti+Npwave+1)+ gg_get( ht->whist[ll] ,i,ti+Npwave+2);
					wnext_wave = wnext;
					wlast_wave=wlast;
				}else{
					wlast=invlaid_wval;wnext=invlaid_wval;wlast_wave=invlaid_wval; wnext_wave=invlaid_wval;
				}
				if( ti>3 && ti< TT-Npwave ){
					wlast_wave = gg_get( ht->whist[ll] ,i,ti-1) + gg_get( ht->whist[ll] ,i,ti-2)+ gg_get( ht->whist[ll] ,i,ti-3);
					wnext_wave = gg_get( ht->whist[ll] ,i,ti+Npwave) + gg_get( ht->whist[ll] ,i,ti+Npwave+1)+ gg_get( ht->whist[ll] ,i,ti+Npwave+2);
				}else{
					wlast_wave=invlaid_wval; wnext_wave=invlaid_wval;
				}

				for(si=0;si<Npwave;si++){
					ti = si + wi * Npwave;
					ri_wv = recIndic[ gsl_vector_int_get(ht->Ahist[ll],ti)]==0? 0 : 1; //0 if expansion, 1 if recession


                    if(ti<TT-1 && ti>0){ //can't get transitions unless ti+1 defined
						if( wlast>invlaid_wval && wnext>invlaid_wval ) valid_w ++;
						else{
							continue;
						}
                        if( ggi_get(ht->uhist[ll],i,ti) == 0 && wlast>invlaid_wval && wnext >invlaid_wval ){

                        	if( ggi_get(ht->uhist[ll],i,ti+1) ==1 ){
                        		int sep = 1;
                        		w_EU = (wnext- wlast)/(wnext + wlast)*2.; //not yet sure if this will be a switch or not, so just store it for now.
                                sw_spell = ggi_get(ht->jhist[ll],i,ti-1);
                                wlast_waveEU = wlast_wave;
                                sw_occwgEU = occwg[ll*JJ+ggi_get(ht->jhist[ll],i,ti-1)];
                            }else{
                                if( ggi_get(ht->J2Jhist[ll],i,ti) ==1 ){

                                    if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti-1)   ){
                                        w_EEsw_wi = (wnext- wlast)/(wnext + wlast)*2;
	                                    wwv_EEsw_wi  = wnext_wave + wlast_wave>0? (wnext_wave - wlast_wave )/(wnext_wave + wlast_wave)*2: 0. ;
	                                    wwv_EEoccsw_wi = (occwg[ll*JJ+ggi_get(ht->jhist[ll],i,ti+1)] -occwg[ll*JJ+ggi_get(ht->jhist[ll],i,ti-1)])/
			                                    (occwg[ll*JJ+ggi_get(ht->jhist[ll],i,ti+1)] +occwg[ll*JJ+ggi_get(ht->jhist[ll],i,ti-1)])*2;
                                        if(!gsl_finite(wwv_EEoccsw_wi)) wwv_EEoccsw_wi=0.;
	                                    I_EEwvsw_wi = wnext_wave + wlast_wave>0 ? 1 : 0;
	                                    I_EEsw_wi = 1;
                                    }else{
                                        w_EEns_wi = (wnext- wlast)/(wnext + wlast)*2;
                                        I_EEns_wi = 1;
                                    }
                                }else{  // not EE
                                    if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti-1) ){
                                        w_stsw_wi = (wnext- wlast)/(wnext + wlast)*2;
                                        I_stsw_wi = 1;

                                    }else{
                                        w_stns_wi = (wnext- wlast)/(wnext + wlast)*2;
                                        I_stns_wi = 1;
                                    }
                                }
                            }
                        // unemployed (and I saw you separate)
                        }else if(ggi_get(ht->uhist[ll],i,ti) ==1 && wlast> invlaid_wval && wnext > invlaid_wval ){ //unemployed
							if( ggi_get(ht->uhist[ll],i,ti+1) ==0 ) valid_EUE ++;
                            // find a job (and I saw you separate)
							if( ggi_get(ht->uhist[ll],i,ti+1) ==0  && sw_spell>-1  ){
								// switch occupations
                                if( ggi_get(ht->jhist[ll],i,ti+1) != sw_spell){
                                    w_UEsw_wi = (wnext- wlast)/(wnext + wlast)*2;
                                    w_EUsw_wi = w_EU;
                                    wwv_UEsw_wi = wnext_wave + wlast_waveEU>0? (wnext_wave - wlast_waveEU)/(wnext_wave + wlast_waveEU)*2: 0.;
	                                wwv_UEoccsw_wi = (occwg[ll*JJ+ggi_get(ht->jhist[ll],i,ti+1)] +sw_occwg) >0 ?
	                                		(occwg[ll*JJ+ggi_get(ht->jhist[ll],i,ti+1)] - sw_occwgEU )/
	                                                 (occwg[ll*JJ+ggi_get(ht->jhist[ll],i,ti+1)] +sw_occwgEU)*2 : 0.;
	                                I_UEwvsw_wi = wnext_wave + wlast_waveEU>0 && gsl_finite(wwv_UEoccsw_wi) ? 1:0;
                                    I_UEsw_wi = 1;
                                    I_EUsw_wi = 1;

                                }else{
                                    w_UEns_wi = (wnext- wlast)/(wnext + wlast)*2;
                                    w_EUns_wi = w_EU;
                                    I_UEns_wi = 1;
                                    I_EUns_wi = 1;
                                }
                            }
                        }
                    }
				}
				if( wi>3 && wi< TT/Npwave-3 && Anan==1 && wnext>0. && wlast>0. && ggi_get(ht->uhist[ll],i,ti) ==0  ){
					double wchng = (wnext - wlast)/(wnext + wlast)*2;
					if(isfinite(wchng)){
                        w_all[idx_all] = wchng;
                        idx_all++;
                        w_all_rec[ri_wv][idx_all_rec[ri_wv]] = wchng;
                        idx_all_rec[ri_wv]++;
                    }
				}
				if( I_EUns_wi==1 || I_EUsw_wi==1 || I_EEns_wi==1 || I_EEsw_wi==1 || sep==1 ){
					// can't change a job and be a stayer
					I_stsw_wi=0;I_stns_wi=0;
				}
				if( I_EEwvsw_wi ==1 && idx_EEwvsw < NswE && isfinite(wwv_EEsw_wi) && isfinite(wwv_EEoccsw_wi) ){
					wwv_EEsw[idx_EEwvsw] = wwv_EEsw_wi;
					wwv_EEoccsw[idx_EEwvsw]= wwv_EEoccsw_wi;
					idx_EEwvsw ++;
				}
				if(I_EEsw_wi ==1 && idx_EEsw < NswE && isfinite(w_EEsw_wi)){
					w_EEsw[idx_EEsw] =w_EEsw_wi;
					idx_EEsw ++;
					//now the cycle-specific dist
					w_MVswrec[ri_wv][idx_MVswrec[ri_wv]] = w_EEsw_wi;
					idx_MVswrec[ri_wv]++;
				}
				if(I_EEns_wi==1 && idx_EEns < NJ2J && isfinite(w_EEns_wi)){
					w_EEns[idx_EEns] = w_EEns_wi;
					idx_EEns ++;
					//now the cylce-specific dist
					w_MVnsrec[ri_wv][idx_MVnsrec[ri_wv]] = w_EEns_wi;
					idx_MVnsrec[ri_wv]++;
				}
				if(I_stsw_wi==1 && idx_stsw < NswSt && isfinite(w_stsw_wi)){
					w_stsw[idx_stsw] = w_stsw_wi;
					idx_stsw ++;
				}
				if(I_stns_wi==1 && idx_stns < Nemp && isfinite(w_stns_wi) ){
					w_stns[idx_stns] = w_stns_wi;
					idx_stns ++;
				}
				if( I_UEwvsw_wi==1 && idx_UEwvsw < NswU &&
				    isfinite(wwv_UEsw_wi) && isfinite(wwv_UEoccsw_wi)){
					wwv_UEsw[idx_UEwvsw] = wwv_UEsw_wi;
					wwv_UEoccsw[idx_UEwvsw] = wwv_UEoccsw_wi;
					idx_UEwvsw ++;
				}
				if(I_UEsw_wi==1 && idx_EUsw < NswU && isfinite(w_UEsw_wi) && isfinite(w_EUsw_wi)){
					w_UEsw[idx_UEsw] = w_UEsw_wi;
					w_EUsw[idx_EUsw] = w_EUsw_wi;

					idx_UEsw ++;
					idx_EUsw ++;

					w_MVswrec[ri_wv][idx_MVswrec[ri_wv]] = w_UEsw_wi;
					idx_MVswrec[ri_wv]++;
					w_MVswrec[ri_wv][idx_MVswrec[ri_wv]] = w_EUsw_wi;
					idx_MVswrec[ri_wv]++;
				}
				if(I_UEns_wi==1 && idx_EUns < Nspell && isfinite(w_UEns_wi) && isfinite(w_EUns_wi)){
					w_UEns[idx_UEns] = w_UEns_wi;
					w_EUns[idx_EUns] = w_EUns_wi;
					idx_UEns ++;
					idx_EUns ++;

					w_MVnsrec[ri_wv][idx_MVnsrec[ri_wv]] = w_UEns_wi;
					idx_MVnsrec[ri_wv]++;
					w_MVnsrec[ri_wv][idx_MVnsrec[ri_wv]] = w_EUns_wi;
					idx_MVnsrec[ri_wv]++;
				}
			} //wi waves
		} // i indiv
	} // ll paths

	if(verbose>1){
		printf("Uneployment rate is %f\n",st->unrate);
		printf("w_stns is allocated %d, needs %d\n",Nemp,idx_stns);
		printf("w_stsw is allocated %d, needs %d\n",NswSt,idx_stsw);
		printf("w_EEsw is allocated %d, needs %d\n",NswE,idx_EEsw);
		printf("w_EEns is allocated %d, needs %d\n",NJ2J,idx_EEns);
		printf("w_EUsw is allocated %d, needs %d\n",NswU,idx_EUsw);
		printf("w_EUns is allocated %d, needs %d\n", Nspell,idx_EUns);
		printf("w_UEsw is allocated %d, needs %d\n",NswU,idx_UEsw);
		printf("w_UEns is allocated %d, needs %d\n",Nspell,idx_UEns);
		printf("The number of recession vs expansion is %d, %d \n", idx_all_rec[1],idx_all_rec[0]);

	}


	/* Doing this sorting in w_qtls
	gsl_sort(w_EUns, 1,idx_EUns);
	gsl_sort(w_EUsw, 1,idx_EUsw);

	gsl_sort(w_stns, 1,idx_stns);
	gsl_sort(w_stsw, 1,idx_stsw);

	gsl_sort(w_UEns, 1,idx_UEns);
	gsl_sort(w_UEsw, 1,idx_UEsw);

	gsl_sort(w_EEns, 1,idx_EEns);
	gsl_sort(w_EEsw, 1,idx_EEsw);

	for(ri=0;ri<2;ri++){
		gsl_sort(w_MVswrec[ri],1,idx_MVswrec[ri]);
		gsl_sort(w_MVnsrec[ri],1,idx_MVnsrec[ri]);
	}*/
	int nonnum=0;
	for(ri=0;ri<idx_EEwvsw;ri++){
		if(!gsl_finite(wwv_EEsw[ri]) ) nonnum ++;
	}
	if( nonnum >0 )
		printf("Not a number on wwv_EEsw %d times \n ", nonnum);
	nonnum =0;
	for(ri=0;ri<idx_EEwvsw;ri++){
		if(!gsl_finite(wwv_EEoccsw[ri]) ) nonnum ++;
	}
	if( nonnum >0 )
		printf("Not a number on wwv_EEoccsw %d times \n ", nonnum);
	nonnum =0;
	for(ri=0;ri<idx_UEwvsw;ri++){
		if(!gsl_finite(wwv_UEoccsw[ri]) ) nonnum ++;
	}
	if( nonnum >0 )
		printf("Not a number on wwv_UEoccsw %d times \n ", nonnum);
	nonnum =0;
	for(ri=0;ri<idx_UEwvsw;ri++){
		if(!gsl_finite(wwv_UEsw[ri]) ) nonnum ++;
	}
	if( nonnum >0 )
		printf("Not a number on wwv_UEsw %d times \n ", nonnum);


	st->corrEE_wgocc  = idx_EEwvsw>2 ? gsl_stats_correlation( wwv_EEsw,1,wwv_EEoccsw ,1,(size_t)idx_EEwvsw) : 0.;
	st->corrEUE_wgocc = idx_UEwvsw>2 ? gsl_stats_correlation( wwv_UEsw,1,wwv_UEoccsw ,1,(size_t)idx_UEwvsw) : 0.;
	if(verbose>1){
	//	printf("The valid w observations are %d, valid EUE observations are %d \n",valid_w, valid_EUE);
		printf("EE  corr of wage and occwg is %f, based on %d obs \n", st->corrEE_wgocc , idx_EEwvsw);
		printf("EUE corr of wage and occwg is %f, based on %d obs \n", st->corrEUE_wgocc, idx_UEwvsw);
	}

	w_qtls(w_stns,1,idx_stns,qtlgrid,st->stns_qtls);
	w_qtls(w_stsw,1,idx_stsw,qtlgrid,st->stsw_qtls);
	w_qtls(w_EEns,1,idx_EEns,qtlgrid,st->EEns_qtls);
	w_qtls(w_EEsw,1,idx_EEsw,qtlgrid,st->EEsw_qtls);
	w_qtls(w_EUns,1,idx_EUns,qtlgrid,st->EUns_qtls);
	w_qtls(w_EUsw,1,idx_EUsw,qtlgrid,st->EUsw_qtls);
	w_qtls(w_UEns,1,idx_UEns,qtlgrid,st->UEns_qtls);
	w_qtls(w_UEsw,1,idx_UEsw,qtlgrid,st->UEsw_qtls);

	if(cf_now==1) {
        w_qtls(w_all, 1, idx_all,qtlgrid, st->all_qtls);
        w_qtls(w_all, 1, idx_all, edgeqtls, st->edge_qtls);
    }
	st->var_wg = gsl_stats_variance(w_all, 1, idx_all);

	for(ri=0;ri<2;ri++){
		st->var_wg_rec[ri] = gsl_stats_variance(w_all_rec[ri],1,idx_all_rec[ri]);
        if(cf_now==1) {
            w_qtls(w_all_rec[ri], 1, idx_all_rec[ri], qtlgrid, st->all_qtls_rec[ri]);
            w_qtls(w_all_rec[ri], 1, idx_all_rec[ri], edgeqtls, st->edge_qtls_rec[ri]);
        }
	}
	double MVsw_qtls[2][Nqtls];
	double MVns_qtls[2][Nqtls];

	for(ri=0;ri<2;ri ++){
		w_qtls(w_MVswrec[ri],1,idx_MVswrec[ri],edgeqtls,MVsw_qtls[ri]);
		w_qtls(w_MVnsrec[ri],1,idx_MVnsrec[ri],edgeqtls,MVns_qtls[ri]);
	}
	for(i=0;i<Nqtls;i++){
		st->MVsw_qtls_ratio[i] = MVsw_qtls[1][i]-MVsw_qtls[0][i];
		st->MVns_qtls_ratio[i] = MVns_qtls[1][i]-MVns_qtls[0][i];
	}

	for(ri=0;ri<2;ri++){
		free(w_MVnsrec[ri]);free(w_MVswrec[ri]);
	}
	free(w_MVnsrec);free(w_MVswrec);
	for(ji=0;ji<JJ;ji++)
		free(occ_gflow[ji]);
	free(occ_gflow);
    for(ji=0;ji<JJ;ji++)
        free(occ_gflowE[ji]);
    free(occ_gflowE);
    for(ji=0;ji<JJ;ji++)
        free(occ_gflowU[ji]);
    free(occ_gflowU);
    for(ji=0;ji<JJ;ji++)
        free(occ_sijflow[ji]);
    free(occ_sijflow);
    for(ji=0;ji<JJ;ji++)
        free(occ_sijflowE[ji]);
    free(occ_sijflowE);
    for(ji=0;ji<JJ;ji++)
        free(occ_sijflowU[ji]);
    free(occ_sijflowU);
	for(ri=0;ri<2;ri++) free(w_all_rec[ri]);free(w_all_rec);
	free(w_stns);
	free(w_stsw);
	free(w_EEns);
	free(w_EEsw);
	free(w_EUns);
	free(w_EUsw);
	free(w_UEns);
	free(w_UEsw);
	free(wwv_EEoccsw);free(wwv_UEoccsw); free(wwv_EEsw); free(wwv_UEsw);
}

void set_xpspace( double *x, struct cal_params *par ){
// this is the flip of set_params that goes from x to params. This takes params and loads an x
    int ii, i,ji;
    ii =0;

    int ci = par->cluster_hr;
    if( ci ==0 || ci == Ncluster){

        x[ii] = par->alpha0     ;ii++;
        x[ii] = par->lambdaUS0  ;ii++;
        x[ii] = par->lambdaUM0  ;ii++;
        x[ii] = par->lambdaES0  ;ii++;
        x[ii] = par->lambdaEM0  ;ii++;
        x[ii] = par->delta_avg  ;ii++;
        x[ii] = par->zloss      ;ii++;
        x[ii] = par->alphaE1    ;ii++;
        x[ii] = par->alphaU1    ;ii++;

        // pair-wise flow rate adjustments
        x[ii] = par->alpha_nf[1][0]; ii++;
        x[ii] = par->alpha_nf[0][1]; ii++;
        x[ii] = par->alpha_nf[0][2]; ii++;
        x[ii] = par->alpha_nf[0][3]; ii++;
        //x[ii] = par->alpha_nf[0][1] ;ii++;; //normalization of the first one? 0.;
        //x[ii] = par->alpha_nf[0][2] ;ii++;
        //x[ii] = par->alpha_nf[0][3] ;ii++;
        //x[ii] = par->alpha_nf[1][2] ;ii++;
        //x[ii] = par->alpha_nf[1][3] ;ii++;
        //x[ii] = par->alpha_nf[2][3] ;ii++;

        ii = Npar_cluster[0] ;

    }
    if(ci==1 || ci == Ncluster){

        x[ii] = par->update_z   ;ii++;
        x[ii] = par->scale_z    ;ii++;
        x[ii] = par->shape_z    ;ii++;

        x[ii] = par->var_pe     ;ii++;
        x[ii] = par->autop      ;ii++;
        x[ii] = par->var_ae     ;ii++;
        x[ii] = x[ii] = par->autoa      = x[ii];ii++;
        x[ii] = par->gdfthr     ;ii++;
        x[ii] = par->stwupdate  ;ii++;
        //par->scale_eps  = x[ii];ii++;
        //par->shape_eps  = x[ii];ii++;
        x[ii] = par->var_eps    ;ii++;
        if(eps_2emg==1){
            x[ii] = par->lshape_eps ;ii++;
            x[ii] = par->rshape_eps ;ii++;
        }
        ii = Npar_cluster[0] + Npar_cluster[1];
    }
    if(ci==2 || ci == Ncluster){
        x[ii] = par->delta_Acoef    ;ii++;
        x[ii] = par->lambdaEM_Acoef ;ii++;
        x[ii] = par->lambdaES_Acoef ;ii++;
        x[ii] = par->lambdaUM_Acoef ;ii++;
        x[ii] = par->lambdaUS_Acoef ;ii++;
        x[ii] = par->zloss_Acoef    ;ii++;
        x[ii] = par->z_Acoef        ;ii++;
        x[ii] = par->z_Amag         ;ii++;
        x[ii] = par->eps_Acoef      ;ii++;
        x[ii] = par->eps_Amag       ;ii++;

        i=10;
        for(ji=0;ji<JJ;ji++){
            if(ji ==0) // norm the first occupation
                par->AloadP->data[ji]= 0.;
            x[ii] = par->AloadP->data[ji];
            if(ji>0) ii++;
        }

    }
}

void read_params(char* name, struct cal_params *par){
	int ii, i, ji;
	int rstatus;
	FILE * parfile;
	parfile = fopen(name,"r");
	double dd;char sd[32];char*endptr;

	int pi =0 ;
    strcpy(sd,"\000");
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->alpha0 = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->lambdaUS0  = dd;
	par->xparopt[pi] = dd; pi++;
    rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
    par->lambdaUM0  = dd;
    par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->lambdaES0 = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->lambdaEM0 = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->delta_avg= dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->zloss=     dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->alphaE1=   dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->alphaU1=   dd;
	par->xparopt[pi] = dd; pi++;
	for(ii=0;ii<JJ;ii++) {
        rstatus = fscanf(parfile, "%s, ", sd);
        dd = strtod(sd, &endptr);
        par->alpha_nf[0][ii] = dd;
        par->xparopt[pi] = dd;
        pi++;
    }
	for(ii=0;ii<JJ;ii++){
	    for(ji=1;ji<JJ;ji++) par->alpha_nf[ji][ii] = ii!=ji? par->alpha_nf[0][ii] : 0.;
	}
    par->alpha_nf[0][0] =0.;
//		for(ji=ii+1;ji<JJ;ji++){
//			rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
//			par->alpha_nf[ii][ji]= dd;
//			par->xparopt[pi] = dd; pi++;
//		}
//	}
	for(ii=0;ii<JJ;ii++){
		for(ji=ii+1;ji<JJ;ji++)
			par->alpha_nf[ji][ii] = -par->alpha_nf[ii][ji];
	}
	// alpha scales, scale on net flow of switching from l to d
	for(ii=0;ii<JJ;ii++) par->alpha_nf[ii][ii] = 0.; //should never be adjusting the diagonal

	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->update_z = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->scale_z  = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->shape_z  = dd;
	par->xparopt[pi] = dd; pi++;

	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->var_pe = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->autop  = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->var_ae = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->autoa  = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->gdfthr = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->stwupdate= dd;
	par->xparopt[pi] = dd; pi++;
	// epsilon :
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	//par->scale_eps =dd;
	//par->xparopt[pi] = dd; pi++;
	//rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,NULL);
	//par->shape_eps =dd;
	//par->xparopt[pi] = dd; pi++;

	par->var_eps= dd;
	par->xparopt[pi] = dd; pi++;

	if(eps_2emg==1){
		rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
		par->lshape_eps= dd;
		par->xparopt[pi] = dd; pi++;
		rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
		par->rshape_eps= dd;
		par->xparopt[pi] = dd; pi++;
	}
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->delta_Acoef   = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->lambdaEM_Acoef= dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->lambdaES_Acoef= dd;
	par->xparopt[pi] = dd; pi++;
    rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
    par->lambdaUM_Acoef = dd;
    par->xparopt[pi] = dd; pi++;
    rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->lambdaUS_Acoef = dd;
	par->xparopt[pi] = dd; pi++;
	rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
	par->zloss_Acoef = dd;
	par->xparopt[pi] = dd; pi++;

    rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
    par->z_Acoef = dd;
    par->xparopt[pi] = dd; pi++;
    rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
    par->z_Amag = dd;
    par->xparopt[pi] = dd; pi++;
    rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
    par->eps_Acoef = dd;
    par->xparopt[pi] = dd; pi++;
    rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
    par->eps_Amag = dd;
    par->xparopt[pi] = dd; pi++;


	for(ji=0;ji<JJ;ji++){
		if(ji>=1){
			rstatus = fscanf(parfile,"%s, ",sd);dd = strtod(sd,&endptr);
			if(rstatus ==0 )
				rstatus = fscanf(parfile,"%s\n",sd);dd = strtod(sd,&endptr);
			par->AloadP->data[ji] = dd;
			par->xparopt[pi] = dd; pi++;
		}

	}

	fclose(parfile);
}

void print_targets(char *name, struct stats * st, struct stats * dat){
    FILE * tgtfile;
    tgtfile = fopen(name,"w+");
    int i;

    fprintf(tgtfile," ,model,data \n ");

    double dat_dur = dat->udur_sw / dat->udur_nosw;
    double mod_dur = gsl_finite(st->udur_nosw) && gsl_finite(st->udur_sw) && st->udur_nosw>0 ?
                     st->udur_sw / st->udur_nosw : 0.;

    // weighting duration and randomness low (0.1)
    fprintf(tgtfile,"J2J,%f,%f \n",st->J2Jprob , dat->J2Jprob);
    fprintf(tgtfile,"UE,%f,%f \n",st->findrate , dat->findrate);
    fprintf(tgtfile,"EU,%f,%f \n",st->seprate , dat->seprate);
    fprintf(tgtfile,"sw Prob|EE,%f,%f \n", st->swProb_EE , dat->swProb_EE);
    fprintf(tgtfile,"sw Prob|EU,%f,%f \n",st->swProb_U , dat->swProb_U);
    fprintf(tgtfile,"sw Prob|st,%f,%f \n",st->swProb_st , dat->swProb_st) ;
    fprintf(tgtfile,"U duration,%f,%f \n",mod_dur , dat_dur) ;
    fprintf(tgtfile,"var flows E,%f,%f \n",st->varGflowsE , dat->varGflowsE);
    fprintf(tgtfile,"var flows U,%f,%f \n",st->varGflowsU , dat->varGflowsU);

    // -------------------------------------
    for(i=0;i<JJ;i++)
        fprintf(tgtfile,"net flow to %d,%f,%f \n",i,
                st->occ_margflow[i],dat->occ_margflow[i]);
    //for(i=0;i<nflows;i++)
    //    fprintf(tgtfile,"net flow %d,%f,%f \n",i,
    //    st->occ_netflow[i] , dat->occ_netflow[i]);

    fprintf(tgtfile,"swProb U ratio,%f,%f \n",st->swProb_U_ratio , dat->swProb_U_ratio);
    fprintf(tgtfile,"swProb E ratio,%f,%f \n",st->swProb_EE_ratio , dat->swProb_EE_ratio) ;
    fprintf(tgtfile,"UE ratio ,%f,%f \n",st->fndrt_ratio , dat->fndrt_ratio) ;
    fprintf(tgtfile,"EU ratio,%f,%f \n",st->seprt_ratio , dat->seprt_ratio) ;
    fprintf(tgtfile,"EE ratio,%f,%f \n",st->EE_ratio , dat->EE_ratio) ;
    //err_vec[ii + 5] = (st->corr_wgocc - dat->corr_wgocc) * 2 / (st->corr_wgocc + dat->corr_wgocc);
    //put them all in one file if we're at the optimal
    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->stns_qtls, Nqtls, 0);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->stns_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->EEns_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->EEns_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->EUns_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->EUns_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->UEns_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->UEns_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->MVns_qtls_ratio, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->MVns_qtls_ratio, Nqtls, 1);

    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->stsw_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->stsw_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->EEsw_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->EEsw_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->EUsw_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->EUsw_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->UEsw_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->UEsw_qtls, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", st ->MVsw_qtls_ratio, Nqtls, 1);
    printarray("qtls_moddatEEEUUEMV_nssw.csv", dat->MVsw_qtls_ratio, Nqtls, 1);

    fclose(tgtfile);
}


void print_params(double *x, int n, struct cal_params * par){
	int ii, i,ji;

	ii =0;
	parhist = fopen(parhi_f,"a+");

	fprintf(parhist,"%8.6f, ", par->alpha0);

	fprintf(parhist,"%8.6f, ", par->lambdaUS0);
    fprintf(parhist,"%8.6f, ", par->lambdaUM0);
	fprintf(parhist,"%8.6f, ", par->lambdaES0);
	fprintf(parhist,"%8.6f, ", par->lambdaEM0);
	fprintf(parhist,"%8.6f, ", par->delta_avg);
	fprintf(parhist,"%8.6f, ", par->zloss    );
	fprintf(parhist,"%8.6f, ", par->alphaE1  );
	fprintf(parhist,"%8.6f, ", par->alphaU1  );
    fprintf(parhist,"%8.6f, ", par->alpha_nf[1][0]);
	for(ji=1;ji<JJ;ji++){
		fprintf(parhist,"%8.6f, ", par->alpha_nf[0][ji]);
	}
	fprintf(parhist,"%8.6f, ", par->update_z );
	fprintf(parhist,"%8.6f, ", par->scale_z  );
	fprintf(parhist,"%8.6f, ", par->shape_z  );

	fprintf(parhist,"%8.6f, ", par->var_pe   );
	fprintf(parhist,"%8.6f, ", par->autop    );
	fprintf(parhist,"%8.6f, ", par->var_ae   );
	fprintf(parhist,"%8.6f, ", par->autoa    );
	fprintf(parhist,"%8.6f, ", par->gdfthr   );
	fprintf(parhist,"%8.6f, ", par->stwupdate);
	fprintf(parhist,"%8.6f, ", par->var_eps  );

	if(eps_2emg==1){
		fprintf(parhist,"%8.6f, ", par->lshape_eps);
		fprintf(parhist,"%8.6f, ", par->rshape_eps);
	}
	//fprintf(parhist,"%8.6f, ", par->scale_eps);
	//fprintf(parhist,"%8.6f, ", par->shape_eps);


	fprintf(parhist,"%8.6f, ", par->delta_Acoef   );
	fprintf(parhist,"%8.6f, ", par->lambdaEM_Acoef);
	fprintf(parhist,"%8.6f, ", par->lambdaES_Acoef);
	fprintf(parhist,"%8.6f, ", par->lambdaUS_Acoef);
    fprintf(parhist,"%8.6f, ", par->lambdaUM_Acoef);
	fprintf(parhist,"%8.6f, ", par->zloss_Acoef   );
    fprintf(parhist,"%8.6f, ", par->z_Acoef       );
    fprintf(parhist,"%8.6f, ", par->z_Amag        );
    fprintf(parhist,"%8.6f, ", par->eps_Acoef     );
    fprintf(parhist,"%8.6f, ", par->eps_Amag      );

    for(ji=0;ji<JJ;ji++){
		if(ji>=1){
			fprintf(parhist,"%8.6f", par->AloadP->data[ji]);
            if(ji<JJ-1) fprintf(parhist,",");
		}
	}
	fprintf(parhist,"\n");

	fclose(parhist);


}


void set_params( double * x, int n, struct cal_params * par,int ci){
	int ii, i,ji;
	ii =0;

	if( ci ==0 || ci == Ncluster){

		par->alpha0     = x[ii];ii++;
		par->lambdaUS0  = x[ii];ii++;
        par->lambdaUM0  = x[ii];ii++;
		par->lambdaES0 = x[ii];ii++;
		par->lambdaEM0 = x[ii];ii++;
		par->delta_avg = x[ii];ii++;
		par->zloss     = x[ii];ii++;
		par->alphaE1    = x[ii];ii++;
		par->alphaU1    = x[ii];ii++;

        // regardless of origin, get this as the destination bonus
        for(i=0;i<JJ;i++) par->alpha_nf[i][0] = x[ii];
        ii++;
        for(i=0;i<JJ;i++) par->alpha_nf[i][1] = x[ii];
        ii++;
        for(i=0;i<JJ;i++) par->alpha_nf[i][2] = x[ii];
        ii++;
        for(i=0;i<JJ;i++) par->alpha_nf[i][3] = x[ii];
        ii++;
		/* pair-wise flow rate adjustments
		par->alpha_nf[0][1]   = x[ii];ii++;; //normalization of the first one? 0.;
		par->alpha_nf[0][2]   = x[ii];ii++;
		par->alpha_nf[0][3]   = x[ii];ii++;
		par->alpha_nf[1][2]   = x[ii];ii++;
		par->alpha_nf[1][3]   = x[ii];ii++;
		par->alpha_nf[2][3]   = x[ii];ii++;
        // make the matrix symmetric
		for(ii=0;ii<JJ;ii++){
			for(ji=ii+1;ji<JJ;ji++)
				par->alpha_nf[ji][ii] = -par->alpha_nf[ii][ji];
		}
        */
		// alpha scales, scale on net flow of switching from l to d
		for(i=0;i<JJ;i++) par->alpha_nf[i][i] = 0.; //should never be adjusting the diagonal
		ii = Npar_cluster[0] ;

	}
	if(ci==1 || ci == Ncluster){

		par->update_z  = x[ii];ii++;
		par->scale_z    = x[ii];ii++;
		par->shape_z    = x[ii];ii++;

		par->var_pe     = x[ii];ii++;
		par->autop      = x[ii];ii++;
		par->var_ae     = x[ii];ii++;
		par->autoa      = x[ii];ii++;
		par->gdfthr     = x[ii];ii++;
		par->stwupdate  = x[ii];ii++;
		//par->scale_eps  = x[ii];ii++;
		//par->shape_eps  = x[ii];ii++;
		par->var_eps    = x[ii];ii++;
		if(eps_2emg==1){
			par->lshape_eps = x[ii];ii++;
			par->rshape_eps = x[ii];ii++;
		}
		ii = Npar_cluster[0] + Npar_cluster[1];
	}
	if(ci==2 || ci == Ncluster){
		par->delta_Acoef     = x[ii];ii++;
		par->lambdaEM_Acoef  = x[ii];ii++;
		par->lambdaES_Acoef  = x[ii];ii++;
		par->lambdaUS_Acoef  = x[ii];ii++;
        par->lambdaUM_Acoef  = x[ii];ii++;
		par->zloss_Acoef     = x[ii];ii++;
		par->z_Acoef         = x[ii];ii++;
        par->z_Amag          = x[ii];ii++;
        par->eps_Acoef       = x[ii];ii++;
        par->eps_Amag        = x[ii];ii++;

		for(ji=0;ji<JJ;ji++){
			if(ji ==0) // norm the first occupation
				par->AloadP->data[ji]= 0.;
			else par->AloadP->data[ji] = x[ii];
			if(ji>0) ii++;
		}
	}
	if(verbose>1){
		printf("Set %d parameters at vector: \n(", ii);
		for(i=0;i<Npar_cluster[ci]-1;i++)
			printf("%8s,", parnames[ci][i]);
		i = Npar_cluster[ci]-1;
		printf("%8s)\n(",parnames[ci][i]);

		for(i=0;i<Npar_cluster[ci]-1;i++)
			printf("%8.3f,", x[i]);
		i = Npar_cluster[ci]-1;
		printf("%8.3f)\n", x[i]);
	}
	// These introduced additional cross-derivatives based on alpha1 that I don't want.
}


void set_dat( struct stats * dat){

	int i;

	dat->J2Jprob = 0.03404871;
	dat->findrate = 0.3946638;
	dat->seprate  =  0.02227534;
	dat->swProb_EE = 0.2685099; // 2 digit: 0.5026808;
	dat->swProb_U = 0.2891623; // 2 digit: 0.5362486;
	dat->swProb_st = 0.01066545; // 2digit: 0.01999708;
	dat->udur_nosw=  4.781546;
	dat->udur_sw=  5.09012;

	dat-> occ_ten = 10.43*12;
	dat->corrEE_wgocc = 0.3872327;
	dat->corrEUE_wgocc = 0.52;

    double stsw[] = {-0.35297299, -0.12816599,  0.01820793,  0.19858560,  0.48457859 };
    memcpy(dat->stsw_qtls ,stsw,Nqtls*sizeof(double));
    double stns[] = {-0.196993647, -0.066443456, -0.001802062,  0.086715528,  0.242205102 };
    memcpy(dat->stns_qtls ,stns,Nqtls*sizeof(double));

    double EUsw[] = {-2.2304440, -1.2123117, -0.5063948, 0.2349499,  1.2261816};
    memcpy(dat->EUsw_qtls ,EUsw,Nqtls*sizeof(double));
    double EUns[] = {-1.9071055, -1.0858343, -0.4433699,  0.1309994,  0.9952223};
    memcpy(dat->EUns_qtls ,EUsw,Nqtls*sizeof(double));

    double UEsw[] = {-1.1040582, -0.4722635,  0.1518163,  1.0317821,  2.1701233};
    memcpy(dat->UEsw_qtls, UEsw,Nqtls*sizeof(double));
    double UEns[] = {-1.060892873, -0.489657380,  0.004472182,  0.742813794,  1.706427430};
    memcpy(dat->UEns_qtls, UEns,Nqtls*sizeof(double));

    double EEsw[] = {-0.6148813, -0.2173941,  0.1263684,  0.6126602,  1.3109891};
	memmove(dat->EEsw_qtls, EEsw,Nqtls*sizeof(double));
    double EEns[] = {-0.51813528, -0.18758653,  0.06988849,  0.41210494,  0.95124483};
    memmove(dat->EEns_qtls, EEns,Nqtls*sizeof(double));


    dat->seprt_ratio = 0.7460479;
    dat->fndrt_ratio =  1.087637;
    dat->EE_ratio =  1.184628;

    dat->swProb_EE_ratio = 1./0.9337642; // these were computed backwards: needs to be higher in expansion and we are computing ratio as expansion/recession
    dat->swProb_U_ratio  =  1./0.9035369;

	dat->doubleswU = 0.2968578;
	dat->doubleswE = 0.2771865;

    double datMVsw_qtls_ratio[] = {-0.5745752 ,-0.4187927 ,-0.1884432 ,-0.4010306 ,-0.3847628 };
    double datMVns_qtls_ratio[] = {-0.1175533 ,-0.1703967 ,-0.1478998, -0.2348585, -0.1054350};

    memcpy(dat->MVns_qtls_ratio,datMVns_qtls_ratio,Nqtls*sizeof(double));
    memcpy(dat->MVsw_qtls_ratio,datMVsw_qtls_ratio,Nqtls*sizeof(double));

    // Occupation numbers 1: NRC, 2: RC, 3: NRM, 4:RM

    // row major- order of net flows matrix. i.e. flow from 1 to 2, 1 to 3, 1 to 4, 2 to 3, 2 to 4, 3 to 4
    double netflows[] = {-0.003475171, -0.004571549, -0.0008167612, -0.010466888, 0.0050908438, 0.0105925984};
    double nrmflows = 0;
    for(i=0;i<nflows;i++) nrmflows = netflows[i] >0 ? netflows[i]+nrmflows : -netflows[i]+nrmflows;
	for(i=0;i<nflows;i++) netflows[i] = netflows[i]/nrmflows;
	memcpy(dat->occ_netflow, netflows, nflows*sizeof(double));
        dat->nrmflows = nrmflows; dat->nrmflowsE = 0.;dat->nrmflowsU = 0.;
	double netflowsU[] = {-0.004461320, -0.003989968, 0.0009252794 , -0.008512734, 0.0044061444, 0.0055699451};
	double netflowsE[] = {-0.003888313, -0.004349818, -0.001930180, -0.013183258,  0.005698791, 0.014134908};
    for(i=0;i<nflows;i++) dat->nrmflowsU = netflowsU[i] >0 ? netflowsU[i]+dat->nrmflowsU : -netflowsU[i]+dat->nrmflowsU;
	for(i=0;i<nflows;i++) dat->nrmflowsE = netflowsE[i] >0 ? netflowsE[i]+dat->nrmflowsE : -netflowsE[i]+dat->nrmflowsE;

	memcpy(dat->occ_netflowU,netflowsU, sizeof(double)*nflows);
	memcpy(dat->occ_netflowE,netflowsE, sizeof(double)*nflows);

	double varnetflowsU = gsl_stats_variance(netflowsU,1,nflows);
	double varnetflowsE = gsl_stats_variance(netflowsE,1,nflows);
    dat->varGflowsE = 0.0223444648;
    dat->varGflowsU = 0.021788583;

    double datoccmarg[] = { 0.1790146, 0.3378744, 0.2322478, 0.2508631};
    memcpy(dat->occ_margflow,datoccmarg,sizeof(double )*JJ);
    if(verbose>=2) {
        printf(" \n +++++++++++++++ \n data varnetflowsU: %f  data varnetflowsE: %f \n",
               varnetflowsU, varnetflowsE);
        printf(" \n +++++++++++++++ \n data absdevnetflowsU: %f  data absdevnetflowsE: %f \n",
               dat->nrmflowsU, dat->nrmflowsE);

        printf(" \n +++++++++++++++ \n data varGflowsU: %f  data varnGflowsE: %f \n",
               dat->varGflowsU, dat->varGflowsE);
    }
}

double param_dist( double * x, struct cal_params *par , int Npar, double * err_vec , int Nerr){

	// Takes in a parameters vector and pumps out an error vector

	int i,ii,ji,ai;
	int success;

	struct valfuns vf;
	struct polfuns pf;
	struct shocks sk;
	struct hists ht;

	struct stats st;
	struct stats dat;

	alloc_qtls(&dat);

	// set data moments and parameter values
    set_dat(&dat);

    set_params( x, Npar, par, par->cluster_hr);
	print_params(x, Npar, par);
	if(par->shape_z ==0.){
	    //this seemed to happen for some reason?
	    printf("+++++++++++++++++++++++++++++++++++++++++++\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++\n");
        printf("Got 0 shape_z (illegal). x[0,1,2] = (%f,%f,%f) \n",x[0],x[1],x[2]);
        printf("+++++++++++++++++++++++++++++++++++++++++++\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++\n");
        for(i=0;i<Nerr;i++) err_vec[i] = ((double)(rand()%100) )/ 100. ;
        return NAN;
	}

	alloc_qtls(&st);
	allocate_mats(&vf,&pf,&ht,&sk);

	init_pf(&pf,par);
	init_vf(&vf,par);
	//set the markov transition matrix for P
	rouwenhorst(par->autop,pow(par->var_pe,0.5),par->Ptrans[0],par->Plev);
	for(i=0;i<JJ;i++){
		gsl_matrix_memcpy(par->Ptrans[i],par->Ptrans[i]);
	}
	// set the markov transition matrix for A
	if(NA>2){
		rouwenhorst(par->autoa,pow(par->var_ae,0.5),par->Atrans,par->Alev);
	}else{
		double probRec = .21; double probExp=(1.-probRec);
		double qAtrans = (probExp/probRec+par->autoa)/(probExp/probRec+1);
		double pAtrans = par->autoa + 1. - qAtrans;
		double zAlev   = pow( par->var_ae ,0.5);
		gg_set(par->Atrans,0,0,   pAtrans);gg_set(par->Atrans,0,1,1.-pAtrans);
		gg_set(par->Atrans,1,0,1.-qAtrans);gg_set(par->Atrans,1,1,   qAtrans);
		gsl_vector_set(par->Alev,0,-zAlev);gsl_vector_set(par->Alev,1,zAlev);
	}

	// setup the z matrix
	double zlev1[NZ]; double zprob1[NZ];
	success = disc_Weibull( zprob1,zlev1,NZ, 0., par->scale_z, par->shape_z );
	double zprobsum =0.;
	for(ii=0;ii<NZ;ii++){
		gsl_vector_set( par->zlev,ii,zlev1[ii] );
		zprobsum += zprob1[ii];
		for(ai=0;ai<NA;ai++)
			gsl_vector_set( par->zprob[ai],ii,zprob1[ii] );
	}
    convprobs(par->zprob[0],par->zlev,par->z_Acoef,par->z_Amag);

	if(eps_2emg==1){
		//success = disc_Weibull( par->epsprob->data,par->epslev->data, NE, 0.,par->scale_eps,par->shape_eps );
        for(ai=0;ai<NA;ai++) success = disc_2emg(par->epsprob[ai]->data,par->epslev->data,(int)par->epsprob[ai]->size,
				0.,par->var_eps,par->lshape_eps,par->rshape_eps);
        // convolve with cyclical adjustment for recession:
        convprobs(par->epsprob[0],par->epslev,par->eps_Acoef,par->eps_Amag);
	}else{
		success =0;
		double epsub[NE-1];

		epsub[0] =gsl_cdf_gaussian_Pinv( 1./(double)NE,par->var_eps );
		for(ii=1;ii<NE-1;ii++)
			epsub[ii] =gsl_cdf_gaussian_Pinv( ((double)ii+1.)/(double)NE,par->var_eps );
		for(ii=1;ii<NE-1;ii++){
			par->epslev->data[ii] =(epsub[ii] + epsub[ii-1])/2;
		}
		par->epslev->data[0] = epsub[0] - (par->epslev->data[1]);
		par->epslev->data[NE-1] = 2*epsub[NE-2] -par->epslev->data[NE-2];
		for(ai=0;ai<NA;ai++) gsl_vector_set_all(par->epsprob[ai],1./(double)NE);
        convprobs(par->epsprob[0],par->epslev,par->eps_Acoef,par->eps_Amag);
	}
	if(success > 0){
		printf(" Did not compute the distribution properly");
	}
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	for(i=0;i<NS;i++) gsl_vector_set(par->xSlev,i, .010984* (double)i/(double) (NS-1)); // really just 0 or 1--- matches 1% increase due to employer tenure in KM2009
	double occ_size_dat[] = {0.2636133, 0.3117827, 0.1493095, 0.2752945};
	memcpy( par->jprob->data, occ_size_dat,sizeof(occ_size_dat));


	// ensure don't voluntarily quit, except for when z is the lowest value
	double w_hr;
	for(ji=0;ji<JJ;ji++){
		for (ii = 0; ii < NN; ii++) {
			int ai = ii / (NP * NS * NZ * NE);
			int pi = (ii - ai * NP * NS * NZ * NE) / (NS * NZ * NE);
			int si = (ii - ai * NP * NS * NZ * NE - pi * NS * NZ * NE) / (NZ * NE);
			int zi =
					(ii - ai * NP * NS * NZ * NE - pi * NS * NZ * NE - si * NZ * NE) / NE;
			int ti = ii - ai * NP * NS * NZ * NE - pi * NS * NZ * NE - si * NZ * NE -
			         zi * NE;

			w_hr = exp(par->AloadP->data[ji] * par->Alev->data[ai] +
			                          par->Plev->data[pi] +
			                          par->epslev->data[ti] +
			                          par->zlev->data[zi] +
			                          par->xSlev->data[si] +
			                          occ_wlevs[ji]) + wage_lev;
			if( w_hr< b && w_hr>0 && zi>0) b = w_hr;
		}
	}
	if(verbose>0){
		printf("Lowest wage (and b) is:%f \n", b);
	}



// print out the grids

	if(print_lev>1 && par->rank==0){
		printvec("Alev.csv", par->Alev,0);printvec("zlev.csv", par->zlev,0);
		printmat("Atrans.csv",par->Atrans,0);printmat("Ptrans.csv",par->Ptrans[0],0);
		printvec("zlev.csv", par->zlev,0);
		printvec("Plev.csv", par->Plev,0);
		printvec("epsprob.csv", par->epsprob[0],0);
        printvec("epsprob.csv", par->epsprob[1],1);
		printvec("epslev.csv", par->epslev,0);
	}

	success = draw_shocks(&sk);
	if(verbose>2 && success != 0) printf("Problem drawing shocks\n");
	if(nosolve==0) {
		success = sol_dyn(par, &vf, &pf, &sk);
		if (verbose > 2 && success != 0) printf("Problem solving model\n");

		success = sim(par, &vf, &pf, &ht, &sk);
		if (verbose > 2 && success != 0) printf("Problem simulating model\n");

		success = sum_stats(par, &vf, &pf, &ht, &sk, &st);

		//form error vector

		ii = 0;
		double dat_dur = dat.udur_sw / dat.udur_nosw;
		double mod_dur = gsl_finite(st.udur_nosw) && gsl_finite(st.udur_sw) && st.udur_nosw>0 ?
				st.udur_sw / st.udur_nosw : 0.;
		if (par->cluster_hr == 0 || par->cluster_hr == Ncluster) {
		    // weighting duration and randomness low (0.1)
			err_vec[ii] = (st.J2Jprob - dat.J2Jprob) * 2 / (st.J2Jprob + dat.J2Jprob);
			ii++;
			err_vec[ii] = (st.findrate - dat.findrate) * 2 / (st.findrate + dat.findrate);
			ii++;
			err_vec[ii] = (st.seprate - dat.seprate) * 2 / (st.seprate + dat.seprate);
			ii++;
			err_vec[ii] = (st.swProb_EE - dat.swProb_EE) * 2 / (st.swProb_EE + dat.swProb_EE);
			ii++;
			err_vec[ii] = (st.swProb_U - dat.swProb_U) * 2 / (st.swProb_U + dat.swProb_U);
			ii++;
			err_vec[ii] = (st.swProb_st - dat.swProb_st) * 2 / (st.swProb_st + dat.swProb_st);
			ii++;
			err_vec[ii] = 0.1*(mod_dur - dat_dur) *2 / (mod_dur + dat_dur);
			ii ++;
			// -------------------------------------
			// Should this go in 1st cluster or second?
            err_vec[ii]   = 0.5*(st.varGflowsE - dat.varGflowsE)/(st.varGflowsE + dat.varGflowsE);
			//err_vec[ii] = .5*(st.corrEE_wgocc  - dat.corrEE_wgocc); //*2 / (st.corr_wgocc + dat.corr_wgocc);
			//err_vec[ii] = (st.doubleswE - dat.doubleswE) /(st.doubleswE + dat.doubleswE);
			ii++;
            err_vec[ii]   = 0.5*(st.varGflowsU - dat.varGflowsU)/(st.varGflowsU + dat.varGflowsU);
			//err_vec[ii] = (st.doubleswU - dat.doubleswU) /(st.doubleswU + dat.doubleswU);
			//err_vec[ii] = .5*(st.corrEUE_wgocc - dat.corrEUE_wgocc); //*2 / (st.corr_wgocc + dat.corr_wgocc);
			ii++;
			// -------------------------------------
			//for(i=0;i<nflows;i++)
			//	err_vec[ii+i] = (st.occ_netflow[i] - dat.occ_netflow[i])
			//			/(double)nflows; // was atan weight all the net flows as one target
		    for(i=0;i<JJ;i++)
		        err_vec[ii+i] = (st.occ_margflow[i]-dat.occ_margflow[i])  ;

			ii = Ntgt_cluster[0];
		}


		if (par->cluster_hr == 1 || par->cluster_hr == Ncluster) {

			for (i = 0; i < Nqtls; i++)
				err_vec[ii + i] =
						(st.stns_qtls[i] - dat.stns_qtls[i]) / (double) Nqtls;// *2/(st.stns_qtls[i]+dat.stns_qtls[i])
			ii += Nqtls;
			for (i = 0; i < Nqtls; i++)
				err_vec[ii + i] =
						(st.stsw_qtls[i] - dat.stsw_qtls[i]) / (double) Nqtls;// *2/(st.stsw_qtls[i]+dat.stsw_qtls[i])
			ii += Nqtls;
			for (i = 0; i < Nqtls; i++)
				err_vec[ii + i] =
						(st.EEns_qtls[i] - dat.EEns_qtls[i]) / (double) Nqtls;// *2/(st.EEns_qtls[i]+dat.EEns_qtls[i])
			ii += Nqtls;
			for (i = 0; i < Nqtls; i++)
				err_vec[ii + i] =
						(st.EEsw_qtls[i] - dat.EEsw_qtls[i]) / (double) Nqtls;// *2/(st.EEsw_qtls[i]+dat.EEsw_qtls[i])
			ii += Nqtls;
			for (i = 0; i < Nqtls; i++)
				err_vec[ii + i] =
						(st.EUns_qtls[i] - dat.EUns_qtls[i]) / (double) Nqtls;// *2/(st.EUns_qtls[i]+dat.EUns_qtls[i])
			ii += Nqtls;
			for (i = 0; i < Nqtls; i++)
				err_vec[ii + i] =
						(st.EUsw_qtls[i] - dat.EUsw_qtls[i]) / (double) Nqtls;// *2/(st.EUsw_qtls[i]+dat.EUsw_qtls[i])
			ii += Nqtls;
			for (i = 0; i < Nqtls; i++)
				err_vec[ii + i] =
						(st.UEns_qtls[i] - dat.UEns_qtls[i]) / (double) Nqtls;// *2/(st.UEns_qtls[i]+dat.UEns_qtls[i])
			ii += Nqtls;
			for (i = 0; i < Nqtls; i++)
				err_vec[ii + i] =
						(st.UEsw_qtls[i] - dat.UEsw_qtls[i]) / (double) Nqtls;// *2/(st.UEsw_qtls[i]+dat.UEsw_qtls[i])
			ii += Nqtls;
		}
		if (par->cluster_hr == 2 || par->cluster_hr == Ncluster) {
			err_vec[ii + 0] = dat.swProb_U_ratio < 1. ?  -(st.swProb_U_ratio - dat.swProb_U_ratio) / (dat.swProb_U_ratio- 1.) :
			        (st.swProb_U_ratio - dat.swProb_U_ratio) / (dat.swProb_U_ratio- 1.);
			err_vec[ii + 1] = dat.swProb_EE_ratio <1. ? -(st.swProb_EE_ratio - dat.swProb_EE_ratio) / (dat.swProb_EE_ratio- 1.) :
                              (st.swProb_EE_ratio - dat.swProb_EE_ratio)  / ( dat.swProb_EE_ratio- 1.) ;
			err_vec[ii + 2] = dat.fndrt_ratio  < 1. ? - (st.fndrt_ratio - dat.fndrt_ratio)  / (dat.fndrt_ratio- 1.) :
			        (st.fndrt_ratio - dat.fndrt_ratio)  / (dat.fndrt_ratio- 1.) ;
			err_vec[ii + 3] = dat.seprt_ratio < 1. ? -(st.seprt_ratio - dat.seprt_ratio)  / (dat.seprt_ratio- 1.) :
			        (st.seprt_ratio - dat.seprt_ratio)  / (dat.seprt_ratio- 1.) ;
			err_vec[ii + 4] = dat.EE_ratio < 1. ?  -(st.EE_ratio - dat.EE_ratio) / ( dat.EE_ratio- 1.) :
                              (st.EE_ratio - dat.EE_ratio) / ( dat.EE_ratio- 1.) ;

			ii +=5;
			for(i=0;i<Nqtls; i++)
				err_vec[ii+i] = (st.MVsw_qtls_ratio[i] - dat.MVsw_qtls_ratio[i])/(double)Nqtls;
			ii += Nqtls;
			for(i=0;i<Nqtls; i++)
				err_vec[ii+i] = (st.MVns_qtls_ratio[i] - dat.MVns_qtls_ratio[i])/(double)Nqtls;
		}

		if (print_lev >= 2) {
			printarray("stns_qtls.csv", st.stns_qtls, Nqtls, 0);
			printarray("stsw_qtls.csv", st.stsw_qtls, Nqtls, 0);
			printarray("EEns_qtls.csv", st.EEns_qtls, Nqtls, 0);
			printarray("EEsw_qtls.csv", st.EEsw_qtls, Nqtls, 0);
			printarray("EUns_qtls.csv", st.EUns_qtls, Nqtls, 0);
			printarray("EUsw_qtls.csv", st.EUsw_qtls, Nqtls, 0);
			printarray("UEns_qtls.csv", st.UEns_qtls, Nqtls, 0);
			printarray("UEsw_qtls.csv", st.UEsw_qtls, Nqtls, 0);

		}
        if(opt_print==1) {
            //put them all in one file if we're at the optimal
            printarray("qtls_stEEEUUEMV_nssw.csv", st.stns_qtls, Nqtls, 0);
            printarray("qtls_stEEEUUEMV_nssw.csv", st.EEns_qtls, Nqtls, 1);
            printarray("qtls_stEEEUUEMV_nssw.csv", st.EUns_qtls, Nqtls, 1);
            printarray("qtls_stEEEUUEMV_nssw.csv", st.UEns_qtls, Nqtls, 1);
            printarray("qtls_stEEEUUEMV_nssw.csv", st.MVns_qtls_ratio, Nqtls, 1);

            printarray("qtls_stEEEUUEMV_nssw.csv", st.stsw_qtls, Nqtls, 1);
            printarray("qtls_stEEEUUEMV_nssw.csv", st.EEsw_qtls, Nqtls, 1);
            printarray("qtls_stEEEUUEMV_nssw.csv", st.EUsw_qtls, Nqtls, 1);
            printarray("qtls_stEEEUUEMV_nssw.csv", st.UEsw_qtls, Nqtls, 1);
            printarray("qtls_stEEEUUEMV_nssw.csv", st.MVsw_qtls_ratio, Nqtls, 1);
            print_targets("targetfit.csv",&st,&dat);
            printarray("alpha_nf.csv",glb_par->alpha_nf[0],JJ,0);
            for(i=1;i<JJ;i++)printarray("alpha_nf.csv",glb_par->alpha_nf[i],JJ,1);

        }

	}else{
		gsl_rng * RNG = gsl_rng_alloc(gsl_rng_default);
		for(i=0;i<Nerr;i++){
			err_vec[i] = gsl_rng_uniform(RNG);
		}
		gsl_rng_free(RNG);

	}
	int nanerr_flag =0;
    for(i=0;i<Nerr;i++){
        if(gsl_finite(err_vec[i])!=1)
            nanerr_flag = 1;
    }
    if(nanerr_flag ==1){
        printf("NaN error at: \n");
        printf(" Target vector \n ");
        //for(i=0;i<Nerr;i++)
        //    printf("%8s,", tgtnames[Ncluster][i]);
        //printf("\n");
        for(i=0;i<Nerr;i++)
            printf("%8.3f,",err_vec[i]);
        printf("parameter vector, length %d: \n(", ii);
        int ci =par->cluster_hr;
        for(i=0;i<Npar_cluster[ci]-1;i++)
            printf("%8s,", parnames[ci][i]);
        i = Npar_cluster[ci]-1;
        printf("%8s)\n(",parnames[ci][i]);

        for(i=0;i<Npar_cluster[ci]-1;i++)
            printf("%8.3f,", x[i]);
        i = Npar_cluster[ci]-1;
        printf("%8.3f)\n", x[i]);
        printf("The other stats are: \n");
        printf("unrate %f, UE rt %f, EU rt %f, EE rt %f", st.unrate,st.findrate,st.seprate, st.J2Jprob);
        printf("swProbU %f, swProbE %f, swProbst %f", st.swProb_U,st.swProb_EE,st.swProb_st);


    }


    double quad_dist =0;
	for(i=0;i<Nerr;i++)
		quad_dist += err_vec[i]*err_vec[i];

	if(print_lev>2){
		char zname[21],epsname[21];
		sprintf(zname,"zlev_%06.3f.csv",quad_dist);
		sprintf(epsname,"epslev_%06.3f.csv",quad_dist);
		printvec(zname, par->zlev,0);
		printvec(epsname, par->zlev,0);
	}

	free_mats(&vf,&pf,&ht,&sk);

	free_qtls(&st); free_qtls(&dat);
	return(quad_dist);

}

void compute_derivs(double**cal_Jacobian, struct cal_params * par0, int printJac ){

    int xi,yi;
    double * x_pspace  = malloc(sizeof(double)*Nparams);
    double * xp_pspace = malloc(sizeof(double)*Nparams); //move up
    double * xm_pspace = malloc(sizeof(double)*Nparams); //move down
    double * x_del = malloc(sizeof(double)*Nparams);

    par0->cluster_hr = Ncluster;
    set_xpspace(x_pspace,par0);
    double err0[Ntargets];
    param_dist( x_pspace,par0, Nparams, err0, Ntargets );
    int pl_old = print_lev;
    print_lev =0;

    for(xi=0;xi<Nparams;xi++){
        if(verbose>1)
            printf("computing derivatives for %d", xi);
        int xxi;
        for(xxi=0;xxi<Nparams;xxi++){ //set the baseline
            xp_pspace[xxi]=x_pspace[xxi];
            xm_pspace[xxi]=x_pspace[xxi];
        }

        xp_pspace[xi] = gsl_min( x_pspace[xi] +
                0.005*( par0->param_ub[xi]- par0->param_lb[xi]),
                par0->param_ub[xi]);
        xm_pspace[xi] = gsl_max( x_pspace[xi] -
                                 0.005*( par0->param_ub[xi]- par0->param_lb[xi]),
                                 par0->param_lb[xi]);
        x_del[xi] = xp_pspace[xi] - xm_pspace[xi];
        double err_p[Ntargets];double err_m[Ntargets];
        param_dist(xp_pspace,par0,Nparams,err_p,Ntargets);
        param_dist(xm_pspace,par0,Nparams,err_m,Ntargets);

        for(yi=0;yi<Ntargets;yi++)
            cal_Jacobian[xi][yi] = (err_p[yi] - err_m[yi])/x_del[xi];
    }
    if(printJac==1){
        FILE* matfile;matfile = fopen("cal_Jacobian.csv", "w");
        for(xi=0;xi<Nparams;xi++){
            fprintf(matfile, "%s, ",parnames_cluall[xi]);
            for(yi=0;yi<Ntargets-1;yi++)
                fprintf(matfile,"%f, ",cal_Jacobian[xi][yi]);
            fprintf(matfile,"%f\n",cal_Jacobian[xi][Ntargets-1]);
        }
        fclose(matfile);
    }
    for(xi=0;xi<Nparams;xi++)
        free(cal_Jacobian[xi]);
    free(cal_Jacobian);
    print_lev = pl_old;

    free(x_pspace);free(xp_pspace);free(xm_pspace);free(x_del);


}

double f_wrapper_nlopt(unsigned n, const double * x, double * grad, void * par){
	// this is going to wrap param_dist

	int i,mv,nv;
	int verbose_old,print_lev_old;
	struct cal_params*   cpar = (struct cal_params*)par;
	double dist;
	double * errvec = malloc(sizeof(double)*Ntgt_cluster[cpar->cluster_hr] );
	int ci = cpar->cluster_hr;
	if(verbose>1) printf("Entering Nelder-Meade evaluation\n");
	verbose_old = verbose;
	print_lev_old = print_lev;

	verbose=0;
	print_lev=0;

	nv = n;
	mv = Ntgt_cluster[cpar->cluster_hr];

	double * x_pspace = malloc(sizeof(double)*n);
	for(i=0;i<n;i++) // convert from [0,1] domain to parameter space
		x_pspace[i] = x[i] * (cpar->cluster_ub[ci][i]-cpar->cluster_lb[ci][i])
				+cpar->cluster_lb[ci][i];

	if(verbose_old>1){
		printf("eval at: ");
		for(i=0;i<n-1;i++)
			printf("%f,",x_pspace[i]);
		printf("%f\n",x_pspace[n-1]);
	}

	dist = param_dist( x_pspace, cpar, (int) n ,  errvec, Ntgt_cluster[cpar->cluster_hr] );

	verbose=verbose_old;
	print_lev = print_lev_old;

    if(print_lev>1){
	    int chr =0;
	    calhist = fopen(calhi_f,"a+");
	    fprintf(calhist,"%f,",dist);
	    if(ci<Ncluster) {
		    for (i = 0; i < ci; i++)chr += Ntgt_cluster[i];
		    for(i=0;i<chr;i++) fprintf(calhist," , ");
	    }

	    for(i=0;i<mv;i++)
		    fprintf(calhist,"%f,",errvec[i]);
	    for(i=(mv+chr);i<Ntargets;i++) fprintf(calhist," , ");

	    chr = 0; for(i=0;i<ci;i++) chr += Npar_cluster[i];
	    if(ci<Ncluster){
		    for(i=0;i<chr;i++) fprintf(calhist," , ");
	    }
	    for(i=0;i<nv-1;i++)
		    fprintf(calhist,"%f,",x_pspace[i]);
	    fprintf(calhist,"%f\n",x_pspace[nv-1]);
	    fclose(calhist);
    }
	if(  dist <solver_state[0]){
		solver_state[0] = dist;
		for(i=0;i<n;i++) solver_state[i+1] = x_pspace[i];
	}

	free(x_pspace);free(errvec);
	return dist;

}

void fix_wagedists(struct stats *dat, struct stats *mod, struct hists * ht ){
    // renorms the wage changes to be equal to the data.

}


void dfovec_iface_(double * errvec, double * x, int * n, int* m){
	// this is going to interface with dfovec.f and call param_dist using the global params
	int nv = *n;
	int mv = *m;
	int i,ci;
	int verbose_old,print_lev_old;
	int cf_old;
	double dist;

	if(verbose>1) printf("Entering DFBOLS evaluation. First param %f\n",  x[0]);

	//for some reason, this may get an all 0:
	double xneq0 =0;
	for(i=0;i<nv;i++) xneq0 += x[i];
	if(xneq0 < 1e-6 && xneq0 >-1e-6 ) {
	    return;
	}
	verbose_old = verbose;
	print_lev_old = print_lev;
    cf_old = cf_now;
    cf_now = 1;
	verbose=0;
	print_lev=0;
	ci = glb_par->cluster_hr;
	double * x_pspace = malloc(sizeof(double)*nv);
	for(i=0;i<nv;i++) // convert from [0,1] domain to parameter space
		x_pspace[i] = x[i] * (glb_par->cluster_ub[ci][i]-glb_par->cluster_lb[ci][i])
				+glb_par->cluster_lb[ci][i];
//	if(verbose_old>1){
//		printf("Evaluating a DFBOLS point at (param space): \n");
//		for(i=0;i<nv-1;i++)
//			printf("%f, ", x_pspace[i]);
//		printf("%f\n", x_pspace[nv-1]);
//	}

	dist = param_dist( x_pspace, glb_par, nv ,errvec, mv );
	//dist =  (double)(rand() % 1000) /1000.;


	verbose=verbose_old;
	print_lev = print_lev_old;
    cf_now = cf_old;
	if(print_lev>1){
		int chr =0;
		calhist = fopen(calhi_f,"a+");
		fprintf(calhist,"%f,",dist);
		if(ci<Ncluster) {
			for (i = 0; i < ci; i++)chr += Ntgt_cluster[i];
			for(i=0;i<chr;i++) fprintf(calhist," , ");
		}

		for(i=0;i<mv;i++)
			fprintf(calhist,"%f,",errvec[i]);
		for(i=(mv+chr);i<Ntargets;i++) fprintf(calhist," , ");

		chr = 0; for(i=0;i<ci;i++) chr += Npar_cluster[i];
		if(ci<Ncluster){
			for(i=0;i<chr;i++) fprintf(calhist," , ");
		}
		for(i=0;i<nv-1;i++)
			fprintf(calhist,"%f,",x_pspace[i]);
		fprintf(calhist,"%f\n",x_pspace[nv-1]);
		fclose(calhist);
	}
	if(verbose>1){
	    printf(" Error is: \n ");
	    for(i=0;i<mv;i++)
            printf("%f,",errvec[i]);

	}
	if(  dist <solver_state[0]){
		solver_state[0] = dist;
		for(i=0;i<nv;i++) solver_state[i+1] = x_pspace[i];
	}

	free(x_pspace);
}
void convprobs(gsl_vector* probs, gsl_vector * levs, double coef,double magn){
    int ii;
    int N = probs->size;
    double probs0[N],lev0,levN;
    for(ii=0;ii<N;ii++) probs0[ii] = probs->data[ii];
    lev0 = levs->data[0];levN = levs->data[N-1];
    double mixlevs[N];
    for(ii=0;ii<N;ii++) mixlevs[ii] = levs->data[ii];
    double mixprobs[N];
    for(ii=0;ii<N;ii++) mixprobs[ii] = coef>= 0 ? pow((mixlevs[ii] - lev0)/(levN-lev0), coef) :
            1.-pow((mixlevs[ii] - lev0)/(levN-lev0), -coef);
    double distint =0.;
    for( ii=0;ii<N;ii++) distint += mixprobs[ii];
    for( ii=0;ii<N;ii++) mixprobs[ii] = mixprobs[ii]/distint;
    for(ii=0;ii<N;ii++)
        probs->data[ii] = probs0[ii]*(1.-magn) + magn * mixprobs[ii];
    //printf("made probs\n");
}


int draw_shocks(struct shocks * sk){
	int ji,i,ti,ll;

	int seed;

	#pragma omp parallel for private(ll, ti,i,ji)
	for(ll=0;ll<Npaths;ll++){
		gsl_rng * rng0 = gsl_rng_alloc( gsl_rng_default );
		seed = 12281951 + 10*omp_get_thread_num();
		gsl_rng_set( rng0, seed );

		for(ti=0;ti<TTT;ti++){
			gsl_vector_set(sk->Asel[ll],ti, gsl_rng_uniform(rng0));
			for(i=0;i<Nsim;i++){
				gg_set( sk->zsel[ll]      ,i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->epssel[ll]     ,i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->lambdaEMsel[ll],i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->lambdaESsel[ll],i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->lambdaUSsel[ll],i,ti, gsl_rng_uniform(rng0) );
                gg_set( sk->lambdaUMsel[ll],i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->xSsel[ll]     ,i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->xGsel[ll]     ,i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->jsel[ll]      ,i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->dsel[ll]      ,i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->msel[ll]      ,i,ti, gsl_rng_uniform(rng0) );
			}
			for(ji=0;ji<JJ;ji++){
				gg_set(sk->Psel[ll] , ti,ji, gsl_rng_uniform(rng0));
			}
		}
		gsl_rng_free(rng0);
	}

	return 0;
}

int w_qtls( double* vec, int stride, int len, double * qtlgrid, double * qtls_out){
    int qi ,ii,li;
    int status;
    status =0;
    gsl_sort(vec,(size_t)stride, (size_t)len);
    if(len<=Nqtls){
        status++;
        if(len>2){
            double obsgrid[len];for(qi=0;qi<len;qi++)obsgrid[qi] = (double)qi/((double)len-1);
            gsl_interp * qtlinterp = gsl_interp_alloc(gsl_interp_linear,len);
            gsl_interp_accel* qtlaccl = gsl_interp_accel_alloc();
            gsl_interp_init(qtlinterp,obsgrid, vec,len);
            for(qi=0;qi<Nqtls;qi++)
                qtls_out[qi] = gsl_interp_eval(qtlinterp, obsgrid,vec,qtlgrid[qi],qtlaccl);
            gsl_interp_free(qtlinterp);gsl_interp_accel_free(qtlaccl);

        }
        else{ //nothing do to
            for(qi=0;qi<Nqtls;qi++)
                qtls_out[qi]=0.;
        }
        return status;
    }else {
        for (qi = 0; qi < Nqtls; qi++) {
            qtls_out[qi] = gsl_stats_quantile_from_sorted_data(vec, (size_t) stride, (size_t) len, qtlgrid[qi]);
        }

        //check the order (i.e. must be increasing):
        gsl_sort(qtls_out, 1, Nqtls);
        // check for NaNs:
        for (qi = 0; qi < Nqtls; qi++) {
            if (!gsl_finite(qtls_out[qi])) {
                if (verbose > 0) printf(" NaN qtl! ");
                status++;
                if (qi == 0)
                    qtls_out[qi] = gsl_finite(qtls_out[qi + 1]) ? qtls_out[qi + 1] : vec[0];
                else if (qi == Nqtls - 1)
                    qtls_out[qi] = gsl_finite(qtls_out[qi - 1]) ? qtls_out[qi - 1] : vec[len - 1];
            }
        }
        //first pass over. Now do it again, just in case:
        for (qi = 0; qi < Nqtls; qi++) {
            if (gsl_finite(qtls_out[qi]) == 0) {
                status++;
                double obsgrid[len];for(qi=0;qi<len;qi++)obsgrid[qi] = (double)qi/((double)len-1);
                gsl_interp * qtlinterp = gsl_interp_alloc(gsl_interp_linear,len);
                gsl_interp_accel* qtlaccl = gsl_interp_accel_alloc();
                gsl_interp_init(qtlinterp,obsgrid, vec,len);
                for(qi=0;qi<Nqtls;qi++)
                    qtls_out[qi] = gsl_interp_eval(qtlinterp, obsgrid,vec,qtlgrid[qi],qtlaccl);
                gsl_interp_free(qtlinterp);gsl_interp_accel_free(qtlaccl);
                break;
            }
        }

        return status;
    }
}


void allocate_pars( struct cal_params * par){
	int ji,ii;
	TTT = TT + burnin;
	par->Plev    = gsl_vector_calloc(NP) ;
	par->Alev    = gsl_vector_calloc(NA) ;
	par->xSlev   = gsl_vector_calloc(NS) ;
	par->zlev    = gsl_vector_calloc(NZ) ;
	par->epslev    = gsl_vector_calloc(NE) ;

	par->AloadP = gsl_vector_calloc(JJ);
	par->Ptrans = malloc(sizeof(gsl_matrix*)*JJ);
	par->param_lb = malloc(sizeof(double)*Nparams);
	par->param_ub = malloc(sizeof(double)*Nparams);
	for(ji=0;ji<JJ;ji++){
		par-> Ptrans[ji] = gsl_matrix_calloc(NP,NP);
	}
	par->Atrans = gsl_matrix_calloc(NA,NA);
	par->xStrans = gsl_matrix_calloc(NS,NS);
//	par->ztrans = gsl_matrix_calloc(NZ,NZ);
    par->epsprob = malloc(sizeof(gsl_vector*)*NA);
    for(ii=0;ii<NA;ii++)
        par->epsprob[ii] = gsl_vector_calloc(NE);
	par->zprob = malloc(sizeof(gsl_vector*)*NA);
	for(ii=0;ii<NA;ii++)
		par->zprob[ii] = gsl_vector_calloc(NZ);
	par->jprob = gsl_vector_calloc(JJ);
	par->alpha_nf = malloc(sizeof(double*)*JJ);
	for(ji=0;ji<JJ;ji++){
		par->alpha_nf[ji] = malloc(sizeof(double)*JJ);
		for(ii=0;ii<JJ;ii++)
			par->alpha_nf[ji][ii] = 0.;
	}

}


void allocate_mats( struct valfuns * vf, struct polfuns * pf, struct hists * ht, struct shocks * sk ){

	int j;

	NN = NA*NP*NS*NZ*NE;
	NUN = NA*NP*NS*NZ;
	TTT = TT + burnin;

	alloc_valfuns(vf);
	alloc_hists(ht);
	alloc_shocks(sk);
	pf -> mE =  gsl_matrix_calloc( NN,JJ );
	pf -> mU =  gsl_matrix_calloc( NUN,JJ );
	pf -> sE = malloc(sizeof(gsl_matrix*)*JJ);
	pf -> sU = malloc(sizeof(gsl_matrix*)*JJ);
	for (j = 0; j <JJ ; j++) {
		pf->sE[j] = gsl_matrix_calloc(NN,JJ);
		pf->sU[j] = gsl_matrix_calloc(NUN,JJ);
	}
}

void alloc_valfuns(struct valfuns *vf ){

	vf -> WE =  gsl_matrix_calloc( NN,JJ );
	vf -> WU =  gsl_matrix_calloc( NUN,JJ );
	vf -> RE =  gsl_matrix_calloc( NN,JJ );
	vf -> RU =  gsl_matrix_calloc( NUN,JJ );

	vf -> WEdist = gsl_matrix_calloc( NN,JJ );

}

void alloc_shocks(struct shocks * sk){
	int j;

	sk->Asel       = malloc(sizeof(gsl_vector*)*Npaths);
	sk->lambdaEMsel = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->lambdaESsel = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->lambdaUSsel = malloc(sizeof(gsl_matrix*)*Npaths);
    sk->lambdaUMsel = malloc(sizeof(gsl_matrix*)*Npaths);
    sk->Psel       = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->epssel      = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->xGsel      = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->xSsel      = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->zsel       = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->jsel       = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->dsel       = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->msel       = malloc(sizeof(gsl_matrix*)*Npaths);

	for(j=0;j<Npaths;j++){
		sk->Asel[j] = gsl_vector_calloc(TTT);
		sk->Psel[j]       = gsl_matrix_calloc(TTT,JJ);
		sk->lambdaEMsel[j] = gsl_matrix_calloc(Nsim,TTT);
		sk->lambdaESsel[j] = gsl_matrix_calloc(Nsim,TTT);
		sk->lambdaUSsel[j] = gsl_matrix_calloc(Nsim,TTT);
        sk->lambdaUMsel[j] = gsl_matrix_calloc(Nsim,TTT);
		sk->epssel[j]     = gsl_matrix_calloc(Nsim,TTT);
		sk->xGsel[j]      = gsl_matrix_calloc(Nsim,TTT);
		sk->xSsel[j]      = gsl_matrix_calloc(Nsim,TTT);
		sk->zsel[j]       = gsl_matrix_calloc(Nsim,TTT);
		sk->jsel[j]       = gsl_matrix_calloc(Nsim,TTT);
		sk->dsel[j]       = gsl_matrix_calloc(Nsim,TTT);
		sk->msel[j]       = gsl_matrix_calloc(Nsim,TTT);
	}
}

void memcpy_shocks(struct shocks * sk_dest , struct shocks * sk_orig){
	int j;

	for(j=0;j<Npaths;j++){
		gsl_vector_memcpy(sk_dest->Asel[j],sk_orig->Asel[j]);
		gsl_matrix_memcpy(sk_dest->Psel[j],sk_orig->Psel[j]);
		gsl_matrix_memcpy(sk_dest->lambdaEMsel[j],sk_orig->lambdaEMsel[j]);
		gsl_matrix_memcpy(sk_dest->lambdaESsel[j],sk_orig->lambdaESsel[j]);
		gsl_matrix_memcpy(sk_dest->lambdaUSsel[j],sk_orig->lambdaUSsel[j]);
        gsl_matrix_memcpy(sk_dest->lambdaUMsel[j],sk_orig->lambdaUMsel[j]);
		gsl_matrix_memcpy(sk_dest->epssel[j],sk_orig->epssel[j]);
		gsl_matrix_memcpy(sk_dest->xGsel[j],sk_orig->xGsel[j]);
		gsl_matrix_memcpy(sk_dest->xSsel[j],sk_orig->xSsel[j]);
		gsl_matrix_memcpy(sk_dest->zsel[j],sk_orig->zsel[j]);
		gsl_matrix_memcpy(sk_dest->jsel[j],sk_orig->jsel[j]);
		gsl_matrix_memcpy(sk_dest->dsel[j],sk_orig->dsel[j]);
		gsl_matrix_memcpy(sk_dest->msel[j],sk_orig->msel[j]);
	}
}

void alloc_hists( struct hists *ht ){
	int ll,ji;

	ht->xShist = malloc(sizeof(gsl_matrix_int *)*Npaths);
    ht->zhist = malloc(sizeof(gsl_matrix_int *)*Npaths);
    ht->epshist = malloc(sizeof(gsl_matrix_int *)*Npaths);

	ht->uhist = malloc(sizeof(gsl_matrix_int *)*Npaths);
	ht->whist = malloc(sizeof(gsl_matrix      *)*Npaths);
	ht->jhist = malloc(sizeof(gsl_matrix_int *)*Npaths);
	ht->Ahist = malloc(sizeof(gsl_vector_int *)   *Npaths);
	ht->Phist = malloc(sizeof(gsl_matrix_int *)   *Npaths);
	ht->J2Jhist = malloc(sizeof(gsl_matrix_int *) *Npaths);
	ht->sijhist =  malloc(sizeof(gsl_matrix**)*Npaths);
	for(ll=0;ll<Npaths;ll++){
		ht->uhist[ll] = gsl_matrix_int_calloc(Nsim,TT);
		ht->whist[ll] = gsl_matrix_calloc(Nsim,TT);
		ht->xShist[ll] = gsl_matrix_int_calloc(Nsim,TT);
        ht->epshist[ll] = gsl_matrix_int_calloc(Nsim,TT);
        ht->zhist[ll] = gsl_matrix_int_calloc(Nsim,TT);
		ht->jhist[ll] = gsl_matrix_int_calloc(Nsim,TT);
		ht->Ahist[ll] = gsl_vector_int_calloc(TT);
		ht->Phist[ll] = gsl_matrix_int_calloc(Nsim,TT);
		ht->J2Jhist[ll] = gsl_matrix_int_calloc(Nsim,TT);
		ht->sijhist[ll] = malloc(sizeof(gsl_matrix*)*JJ);
		for(ji=0;ji<JJ;ji++)
			ht->sijhist[ll][ji] = gsl_matrix_calloc(Nsim,TT);
	}
}

void alloc_qtls( struct stats *st ){

	int ri;
	st->EEns_qtls = malloc(Nqtls*sizeof(double));
	st->EEsw_qtls = malloc(Nqtls*sizeof(double));
	st->EUns_qtls = malloc(Nqtls*sizeof(double));
	st->EUsw_qtls = malloc(Nqtls*sizeof(double));
	st->UEns_qtls = malloc(Nqtls*sizeof(double));
	st->UEsw_qtls = malloc(Nqtls*sizeof(double));
	st->stns_qtls = malloc(Nqtls*sizeof(double));
	st->stsw_qtls = malloc(Nqtls*sizeof(double));
	st->all_qtls = malloc(Nqtls*sizeof(double));
	st->all_qtls_rec = malloc(sizeof(double*)*2);
	for(ri=0;ri<2;ri++)st->all_qtls_rec[ri] = malloc(Nqtls*sizeof(double));
	st->MVns_qtls_ratio = malloc(Nqtls*sizeof(double));
	st->MVsw_qtls_ratio = malloc(Nqtls*sizeof(double));
	nflows = fact_int(JJ)/(fact_int(2)*fact_int(JJ-2));
	st->occ_netflow = malloc( nflows*sizeof(double) );

	st->occ_netflowE = malloc( nflows*sizeof(double) );
	st->occ_netflowU = malloc( nflows*sizeof(double) );
    st->occ_margflow = malloc( JJ*sizeof(double));

	st->edge_qtls = malloc(Nqtls*sizeof(double));
	st->edge_qtls_rec = malloc(2*sizeof(double*));
	for(ri=0;ri<2;ri++)st->edge_qtls_rec[ri] = malloc(Nqtls*sizeof(double));

}

void free_qtls( struct stats *st ){

	int ri;
	free(st->EEns_qtls);free(st->EEsw_qtls);
	free(st->EUns_qtls);free(st->EUsw_qtls);
	free(st->UEns_qtls);free(st->UEsw_qtls);
	free(st->stns_qtls);free(st->stsw_qtls);
	free(st->all_qtls);
	for(ri=0;ri<2;ri++)free(st->all_qtls_rec[ri]); free(st->all_qtls_rec);
	free(st->MVns_qtls_ratio);free(st->MVsw_qtls_ratio);
	free(st->occ_netflow);
	free(st->occ_netflowU);free(st->occ_netflowE);
	free(st->occ_margflow);
	free(st->edge_qtls);
	for(ri=0;ri<2;ri++)free(st->edge_qtls_rec[ri]);
	free(st->edge_qtls_rec);

}

void free_mats(struct valfuns * vf, struct polfuns * pf, struct hists *ht,struct shocks * sk){

	int j;

	gsl_matrix_free(pf->mU);
	gsl_matrix_free(pf->mE);
	for (j = 0; j <JJ ; j++) {
		gsl_matrix_free(pf->sE[j]);
		gsl_matrix_free(pf->sU[j]);
	}
	free( pf->sE );
	free( pf->sU );
	free_shocks(sk);

	free_valfuns(vf);
	free_hists(ht);
}
void free_valfuns(struct valfuns *vf){
	gsl_matrix_free(vf->WE);
	gsl_matrix_free(vf->RE);
	gsl_matrix_free(vf->WU);
	gsl_matrix_free(vf->RU);
	gsl_matrix_free(vf->WEdist);
}

void free_pars( struct cal_params * par){
	int ji,ii;

	gsl_vector_free(par->Plev);
	gsl_vector_free(par->Alev);
	gsl_vector_free(par->xSlev);
	gsl_vector_free(par->zlev);
	gsl_vector_free(par->epslev);
	gsl_vector_free(par->AloadP);
	for(ji=0;ji<JJ;ji++){
		gsl_matrix_free(par-> Ptrans[ji]);
	}
	free(par->Ptrans);
	free(par->param_lb);
	free(par->param_ub);
	gsl_matrix_free(par->Atrans);
	gsl_matrix_free(par->xStrans);
//	gsl_matrix_free(par->ztrans);
    for(ii=0;ii<NA;ii++) gsl_vector_free(par->epsprob[ii]);
	for(ii=0;ii<NA;ii++) gsl_vector_free(par->zprob[ii]);
	free(par->zprob);
	gsl_vector_free(par->jprob);
	for(ji=0;ji<JJ;ji++)
		free(par->alpha_nf[ji]);
	free(par->alpha_nf);

}

void free_hists( struct hists *ht ){
	int ll,ji;
	for(ll=0;ll<Npaths;ll++){
		gsl_matrix_int_free(ht->uhist[ll]);
		gsl_matrix_free(ht->whist[ll]);
		gsl_matrix_int_free(ht->xShist[ll]);
        gsl_matrix_int_free(ht->epshist[ll]);
        gsl_matrix_int_free(ht->zhist[ll]);

        gsl_matrix_int_free(ht->jhist[ll]);
		gsl_vector_int_free(ht->Ahist[ll]);
		gsl_matrix_int_free(ht->Phist[ll]);
		gsl_matrix_int_free(ht->J2Jhist[ll]);
		for(ji=0;ji<JJ;ji++)
			gsl_matrix_free(ht->sijhist[ll][ji]);
	}

	free(ht->uhist);
	free(ht->whist);
	for(ll=0;ll<Npaths;ll++){
		free(ht->sijhist[ll]);
	}
	free(ht->xShist);
    free(ht->zhist);
    free(ht->epshist);
	free(ht->sijhist);
	free(ht->jhist);
	free(ht->Ahist);
	free(ht->Phist);
	free(ht->J2Jhist);
}
void free_shocks(struct shocks * sk){
	int j;
	for(j=0;j<Npaths;j++){
		gsl_vector_free(sk->Asel[j]);
		gsl_matrix_free(sk->Psel[j]);
		gsl_matrix_free(sk->lambdaEMsel[j]);
		gsl_matrix_free(sk->lambdaESsel[j]);
		gsl_matrix_free(sk->lambdaUSsel[j]);
        gsl_matrix_free(sk->lambdaUMsel[j]);
		gsl_matrix_free(sk->epssel[j]);
		gsl_matrix_free(sk->xGsel[j]);
		gsl_matrix_free(sk->xSsel[j]);
		gsl_matrix_free(sk->zsel[j]);
		gsl_matrix_free(sk->jsel[j]);
		gsl_matrix_free(sk->dsel[j]);
		gsl_matrix_free(sk->msel[j]);
	}
	free(sk->Asel);
	free(sk->Psel);
	free(sk->lambdaESsel);
	free(sk->lambdaUSsel);
    free(sk->lambdaUMsel);
	free(sk->lambdaEMsel);

	free(sk->epssel);
	free(sk->xGsel);
	free(sk->xSsel);
	free(sk->zsel);
	free(sk->jsel);
	free(sk->dsel);
	free(sk->msel);
}

void init_pf( struct polfuns *pf ,struct cal_params * par){

	int ji;
	gsl_matrix_set_all(pf->mE,0.5);
	gsl_matrix_set_all(pf->mU,0.5);
	for(ji=0;ji<JJ;ji++) gsl_matrix_set_all(pf->sE[ji],1./(double)(JJ-1));
	for(ji=0;ji<JJ;ji++) gsl_matrix_set_all(pf->sU[ji],1./(double)(JJ-1));

}

void memcpy_pf(struct polfuns *pf_dest, struct polfuns * pf_orig ){
	int ji,ii;

	gsl_matrix_memcpy(pf_dest->mU,pf_orig->mU);
	gsl_matrix_memcpy(pf_dest->mE,pf_orig->mE);
	for (ji = 0; ji <JJ ; ji++) {
		gsl_matrix_memcpy(pf_dest->sE[ji],pf_orig->sE[ji]);
		gsl_matrix_memcpy(pf_dest->sU[ji],pf_orig->sU[ji]);
	}
}

void init_vf( struct valfuns *vf ,struct cal_params * par){
	int ji;

	gsl_matrix_set_all( vf->WE, 0.);
	gsl_matrix_set_all( vf->WU, 0.);
	gsl_matrix_set_all( vf->RE, 0.);
	gsl_matrix_set_all( vf->RU, 0.);
}

void memcpy_vf(struct valfuns *vf_dest, struct valfuns * vf_orig){
	gsl_matrix_memcpy(vf_dest->WE,vf_orig->WE);
	gsl_matrix_memcpy(vf_dest->RE,vf_orig->RE);

	gsl_matrix_memcpy(vf_dest->WU,vf_orig->WU);
	gsl_matrix_memcpy(vf_dest->RU,vf_orig->RU);
	gsl_matrix_memcpy(vf_dest->WEdist,vf_orig->WEdist);

}
