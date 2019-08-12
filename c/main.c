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


// No solve to debug the optimzation routine
int nosolve = 0;


// declare some parameters that will be global scope
int const JJ = 3;     // number of occupations
int const NG = 2;     // number of human capital types
int const NS = 2;     // number of occ tenure types
int const NZ = 7;     // 7 number of occ match-quality types
int const NE = 7;     // 7 number of firm epsilon match-quality types
int const NP = 5;     // number of occupation-specific productivies
int const NA = 5;     // number of aggregate productivities

int NN ,NUN;


int const TT      = 12*15;    // periods per simulation path
int const burnin  = 12;       // number of periods to throw away each time
int        TTT ;
int const Npaths  = 32;//128      // number of simulation paths to draw
int const Nsim    = 1000;

int const Npwave  = 4;
int const Anan    = 1;
int const Nqtls   = 5; // number of quantiles that will use to compare distributions
double     qtlgrid[] = {0.1,0.25,0.5,0.75,0.9};
int        Nparams = 16;
int        Ntargets= 8;

int        Ncluster = 3;
int        Npar_cluster[4] ={8,0,0,0}; // first the flows parameters, then the dist params then the cyclcical parameters
int        Ntgt_cluster[4] ={8,0,0,0};



int eps_2emg = 1; //should we use a double-exponentially modified gaussian, or just normal
int const nstarts = 1;  // how many starts per node

int verbose = 3;
int print_lev = 3;

int maxiter = 4000; //2000
double vftol = 1e-4;
double rhotightening = .5;
double caltol = 1e-3;

double beta	= 0.997;		// discount factor
double b 	= 0.10; 		// unemployment benefit
double wage_lev = 1;        // will be a shifter so the average wage is >0

double urt_avg = .055;     // average separation rate

// specifics for the calibration routine.
FILE * calhist;
char calhi_f[] = "calhistX.csv";
FILE * parhist;
char parhi_f[] = "paramhistX.csv";
double * solver_state;
int solX =0;

char * parnames_clu0[] = {"alphaE0","alphaU0","alpha1","lambdaU0","lambdaES","lambdaEM","delta","zloss"};
char * parnames_clu1[] = {"var_ze","autoz","var_pe","autop","var_a","autoa","gdfather","stwupdt","var_eps","ltl_eps","rtl_eps"};
char * parnames_clu2[] = {"Delta_A","lamEM_A","lamES_A","lamU_A","occ1","occ2"};
char * parnames_cluall[] = {"alphaE0","alphaU0","alpha1","lambdaU0","lambdaES","lambdaEM","delta","zloss",
                            "var_ze","autoz","var_pe","autop","var_a","autoa","gdfather","stwupdt","var_eps","ltl_eps","rtl_eps",
                            "Delta_A","lamEM_A","lamES_A","lamU_A","occ1","occ2"};
char ** parnames[] = {parnames_clu0,parnames_clu1,parnames_clu2,parnames_cluall};
char * tgtnames_clu0[] = {"J2J_err","fnd_err","sep_err","swEE_err","swU_err","swSt_err","dur_ratio"};
char ** tgtnames_clu1;
char * tgtnames_clu2[] = {"swPrUratio","swPrEEratio","frtratio","seprtratio","EEratio","cor(occ w)"};


struct cal_params{
	int cluster_hr, rank;
	double gdfthr, lambdaEM0, lambdaES0, lambdaU0;
	double alphaU0; 	// scale of alpha function
	double alphaU1;		// concavity of alpha function
	double alphaE0; 	// scale of alpha function
	double alphaE1;		// concavity of alpha function
	double kappa;		// cost of switching
	double autoa;		// persistence of aggregate shock
	double autop;		// persistence of occ-specific shock
	double autoz; 		// persistence of match-quality shock
	double var_ae;		// innovations to aggregate shock
	double var_pe;		// innovations to occ-specific shock
	double var_ze;		// innovations to match-quality shock
	double var_eps;     // variance of epsilon
	double lshape_eps;  // left-tail skewness of epsilon (exponential parameter)
	double rshape_eps;  // right-tail skewness of epsilon
	double wage_curve;  // curviness of wage function
    double delta_Acoef;
	double delta_avg;     // average separation rate
	double lambdaU_Acoef,lambdaES_Acoef,lambdaEM_Acoef;
    double zloss;
	double stwupdate; // update rate for wages of stayers
    double * xparopt, *xsolopt;

    double * param_lb; double * param_ub;  // the upper and lower bounds for all parameters
    double ** cluster_lb; double** cluster_ub;

	gsl_vector * AloadP; //loading on A for each P
	gsl_vector * Plev;
	gsl_vector * Alev;
	gsl_vector * xGlev;
	gsl_vector * xSlev;
	gsl_vector * zlev;
	gsl_vector * epslev;

	gsl_matrix ** Ptrans; // markov transition matrix for P, for each J
	gsl_matrix * Atrans;
	gsl_matrix * xGtrans;
	gsl_matrix * xStrans;
	gsl_matrix * ztrans;
	gsl_vector * zprob; // distribution from which to draw when z first drawn
	gsl_vector * epsprob; // iid, so just a vector of probabilities
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
	gsl_matrix ** lambdaUsel;
	gsl_matrix ** lambdaMsel;
	gsl_matrix ** lambdaSsel;
	gsl_matrix ** Psel;
	gsl_vector ** Asel;


};

struct hists{
	gsl_matrix_int ** uhist;
	gsl_matrix     ** whist;
	gsl_matrix_int*** xhist;
	gsl_matrix_int ** jhist;
	gsl_vector_int ** Ahist;
	gsl_matrix_int ** Phist;
	gsl_matrix_int ** J2Jhist;

};

struct stats{
	double swProb_U;
	double swProb_st;
	double swProb_EE;
	double J2Jprob;
	double unrate;
	double findrate;
	double seprate;
	double udur_nosw;
	double udur_sw;
	double corr_wgocc;

	double seprt_ratio;
	double fndrt_ratio;
	double EE_ratio;
	double swProb_EE_ratio;
	double swProb_U_ratio;

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
};

void allocate_pars( struct cal_params * par );
void free_pars( struct cal_params * par);
void allocate_mats(struct valfuns * vf, struct polfuns * pf, struct hists * ht, struct shocks * sk);
void alloc_valfuns(struct valfuns *vf );
void alloc_hists( struct hists *ht );

void free_mats(struct valfuns * vf, struct polfuns * pf, struct hists *ht, struct shocks * sk);
void free_valfuns(struct valfuns *vf);
void free_hists( struct hists *ht);

void init_pf( struct polfuns *pf ,struct cal_params * par);
void init_vf( struct valfuns *vf ,struct cal_params * par);

int draw_shocks(struct shocks * sk);

// solve, simulate and compute summary statistics
int sol_dyn( struct cal_params * par, struct valfuns * vf, struct polfuns * pf, struct shocks * sk);
int sim( struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk );
int sum_stats(   struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk, struct stats *st  );

double param_dist( double * x, struct cal_params *par, int Npar, double * err_vec , int Nerr);

void set_dat( struct stats * );
void set_params( double * x, int n, struct cal_params * par,int ci);
void print_params(double *x, int n, struct cal_params * par);
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
void w_qtls( double* vec, int stride, int len, double * qtls_out ){
	int qi ;
	if(Nqtls ==5){
		qtls_out[0] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.1);
		qtls_out[1] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.25);
		qtls_out[2] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.5);
		qtls_out[3] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.75);
		qtls_out[4] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.9);
	}
	else{

		for(qi = 0; qi<Nqtls;qi++)
			qtls_out[qi] = gsl_stats_quantile_from_sorted_data(vec,(size_t)stride,(size_t)len, (double)(qi+1)/(double)(Nqtls+1)  );
	}
	//check the order (i.e. must be increasing):
	gsl_sort( qtls_out, 1, Nqtls );
	// check for NaNs:
	for( qi = 0 ; qi<Nqtls;qi++){
		if(isnan(qtls_out[qi])){
			if(verbose>0) printf(" NaN qtl! ");
			if(qi>0 && qi<Nqtls-1){
				qtls_out[qi] = (qtls_out[qi+1]-qtls_out[qi-1])/(qtlgrid[qi+1]-qtlgrid[qi-1])*(qtlgrid[qi]-qtlgrid[qi-1]) + qtls_out[qi-1];
			}else if(qi==0)
				qtls_out[qi] =qtls_out[qi+1];
			else if(qi==Nqtls-1)
				qtls_out[qi] =qtls_out[qi-1];
		}
	}

}

int main(int argc,char *argv[] ) {

	int i,ii,ji,j, ci;
	int success, rank,nnodes;

	int cal_now;

	struct cal_params par;
	struct valfuns vf;
	struct polfuns pf;
	struct shocks sk;
	struct hists ht;

	struct stats st;


	// actually take this from argc
	cal_now = 1;

	Npar_cluster[0] = 8;
	Npar_cluster[1] = 9;
	if(eps_2emg==1) Npar_cluster[1] +=2;
	Npar_cluster[2] = 4+JJ-JJ/3; //cyclical parameters
	Nparams = 0;
	for(i=0;i<Ncluster;i++)Nparams += Npar_cluster[i];

	Ntgt_cluster[0] = Npar_cluster[0];
	Ntgt_cluster[1] = Nqtls*8;
	Ntgt_cluster[2] = 5+Nqtls*2;
	Ntargets =0;
	for(i=0;i<Ncluster;i++) Ntargets += Ntgt_cluster[i];

	Npar_cluster[Ncluster] = 0;
	for(ci=0;ci<Ncluster;ci++) Npar_cluster[Ncluster] +=Npar_cluster[ci];
	Ntgt_cluster[Ncluster] = 0;
	for(ci=0;ci<Ncluster;ci++) Ntgt_cluster[Ncluster] +=Ntgt_cluster[ci];

	if(verbose>1){
		printf("Solving for 3 clusters, sized %d,%d,%d\n", Npar_cluster[0],Npar_cluster[1],Npar_cluster[2]);
		printf("Targeting 3 clusters, sized %d,%d,%d\n", Ntgt_cluster[0],Ntgt_cluster[1],Ntgt_cluster[2]);
	}

	st.EEns_qtls = malloc(Nqtls*sizeof(double));st.EEsw_qtls = malloc(Nqtls*sizeof(double));
	st.EUns_qtls = malloc(Nqtls*sizeof(double));st.EUsw_qtls = malloc(Nqtls*sizeof(double));
	st.UEns_qtls = malloc(Nqtls*sizeof(double));st.UEsw_qtls = malloc(Nqtls*sizeof(double));
	st.stns_qtls = malloc(Nqtls*sizeof(double));st.stsw_qtls = malloc(Nqtls*sizeof(double));

	tgtnames_clu1  = malloc(sizeof(char)*8*Nqtls*11);
	char tgthr[11];

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
	gsl_matrix_set_all(par.Atrans, 0.025/( (double)NA-1. ));
	for(ii=0;ii<NA;ii++) gg_set(par.Atrans,ii,ii,0.975);
	gsl_matrix_set_all(par.xGtrans, 0.025/( (double)NG-1. ));
	for(ii=0;ii<NG;ii++) gg_set(par.xGtrans,ii,ii,0.975);

	par.var_ze = 0.025*0.025; par.autoz = 0.975;
	par.zloss  = 0.025;

	gsl_matrix_view ztransLower =gsl_matrix_submatrix(par.ztrans,1,1,NZ-1,NZ-1);
	gsl_vector_view zlevLower = gsl_vector_subvector(par.zlev,1,NZ-1);
	gsl_vector_view zprobLower = gsl_vector_subvector(par.zprob,1,NZ-1);
	rouwenhorst(par.autoz ,pow(par.var_ze,0.5),& (ztransLower.matrix), &(zlevLower.vector) );
	ergod_dist( &(ztransLower.matrix) , &(zprobLower.vector));
	gsl_vector_set(par.zlev,0,5*par.zlev->data[1]);
	gsl_vector_set(par.zprob,0,0.);
	for(ii=0;ii<NZ;ii++){
		gg_set(par.ztrans,ii,0,par.zloss);
		gg_set(par.ztrans,0,ii,gsl_vector_get(par.zprob,ii));
		gsl_vector_view zrow = gsl_matrix_row(par.ztrans,ii);
		if(ii>0){
			double zrowsum = 0;
			int ci;
			for( ci=0;ci<NZ;ci++ ) zrowsum += gsl_vector_get( &(zrow.vector),ci );
			gsl_vector_scale(  & (zrow.vector) ,1./zrowsum);
		}
	}

	gsl_matrix_set_all(par.xStrans,.025/( (double)NS-1.) );
	for(ii=0;ii<NS;ii++) gg_set(par.xStrans,ii,ii,.975);
	gsl_vector_set_all(par.epsprob, 1./(double) NE);
	gsl_vector_set_all(par.AloadP,1.0);
	par.autoa = 0.95;par.var_ae = 0.01*0.01;
	rouwenhorst(par.autoa,pow(par.var_ae,0.5),par.Atrans,par.Alev);

	par.autop = 0.95; par.var_pe = 0.02*0.2;
	rouwenhorst(par.autop,pow(par.var_pe,0.5),par.Ptrans[0],par.Plev);
	for(i=1;i<JJ;i++){
		gsl_matrix_memcpy(par.Ptrans[i],par.Ptrans[0]);
	}

	int ai = 0;int pi = 0;int gi = 0;int si = 0;int zi = 1;	int thi = 0; ji=0;



	par.alphaE1 = 0.5;
	par.alphaE0 = 0.05*pow((double)(JJ-1),-par.alphaE1);
	par.alphaU1 = 0.5;
	par.alphaU0 = 1.*pow((double)(JJ-1),-par.alphaU1);
	par.lambdaU0  = 0.2 / 0.5;
	par.lambdaES0 = 0.01;
	par.lambdaEM0 = 0.8;
	par.kappa     = 0.00 ;
	par.gdfthr    = 0.5 ;
	par.wage_curve= 0.0 ;
	par.delta_avg = 0.01;
	par.delta_Acoef = 0.;
	par.lambdaU_Acoef = 0.0;
    par.lambdaES_Acoef = 0.0;
	par.lambdaEM_Acoef = 0.0;

	par.var_eps = 0.1;
	par.lshape_eps = 1.0;par.rshape_eps = 1.0;


    double wage_lev0 = exp(par.AloadP->data[ji] * par.Alev->data[ai] +
                               par.Plev->data[pi] +
                               par.epslev->data[thi] +
                               par.zlev->data[zi] +
                               par.xSlev->data[si] +
                               par.xGlev->data[gi]);
    wage_lev = 0.;


	// parameter space:
	// alphaE , alphaU, alphaU_1, lambdaU,lambdaES, lambdaEM, delta_avg, zloss_prob
	par.param_lb[0] = 0.010; par.param_ub[0] = 0.75;
	par.param_lb[1] = 0.010; par.param_ub[1] = 1.00;
	par.param_lb[2] = 0.350; par.param_ub[2] = 0.75;

	par.param_lb[3] = 0.002; par.param_ub[3] = 0.50;
	par.param_lb[4] = 0.002; par.param_ub[4] = 0.15;
	par.param_lb[5] = 0.002; par.param_ub[5] = 1.00;
	par.param_lb[6] = 0.001; par.param_ub[6] = 0.05;
	par.param_lb[7] = 0.001; par.param_ub[7] = 0.10;

	if(Ncluster>1){
		ii = Npar_cluster[0];
		// var_ze, autooz,
		par.param_lb[0+ii] = 0.001; par.param_ub[0+ii] = 0.050;
		par.param_lb[1+ii] = 0.600; par.param_ub[1+ii] = 0.999;
		//var_pe, autop, var_ae,autoa, gdfather, stwupdate
		par.param_lb[2+ii]= 0.001;  par.param_ub[2+ii]= 0.050;
		par.param_lb[3+ii]= 0.600;  par.param_ub[3+ii]= 0.999;
		par.param_lb[4+ii]= 0.001;  par.param_ub[4+ii]= 0.050;
		par.param_lb[5+ii]= 0.600;  par.param_ub[5+ii]= 0.999;
		par.param_lb[6+ii]= 0.01;   par.param_ub[6+ii]= 0.99;
		par.param_lb[7+ii]= 0.01;   par.param_ub[7+ii]= 0.99;
		//var_eps, lshape_eps, rshape_eps
		par.param_lb[8+ii]= 0.001;  par.param_ub[8+ii]= 0.10; //std of 0.5 as upper limit
		if(eps_2emg == 1){
			par.param_lb[9 +ii]= 0.500;  par.param_ub[9 +ii]= 2.500;
			par.param_lb[10+ii]= 0.500;  par.param_ub[10+ii]= 2.500;
		}

	}
	if(Ncluster>2){
		ii = Npar_cluster[1]+Npar_cluster[0];
		//delta_Acoef, lambdaEM, lambdaES, lambdaU (can up-to double)
		par.param_lb[0+ii]=0.001;   par.param_ub[0+ii] = 0.999;
		par.param_lb[1+ii]=0.001;   par.param_ub[1+ii] = 0.999;
		par.param_lb[2+ii]=0.001;   par.param_ub[2+ii] = 0.999;
		par.param_lb[3+ii]=0.001;   par.param_ub[3+ii] = 0.999;
		ii +=4;
		for(ji=0;ji<JJ;ji++){
			if(ji%3==1){
				par.param_lb[ii]=0.001;
				par.param_ub[ii] = 1.999;
				ii++;
			}
			if(ji%3==1){
				par.param_lb[ii]=0.001;
				par.param_ub[ii] = 0.999;
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


	int DFBOLS_ALG;
#ifdef _DFBOLS_USE
	DFBOLS_ALG =1;
#else
	DFBOLS_ALG =0;
#endif


	// branch for MPI here
#ifdef _MPI_USE
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status mpistatus;

#else
	rank = 0;
	nnodes = 1;

#endif
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

	sprintf(calhi_f,"calhist%d.csv",rank);
	calhist = fopen(calhi_f,"w+");
	fprintf(calhist,"dist,J2J_err,fnd_err,sep_err,swEE_err,swU_err,swSt_err,dur_ratio,corr_w_occw,");
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
		for(ii=0;ii<Nqtls;ii++) fprintf(calhist," MVsw%f, ",qtlgrid[ii] );
		for(ii=0;ii<Nqtls;ii++) fprintf(calhist," MVns%f, ",qtlgrid[ii] );
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
	if(rank==0){
		MPI_Scatter( x0starts , nsend,MPI_DOUBLE,x0starts_j,nsend,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
#else
	for(i=0;i<nstarts*Nparams;i++){
		x0starts_j[i] = x0starts[i];
	}
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


		if(i%2 == rank%2 ){

			if(verbose>0)
				printf("Beginning to evaluate a DFBOLS start point \n");

			glb_par = &par;


			// loop over clusters
			for(ci=0;ci<Ncluster+1;ci++){

				set_params(par.xparopt,Nparams,&par,Ncluster); //sets all of the parameters to the initial guess
				par.cluster_hr = ci;
				if(cal_now==0){
					ci = Ncluster;
					par.cluster_hr = ci;
				}
				solver_state = malloc(sizeof(double)*(Npar_cluster[ci]+1));
				solver_state[0] = 1e4;
				// allocations for DFBOLS
				x0_clu = malloc(Npar_cluster[ci]*sizeof(double));
				int npt = 2*Npar_cluster[ci] +1;
				double rhobeg = 0.5*pow((double)(nnodes*nstarts),-1./(double)Npar_cluster[ci]);
				double rhoend = 1e-3;
				if(ci==Ncluster) rhobeg /= 5;

				int maxfun = 400*(Npar_cluster[ci]+1);
			//	int	maxfun = 2*Npar_cluster[ci]+2;

				double *wspace = calloc( (npt+5)*(npt+Npar_cluster[ci])+3*Npar_cluster[ci]*(Npar_cluster[ci]+5)/2 ,sizeof(double) );

				double*dfbols_lb,*dfbols_ub;
				dfbols_lb = calloc(Npar_cluster[ci],sizeof(double));
				dfbols_ub = calloc(Npar_cluster[ci],sizeof(double));
				par.cluster_lb[ci] = malloc(Npar_cluster[ci]*sizeof(double));
				par.cluster_ub[ci] = malloc(Npar_cluster[ci]*sizeof(double));
				for(ii=0;ii<Npar_cluster[ci];ii++){ dfbols_lb[ii]=0.;dfbols_ub[ii]=1.; }
				int par_hr =0;
				// set up the bounds as a subset of the whole set of bounds
				if(ci< Ncluster){
					for(ii=0;ii<ci;ii++ ){ par_hr += Npar_cluster[ii];}
				}
				for(ii=0;ii<Npar_cluster[ci];ii++){
					x0_clu[ii] = x0[ii+par_hr];
					par.cluster_lb[ci][ii] = par.param_lb[ par_hr + ii];
					par.cluster_ub[ci][ii] = par.param_ub[ par_hr + ii];
				}
				if(verbose>1){
					printf("Bounds are: \n    ");
					for(ii=0;ii<Npar_cluster[ci];ii++){ printf("%8s,",parnames[ci][ii]); }
					printf("\n");
					printf("lb: "); for(ii=0;ii<Npar_cluster[ci];ii++){ printf("%8.6f,",par.cluster_lb[ci][ii]); }
					printf("\n");
					printf("ub: "); for(ii=0;ii<Npar_cluster[ci];ii++){ printf("%8.6f,",par.cluster_ub[ci][ii] ); }
					printf("\n");
				}

				int dfbols_printlev = print_lev > 3? 3:print_lev; dfbols_printlev = print_lev <0 ? 0:print_lev;
				if( cal_now ==1 ){
					#ifdef _DFBOLS_USE
				    bobyqa_h_(&(Npar_cluster[ci]),&npt,x0_clu,dfbols_lb,dfbols_ub,&rhobeg,&rhoend,&dfbols_printlev ,&maxfun,wspace,&(Ntgt_cluster[ci]));
					#endif

					for(ii=0;ii<Npar_cluster[ci];ii++) par.xparopt[ii+par_hr] = solver_state[ii+1];
					for(ii=0;ii<Npar_cluster[ci];ii++){
						x0[ii+par_hr] = (par.xparopt[ii+par_hr] - par.cluster_lb[ci][ii])/(par.cluster_ub[ci][ii]-par.cluster_lb[ci][ii]);
					}
				}else{
					ci = Ncluster;
					par.cluster_hr = ci;
					for(ii=0;ii<Nparams;ii++)
						par.xparopt[ii] = x0[ii] * (par.param_ub[ii]-par.param_lb[ii])
						                  +par.param_lb[ii];

				}
				par.cluster_hr= Ncluster;
				dist = param_dist(par.xparopt, & par ,Nparams,err,Ntargets);


				if(verbose>1){
					printf("error is %f at vector:  (", dist);
					for(ii=0;ii<Ntgt_cluster[ci]-1;ii++)
						printf("%f,", err[ii]);
					printf("%f)\n", err[Ntgt_cluster[ci]-1]);
					printf("evaluated at (");
					for(ii=0;ii<Npar_cluster[ci]-1;ii++)
						printf("%f,", par.xparopt[par_hr+ii]);
					printf("%f)\n", par.xparopt[par_hr+Npar_cluster[ci]-1]);
				}

				free(wspace);free(dfbols_lb);free(dfbols_ub);
				free(par.cluster_ub[ci]);free(par.cluster_lb[ci]);
				free(x0_clu);
			}

			free(solver_state);
			if(verbose>0)
				printf("Evaluated a DFBOLS start point \n");
		}else{


			if(verbose>0)
				printf("Beginning to evaluate a NM start point \n");
			for(ci=0;ci<Ncluster+1;ci++){
				par.cluster_hr = ci;

				x0_clu = malloc(Npar_cluster[ci]*sizeof(double));
				solver_state = malloc(sizeof(double)*(Npar_cluster[ci] +1));
				solver_state[0] = 1e4;
				nlopt_opt opt  = nlopt_create(NLOPT_LN_SBPLX,(unsigned)Npar_cluster[ci] );

				double*nlopt_lb,*nlopt_ub;
				nlopt_lb = calloc(Npar_cluster[ci],sizeof(double));
				nlopt_ub = calloc(Npar_cluster[ci],sizeof(double));
				par.cluster_ub[ci] = malloc(Npar_cluster[ci]*sizeof(double));
				par.cluster_lb[ci] = malloc(Npar_cluster[ci]*sizeof(double));

				for(ii=0;ii<Npar_cluster[ci];ii++){ nlopt_lb[ii]=0.;nlopt_ub[ii]=1.; }
				nlopt_set_lower_bounds(opt,nlopt_lb);
				nlopt_set_upper_bounds(opt,nlopt_ub);

				// set up the bounds as a subset of the whole set of bounds
				int par_hr =0;
				if(ci< Ncluster){
					for(ii=0;ii<ci;ii++ ){ par_hr += Npar_cluster[ii];} // offset if we're not using all of the parameters
				}
				for(ii=0;ii<Npar_cluster[ci];ii++){
					x0_clu[ii] = x0[ii+par_hr];
					par.cluster_lb[ci][ii] = par.param_lb[ par_hr + ii];
					par.cluster_ub[ci][ ii] = par.param_ub[ par_hr + ii];
				}
				if(verbose>1){
					printf("Bounds are: \n    ");
					for(ii=0;ii<Npar_cluster[ci];ii++){ printf("%8s,",parnames[ci][ii]); }
					printf("\n");
					printf("lb: "); for(ii=0;ii<Npar_cluster[ci];ii++){ printf("%f,",par.cluster_lb[ci][ii]); }
					printf("\n");
					printf("ub: "); for(ii=0;ii<Npar_cluster[ci];ii++){ printf("%f,",par.cluster_ub[ci][ii] ); }
					printf("\n");
				}

				double tol= 1e-4;
				nlopt_set_xtol_rel(opt,tol);
				nlopt_set_ftol_abs(opt,tol);
				nlopt_set_ftol_rel(opt,tol);
			//	nlopt_set_initial_step(opt, .4);//0.5*pow((double)(nnodes*nstarts),-1./(double)Npar_cluster[ci]) );
				nlopt_set_maxeval(opt,100*(Npar_cluster[ci] +1));

				struct cal_params * par_pt = &par;
				//glb_par = &par;
				nlopt_set_min_objective(opt,f_wrapper_nlopt, (void*) par_pt);
				nlopt_optimize(opt,x0_clu,&dist);

				for(ii=0;ii<Npar_cluster[ci];ii++) par.xparopt[ii+par_hr] = solver_state[ii+1];
				for(ii=0;ii<Npar_cluster[ci];ii++){
					x0[ii+par_hr] = (par.xparopt[ii] - par.cluster_lb[ci][ii])/(par.cluster_ub[ci][ii]-par.cluster_lb[ci][ii]);
				}
				dist = param_dist(par.xparopt, & par ,Nparams,err,Ntargets);

				if(verbose>1){
					printf("error is %f at vector:  (", dist);
					for(ii=0;ii<Ntgt_cluster[ci]-1;ii++)
						printf("%f,", err[ii]);
					printf("%f)\n", err[Ntgt_cluster[ci]-1]);
					printf("evaluated at (");
					for(ii=0;ii<Npar_cluster[ci]-1;ii++)
						printf("%f,", par.xparopt[ii]);
					printf("%f)\n", par.xparopt[Npar_cluster[ci]-1]);
				}

				free(solver_state);
				nlopt_destroy(opt);
				free(nlopt_lb);free(nlopt_ub);
				free(x0_clu);
				free(par.cluster_ub[ci]);free(par.cluster_lb[ci]);
				if(verbose>0)
					printf("Evaluated a NM start point \n");
			}
		}
		caldist_j[i] = dist;
		for(ii=0;ii<Nparams;ii++) calx_j[i*Nparams + ii] = x0[ii];

#ifdef _MPI_USE
		int ngather = nstarts*Nparams;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(caldist_j,nstarts,MPI_DOUBLE,caldist,nstarts,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(calx_j,ngather ,MPI_DOUBLE,calx,ngather ,MPI_DOUBLE,0,MPI_COMM_WORLD);
#else
		for(ii=0;ii<nstarts;ii++) caldist[ii] = caldist_j[ii];
		for(ii=0;ii<nstarts*Nparams;ii++) calx[ii] = calx_j[ii];
#endif
		if(rank==0){
			if(verbose>0)
				printf("Looking for best point \n");
			for(j=0;j<nnodes;j++){
				dist = caldist[j+i*nnodes];

				if (dist<mindist){
					mindist = dist;
					for(ii=0;ii<Nparams;ii++) par.xparopt[ii] = calx[j*Nparams+ ii];
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


	// gather all of the mindist back to node 1

	free(x0starts_j);
	free(caldist_j);free(calx_j);

	free(x0starts);
	free(caldist);free(calx);

	free(err); free(x0);
	free(par.xparopt);

	free_mats(&vf,&pf,&ht,&sk);
	free_pars(&par);
	free(st.EEns_qtls);free(st.EEsw_qtls);
	free(st.EUns_qtls);free(st.EUsw_qtls);
	free(st.UEns_qtls);free(st.UEsw_qtls);
	free(st.stns_qtls);free(st.stsw_qtls);


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
			int ai = ii / (NP * NG * NS * NZ * NE);
			int pi = (ii - ai * NP * NG * NS * NZ * NE) / (NG * NS * NZ * NE);
			int gi = (ii - ai * NP * NG * NS * NZ * NE - pi * NG * NS * NZ * NE) / (NS * NZ * NE);
			int si = (ii - ai * NP * NG * NS * NZ * NE - pi * NG * NS * NZ * NE - gi * NS * NZ * NE) / (NZ * NE);
			int zi =
					(ii - ai * NP * NG * NS * NZ * NE - pi * NG * NS * NZ * NE - gi * NS * NZ * NE - si * NZ * NE) / NE;
			int ti = ii - ai * NP * NG * NS * NZ * NE - pi * NG * NS * NZ * NE - gi * NS * NZ * NE - si * NZ * NE -
			         zi * NE;

			wagevec[ii][ji] = exp(par->AloadP->data[ji] * par->Alev->data[ai] +
			                  par->Plev->data[pi] +
			                  par->epslev->data[ti] +
			                  par->zlev->data[zi] +
			                  par->xSlev->data[si] +
			                  par->xGlev->data[gi]) + wage_lev;
		}
	}

	for(ji=0;ji<JJ;ji++) {
		for (ii = 0; ii < NN; ii++) {
			gg_set(vf0.WE, ii, ji,
			               pow(wagevec[ii][ji], 1.-par->wage_curve) / (1. - par->wage_curve) / (1. - beta));
		}
	}
	for(ii=0;ii<NUN;ii++){
		gg_set(vf0.WU,ii,0,
		               (1.-par->lambdaU0)*pow(b, 1.-par->wage_curve) / (1. - par->wage_curve)
		               /(1.-beta) +
		               par->lambdaU0*pow(wagevec[ii*NE+NE/2][0], 1.-par->wage_curve)/(1.-par->wage_curve)
		               /(1.-beta));
		for(ji=1;ji<JJ;ji++)  gg_set(vf0.WU,ii,ji, gg_get(vf0.WU,ii,0));
	}


	for(viter = 0;viter<maxiter;viter++){

		for(ji=0;ji<JJ;ji++){
			#pragma omp parallel for private(ii) firstprivate(ji)
			for(ii=0;ii<NN;ii++){
				int ai = ii/(NP*NG*NS*NZ*NE);
				int pi = (ii - ai*NP*NG*NS*NZ*NE)/(NG*NS*NZ*NE);
				int gi = (ii - ai*NP*NG*NS*NZ*NE - pi*NG*NS*NZ*NE)/(NS*NZ*NE);
				int si = (ii - ai*NP*NG*NS*NZ*NE - pi*NG*NS*NZ*NE - gi*NS*NZ*NE)/(NZ*NE) ;
				int zi = (ii - ai*NP*NG*NS*NZ*NE - pi*NG*NS*NZ*NE - gi*NS*NZ*NE -  si*NZ*NE)/NE ;
				int ti = ii -  ai*NP*NG*NS*NZ*NE - pi*NG*NS*NZ*NE - gi*NS*NZ*NE - si*NZ*NE - zi*NE;

				int iU = ai*NP*NG*NS*NZ + pi*NG*NS*NZ + gi*NS*NZ + si*NZ +zi;

                double lambdaEMhr = par->lambdaEM0 * (1. + par->lambdaEM_Acoef*gsl_vector_get(par->Alev,ai));
                double lambdaEShr = par->lambdaES0 * (1. + par->lambdaES_Acoef*gsl_vector_get(par->Alev,ai));

				double delta_hr = par->delta_avg * (1. + par->delta_Acoef * gsl_vector_get(par->Alev,ai) );
				//compute expectations over A, Pt
				double EAPWE = 0.;
				int aai, ppi;
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++)  EAPWE += gsl_max(gg_get( vf0.WE,  aai*NP*NG*NS*NZ*NE + ppi*NG*NS*NZ*NE+ gi*NS*NZ*NE + si*NZ*NE +zi*NE+ ti ,ji),
					                                          gg_get( vf0.WU,  aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ   + gi*NS*NZ    + si*NZ    +zi        ,ji))*
							gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
				}

				double EtTWE = 0.;
				int tti,zzi;
				double ttiprod = 0.;
				for(tti=ti;tti<NE;tti++) ttiprod+=gsl_vector_get(par->epsprob,tti); //make sure tpobs sum to 1
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++){
						for(tti=ti;tti<NE;tti++)
							EtTWE +=  gsl_max(gg_get(vf0.WE, aai*NP*NG*NS*NZ*NE + ppi*NG*NS*NZ*NE+ gi*NS*NZ*NE + si*NZ*NE +zi*NE+ tti ,ji) ,
							                  gg_get(vf0.WU, aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ   + gi*NS*NZ    + si*NZ    +zi         ,ji))*
									gsl_vector_get(par->epsprob,tti)/ttiprod *gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
					}
				}
				double EtWE = 0.;
				for(aai=0;aai<NA;aai++) {
					for (ppi = 0; ppi < NP; ppi++) {

						for (tti = 0; tti < NE; tti++)
							EtWE += gsl_max(gg_get(vf0.WE, aai*NP*NG*NS*NZ*NE + ppi*NG*NS*NZ*NE + gi*NS*NZ*NE + si*NZ*NE + zi*NE + tti, ji) ,
									        gg_get(vf0.WU, aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ    + gi*NS*NZ    + si*NZ    + zi         , ji))*
							        par->epsprob->data[tti]*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
					}
				}
				int jji;
				double REhr = -par->kappa;
				double *EzWE =malloc(JJ*sizeof(double));
				double *EztWE=malloc(JJ*sizeof(double));
				for(jji=0;jji<JJ;jji++){
					EzWE[jji]  = 0.;
					EztWE[jji] = 0.;
					if(jji!=ji){

						for(aai=0;aai<NA;aai++){
							for(ppi=0;ppi<NP;ppi++) {
								for (zzi = 0; zzi < NZ; zzi++)
									EzWE[jji] += gsl_max(gg_get(vf0.WE,
									                            aai*NP*NG*NS*NZ*NE + ppi*NG*NS*NZ*NE +gi*NS*NZ*NE + zzi*NE + ti, jji) ,
									                     gg_get(vf0.WU,
									                            aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ    +gi*NS*NZ    + zzi        , jji)  )*
									             gsl_vector_get(par->zprob, zzi)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);

							}
						}
						for(aai=0;aai<NA;aai++){
							for(ppi=0;ppi<NP;ppi++) {
								for(tti=0;tti<NE;tti++){
									for(zzi=0;zzi<NZ;zzi++)
										EztWE[jji] += gsl_max(gg_get(vf0.WE,aai*NP*NG*NS*NZ*NE + ppi*NG*NS*NZ*NE+ gi*NS*NZ*NE + zzi*NE + tti,jji) ,
												              gg_get(vf0.WU,aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ   + gi*NS*NZ    + zzi         ,jji))*
												(par->zprob->data[zzi])*(par->epsprob->data[tti])*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
								}
							}
						}
						// constructing RE
						REhr += par->alphaE0* pow( gg_get( pf->sE[jji],ii,ji),1.-par->alphaE1)*(lambdaEMhr*EztWE[jji] +(1.- lambdaEMhr)*EzWE[jji] ) ;
					}
				}
				double totalphaS = 0;
				for(jji=0;jji<JJ;jji++){
					if(jji!=ji) totalphaS += par->alphaE0*pow(gg_get(pf->sE[jji],ii,ji),1.-par->alphaE1);
				}
				REhr +=   (1. - gsl_min( totalphaS,1. ))* EAPWE;
				if(gsl_finite(REhr)==0){
					printf("Uhoh. Bad REhr");
				}

				gg_set( vf->RE,ii,ji,REhr);
                double mhr = exp((REhr-EAPWE)/rhotightening)/
                             ( exp((REhr-EAPWE)/rhotightening)+ 1.);
				if( isinf(mhr) | isnan(mhr) ){
					mhr = REhr >= EAPWE ? 1. : 0. ;
				}

				gg_set( pf->mE,ii,ji, mhr );

				//set search direction for next iteration:
				double sEjiDenom = 0.;
				for(jji=0;jji<JJ;jji++){
					double dirreturn = -par->kappa+ lambdaEMhr*EztWE[jji] + (1.-lambdaEMhr)*EzWE[jji]- EAPWE;
					dirreturn = dirreturn < 0. ? 0. : dirreturn;
					if(jji!=ji) sEjiDenom += pow(dirreturn , 1./par->alphaE1);
				}
				if(sEjiDenom >0){
					for(jji=0;jji<JJ;jji++){
						if(jji!=ji){
							double dirreturn =-par->kappa+ lambdaEMhr*EztWE[jji] + (1.-lambdaEMhr)*EzWE[jji]- EAPWE;
							dirreturn = dirreturn< 0 ? 0. : dirreturn;
							gg_set(pf->sE[jji] ,ii,ji, pow(dirreturn , 1./par->alphaE1)/sEjiDenom);
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


				// update the value function
				double WEhr = pow(wagevec[ii][ji], 1.-par->wage_curve)/(1.-par->wage_curve)+
							  beta*delta_hr*gg_get(vf0.WU,iU,ji) +
				              beta*(1.- delta_hr )*(
				              		gg_get(pf->mE,ii,ji)*gg_get( vf->RE,ii,ji) +
				              		(1.-gg_get(pf->mE,ii,ji))*(par->gdfthr*lambdaEShr*EtWE + (1.-par->gdfthr)*lambdaEShr*EtTWE+
				                                                (1.-lambdaEShr)*EAPWE )  );

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
				int ai =  ii/(NP*NG*NS*NZ);
				int pi = (ii - ai*NP*NG*NS*NZ)/(NG*NS*NZ);
				int gi = (ii - ai*NP*NG*NS*NZ - pi*NG*NS*NZ)/(NS*NZ);
				int si = (ii - ai*NP*NG*NS*NZ - pi*NG*NS*NZ - gi*NS*NZ)/(NZ) ;
				int zi =  ii - ai*NP*NG*NS*NZ - pi*NG*NS*NZ - gi*NS*NZ -  si*NZ ;

				int jji,zzi,tti,aai,ppi;

				double lambdaUhr = par->lambdaU0*(1 + par->lambdaU_Acoef*par->Alev->data[ai]);

				double EAPWU = 0.;
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++) {
						EAPWU += gg_get(vf0.WU, aai*NP*NG*NS*NZ+ppi*NG*NS*NZ+gi*NS*NZ + si*NZ + zi,ji)*
								gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
					}
				}
				double EzWU[JJ];
				double RUhr = -par->kappa;
				for(jji=0;jji<JJ;jji++){
					EzWU[jji]=0;
					double EtWE_jji=0;
					if(jji!=ji){
						// switchers might find a job right away:
						for(aai=0;aai<NA;aai++){
							for(ppi=0;ppi<NP;ppi++) {
								for(tti=0;tti<NE;tti++){
									int iihr = aai*NP*NG*NS*NZ*NE+ppi*NG*NS*NZ*NE+gi*NS*NZ*NE+si*NZ*NE+zi*NE+tti;
									EtWE_jji += gsl_max(gg_get(vf0.WE, iihr, jji), gg_get(vf0.WU,ii,jji) )*
									        gsl_vector_get(par->epsprob,tti)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[jji], pi, ppi);}
							}
						}
						for(aai=0;aai<NA;aai++){
							for(ppi=0;ppi<NP;ppi++) {
								for( zzi=0;zzi<NZ;zzi++ )
									EzWU[jji] += gg_get(vf0.WU, aai*NP*NG*NS*NZ+ppi*NG*NS*NZ+gi*NS*NZ + zzi,jji)*
											gsl_vector_get(par->zprob,zzi)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[jji],pi,ppi);
							}
						}
						RUhr += par-> alphaU0*pow(gg_get(pf->sU[jji],ii,ji),1.-par->alphaU1 )*(EzWU[jji]*(1.-lambdaUhr)+
						                                                                       lambdaUhr*EtWE_jji);
					}
				}
				double totalalphaS = 0;
				for(jji=0;jji<JJ;jji++)
					totalalphaS += par-> alphaU0*pow(gg_get(pf->sU[jji],ii,ji),1.-par->alphaU1 );
				RUhr += (1.-totalalphaS )*EAPWU;
				gg_set(vf->RU, ii,ji, RUhr);

				double EtWE=0;
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++) {
						for(tti=0;tti<NE;tti++){
							int iihr = aai*NP*NG*NS*NZ*NE+ppi*NG*NS*NZ*NE+gi*NS*NZ*NE+si*NZ*NE+zi*NE+tti;
							EtWE += gsl_max(gg_get(vf0.WE, iihr, ji), vf0.WU->data[ii*vf0.WU->tda+ji] )*
							        gsl_vector_get(par->epsprob,tti)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji], pi, ppi);}
					}
				}
				double mhr = exp(RUhr/rhotightening-((1.-lambdaUhr)*EAPWU+
                        lambdaUhr*gsl_max(EtWE,EAPWU)  )/rhotightening)/
				               (exp(RUhr/rhotightening-((1.-lambdaUhr)*EAPWU+
                                       lambdaUhr*gsl_max(EtWE,EAPWU)  )/rhotightening)+1.);
				if( isinf(mhr) | isnan(mhr) ){
					mhr = RUhr > (1.-lambdaUhr)*EAPWU+ lambdaUhr*gsl_max(EtWE,EAPWU) ?
							1. : 0. ;
				}
				gg_set(pf->mU,ii,ji, mhr );
				//search dir for next iterations
				double sUdenom = 0.;
				for(jji=0;jji<JJ;jji++){
					double dirreturn = -par->kappa + EzWU[jji]- EAPWU;
					dirreturn =dirreturn<0 ? 0. : dirreturn;
					if(jji!=ji) sUdenom += pow(dirreturn ,1./par->alphaU1);
				}
				if(sUdenom > 0 ) {
					for (jji = 0; jji < JJ; jji++) {
						if (jji != ji) {
							double dirreturn = -par->kappa + EzWU[jji] - EAPWU;
							dirreturn = dirreturn < 0 ? 0. : dirreturn;
							gg_set(pf->sU[jji], ii, ji, pow(dirreturn, 1. / par->alphaU1) / sUdenom);
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

				gg_set( vf->WU, ii, ji, pow(b,1.-par->wage_curve)/(1.-par->wage_curve) +
						beta*( pf->mU->data[ii*pf->mU->tda + ji]*vf->RU->data[ii*vf->RU->tda+ji] +
						(1.-pf->mU->data[ii*pf->mU->tda + ji])*( (1.-lambdaUhr)*EAPWU + lambdaUhr*EtWE ) ) );

			}
		}

		gsl_matrix_memcpy(vf0.WU,vf->WU);
		double maxdist = gsl_max( gsl_matrix_max(vf->WEdist),-gsl_matrix_min(vf->WEdist));
		if(viter % 200 == 0 && verbose >1 )  printf("Max distance is %f on iteration %d \n", maxdist,viter);
		double tolhr = pow(viter+1,-1)*1e-2+vftol;
		if( maxdist < tolhr ){
			success = 0;
			break;
		}else{
			success = (int) maxdist;
		}

	}
	for(ii=0;ii<NN;ii++)free(wagevec[ii]);
	free(wagevec);
	free_valfuns(&vf0);


    return success;

}

int sim( struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk ){

	int i,ti,ll,ji;
	int distiter;

	double cmzprob[NZ];
	double cmepsprob[NE],cmepsprob0[NE];
	double cmjprob[JJ];
	double cmAtrans[NA][NA];
	double cmPtrans[JJ][NP][NP];
	double cmxGtrans[NG][NG];
	double cmxStrans[NS][NS];
	double cmztrans[NZ][NZ];
	double cmxGprob[NG];
	double cmxSprob[NS];

	gsl_matrix ** swprob_hist = malloc(sizeof(gsl_matrix *)*Npaths);
	gsl_matrix ** WE_hist = malloc(sizeof(gsl_matrix *)*Npaths);
	gsl_matrix ** WU_hist = malloc(sizeof(gsl_matrix *)*Npaths);
	gsl_matrix ** epssel_hist = malloc(sizeof(gsl_matrix *)*Npaths);

	for(ll=0;ll<Npaths;ll++){
		swprob_hist[ll] = gsl_matrix_calloc(Nsim,TT);
		WE_hist[ll] = gsl_matrix_calloc(Nsim,TT);
		WU_hist[ll] = gsl_matrix_calloc(Nsim,TT);
		epssel_hist[ll] = gsl_matrix_calloc(Nsim,TTT);
	}

	cmjprob[0] = par->jprob->data[0];
	for(ji=0;ji<JJ-1;ji++) cmjprob[ji+1] = par->jprob->data[ji+1] + cmjprob[ji];
	cmepsprob[0] = par->epsprob->data[0];
	for(i=0;i<NE-1;i++) cmepsprob[i+1] = par->epsprob->data[i+1] + cmepsprob[i];
	for(i=0;i<NE;i++) cmepsprob0[i]=cmepsprob[i];
	cmzprob[0] = par->zprob->data[0];
	for(i=0;i<NZ-1;i++) cmzprob[i+1] = par->zprob->data[i+1] + cmzprob[i];
	gsl_matrix * xtrans2 = gsl_matrix_calloc(NG,NG);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,par->xGtrans,par->xGtrans,0.,xtrans2);
	gsl_matrix * xGprob = gsl_matrix_calloc(NG,NG);
	for(i=0;i<NG*NG*NG*NG*NG*NG*NG*NG;i++) {
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,par->xGtrans,xtrans2,0.,xGprob);
		gsl_matrix_memcpy(xtrans2,xGprob);
	}
	cmxGprob[0]=gg_get(xGprob,0,0);
	for(i=0;i<NG-1;i++) cmxGprob[i+1] = gg_get(xGprob,0,i+1) + cmxGprob[i];
	gsl_matrix_free(xtrans2);gsl_matrix_free(xGprob);

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
	for(i=0;i<NZ;i++){
		for(ti=0;ti<NZ;ti++) cmztrans[i][ti] = ti>0 ? gg_get(par->ztrans,i,ti) + cmztrans[i][ti-1] : gg_get(par->ztrans,i,ti);
	}
	for(i=0;i<NG;i++){
		for(ti=0;ti<NG;ti++) cmxGtrans[i][ti] = ti>0 ? gg_get(par->xGtrans,i,ti) + cmxGtrans[i][ti-1] : gg_get(par->xGtrans,i,ti);
	}
	for(i=0;i<NS;i++){
		for(ti=0;ti<NS;ti++) cmxStrans[i][ti] = ti>0 ? gg_get(par->xStrans,i,ti) + cmxStrans[i][ti-1] : gg_get(par->xStrans,i,ti);
	}
	for(distiter = 0;distiter<maxiter;distiter++){


		#pragma omp parallel for private(i,ti,ll,ji) firstprivate(cmzprob,cmepsprob,cmjprob,cmAtrans,cmPtrans,cmxGtrans,cmxStrans,cmztrans)
		for(ll=0;ll<Npaths;ll++){
			int ** xt,**xtm1; //xG.xS.z.eps
			int At,Atm1;
			int Pt[JJ],Ptm1[JJ];
			int * jt,*jtm1;
			int * ut,*utm1;

			xt = malloc(sizeof(int*)*Nsim);
			for(i=0;i<Nsim;i++) xt[i] = malloc(sizeof(int)*4);
			jt = malloc(sizeof(int)*Nsim);
			xtm1 = malloc(sizeof(int*)*Nsim);
			for(i=0;i<Nsim;i++) xtm1[i] = malloc(sizeof(int)*4);
			jtm1 = malloc(sizeof(int)*Nsim);
			ut = malloc(sizeof(int)*Nsim);
			utm1 = malloc(sizeof(int)*Nsim);

			int ai,gi,si,zi,thi,xi,ii,iU,jji;

			// initial productivity:
			At =NA/2;
			for(ji=0;ji<JJ;ji++) Pt[ji] = NP/2;

			// initial allocation for all these guys
			for(i=0;i<Nsim;i++) {
				jt[i] = 0;
				for (ji = 0; ji < JJ; ji++) if (gg_get(sk->jsel[ll], i, 0) > cmjprob[ji]) ++ jt[i] ;
				xt[i][0] = 0;
				for (gi = 0; gi < NG; gi++) if (gg_get(sk->xGsel[ll], i, 0) > cmxGprob[gi]) ++ xt[i][0] ;
				xt[i][1] = 0;
				for (si = 0; si < NS; si++) if (gg_get(sk->xSsel[ll], i, 0) > cmxSprob[si]) ++ xt[i][1] ;
				xt[i][2] =0;
				for(zi=0;zi<NZ;zi++) if( gg_get(sk->zsel[ll],i,0) > cmzprob[zi] ) ++ xt[i][2] ;
				xt[i][3] =0;
				for(thi=0;thi<NE;thi++) if( gg_get(sk->epssel[ll],i,0) > cmepsprob0[thi] ) ++ xt[i][3] ;

				if( gg_get( sk->lambdaUsel[ll],i,0) <  urt_avg ){
					ut[i] =1;
				}else{
					ut[i] =0;
				}
			}
			for(ti=0;ti<TTT;ti++){
				// increment shocks
				Atm1 = At;
				At=0;
				for(ai=0;ai<NA;ai++) if(gsl_vector_get(sk->Asel[ll],ti) > cmAtrans[Atm1][ai]) ++At;
				if( ti>=burnin ) gsl_vector_int_set(ht->Ahist[ll],ti-burnin, At);
				for(ji=0;ji<JJ;ji++) {
					Ptm1[ji] = Pt[ji];
					Pt[ji]=0;
					for(ai=0;ai<NP;ai++) if(gg_get(sk->Psel[ll],ti,ji) > cmPtrans[ji][Ptm1[ji]][ai]) ++Pt[ji];
				}
				double lambdaEMhr[JJ];double lambdaEShr[JJ]; double lambdaUhr[JJ];
				for(i=0;i<JJ;i++){
					lambdaEMhr[i] = (par->lambdaEM0*(1. + par->lambdaEM_Acoef*At));
					lambdaEShr[i] = (par->lambdaES0*(1. + par->lambdaES_Acoef*At));
					lambdaUhr[i]  = (par->lambdaU0 *(1. + par->lambdaU_Acoef *At));
				}

				for(i=0;i<Nsim;i++){
					// save wages using values carried from last period
					if(ti>=burnin) {
						double wagehr = exp(par->Alev->data[At]) + exp(par->Plev->data[Pt[jt[i]]]) *
						                                               exp(par->xGlev->data[xt[i][0]]) *
						                                               exp(par->xSlev->data[xt[i][1]]) *
						                                               exp(par->zlev->data[xt[i][2]]) *
						                                               exp(par->epslev->data[xt[i][3]]) + wage_lev;
						gg_set(ht->whist[ll], i, ti-burnin,wagehr);
					}

					// save this guy's state
					jtm1[i] = jt[i];
					utm1[i] = ut[i];
					for(xi=0;xi<4;xi++) xtm1[i][xi] = xt[i][xi];

					// increment individual specific shocks:
					xt[i][0] = 0;
					for (gi = 0; gi < NG; gi++) if (gg_get(sk->xGsel[ll], i, ti) > cmxGtrans[xtm1[i][0]][gi]) ++ xt[i][0] ;
					xt[i][1] = 0;
					for (si = 0; si < NS; si++) if (gg_get(sk->xSsel[ll], i, ti) > cmxStrans[xtm1[i][1]][si]) ++ xt[i][1] ;
					xt[i][2] =0;
					for(zi=0;zi<NZ;zi++) if( gg_get(sk->zsel[ll],i,ti) > cmztrans[xtm1[i][2]][zi] ) ++ xt[i][2] ;
					// xt[i][3]? no need to redraw theta unless later there's a job switch
					xt[i][3] = xtm1[i][3];

					ii = At*NP*NG*NS*NZ*NE + Pt[jt[i]]*NG*NS*NZ*NE + xt[i][0]*NS*NZ*NE+xt[i][1]*NZ*NE+xt[i][2]*NE+xt[i][3];
					iU = At*NP*NG*NS*NZ    + Pt[jt[i]]*NG*NS*NZ    + xt[i][0]*NS*NZ   +xt[i][1]*NZ   +xt[i][2];

					//put the shocks in the history matrix
					if(ti>=burnin){
						for(xi=0;xi<4;xi++)	ggi_set(ht->xhist[ll][xi],i, ti-burnin,xtm1[i][xi]);
						ggi_set(ht->uhist[ll],i,ti-burnin,ut[i]);
						ggi_set(ht->jhist[ll],i,ti-burnin,jt[i]);
						ggi_set(ht->Phist[ll],i,ti-burnin, Pt[jt[i]]);

						gg_set(WE_hist[ll],i,ti-burnin, gg_get(vf->WE,ii,jt[i]));
						gg_set(WU_hist[ll],i,ti-burnin, gg_get(vf->WU,iU,jt[i]));
					}

					//evaluate decision rules for the worker:
					if(utm1[i]==0){

						// employed workers' choices
						double delta_hr =par->delta_avg * (1. + par->delta_Acoef * par->Alev->data[At]);
						// separate if zi =0 or unemployment shock or want to separate
						if( delta_hr >gg_get(sk->dsel[ll],i,ti) || xt[i][2]==0 || gg_get(vf->WU,iU,jt[i]) > gg_get(vf->WE,ii,jt[i])  ){
//						if( delta_hr >gg_get(sk->dsel[ll],i,ti) || xt[i][2]==0 ){
							ut[i] = 1;
							xt[i][2] =0;
							if(ti<TTT-1) for(zi=0;zi<NZ;zi++) if( gg_get(sk->zsel[ll],i,ti+1) > cmzprob[zi] ) ++ xt[i][2] ;
						}else{
							ut[i] = 0;
							// stay or go?
							if( ti>=burnin ) gg_set(swprob_hist[ll],i,ti-burnin,  gg_get(pf->mE,ii,jt[i]) );
							if( gg_get(sk->msel[ll],i,ti) < gg_get(pf->mE,ii,jt[i])  ){
								//RE
								double sij[JJ]; sij[0] =0.;

								for(jji=0;jji<JJ;jji++) sij[jji] = jji > 0 ? sij[jji-1] + par->alphaE0*pow(gg_get(pf->sE[jji],ii,jt[i]),1.-par->alphaE1):
										par->alphaE0*pow(gg_get(pf->sE[jji],ii,jt[i]),1.-par->alphaE1);
								if(gg_get(sk->jsel[ll], i, ti) > sij[JJ-1] ){
									jt[i] = jtm1[i];
								}else{
									// successfully switched: if( gg_get(sk->jsel[ll],i,ti) <  alpha(sij[JJ]))
									jt[i] = 0;
									for(jji=0;jji<JJ;jji++){
										if (gg_get(sk->jsel[ll], i, ti) > sij[jji]){
											jt[i] ++;
										}
									}
								}
								if( jt[i] != jtm1[i] ){ // switchers
									if( gg_get(sk->lambdaMsel[ll],i,ti)<  lambdaEMhr[jt[i]] ) {
										xt[i][1] = 0; // lose specific skill
										// draw a new z:
										xt[i][2] =0;
										for(zi=0;zi<NZ;zi++){
											if( gg_get(sk->zsel[ll],i,ti) > cmzprob[zi] ) ++ xt[i][2] ;
										}

										// draw a new epsilon
										gg_set(epssel_hist[ll],i,ti,gg_get(sk->epssel[ll],i,ti));
										xt[i][3] =0;
										for(thi =0; thi<NE;thi++){ if( gg_get( sk->epssel[ll],i,ti ) > cmepsprob[thi] ) ++xt[i][3];}
										xt[i][3] = xt[i][3]>NE-1 ? NE-1: xt[i][3];
										int ip =At*NP*NG*NS*NZ*NE + Pt[jt[i]]*NG*NS*NZ*NE + xt[i][0]*NS*NZ*NE+xt[i][1]*NZ*NE+xt[i][2]*NE+xt[i][3];

										if( gg_get(vf->WE, ip , jt[i]) >= gg_get(vf->WE, ii,jtm1[i] ) ){
											if(ti>=burnin) ggi_set(ht->J2Jhist[ll], i, ti-burnin, 1);
										}else{//don't switch
											xt[i][1] = xtm1[i][1];
											xt[i][2] = xtm1[i][2];
											xt[i][3] = xtm1[i][3];
											jt[i] = jtm1[i];
										}
									}
								}// else nothing happens (except paid kappa)


							}else{
								//WEs
								if( gg_get(sk->lambdaSsel[ll],i,ti) < lambdaEShr[jt[i]] ) {
									if(ti>=burnin) ggi_set(ht->J2Jhist[ll], i, ti-burnin, 1);
									gg_set(epssel_hist[ll],i,ti,gg_get(sk->epssel[ll],i,ti));
									if (gg_get(sk->lambdaUsel[ll], i, ti) < par->gdfthr) { //godfather (gamma) shock?
										xt[i][3] = 0;
										for (thi = 0; thi < NE; thi++){
											if (gg_get(sk->epssel[ll], i, ti) > cmepsprob[thi])++xt[i][3];
										}
										xt[i][3] = xt[i][3] > NE - 1 ? NE - 1 : xt[i][3];
									}else { //climbing the ladder!
										if( (gg_get(sk->epssel[ll], i, ti) > cmepsprob[xt[i][3]]) ){
											xt[i][3] = 0;
											for (thi = 0; thi < NE; thi++){
												if (gg_get(sk->epssel[ll], i, ti) > cmepsprob[thi])++xt[i][3];
											}
											xt[i][3] = xt[i][3] > NE - 1 ? NE - 1 : xt[i][3];
										}
									}

								}
							}
							// impose wage stickiness for stayers. Doing ti> burnin to start at period 1, not 0
							if( ti>burnin){
								if( ggi_get(ht->J2Jhist[ll], i, ti-burnin) ==0 && jt[i] == jtm1[i]){
									if(  gg_get(sk->lambdaUsel[ll],i,ti) < 1.-par->stwupdate ) //lambdaUsel not ever used before here
									gg_set(ht->whist[ll], i, ti-burnin,
									       gg_get(ht->whist[ll], i, ti-burnin-1) );
								}
							}
						}
					// ut ==1 unemployed
					}else{
						ii = At*NP*NG*NS*NZ + Pt[jt[i]]*NG*NS*NZ + xt[i][0]*NS*NZ +xt[i][1]*NZ +xt[i][2];
						if( ti>=burnin){ gg_set(ht->whist[ll],i,ti-burnin,0.);}
						// stay or go?
						if( ti>=burnin ) gg_set(swprob_hist[ll],i,ti-burnin,  gg_get(pf->mU,ii,jt[i]) );
						if( gg_get(pf->mU,ii,jt[i])>gg_get(sk->msel[ll],i,ti) ){
							//RU
							double sij[JJ]; sij[0]=par->alphaU0*pow(gg_get(pf->sU[0],ii,jt[i]), 1.-par->alphaU1);
							for(jji=1;jji<JJ;jji++) sij[jji] = sij[jji-1] + par->alphaU0*pow(gg_get(pf->sU[jji],ii,jt[i]),1.-par->alphaU1);

							if(gg_get(sk->jsel[ll],i,ti) > sij[JJ-1]){
								jt[i] = jtm1[i];
								ut[i] = 1;

							}else{
								// successfully switched: if( gg_get(sk->jsel[ll],i,ti) <  sij[JJ])
								jt[i] = 0;
								for(jji=0;jji<JJ;jji++)
									if (gg_get(sk->jsel[ll], i, ti) > sij[jji]) jt[i]++;
							}
							if( jt[i] != jtm1[i] ){ // switchers
								xt[i][1] = 0; // lose specific skill
								// draw a new z:
								xt[i][2] =0;
								for(zi=0;zi<NZ;zi++){
									if( gg_get(sk->zsel[ll],i,ti) > cmzprob[zi] ) ++ xt[i][2] ;
								}
								if( gg_get(sk->lambdaUsel[ll],i,ti)<  lambdaUhr[jt[i]] ){
									// draw a new epsilon
									xt[i][3] = 0;
									gg_set(epssel_hist[ll],i,ti,gg_get(sk->epssel[ll],i,ti));
									for(thi =0;thi<NE;thi++) if( gg_get( sk->epssel[ll],i,ti ) > cmepsprob[thi] ) ++xt[i][3];
									int iM = At*NP*NG*NS*NZ*NE + Pt[jt[i]]*NG*NS*NZ*NE + xt[i][0]*NS*NZ*NE+xt[i][1]*NZ*NE+xt[i][2]*NE+xt[i][3];
									if( gg_get(vf->WE,iM,jt[i]) < gg_get(vf->WU,ii,jt[i])){
										ut[i] = 1;
										xt[i][3] = xtm1[i][3];
									}else{
										ut[i] = 0;
									}
								}
							}// else nothing happens (except paid kappa)
						// not moving
						}else{
							if( gg_get(sk->lambdaUsel[ll],i,ti)< lambdaUhr[jt[i]] ){ // found a job??
								xt[i][3] = 0;
								gg_set(epssel_hist[ll],i,ti,gg_get(sk->epssel[ll],i,ti));
								for(thi =0;thi<NE;thi++) if( gg_get( sk->epssel[ll],i,ti ) > cmepsprob[thi] ) ++xt[i][3];
								int iM = At*NP*NG*NS*NZ*NE + Pt[jt[i]]*NG*NS*NZ*NE + xt[i][0]*NS*NZ*NE+xt[i][1]*NZ*NE+xt[i][2]*NE+xt[i][3];
								if( gg_get(vf->WE,iM,jt[i]) < gg_get(vf->WU,ii,jt[i])){
									ut[i]  = 1;
									xt[i][3] = xtm1[i][3];
								}else{
									ut[i] = 0;
								}
							}
						}
					}
				} // i=1:NSim

			}


			for(i=0;i<Nsim;i++) free(xt[i]);
			for(i=0;i<Nsim;i++) free(xtm1[i]);
			free(xt);free(xtm1);free(jt);free(jtm1);free(ut);free(utm1);
		} // end omp loop over ll
		double nemp =0.;
		double topepsprob = cmepsprob0[NE-1] - cmepsprob0[NE-2];
		for(i=0;i<NE;i++)cmepsprob0[i]=0.;
		for(i=0;i<NZ;i++)cmzprob[i]=0.;
		for(ll=0;ll<Npaths;ll++){
			for(ti=0;ti<TT;ti++){
				for(i=0;i<Nsim;i++){
					if( gsl_matrix_int_get(ht->uhist[ll],i,ti)==0 ){
						int epshr = gsl_matrix_int_get( ht->xhist[ll][3],i,ti);
						cmepsprob0[ epshr ] += 1.;
						int zhr =  gsl_matrix_int_get( ht->xhist[ll][2],i,ti);
						cmzprob[zhr] += 1.;
						nemp += 1.;
					}
				}
			}
		}
		for(i=0;i<NE;i++)cmepsprob0[i] *= 1./nemp;
		for(i=0;i<NZ;i++)cmzprob[i] *= 1./nemp;

		for(i=1;i<NE;i++) cmepsprob0[i] += cmepsprob0[i-1];
		for(i=1;i<NZ;i++) cmzprob[i] += cmzprob[i-1];
		topepsprob -= (cmepsprob0[NE-1] - cmepsprob0[NE-2] );
		if( gsl_max(topepsprob,-topepsprob)<1e-6 && distiter >2)
			break;

	}

	if(print_lev >=2){
		int append =0;
		for(ll=0;ll<Npaths;ll++){
			printmat("whist.csv",ht->whist[ll],append );
			printmat_int("uhist.csv",ht->uhist[ll],append );
			printvec_int("Ahist.csv",ht->Ahist[ll],append );
			printmat_int("Phist.csv",ht->Phist[ll],append );
			printmat_int("jhist.csv",ht->jhist[ll],append );
			printmat_int("xGhist.csv",ht->xhist[ll][0],append );
			printmat_int("xShist.csv",ht->xhist[ll][1],append );
			printmat_int("zhist.csv",ht->xhist[ll][2],append );
			printmat_int("epshist.csv",ht->xhist[ll][3],append );

			printmat("swprob_hist.csv", swprob_hist[ll],append );
			printmat("WUhist.csv", WU_hist[ll],append );
			printmat("WEhist.csv", WE_hist[ll],append );
			printmat("epsselhist.csv", epssel_hist[ll],append );
			printmat_int("J2Jhist.csv", ht->J2Jhist[ll],append );
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

	int ll, i,ti,wi,si,ji;
	int ri;

	int Nemp=0, Nunemp = 0, Nnosep=0, Nfnd =0, NswU = 0, NswE=0, NJ2J=0, Nspell=0, NswSt=0,Ndur_sw=0,Ndur_nosw=0,Nsep=0,NErisksep=0;
	int Nemp_rec[2]={0,0}, Nunemp_rec[2] = {0,0}, Nnosep_rec[2]={0,0}, Nfnd_rec[2]={0,0}, NswU_rec[2]= {0,0}, NswE_rec[2]={0,0}, NJ2J_rec[2]={0,0}, Nspell_rec[2]={0,0}, NswSt_rec[2]={0,0},Ndur_sw_rec[2]={0,0},Ndur_nosw_rec[2]={0,0};
	int Nsep_rec[2]={0,0},NErisksep_rec[2]={0,0};
	double * occwg = calloc(JJ*TT*Npaths,sizeof(double));
	int    * occsz = calloc(JJ*TT*Npaths,sizeof(int)   );

	double invlaid_wval = -1e6;
	gsl_vector * Atrans_ergod = gsl_vector_calloc(NA);
	int recIndic[NA];double cumAtrans_ergod[NA];
	ergod_dist(par->Atrans,Atrans_ergod);
	cumAtrans_ergod[0] = Atrans_ergod->data[0];
	for(i=1;i<NA;i++) cumAtrans_ergod[i] =Atrans_ergod->data[i] + cumAtrans_ergod[i-1];
	for(i=0;i<NA;i++) recIndic[i] = cumAtrans_ergod[i]<0.25 ? 1 : 0;

	if(verbose>1){
		int recFrac =0;
		for( ll=0;ll<Npaths;ll++){
			for(ti=0;ti<TT;ti++)
				recFrac += recIndic[gsl_vector_int_get(ht->Ahist[ll],ti)];
		}
		printf(" The fraction of periods in recession is %f, with %d periods ", (double)recFrac /(double)(TT*Npaths) , recFrac);
	}

	for(ri=0;ri<3;ri++){
		#pragma omp parallel for private(ll,ti,i,wi,si) reduction( +: Nemp ) reduction( +:Nunemp) reduction( +:Nnosep) reduction( +: Nfnd) reduction( +: NswU) reduction( +: NswE) reduction( +: Nspell) reduction( +: NswSt) reduction( +: Ndur_sw) reduction( +: Ndur_nosw) reduction(+: Nsep) reduction(+:NErisksep)
		for(ll=0;ll<Npaths;ll++){

			int count_wi;

			int Nemp_ll = 0, Nunemp_ll = 0,Nnosep_ll=0, Nfnd_ll =0, NswU_ll = 0, NswE_ll=0, NJ2J_ll=0,Nspell_ll=0, NswSt_ll =0, Ndur_sw_ll=0, Ndur_nosw_ll=0,Nsep_ll = 0,NErisksep_ll=0;
			int sw_spell=-1,tsep=-1;
			for(i=0;i<Nsim;i++){
				sw_spell = -1;tsep=-1;
				for(wi=0;wi<(TT/Npwave);wi++){//first loop over waves (wi), then loop over reference month (si)
					count_wi=0;
	                int Nemp_wi = 0,Nspell_wi=0,Nfnd_wi=0,NJ2J_wi=0,NswE_wi=0,NswSt_wi=0,Nnosep_wi=0,Nunemp_wi=0,NswU_wi=0,Ndur_sw_wi=0,Ndur_nosw_wi=0;
	                int Nsep_wi = 0,NErisksep_wi=0;
				    for(si=0;si<Npwave;si++){
	                    ti = si + wi * Npwave;
						// should I count this observation? for ri=2, always true. For ri=0 only true when expansion and for ri=1 only true when recession
					    if(ri>= recIndic[ gsl_vector_int_get(ht->Ahist[ll],ti)] && count_wi==0)
					        count_wi = 1;
	                    if(ti<TT-1 && ti>0){
	                        if( ggi_get(ht->uhist[ll],i,ti) ==0 ){
	                            Nemp_wi=1;
		                        NErisksep_wi=1;
	                            occwg[ ll*JJ*TT + ti*JJ + ggi_get(ht->jhist[ll],i,ti)] += gg_get(ht->whist[ll],i ,ti);
		                        occsz[ ll*JJ*TT + ti*JJ + ggi_get(ht->jhist[ll],i,ti)] ++;
	                            if( ggi_get(ht->uhist[ll],i,ti+1) ==1 ){
	                                sw_spell = ggi_get(ht->jhist[ll],i,ti-1);
		                            tsep = ti;
		                            Nnosep_wi =0;
		                            Nsep_wi = 1;
	                            }else{
	                                Nnosep_wi =1;
		                            if( ggi_get(ht->J2Jhist[ll],i,ti) ==1 ){
		                                NJ2J_wi = 1;
		                                if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti-1) ) NswE_wi =1;
		                            }else{  // not EE
		                                if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti-1) ) NswSt_wi = 1;
		                            }
	                            }
	                        }else{
	                            if(sw_spell>-1 && tsep>-1) Nunemp_wi =1;

	                            if( ggi_get(ht->uhist[ll],i,ti+1) ==0 && sw_spell>-1 && tsep>-1 ){
	                                Nfnd_wi =1;
	                                Nspell_wi = 1;
		                            if( ggi_get(ht->jhist[ll],i,ti+1) != sw_spell ){
	                                    NswU_wi = 1; //only count switches at the end of the spell.
	                                    Ndur_sw_wi = ti - tsep+1;
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
				    }
				}
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
    st->J2Jprob  = (double) NJ2J / (double) Nnosep ;
    st->findrate = Nunemp>0 ? (double) Nfnd / (double) Nunemp : 1.;
    st->seprate  = Nemp >0 ? (double) Nsep / (double) NErisksep : 1.; // can use Nfnd because only count separations that find again
    st->swProb_EE = NJ2J>0 ? (double) NswE / (double) NJ2J : 0.;
    st->swProb_U = Nspell>0 ? (double) NswU / (double) Nspell : 0.;
	st->swProb_st = Nnosep - NJ2J > 0 ? (double) NswSt / (double) (Nnosep- NJ2J): 0.;
    st->unrate =  (double) Nunemp / (double) ( Nunemp + Nemp );
    st->udur_nosw = Nspell - NswU>0 ? (double) Ndur_nosw/ (double)(Nspell - NswU): 0.;
	st->udur_sw = NswU>0 ? (double) Ndur_sw/ (double)NswU : 0.;

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



	for(ll=0;ll<Npaths;ll++){
		for( ti=0;ti<TT;ti++ ){
			for(ji=0;ji<JJ;ji++)
				occwg[ll*TT*JJ +ti*JJ + ji ] /= (double)occsz[ll*TT*JJ +ti*JJ + ji ];
		}
	}

    double * w_stns = Nemp>0   ? malloc(sizeof(double) * Nemp  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave );
	double * w_stsw = NswSt>0  ? malloc(sizeof(double) * NswSt ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_EEsw = NswE>0   ? malloc(sizeof(double) * NswE  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );

	double * w_EEns = NJ2J>0   ? malloc(sizeof(double) * NJ2J  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_EUsw = NswU>0   ? malloc(sizeof(double) * NswU  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_EUns = Nspell>0 ? malloc(sizeof(double) * Nspell): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_UEsw = NswU>0   ? malloc(sizeof(double) * NswU  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * w_UEns = Nspell>0 ? malloc(sizeof(double) * Nspell): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );

	double * wwv_EEsw = NswE>0   ? malloc(sizeof(double) * NswE  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );
	double * wwv_occsw = NswE>0   ? malloc(sizeof(double) * NswE  ): malloc(sizeof(double) * Npaths*Nsim *TT/Npwave  );

	double ** w_MVswrec = malloc(sizeof(double*)*2);
	double ** w_MVnsrec = malloc(sizeof(double*)*2);
	w_MVswrec[0] = malloc((NJ2J + Nspell*2)*sizeof(double) );
	w_MVswrec[1] = malloc((NJ2J + Nspell*2)*sizeof(double) );

	w_MVnsrec[0] = malloc((NJ2J + Nspell*2)*sizeof(double) );
	w_MVnsrec[1] = malloc((NJ2J + Nspell*2)*sizeof(double) );

	int idx_EUns =0, idx_UEns =0, idx_EEns=0, idx_EUsw=0,idx_UEsw=0,idx_EEsw=0;
	int idx_stns =0, idx_stsw =0;int valid_w =0;int valid_EUE=0;

	int idx_MVswrec[2]; idx_MVswrec[0]=0; idx_MVswrec[1]=0;
	int idx_MVnsrec[2]; idx_MVnsrec[0]=0; idx_MVnsrec[1]=0;
	for(ll=0;ll<Npaths;ll++){

		int sw_spell=-1;
		double wlast =invlaid_wval, wnext=invlaid_wval, wlast_wave=invlaid_wval, wnext_wave=invlaid_wval;
		double w_EU=0.;
		for(i=0;i<Nsim;i++){
			for(wi=0;wi<TT/Npwave;wi++){
				double w_EEsw_wi =0., w_EEns_wi=0., w_stsw_wi=0.,w_stns_wi=0.,w_UEsw_wi=0.,w_UEns_wi=0.,w_EUsw_wi=0.,w_EUns_wi=0.;
				double  wwv_EEsw_wi=0.,wwv_occsw_wi =0.;
				int    I_EEsw_wi =0,  I_EEns_wi=0 , I_stsw_wi=0 ,I_stns_wi=0 ,I_UEsw_wi=0 ,I_UEns_wi=0 ,I_EUsw_wi=0 ,I_EUns_wi=0 ;
				int    sep=0;
				int    ri_wv = 0;
				wlast=0;wnext=0;
				ti = wi*Npwave;

				ri_wv = ri_wv == 1? recIndic[ gsl_vector_int_get(ht->Ahist[ll],ti)]: 1;
				if( ti>3 && ti< TT-Npwave && Anan==0 ){
					wlast = gg_get( ht->whist[ll] ,i,ti-1) + gg_get( ht->whist[ll] ,i,ti-2)+ gg_get( ht->whist[ll] ,i,ti-3);
					wnext = gg_get( ht->whist[ll] ,i,ti+Npwave) + gg_get( ht->whist[ll] ,i,ti+Npwave+1)+ gg_get( ht->whist[ll] ,i,ti+Npwave+2);
					wnext_wave = wnext;wlast_wave=wlast;
				}else if( wi>3 && wi< TT/Npwave-3 && Anan==1  ){
					int ri=0;
					for(ri=1;ri<Npwave*3+1;ri++)
						wlast += gg_get(ht->whist[ll], i, ti - ri) ;
					for(ri=0;ri<Npwave*3;ri++)
						wnext += gg_get( ht->whist[ll] ,i,ti+ri) ;
					wlast_wave = gg_get( ht->whist[ll] ,i,ti-1) + gg_get( ht->whist[ll] ,i,ti-2)+ gg_get( ht->whist[ll] ,i,ti-3);
					wnext_wave = gg_get( ht->whist[ll] ,i,ti+Npwave) + gg_get( ht->whist[ll] ,i,ti+Npwave+1)+ gg_get( ht->whist[ll] ,i,ti+Npwave+2);

				}else{
					wlast=invlaid_wval;wnext=invlaid_wval;wlast_wave=invlaid_wval; wnext_wave=invlaid_wval;
				}
				for(si=0;si<Npwave;si++){
					ti = si + wi * Npwave;
					ri_wv = ri_wv == 1? recIndic[ gsl_vector_int_get(ht->Ahist[ll],ti)]: 1;

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
                            }else{
                                if( ggi_get(ht->J2Jhist[ll],i,ti) ==1 ){

                                    if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti-1) ){
                                        w_EEsw_wi = (wnext- wlast)/(wnext + wlast)*2;
	                                    wwv_EEsw_wi  = (wnext_wave - wlast_wave )/(wnext+wlast)*2 ;
	                                    wwv_occsw_wi = (occwg[ll*JJ*TT+ti*JJ+ggi_get(ht->jhist[ll],i,ti+1)] -occwg[ll*JJ*TT+ti*JJ+ggi_get(ht->jhist[ll],i,ti-1)])/
			                                    (occwg[ll*JJ*TT+ti*JJ+ggi_get(ht->jhist[ll],i,ti+1)] +occwg[ll*JJ*TT+ti*JJ+ggi_get(ht->jhist[ll],i,ti-1)])*2;
                                        if(! gsl_finite(wwv_occsw_wi)) wwv_occsw_wi=0.;
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
                        }else if(ggi_get(ht->uhist[ll],i,ti) ==1 && wlast> invlaid_wval && wnext > invlaid_wval ){ //unemployed
							if( ggi_get(ht->uhist[ll],i,ti+1) ==0 ) valid_EUE ++;
                            if( ggi_get(ht->uhist[ll],i,ti+1) ==0  && sw_spell>-1  ){
                                if( ggi_get(ht->jhist[ll],i,ti+1) != sw_spell){
                                    w_UEsw_wi = (wnext- wlast)/(wnext + wlast)*2;
                                    w_EUsw_wi = w_EU;
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
				if( I_EUns_wi==1 || I_EUsw_wi==1 || I_EEns_wi==1 || I_EEsw_wi==1 || sep==1 ){
					// can't change a job and be a stayer
					I_stsw_wi=0;I_stns_wi=0;
				}

				if(I_EEsw_wi ==1 && idx_EEsw < NswE){
					w_EEsw[idx_EEsw] =w_EEsw_wi;
					wwv_EEsw[idx_EEsw] = wwv_EEsw_wi;
					wwv_occsw[idx_EEsw]= wwv_occsw_wi;
					idx_EEsw ++;

					//now the cycle-specific dist
					w_MVswrec[ri_wv][idx_MVswrec[ri_wv]] = w_EEsw_wi;
					idx_MVswrec[ri_wv]++;
				}
				if(I_EEns_wi==1 && idx_EEns < NJ2J){
					w_EEns[idx_EEns] = w_EEns_wi;
					idx_EEns ++;
					//now the cylce-specific dist
					w_MVnsrec[ri_wv][idx_MVnsrec[ri_wv]] = w_EEns_wi;
					idx_MVnsrec[ri_wv]++;
				}
				if(I_stsw_wi==1 && idx_stsw < NswSt){
					w_stsw[idx_stsw] = w_stsw_wi;
					idx_stsw ++;
				}
				if(I_stns_wi==1 && idx_stns < Nemp ){
					w_stns[idx_stns] = w_stns_wi;
					idx_stns ++;
				}
				if(I_UEsw_wi==1 && idx_EUsw < NswU && !isnan(w_UEsw_wi)){
					w_UEsw[idx_UEsw] = w_UEsw_wi;
					w_EUsw[idx_EUsw] = w_EUsw_wi;
					idx_UEsw ++;
					idx_EUsw ++;

					w_MVswrec[ri_wv][idx_MVswrec[ri_wv]] = w_UEsw_wi;
					idx_MVswrec[ri_wv]++;
					w_MVswrec[ri_wv][idx_MVswrec[ri_wv]] = w_EUsw_wi;
					idx_MVswrec[ri_wv]++;
				}
				if(I_UEns_wi==1 && idx_EUns < Nspell && !isnan(w_UEns_wi)){
					w_UEns[idx_UEns] = w_UEns_wi;
					w_EUns[idx_EUns] = w_EUns_wi;
					idx_UEns ++;
					idx_EUns ++;

					w_MVnsrec[ri_wv][idx_MVnsrec[ri_wv]] = w_UEns_wi;
					idx_MVnsrec[ri_wv]++;
					w_MVnsrec[ri_wv][idx_MVnsrec[ri_wv]] = w_EUns_wi;
					idx_MVnsrec[ri_wv]++;
				}
			}
		}
	}

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
	}


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
	}
	int nonnum=0;
	for(ri=0;ri<idx_EEsw;ri++){
		if(!gsl_finite(wwv_EEsw[ri]) ) nonnum ++;
	}
	if( nonnum >0 )
		printf("Not a number on wwv_EEsw %d times \n ", nonnum);
	nonnum =0;
	for(ri=0;ri<idx_EEsw;ri++){
		if(!gsl_finite(wwv_occsw[ri]) ) nonnum ++;
	}
	if( nonnum >0 )
		printf("Not a number on wwv_occsw %d times \n ", nonnum);


	st->corr_wgocc= idx_EEsw>0 ? gsl_stats_correlation( wwv_EEsw,1,wwv_occsw ,1,(size_t)idx_EEsw) : 0.;
	//if(verbose>1){
	//	printf("The valid w observations are %d, valid EUE observations are %d \n",valid_w, valid_EUE);
		printf("corr of wage and occwg is %f, based on %d obs \n", st->corr_wgocc, idx_EEsw);
	//}

	w_qtls(w_stns,1,idx_stns,st->stns_qtls);
	w_qtls(w_stsw,1,idx_stsw,st->stsw_qtls);
	w_qtls(w_EEns,1,idx_EEns,st->EEns_qtls);
	w_qtls(w_EEsw,1,idx_EEsw,st->EEsw_qtls);
	w_qtls(w_EUns,1,idx_EUns,st->EUns_qtls);
	w_qtls(w_EUsw,1,idx_EUsw,st->EUsw_qtls);
	w_qtls(w_UEns,1,idx_UEns,st->UEns_qtls);
	w_qtls(w_UEsw,1,idx_UEsw,st->UEsw_qtls);

	double MVsw_qtls[2][Nqtls];
	double MVns_qtls[2][Nqtls];

	for(ri=0;ri<2;ri ++){
		w_qtls(w_MVswrec[ri],1,idx_MVswrec[ri],MVsw_qtls[ri]);
		w_qtls(w_MVnsrec[ri],1,idx_MVnsrec[ri],MVns_qtls[ri]);
	}
	for(ri=0;ri<Nqtls;ri++){
		st->MVsw_qtls_ratio[ri] = MVsw_qtls[1][ri]/MVsw_qtls[0][ri];
		st->MVns_qtls_ratio[ri] = MVns_qtls[1][ri]/MVns_qtls[0][ri];
	}

	for(ri=0;ri<2;ri++){
		free(w_MVnsrec[ri]);free(w_MVswrec[ri]);
	}
	free(w_MVnsrec);free(w_MVswrec);

	free(w_stns);
	free(w_stsw);
	free(w_EEns);
	free(w_EEsw);
	free(w_EUns);
	free(w_EUsw);
	free(w_UEns);
	free(w_UEsw);
}
void print_params(double *x, int n, struct cal_params * par){
	int ii, i,ji;

	ii =0;
	parhist = fopen(parhi_f,"a+");

	fprintf(parhist,"%8.3f,", par->alphaE0);
	fprintf(parhist,"%8.3f,", par->alphaU0);
	fprintf(parhist,"%8.3f,", par->alphaU1);
	fprintf(parhist,"%8.3f,", par->lambdaU0  );
	fprintf(parhist,"%8.3f,", par->lambdaES0);
	fprintf(parhist,"%8.3f,", par->lambdaEM0);
	fprintf(parhist,"%8.3f,", par->delta_avg);
	fprintf(parhist,"%8.3f,", par->zloss);

	fprintf(parhist,"%8.3f,", par->var_ze);
	fprintf(parhist,"%8.3f,", par->autoz);
	fprintf(parhist,"%8.3f,", par->var_pe);
	fprintf(parhist,"%8.3f,", par->autop);
	fprintf(parhist,"%8.3f,", par->var_ae);
	fprintf(parhist,"%8.3f,", par->autoa);
	fprintf(parhist,"%8.3f,", par->gdfthr);
	fprintf(parhist,"%8.3f,", par->stwupdate);
	fprintf(parhist,"%8.3f,", par->var_eps);

	if(eps_2emg==1){
		fprintf(parhist,"%8.3f,", par->lshape_eps);
		fprintf(parhist,"%8.3f,", par->rshape_eps);
	}

	fprintf(parhist,"%8.3f,", par->delta_Acoef);
	fprintf(parhist,"%8.3f,", par->lambdaEM_Acoef);
	fprintf(parhist,"%8.3f,", par->lambdaES_Acoef);
	fprintf(parhist,"%8.3f,", par->lambdaU_Acoef);

	for(ji=0;ji<JJ;ji++){
		if(ji%3==1){
			fprintf(parhist,"%8.3f,", par->AloadP->data[ji]);

		}
		if(ji%3==2){
			fprintf(parhist,"%8.3f", par->AloadP->data[ji]);
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


		par->alphaE0   = x[0];
		par->alphaU0   = x[1];
		par->alphaE1   = x[2];
		par->alphaU1   = x[2];
		par->lambdaU0  = x[3];
		par->lambdaES0 = x[4];
		par->lambdaEM0 = x[5];
		par->delta_avg = x[6];
		par->zloss     = x[7];
		ii = Npar_cluster[0] ;
	}
	if(ci==1 || ci == Ncluster){


		par->var_ze     = x[0+ii];
		par->autoz      = x[1+ii];
		par->var_pe     = x[2+ii];
		par->autop      = x[3+ii];
		par->var_ae     = x[4+ii];
		par->autoa      = x[5+ii];
		par->gdfthr     = x[6+ii];
		par->stwupdate  = x[7+ii];
		par->var_eps    = x[8+ii];
		if(eps_2emg==1){
			par->lshape_eps = x[9 +ii];
			par->rshape_eps = x[10+ii];
		}
		ii = Npar_cluster[0] + Npar_cluster[1];

	}
	if(ci==2 || ci == Ncluster){


		par->delta_Acoef    = x[0+ii];
		par->lambdaEM_Acoef = x[1+ii];
		par->lambdaES_Acoef = x[2+ii];
		par->lambdaU_Acoef  = x[3+ii];
		i=4;
		for(ji=0;ji<JJ;ji++){
			if(ji%3==0)
				par->AloadP->data[ji]= 0.;
			if(ji%3==1){
				par->AloadP->data[ji] = x[i+ii];
				i++;
			}
			if(ji%3==2){
				par->AloadP->data[ji] = x[i+ii-1]*x[i+ii];
				i++;
			}
		}
		ii += i;
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
}


void set_dat( struct stats * dat){

	int i;

	dat->J2Jprob = 0.03404871;
	dat->findrate = 0.3946638;
	dat->seprate  =  0.02227534;
	dat->swProb_EE = 0.5026808;
	dat->swProb_U = 0.5362486;
	dat->swProb_st = 0.01999708;
	dat->udur_nosw= 5.448572;
	dat->udur_sw=  6.397636;
	dat->corr_wgocc = 0.3872327;

    double stsw[] = {-0.33831643,-0.12185895 ,0.01786183,0.18345242,0.44689113};
    memcpy(dat->stsw_qtls ,stsw,Nqtls*sizeof(double));
    double stns[] = {-0.191465375, -0.065237539, -0.001523193, 0.086168956, 0.239770718};
    memcpy(dat->stns_qtls ,stns,Nqtls*sizeof(double));

    double EUsw[] = {-2.25325384, -1.23791693, -0.55627011, 0.07654097, 0.82149574 };
    memcpy(dat->EUsw_qtls ,EUsw,Nqtls*sizeof(double));
    double EUns[] = {-1.85892018, -1.06531582, -0.45595244, 0.02397662, 0.69029482 };
    memcpy(dat->EUns_qtls ,EUsw,Nqtls*sizeof(double));

    double UEsw[] = {-1.14731344, -0.52494189, 0.05439134, 0.76429646, 1.52761627};
    memcpy(dat->UEsw_qtls, UEsw,Nqtls*sizeof(double));
    double UEns[] = {-1.05169527, -0.50902791, -0.08222744, 0.50908986, 1.21678663};
    memcpy(dat->UEns_qtls, UEns,Nqtls*sizeof(double));

    double EEsw[] = {-0.6016007, -0.2139596, 0.1109101, 0.5378298, 1.1205155};
	memmove(dat->EEsw_qtls, EEsw,Nqtls*sizeof(double));
    double EEns[] = {-0.48119935, -0.17915801, 0.05450329, 0.34849799, 0.82659345};
    memmove(dat->EEns_qtls, EEns,Nqtls*sizeof(double));

    dat->seprt_ratio = 0.7460479;
    dat->fndrt_ratio =  1.087637;
    dat->EE_ratio =  1.184628;

    dat->swProb_EE_ratio = 0.973819;
    dat->swProb_U_ratio  = 0.957414;



    double datMVsw_qtls_ratio[] = {1.1714109, 1.2188297, 2.8478620, 0.5859163, 0.7557770};
    double datMVns_qtls_ratio[] = {1.2007256, 1.3091474, 2.4762342, 0.3698798, 0.6930928};

    memcpy(dat->MVns_qtls_ratio,datMVns_qtls_ratio,Nqtls*sizeof(double));
    memcpy(dat->MVsw_qtls_ratio,datMVsw_qtls_ratio,Nqtls*sizeof(double));
}

double param_dist( double * x, struct cal_params *par , int Npar, double * err_vec , int Nerr){

	// Takes in a parameters vector and pumps out an error vector

	int i,ii,ji;
	int success;

	struct valfuns vf;
	struct polfuns pf;
	struct shocks sk;
	struct hists ht;

	struct stats st;
	struct stats dat;

    dat.EEns_qtls = malloc(Nqtls*sizeof(double));
    dat.EEsw_qtls = malloc(Nqtls*sizeof(double));
    dat.EUns_qtls = malloc(Nqtls*sizeof(double));
    dat.EUsw_qtls = malloc(Nqtls*sizeof(double));
    dat.UEns_qtls = malloc(Nqtls*sizeof(double));dat.UEsw_qtls = malloc(Nqtls*sizeof(double));
    dat.stns_qtls = malloc(Nqtls*sizeof(double));dat.stsw_qtls = malloc(Nqtls*sizeof(double));
	dat.MVns_qtls_ratio = malloc(Nqtls*sizeof(double));
	dat.MVsw_qtls_ratio = malloc(Nqtls*sizeof(double));


    // set data moments and parameter values
    set_dat(&dat);
    set_params( x, Npar, par, par->cluster_hr);
	print_params(x, Npar, par);

	st.EEns_qtls = malloc(Nqtls*sizeof(double));st.EEsw_qtls = malloc(Nqtls*sizeof(double));
	st.EUns_qtls = malloc(Nqtls*sizeof(double));st.EUsw_qtls = malloc(Nqtls*sizeof(double));
	st.UEns_qtls = malloc(Nqtls*sizeof(double));st.UEsw_qtls = malloc(Nqtls*sizeof(double));
	st.stns_qtls = malloc(Nqtls*sizeof(double));st.stsw_qtls = malloc(Nqtls*sizeof(double));
	st.MVns_qtls_ratio = malloc(Nqtls*sizeof(double));
	st.MVsw_qtls_ratio = malloc(Nqtls*sizeof(double));

	allocate_mats(&vf,&pf,&ht,&sk);

	init_pf(&pf,par);
	init_vf(&vf,par);
	//set the markov transition matrix for P
	rouwenhorst(par->autop,pow(par->var_pe,0.5),par->Ptrans[0],par->Plev);
	for(i=0;i<JJ;i++){
		gsl_matrix_memcpy(par->Ptrans[i],par->Ptrans[i]);
	}
	// set the markov transition matrix for A
	rouwenhorst(par->autoa,pow(par->var_ae,0.5),par->Atrans,par->Alev);

	// setup the z matrix
	gsl_matrix_view ztransLower =gsl_matrix_submatrix(par->ztrans,1,1,NZ-1,NZ-1);
	gsl_vector_view zlevLower = gsl_vector_subvector(par->zlev,1,NZ-1);
	gsl_vector_view zprobLower = gsl_vector_subvector(par->zprob,1,NZ-1);

	rouwenhorst(par->autoz ,pow(par->var_ze,0.5),& (ztransLower.matrix), &(zlevLower.vector) );
	ergod_dist( &(ztransLower.matrix) , &(zprobLower.vector));
	gsl_vector_set(par->zlev,0,5*par->zlev->data[1]);
	gsl_vector_set(par->zprob,0,0.);
	for(ii=0;ii<NZ;ii++){
		gg_set(par->ztrans,ii,0,par->zloss);
		gg_set(par->ztrans,0,ii,gsl_vector_get(par->zprob,ii));
		gsl_vector_view zrow = gsl_matrix_row(par->ztrans,ii);
		if(ii>0){
			double zrowsum = 0;
			int ci;
			for( ci=0;ci<NZ;ci++ ) zrowsum += gsl_vector_get( &(zrow.vector),ci );
			gsl_vector_scale(  & (zrow.vector) ,1./zrowsum);
		}
	}

	if(eps_2emg==1){
		success = disc_2emg(par->epsprob->data,par->epslev->data,(int)par->epsprob->size,
				0.,par->var_eps,par->lshape_eps,par->rshape_eps);
	}else{
		success =0;
		double epsub[NE-1];

		epsub[0] =gsl_cdf_gaussian_Pinv( 1./(double)NE,par->var_eps );
		for(ii=1;ii<NE-1;ii++)
			epsub[ii] =gsl_cdf_gaussian_Pinv( (double)(ii+1)/(double)NE,par->var_eps );
		for(ii=1;ii<NE-1;ii++){
			par->epslev->data[ii] =(epsub[ii] + epsub[ii-1])/2;
		}
		par->epslev->data[0] = 2*epsub[0] -par->epslev->data[1];
		par->epslev->data[NE-1] = 2*epsub[NE-2] -par->epslev->data[NE-2];
		gsl_vector_set_all(par->epsprob,1./(double)NE);
	}
	if(success > 0){
		printf(" Did not compute the distribution properly");
	}
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// come back and figure out how to do xS and xG !!!!
	for(i=0;i<NG;i++) gsl_vector_set(par->xGlev,i,-0.2  + .4* (double)i/(double) (NG-1));
	for(i=0;i<NS;i++) gsl_vector_set(par->xSlev,i,-0.2  + .4* (double)i/(double) (NS-1));
	for(ji=0;ji<JJ;ji++) gsl_vector_set(par->jprob,ji,1./(double)JJ);


	// ensure don't voluntarily quit, except for when z is the lowest value
	double w_hr;
	for(ji=0;ji<JJ;ji++){
		for (ii = 0; ii < NN; ii++) {
			int ai = ii / (NP * NG * NS * NZ * NE);
			int pi = (ii - ai * NP * NG * NS * NZ * NE) / (NG * NS * NZ * NE);
			int gi = (ii - ai * NP * NG * NS * NZ * NE - pi * NG * NS * NZ * NE) / (NS * NZ * NE);
			int si = (ii - ai * NP * NG * NS * NZ * NE - pi * NG * NS * NZ * NE - gi * NS * NZ * NE) / (NZ * NE);
			int zi =
					(ii - ai * NP * NG * NS * NZ * NE - pi * NG * NS * NZ * NE - gi * NS * NZ * NE - si * NZ * NE) / NE;
			int ti = ii - ai * NP * NG * NS * NZ * NE - pi * NG * NS * NZ * NE - gi * NS * NZ * NE - si * NZ * NE -
			         zi * NE;

			w_hr = exp(par->AloadP->data[ji] * par->Alev->data[ai] +
			                          par->Plev->data[pi] +
			                          par->epslev->data[ti] +
			                          par->zlev->data[zi] +
			                          par->xSlev->data[si] +
			                          par->xGlev->data[gi]) + wage_lev;
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
		printmat("ztrans.csv",par->ztrans,0);
		printvec("Plev.csv", par->Plev,0);
		printvec("epsprob.csv", par->epsprob,0);
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
		double mod_dur = st.udur_sw / st.udur_nosw;
		if (par->cluster_hr == 0 || par->cluster_hr == Ncluster) {
			err_vec[0] = (st.J2Jprob - dat.J2Jprob) * 2 / (st.J2Jprob + dat.J2Jprob);
			err_vec[1] = (st.findrate - dat.findrate) * 2 / (st.findrate + dat.findrate);
			err_vec[2] = (st.seprate - dat.seprate) * 2 / (st.seprate + dat.seprate);
			err_vec[3] = (st.swProb_EE - dat.swProb_EE) * 2 / (st.swProb_EE + dat.swProb_EE);
			err_vec[4] = (st.swProb_U - dat.swProb_U) * 2 / (st.swProb_U + dat.swProb_U);
			err_vec[5] = (st.swProb_st - dat.swProb_st) * 2 / (st.swProb_st + dat.swProb_st);

			err_vec[6] = (mod_dur - dat_dur) * 2 / (mod_dur + dat_dur);

			err_vec[7] = (st.corr_wgocc - dat.corr_wgocc); //*2 / (st.corr_wgocc + dat.corr_wgocc);
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
			err_vec[ii + 0] = (st.swProb_U_ratio - dat.swProb_U_ratio) * 2 / (st.swProb_U_ratio + dat.swProb_U_ratio);
			err_vec[ii + 1] =
					(st.swProb_EE_ratio - dat.swProb_EE_ratio) * 2 / (st.swProb_EE_ratio + dat.swProb_EE_ratio);
			err_vec[ii + 2] = (st.fndrt_ratio - dat.fndrt_ratio) * 2 / (st.fndrt_ratio + dat.fndrt_ratio);
			err_vec[ii + 3] = (st.seprt_ratio - dat.seprt_ratio) * 2 / (st.seprt_ratio + dat.seprt_ratio);
			err_vec[ii + 4] = (st.EE_ratio - dat.EE_ratio) * 2 / (st.EE_ratio + dat.EE_ratio);
			//err_vec[ii + 5] = (st.corr_wgocc - dat.corr_wgocc) * 2 / (st.corr_wgocc + dat.corr_wgocc);
			ii +=5;
			for(i=0;i<Nqtls; i++)
				err_vec[ii+i] = (st.MVsw_qtls_ratio[i] - dat.MVsw_qtls_ratio[i])*2/(st.MVsw_qtls_ratio[i] + dat.MVsw_qtls_ratio[i])/(double)Nqtls;
			ii += Nqtls;
			for(i=0;i<Nqtls; i++)
				err_vec[ii+i] = (st.MVns_qtls_ratio[i] -dat.MVns_qtls_ratio[i])*2/(st.MVns_qtls_ratio[i] + dat.MVns_qtls_ratio[i]) /(double)Nqtls;

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
	}else{
		gsl_rng * RNG = gsl_rng_alloc(gsl_rng_default);
		for(i=0;i<Nerr;i++){
			err_vec[i] = gsl_rng_uniform(RNG);
		}
		gsl_rng_free(RNG);

	}
	double quad_dist =0;
	for(i=0;i<Nerr;i++)
		quad_dist += err_vec[i]*err_vec[i];

	free_mats(&vf,&pf,&ht,&sk);

	free(dat.EEns_qtls);free(dat.EEsw_qtls);
	free(dat.EUns_qtls);free(dat.EUsw_qtls);
	free(dat.UEns_qtls);free(dat.UEsw_qtls);
	free(dat.stns_qtls);free(dat.stsw_qtls);
	free(dat.MVns_qtls_ratio); free(dat.MVsw_qtls_ratio);

	free(st.EEns_qtls);free(st.EEsw_qtls);
	free(st.EUns_qtls);free(st.EUsw_qtls);
	free(st.UEns_qtls);free(st.UEsw_qtls);
	free(st.stns_qtls);free(st.stsw_qtls);

	free(st.MVns_qtls_ratio); free(st.MVsw_qtls_ratio);
	return(quad_dist);

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

void dfovec_iface_(double * f, double * x, int * n, int* m){
	// this is going to interface with dfovec.f and call param_dist using the global params
	int nv = *n;
	int mv = *m;
	int i,ci;
	int verbose_old,print_lev_old;
	double dist;

	if(verbose>1) printf("Entering DFBOLS evaluation\n");

	verbose_old = verbose;
	print_lev_old = print_lev;

	verbose=0;
	print_lev=0;
	ci = glb_par->cluster_hr;
	double * x_pspace = malloc(sizeof(double)*nv);
	for(i=0;i<nv;i++) // convert from [0,1] domain to parameter space
		x_pspace[i] = x[i] * (glb_par->cluster_ub[ci][i]-glb_par->cluster_lb[ci][i])
				+glb_par->cluster_lb[ci][i];
	if(verbose_old>1){
		printf("Evaluating at (in parameter space): ");
		for(i=0;i<nv-1;i++)
			printf("%f, ", x_pspace[i]);
		printf("%f\n", x_pspace[nv-1]);
	}

	dist = param_dist( x_pspace, glb_par, nv ,f, mv );

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
			fprintf(calhist,"%f,",f[i]);
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
		for(i=0;i<nv;i++) solver_state[i+1] = x_pspace[i];
	}

	free(x_pspace);
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
				gg_set( sk->lambdaMsel[ll],i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->lambdaSsel[ll],i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->lambdaUsel[ll],i,ti, gsl_rng_uniform(rng0) );
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

void allocate_pars( struct cal_params * par){
	int ji;
	TTT = TT + burnin;

	par->Plev    = gsl_vector_calloc(NP) ;
	par->Alev    = gsl_vector_calloc(NA) ;
	par->xGlev   = gsl_vector_calloc(NG) ;
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
	par->xGtrans = gsl_matrix_calloc(NG,NG);
	par->xStrans = gsl_matrix_calloc(NS,NS);
	par->ztrans = gsl_matrix_calloc(NZ,NZ);
	par->epsprob = gsl_vector_calloc(NE);
	par->zprob = gsl_vector_calloc(NZ);
	par->jprob = gsl_vector_calloc(JJ);

}


void allocate_mats( struct valfuns * vf, struct polfuns * pf, struct hists * ht, struct shocks * sk ){

	int j;

	NN = NA*NP*NG*NS*NZ*NE;
	NUN = NA*NP*NG*NS*NZ;
	TTT = TT + burnin;

	alloc_valfuns(vf);
	alloc_hists(ht);
	pf -> mE =  gsl_matrix_calloc( NN,JJ );
	pf -> mU =  gsl_matrix_calloc( NUN,JJ );
	pf -> sE = malloc(sizeof(gsl_matrix*)*JJ);
	pf -> sU = malloc(sizeof(gsl_matrix*)*JJ);
	for (j = 0; j <JJ ; j++) {
		pf->sE[j] = gsl_matrix_calloc(NN,JJ);
		pf->sU[j] = gsl_matrix_calloc(NUN,JJ);
	}
	sk->Asel       = malloc(sizeof(gsl_vector*)*Npaths);
	sk->lambdaMsel = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->lambdaSsel = malloc(sizeof(gsl_matrix*)*Npaths);
	sk->lambdaUsel = malloc(sizeof(gsl_matrix*)*Npaths);
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
		sk->lambdaMsel[j] = gsl_matrix_calloc(Nsim,TTT);
		sk->lambdaSsel[j] = gsl_matrix_calloc(Nsim,TTT);
		sk->lambdaUsel[j] = gsl_matrix_calloc(Nsim,TTT);
		sk->epssel[j]     = gsl_matrix_calloc(Nsim,TTT);
		sk->xGsel[j]      = gsl_matrix_calloc(Nsim,TTT);
		sk->xSsel[j]      = gsl_matrix_calloc(Nsim,TTT);
		sk->zsel[j]       = gsl_matrix_calloc(Nsim,TTT);
		sk->jsel[j]       = gsl_matrix_calloc(Nsim,TTT);
		sk->dsel[j]       = gsl_matrix_calloc(Nsim,TTT);
		sk->msel[j]       = gsl_matrix_calloc(Nsim,TTT);
	}
}

void alloc_valfuns(struct valfuns *vf ){

	vf -> WE =  gsl_matrix_calloc( NN,JJ );
	vf -> WU =  gsl_matrix_calloc( NUN,JJ );
	vf -> RE =  gsl_matrix_calloc( NN,JJ );
	vf -> RU =  gsl_matrix_calloc( NUN,JJ );

	vf -> WEdist = gsl_matrix_calloc( NN,JJ );

}

void alloc_hists( struct hists *ht ){
	int ll,ji;

	ht->xhist = malloc(sizeof(gsl_matrix_int **)*Npaths);
	for(ll=0;ll<Npaths;ll++)
		ht->xhist[ll] = malloc(sizeof(gsl_matrix_int *)*4);

	ht->uhist = malloc(sizeof(gsl_matrix_int *)*Npaths);
	ht->whist = malloc(sizeof(gsl_matrix      *)*Npaths);
	ht->jhist = malloc(sizeof(gsl_matrix_int *)*Npaths);
	ht->Ahist = malloc(sizeof(gsl_vector_int *)   *Npaths);
	ht->Phist = malloc(sizeof(gsl_matrix_int *)   *Npaths);
	ht->J2Jhist = malloc(sizeof(gsl_matrix_int *) *Npaths);

	for(ll=0;ll<Npaths;ll++){
		ht->uhist[ll] = gsl_matrix_int_calloc(Nsim,TT);
		ht->whist[ll] = gsl_matrix_calloc(Nsim,TT);
		for(ji=0;ji<4;ji++)
			ht->xhist[ll][ji] = gsl_matrix_int_calloc(Nsim,TT);
		ht->jhist[ll] = gsl_matrix_int_calloc(Nsim,TT);
		ht->Ahist[ll] = gsl_vector_int_calloc(TT);
		ht->Phist[ll] = gsl_matrix_int_calloc(Nsim,TT);
		ht->J2Jhist[ll] = gsl_matrix_int_calloc(Nsim,TT);
	}
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
	for(j=0;j<Npaths;j++){
		gsl_vector_free(sk->Asel[j]);
		gsl_matrix_free(sk->Psel[j]);
		gsl_matrix_free(sk->lambdaMsel[j]);
		gsl_matrix_free(sk->lambdaSsel[j]);
		gsl_matrix_free(sk->lambdaUsel[j]);
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
	free(sk->lambdaSsel);
	free(sk->lambdaUsel);
	free(sk->lambdaMsel);

	free(sk->epssel);
	free(sk->xGsel);
	free(sk->xSsel);
	free(sk->zsel);
	free(sk->jsel);
	free(sk->dsel);
	free(sk->msel);

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
	int ji;

	gsl_vector_free(par->Plev);
	gsl_vector_free(par->Alev);
	gsl_vector_free(par->xGlev);
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
	gsl_matrix_free(par->xGtrans);
	gsl_matrix_free(par->xStrans);
	gsl_matrix_free(par->ztrans);
	gsl_vector_free(par->epsprob);
	gsl_vector_free(par->zprob);
	gsl_vector_free(par->jprob);

}

void free_hists( struct hists *ht ){
	int ll,ji;
	for(ll=0;ll<Npaths;ll++){
		gsl_matrix_int_free(ht->uhist[ll]);
		gsl_matrix_free(ht->whist[ll]);
		for(ji=0;ji<4;ji++)
			gsl_matrix_int_free(ht->xhist[ll][ji]);
		gsl_matrix_int_free(ht->jhist[ll]);
		gsl_vector_int_free(ht->Ahist[ll]);
		gsl_matrix_int_free(ht->Phist[ll]);
		gsl_matrix_int_free(ht->J2Jhist[ll]);
	}

	free(ht->uhist);
	free(ht->whist);
	for(ll=0;ll<Npaths;ll++)
		free(ht->xhist[ll]);
	free(ht->xhist);
	free(ht->jhist);
	free(ht->Ahist);
	free(ht->Phist);
	free(ht->J2Jhist);



}
void init_pf( struct polfuns *pf ,struct cal_params * par){

	int ji;
	gsl_matrix_set_all(pf->mE,0.5);
	gsl_matrix_set_all(pf->mU,0.5);
	for(ji=0;ji<JJ;ji++) gsl_matrix_set_all(pf->sE[ji],1./(double)(JJ-1));
	for(ji=0;ji<JJ;ji++) gsl_matrix_set_all(pf->sU[ji],1./(double)(JJ-1));

}
void init_vf( struct valfuns *vf ,struct cal_params * par){
	int ji;

	gsl_matrix_set_all( vf->WE, 0.);
	gsl_matrix_set_all( vf->WU, 0.);
	gsl_matrix_set_all( vf->RE, 0.);
	gsl_matrix_set_all( vf->RU, 0.);



}
