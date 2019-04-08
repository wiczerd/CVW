//gcc -fopenmp main.c -lgsl -lgslcblas -lm -o CVW.out
// valgrind --leak-check=yes --error-limit=no --track-origins=yes --log-file=CVWvalgrind.log ./CVW_dbg.out &
#ifdef _MKL_USE
#include "mkl_lapacke.h"
#endif
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
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_statistics.h>
// declare some parameters that will be global scope
int static JJ = 3;     // number of occupations
int static NG = 2;     //number of human capital types
int static NS = 2;     //number of occ tenure types
int static NZ = 7;     //number of occ match-quality types
int static NT = 3;     //number of firm theta match-quality types
int static NP = 5;     //number of occupation-specific productivies
int static NA = 5;     //number of aggregate productivities

int NN ,NUN;

int static TT      = 12*15;    // periods per simulation path
int static burnin  = 48;       // number of periods to throw away each time
int static TTT ;
int static Npaths  = 50;      // number of simulation paths to draw
int static Nsim    = 2000;

int static Nqtls   = 5; // number of quantiles that will use to compare distributions

int verbose = 3;
int print_lev = 3;

int maxiter = 5000;
double vftol = 1e-3;
double rhotightening = .5;

double beta	= 0.997;		// discount factor
double b 	= 0.0; 			// unemployment benefit
double wage_lev = 1;        // will be a shifter so the average wage is >0

double delta_avg = .01;     // average separation rate
double urt_avg = .055;     // average separation rate

struct cal_params{
	double gdfthr, lambdaEM0, lambdaES0, lambdaU0, var_les, var_lem, var_lu ;
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
	double wage_curve;  // curviness of wage function
    double delta_Acoef;
    double lambdaU_Acoef,lambdaES_Acoef,lambdaEM_Acoef;
    double zloss;

	gsl_vector * AloadP; //loading on A for each P
	gsl_vector * Plev;
	gsl_vector * Alev;
	gsl_vector * xGlev;
	gsl_vector * xSlev;
	gsl_vector * zlev;
	gsl_vector * tlev;

	gsl_matrix ** Ptrans; // markov transition matrix for P, for each J
	gsl_matrix * Atrans;
	gsl_matrix * xGtrans;
	gsl_matrix * xStrans;
	gsl_matrix * ztrans;
	gsl_vector * zprob; // distribution from which to draw when z first drawn
	gsl_vector * tprob; // iid, so just a vector of probabilities
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
	gsl_matrix ** thsel;
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
	double swProb_EE;
	double J2Jprob;
	double unrate;
	double findrate;
	double seprate;

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

double param_dist( struct cal_params *par, int Npar, double * err_vec , int Nerr);

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
	if(Nqtls ==5){
		qtls_out[0] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.1);
		qtls_out[1] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.25);
		qtls_out[2] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.5);
		qtls_out[3] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.75);
		qtls_out[4] = gsl_stats_quantile_from_sorted_data(vec, (size_t)stride ,(size_t)len, 0.9);
	}
	else{
		int qi ;
		for(qi = 0; qi<Nqtls;qi++)
			qtls_out[qi] = gsl_stats_quantile_from_sorted_data(vec,(size_t)stride,(size_t)len, (double)(qi+1)/(double)(Nqtls+1)  );
	}
}

int main() {

	int i,ii,ji;
	int success;

	struct cal_params par;
	struct valfuns vf;
	struct polfuns pf;
	struct shocks sk;
	struct hists ht;

	struct stats st;

	st.EEns_qtls = malloc(Nqtls*sizeof(double));st.EEsw_qtls = malloc(Nqtls*sizeof(double));
	st.EUns_qtls = malloc(Nqtls*sizeof(double));st.EUsw_qtls = malloc(Nqtls*sizeof(double));
	st.UEns_qtls = malloc(Nqtls*sizeof(double));st.UEsw_qtls = malloc(Nqtls*sizeof(double));
	st.stns_qtls = malloc(Nqtls*sizeof(double));st.stsw_qtls = malloc(Nqtls*sizeof(double));

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
	gsl_matrix_view ztransLower =gsl_matrix_submatrix(par.ztrans,1,1,NZ-1,NZ-1);
	gsl_vector_view zlevLower = gsl_vector_subvector(par.zlev,1,NZ-1);
	gsl_vector_view zprobLower = gsl_vector_subvector(par.zprob,1,NZ-1);
	par.var_ze = 0.025*0.025; par.autoz = 0.975;
	par.zloss  = 0.025;
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
			for( int ci=0;ci<NZ;ci++ ) zrowsum += gsl_vector_get( &(zrow.vector),ci );
			gsl_vector_scale(  & (zrow.vector) ,1./zrowsum);
		}
	}


	//gsl_matrix_set_all(par.ztrans, 0.025/( (double)NZ-1. ));
	//for(ii=0;ii<NZ;ii++) gsl_matrix_set(par.ztrans,ii,ii,0.975);

	gsl_matrix_set_all(par.xStrans,.025/( (double)NS-1.) );
	for(ii=0;ii<NS;ii++) gg_set(par.xStrans,ii,ii,.975);
	gsl_vector_set_all(par.tprob, 1./(double) NT);
	gsl_vector_set_all(par.zprob, 1./((double) NZ));
	gsl_vector_set_all(par.AloadP,1.0);
	par.autoa = 0.95;par.var_ae = 0.01*0.01;
	rouwenhorst(par.autoa,pow(par.var_ae,0.5),par.Atrans,par.Alev);

	par.autop = 0.95; par.var_pe = 0.02*0.2;
	rouwenhorst(par.autop,pow(par.var_pe,0.5),par.Ptrans[0],par.Plev);
	for(i=1;i<JJ;i++){
		gsl_matrix_memcpy(par.Ptrans[i],par.Ptrans[0]);
	}

//	for(i=0;i<NP;i++) gsl_vector_set(par.Plev ,i,-0.05 + .1* (double)i/(double) (NP-1));
	for(i=0;i<NG;i++) gsl_vector_set(par.xGlev,i,-0.2  + .4* (double)i/(double) (NG-1));
	for(i=0;i<NS;i++) gsl_vector_set(par.xSlev,i,-0.2  + .4* (double)i/(double) (NS-1));
	for(i=0;i<NT;i++) gsl_vector_set(par.tlev ,i,-0.1  + .2* (double)i/(double) (NT-1));
	for(ji=0;ji<JJ;ji++) gsl_vector_set(par.jprob,ji,1./(double)JJ);
	par.alphaE1 = 0.5;
	par.alphaE0 = 0.05*pow((double)(JJ-1),-par.alphaE1);
	par.alphaU1 = 0.5;
	par.alphaU0 = 1.*pow((double)(JJ-1),-par.alphaU1);
	par.lambdaU0  = 0.2 / 0.5;
	par.lambdaES0 = 0.01;
	par.lambdaEM0 = 0.8;
	par.kappa     = 0.01 ;
	par.gdfthr    = 0.5 ;
	par.wage_curve= 1.5 ;
	par.delta_Acoef = 0.;
	par.lambdaU_Acoef = 0.0;
    par.lambdaES_Acoef = 0.0;
	par.lambdaEM_Acoef = 0.0;


	int ai = 0;int pi = 0;int gi = 0;int si = 0;int zi = 1;	int thi = 0; ji=0;

	double wage_lev0 = pow(exp(par.AloadP->data[ji] * par.Alev->data[ai]) *
	                      exp(par.Plev->data[pi]) *
	                      exp(par.tlev->data[thi]) *
	                      exp(par.zlev->data[zi]) *
	                      exp(par.xSlev->data[si]) *
	                      exp(par.xGlev->data[gi]), 1.-par.wage_curve) / (1. - par.wage_curve);
	wage_lev = wage_lev - wage_lev0 ;



	// print out the grids
	if(print_lev>1){
		printvec("Alev.csv", par.Alev,0);printvec("zlev.csv", par.zlev,0);
		printmat("Atrans.csv",par.Atrans,0);printmat("Ptrans.csv",par.Ptrans[0],0);
		printmat("ztrans.csv",par.ztrans,0);
		printvec("Thlev.csv", par.tlev,0);printvec("Plev.csv", par.Plev,0);
	}

	success = draw_shocks(&sk);
	if(verbose>2 && success != 0) printf("Problem drawing shocks\n");

	success = sol_dyn( &par, &vf, &pf, &sk);
	if(verbose>2 && success != 0) printf("Problem solving model\n");

	success = sim(&par, &vf, &pf, &ht, &sk);
	if(verbose>2 && success != 0) printf("Problem simulating model\n");

	success = sum_stats(&par,&vf,&pf,&ht,&sk, &st);

	if(verbose>=2){
		printf(" unemployment rate is %f \n", st.unrate);
		printf(" J2J rate is %f  \n", st.J2Jprob);
		printf(" finding and separation rates are %f, %f  \n", st.findrate, st.seprate);
		printf(" occupation switching rates conditional on EE or EUE are %f, %f  \n", st.swProb_EE, st.swProb_U);


	}

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
			int ai = ii / (NP * NG * NS * NZ * NT);
			int pi = (ii - ai * NP * NG * NS * NZ * NT) / (NG * NS * NZ * NT);
			int gi = (ii - ai * NP * NG * NS * NZ * NT - pi * NG * NS * NZ * NT) / (NS * NZ * NT);
			int si = (ii - ai * NP * NG * NS * NZ * NT - pi * NG * NS * NZ * NT - gi * NS * NZ * NT) / (NZ * NT);
			int zi =
					(ii - ai * NP * NG * NS * NZ * NT - pi * NG * NS * NZ * NT - gi * NS * NZ * NT - si * NZ * NT) / NT;
			int ti = ii - ai * NP * NG * NS * NZ * NT - pi * NG * NS * NZ * NT - gi * NS * NZ * NT - si * NZ * NT -
			         zi * NT;

			wagevec[ii][ji] = pow(exp(par->AloadP->data[ji] * par->Alev->data[ai]) *
			                  exp(par->Plev->data[pi]) *
			                  exp(par->tlev->data[ti]) *
			                  exp(par->zlev->data[zi]) *
			                  exp(par->xSlev->data[si]) *
			                  exp(par->xGlev->data[gi]), 1.-par->wage_curve) / (1. - par->wage_curve) + wage_lev;
		}
	}

	for(ji=0;ji<JJ;ji++) {
		for (ii = 0; ii < NN; ii++) {
			gg_set(vf0.WE, ii, ji,
			               wagevec[ii][ji] / (1. - beta));
		}
	}
	for(ii=0;ii<NUN;ii++){
		gg_set(vf0.WU,ii,0,
		               (1.-par->lambdaU0)*b/(1.-beta) + par->lambdaU0*wagevec[ii*NT+NT/2][0]/(1.-beta));
		for(ji=1;ji<JJ;ji++)  gg_set(vf0.WU,ii,ji, gg_get(vf0.WU,ii,0));
	}


	for(viter = 0;viter<maxiter;viter++){

		for(ji=0;ji<JJ;ji++){
			#pragma omp parallel for private(ii) firstprivate(ji)
			for(ii=0;ii<NN;ii++){
				int ai = ii/(NP*NG*NS*NZ*NT);
				int pi = (ii - ai*NP*NG*NS*NZ*NT)/(NG*NS*NZ*NT);
				int gi = (ii - ai*NP*NG*NS*NZ*NT - pi*NG*NS*NZ*NT)/(NS*NZ*NT);
				int si = (ii - ai*NP*NG*NS*NZ*NT - pi*NG*NS*NZ*NT - gi*NS*NZ*NT)/(NZ*NT) ;
				int zi = (ii - ai*NP*NG*NS*NZ*NT - pi*NG*NS*NZ*NT - gi*NS*NZ*NT -  si*NZ*NT)/NT ;
				int ti = ii -  ai*NP*NG*NS*NZ*NT - pi*NG*NS*NZ*NT - gi*NS*NZ*NT - si*NZ*NT - zi*NT;

				int iU = ai*NP*NG*NS*NZ + pi*NG*NS*NZ + gi*NS*NZ + si*NZ +zi;

                double lambdaEMhr = par->lambdaEM0 + par->lambdaEM_Acoef*par->Alev->data[ai];
                double lambdaEShr = par->lambdaES0 + par->lambdaES_Acoef*par->Alev->data[ai];

				double delta_hr = delta_avg + par->delta_Acoef * gsl_vector_get(par->Alev,ai) ;
				//compute expectations over A, Pt
				double EAPWE = 0.;
				int aai, ppi;
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++)  EAPWE += gsl_max(gg_get( vf0.WE,  aai*NP*NG*NS*NZ*NT + ppi*NG*NS*NZ*NT+ gi*NS*NZ*NT + si*NZ*NT +zi*NT+ ti ,ji),
					                                          gg_get( vf0.WU,  aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ   + gi*NS*NZ    + si*NZ    +zi        ,ji))*
							gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
				}

				double EtTWE = 0.;
				int tti,zzi;
				double ttiprod = 0.;
				for(tti=ti;tti<NT;tti++) ttiprod+=gsl_vector_get(par->tprob,tti); //make sure tpobs sum to 1
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++){
						for(tti=ti;tti<NT;tti++)
							EtTWE +=  gsl_max(gg_get(vf0.WE, aai*NP*NG*NS*NZ*NT + ppi*NG*NS*NZ*NT+ gi*NS*NZ*NT + si*NZ*NT +zi*NT+ tti ,ji) ,
							                  gg_get(vf0.WU, aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ   + gi*NS*NZ    + si*NZ    +zi         ,ji))*
									gsl_vector_get(par->tprob,tti)/ttiprod *gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
					}
				}
				double EtWE = 0.;
				for(aai=0;aai<NA;aai++) {
					for (ppi = 0; ppi < NP; ppi++) {

						for (tti = 0; tti < NT; tti++)
							EtWE += gsl_max(gg_get(vf0.WE, aai*NP*NG*NS*NZ*NT + ppi*NG*NS*NZ*NT + gi*NS*NZ*NT + si*NZ*NT + zi*NT + tti, ji) ,
									        gg_get(vf0.WU, aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ    + gi*NS*NZ    + si*NZ    + zi         , ji))*
							        par->tprob->data[tti]*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
					}
				}

				int jji;
				double REhr = -par->kappa;
				double EzWE[JJ];
				double EztWE[JJ];
				for(jji=0;jji<JJ;jji++){
					EzWE[jji]  = 0.;
					EztWE[jji] = 0.;
					if(jji!=ji){

						for(aai=0;aai<NA;aai++){
							for(ppi=0;ppi<NP;ppi++) {
								for (zzi = 0; zzi < NZ; zzi++)
									EzWE[jji] += gsl_max(gg_get(vf0.WE,
									                            aai*NP*NG*NS*NZ*NT + ppi*NG*NS*NZ*NT +gi*NS*NZ*NT + zzi*NT + ti, jji) ,
									                     gg_get(vf0.WU,
									                            aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ    +gi*NS*NZ    + zzi        , jji)  )*
									             gsl_vector_get(par->zprob, zzi)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);

							}
						}
						for(aai=0;aai<NA;aai++){
							for(ppi=0;ppi<NP;ppi++) {
								for(tti=0;tti<NT;tti++){
									for(zzi=0;zzi<NZ;zzi++)
										EztWE[jji] += gsl_max(gg_get(vf0.WE,aai*NP*NG*NS*NZ*NT + ppi*NG*NS*NZ*NT+ gi*NS*NZ*NT + zzi*NT + tti,jji) ,
												              gg_get(vf0.WU,aai*NP*NG*NS*NZ    + ppi*NG*NS*NZ   + gi*NS*NZ    + zzi         ,jji))*
												(par->zprob->data[zzi])*(par->tprob->data[tti])*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji],pi,ppi);
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
					if(jji!=ji) sEjiDenom += pow(-par->kappa+ lambdaEMhr*EztWE[jji] + (1.-lambdaEMhr)*EzWE[jji]- EAPWE , 1/par->alphaE1);
				}
				for(jji=0;jji<JJ;jji++){
					if(jji!=ji){
						gg_set(pf->sE[jji] ,ii,ji, pow(-par->kappa+ lambdaEMhr*EztWE[jji] + (1.-lambdaEMhr)*EzWE[jji]- EAPWE , 1/par->alphaE1)/sEjiDenom);
					}else{ // this is kind of redundant because it should have initialized to 0
						gg_set(pf->sE[jji] ,ii,ji, 0.);
					}
				}


				// update the value function
				double WEhr = wagevec[ii][ji] + beta*delta_hr*gg_get(vf0.WU,iU,ji) +
				              beta*(1.- delta_hr )*(
				              		gg_get(pf->mE,ii,ji)*gg_get( vf->RE,ii,ji) +
				              		(1.-gg_get(pf->mE,ii,ji))*(par->gdfthr*lambdaEShr*EtWE + (1.-par->gdfthr)*lambdaEShr*EtTWE+
				                                                (1.-lambdaEShr)*EAPWE )  );

				gg_set( vf->WE, ii,ji,WEhr);

				gg_set( vf-> WEdist, ii,ji, gg_get(vf->WE,ii,ji) - gg_get(vf0.WE,ii,ji) );


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

				double lambdaUhr = par->lambdaU0 + par->lambdaU_Acoef*par->Alev->data[ai];

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
					if(jji!=ji){
						for(aai=0;aai<NA;aai++){
							for(ppi=0;ppi<NP;ppi++) {
								for( zzi=0;zzi<NZ;zzi++ )
									EzWU[jji] += gg_get(vf0.WU, aai*NP*NG*NS*NZ+ppi*NG*NS*NZ+gi*NS*NZ + zzi,jji)*
											gsl_vector_get(par->zprob,zzi)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[jji],pi,ppi);
							}
						}
					}
					RUhr += par-> alphaU0*pow(gg_get(pf->sU[jji],ii,ji),1.-par->alphaU1 )*EzWU[jji];
				}
				double totalalphaS = 0;
				for(jji=0;jji<JJ;jji++)
					totalalphaS += par-> alphaU0*pow(gg_get(pf->sU[jji],ii,ji),1.-par->alphaU1 );
				RUhr += (1.-totalalphaS )*EAPWU;
				gg_set(vf->RU, ii,ji, RUhr);

				double EtWE=0;
				for(aai=0;aai<NA;aai++){
					for(ppi=0;ppi<NP;ppi++) {
						for(tti=0;tti<NT;tti++){
							int iihr = aai*NP*NG*NS*NZ*NT+ppi*NG*NS*NZ*NT+gi*NS*NZ*NT+si*NZ*NT+zi*NT+tti;
							EtWE += gsl_max(gg_get(vf0.WE, iihr, ji), vf0.WU->data[ii*vf0.WU->tda+ji] )*
							        gsl_vector_get(par->tprob,tti)*gg_get(par->Atrans,ai,aai)*gg_get(par->Ptrans[ji], pi, ppi);}
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
					if(jji!=ji) sUdenom += pow(-par->kappa + EzWU[jji]- EAPWU ,1./par->alphaU1);
				}
				for(jji=0;jji<JJ;jji++){
					if(jji!=ji){
						gg_set(pf->sU[jji],ii,ji,pow(-par->kappa + EzWU[jji] - EAPWU ,1./par->alphaU1)/sUdenom );
					}else{//this is kind of redundant, because it should have initialized to 0
						gg_set(pf->sU[jji],ii,ji,0.);
					}
				}

				gg_set( vf->WU, ii, ji, b + beta*( pf->mU->data[ii*pf->mU->tda + ji]*vf->RU->data[ii*vf->RU->tda+ji] +
						(1.-pf->mU->data[ii*pf->mU->tda + ji])*( (1.-lambdaUhr)*EAPWU + lambdaUhr*EtWE ) ) );

			}
		}

		gsl_matrix_memcpy(vf0.WU,vf->WU);
		double maxdist = gsl_max( gsl_matrix_max(vf->WEdist),-gsl_matrix_min(vf->WEdist));
		if(viter % 200 == 0 && verbose >1 )  printf("Max distance is %f on iteration %d \n", maxdist,viter);

		if( maxdist < vftol ){
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

	double cmzprob[NZ];
	double cmtprob[NT];
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

	for(ll=0;ll<Npaths;ll++){
		swprob_hist[ll] = gsl_matrix_calloc(Nsim,TT);
		WE_hist[ll] = gsl_matrix_calloc(Nsim,TT);
		WU_hist[ll] = gsl_matrix_calloc(Nsim,TT);
	}

	cmjprob[0] = par->jprob->data[0];
	for(ji=0;ji<JJ-1;ji++) cmjprob[ji+1] = par->jprob->data[ji+1] + cmjprob[ji];
	cmtprob[0] = par->tprob->data[0];
	for(i=0;i<NT-1;i++) cmtprob[i+1] = par->tprob->data[i+1] + cmtprob[i];
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

	#pragma omp parallel for private(i,ti,ll,ji) firstprivate(cmzprob,cmtprob,cmjprob,cmAtrans,cmPtrans,cmxGtrans,cmxStrans,cmztrans)
	for(ll=0;ll<Npaths;ll++){
		int ** xt,**xtm1; //xG.xS.z.th
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
			for(thi=0;thi<NT;thi++) if( gg_get(sk->thsel[ll],i,0) > cmtprob[thi] ) ++ xt[i][3] ;

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
				lambdaEMhr[i] = (par->lambdaEM0 + par->lambdaEM_Acoef*At);
				lambdaEShr[i] = (par->lambdaES0 + par->lambdaES_Acoef*At);
				lambdaUhr[i]  = (par->lambdaU0  + par->lambdaU_Acoef *At);
			}

			for(i=0;i<Nsim;i++){
				jtm1[i] = jt[i];
				utm1[i] = ut[i];
				// increment individual specific shocks:
				for(xi=0;xi<4;xi++) xtm1[i][xi] = xt[i][xi];
				xt[i][0] = 0;
				for (gi = 0; gi < NG; gi++) if (gg_get(sk->xGsel[ll], i, 0) > cmxGtrans[xtm1[i][0]][gi]) ++ xt[i][0] ;
				xt[i][1] = 0;
				for (si = 0; si < NS; si++) if (gg_get(sk->xSsel[ll], i, 0) > cmxStrans[xtm1[i][1]][si]) ++ xt[i][1] ;
				xt[i][2] =0;
				for(zi=0;zi<NZ;zi++) if( gg_get(sk->zsel[ll],i,0) > cmztrans[xtm1[i][2]][zi] ) ++ xt[i][2] ;
				// xt[i][3]? no need to redraw theta unless later there's a job switch
				xt[i][3] = xtm1[i][3];

				ii = At*NP*NG*NS*NZ*NT + Pt[jt[i]]*NG*NS*NZ*NT + xt[i][0]*NS*NZ*NT+xt[i][1]*NZ*NT+xt[i][2]*NT+xt[i][3];
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

					if(ti>=burnin) {
						double wagehr = pow(exp(par->Alev->data[At]) + exp(par->Plev->data[Pt[jt[i]]]) *
						                    exp(par->xGlev->data[xt[i][0]]) *
						                    exp(par->xSlev->data[xt[i][1]]) *
						                    exp(par->zlev->data[xt[i][2]]) *
						                    exp(par->tlev->data[xt[i][3]]), 1. - par->wage_curve) /
						                (1. - par->wage_curve) + wage_lev;
						gg_set(ht->whist[ll], i, ti-burnin,wagehr);
					}
					// employed workers' choices
					double delta_hr =delta_avg + par->delta_Acoef * par->Alev->data[At];
					if( delta_hr >gg_get(sk->dsel[ll],i,ti) || gg_get(vf->WU,iU,jt[i]) > gg_get(vf->WE,ii,jt[i]) ){
						ut[i] = 1;
					}else{
						// stay or go?
						if( ti>=burnin ) gg_set(swprob_hist[ll],i,ti-burnin,  gg_get(pf->mE,ii,jt[i]) );
						if( gg_get(pf->mE,ii,jt[i])>gg_get(sk->msel[ll],i,ti) ){
							//RE
							double sij[JJ]; sij[0] =0.;

							for(jji=0;jji<JJ;jji++) sij[jji] = jji > 0 ? sij[jji-1] + gg_get(pf->sE[jji],ii,jt[i]) : gg_get(pf->sE[jji],ii,jt[i]);
							// successfully switched: if( gg_get(sk->jsel[ll],i,ti) <  sij[JJ])
							jt[i] = 0;
							for(jji=0;jji<JJ;jji++){
								if (gg_get(sk->jsel[ll], i, ti) > sij[jji]){
									jt[i] ++;
								}
							}
							if( jt[i] != jtm1[i] ){ // switchers
								xt[i][1] = 0; // lose specific skill
								// draw a new z:
								xt[i][2] =0;
								for(zi=0;zi<NZ;zi++){
									if( gg_get(sk->zsel[ll],i,ti) > cmzprob[zi] ) ++ xt[i][2] ;
								}
								if( gg_get(sk->lambdaMsel[ll],i,ti)<  lambdaEMhr[jt[i]] ) {
									if(ti>=burnin) ggi_set(ht->J2Jhist[ll], i, ti-burnin, 1);
									// draw a new theta
									for(thi =xtm1[i][3];thi<NT;thi++) if( gg_get( sk->thsel[ll],i,ti ) > cmtprob[thi] ) ++xt[i][3];
									xt[i][3] = xt[i][3]>NT-1 ? NT-1: xt[i][3];
								}
							}// else nothing happens (except paid kappa)


						}else{
							//WEm
							if( gg_get(sk->lambdaSsel[ll],i,ti) < lambdaEShr[jt[i]] ) {
								if (gg_get(sk->lambdaUsel[ll], i, ti) < par->gdfthr) { //godfather (gamma) shock?
									if(ti>=burnin) ggi_set(ht->J2Jhist[ll], i, ti-burnin, 1);
									xt[i][3] = 0;
									for (thi = 0; thi < NT; thi++)
										if (gg_get(sk->thsel[ll], i, ti) > cmtprob[thi])++xt[i][3];
									xt[i][3] = xt[i][3] > NT - 1 ? NT - 1 : xt[i][3];
								}else {
									if(ti>=burnin) ggi_set(ht->J2Jhist[ll], i, ti-burnin, 1);
									for (thi = xtm1[i][3]; thi < NT; thi++)
										if (gg_get(sk->thsel[ll], i, ti) > cmtprob[thi])++xt[i][3];
									xt[i][3] = xt[i][3] > NT - 1 ? NT - 1 : xt[i][3];
								}
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
						double sij[JJ];
						for(jji=0;jji<JJ;jji++) sij[jji] = jji>0 ? sij[jji-1] + gg_get(pf->sU[jji],ii,jt[i]): gg_get(pf->sU[jji],ii,jt[i]);
						// successfully switched: if( gg_get(sk->jsel[ll],i,ti) <  sij[JJ])
						jt[i] = 0;
						for(jji=0;jji<JJ;jji++)
							if (gg_get(sk->jsel[ll], i, ti) > sij[jji]) ++jt[i];

						if( jt[i] != jtm1[i] ){ // switchers
							xt[i][1] = 0; // lose specific skill
							// draw a new z:
							xt[i][2] =0;
							for(zi=0;zi<NZ;zi++){
								if( gg_get(sk->zsel[ll],i,ti) > cmzprob[zi] ) ++ xt[i][2] ;
							}
							if( gg_get(sk->lambdaUsel[ll],i,ti)<  lambdaUhr[jt[i]] ){
								// draw a new theta
								xt[i][3] = 0;
								for(thi =0;thi<NT;thi++) if( gg_get( sk->thsel[ll],i,ti ) > cmtprob[thi] ) ++xt[i][3];
								int iM = At*NP*NG*NS*NZ*NT + Pt[jt[i]]*NG*NS*NZ*NT + xt[i][0]*NS*NZ*NT+xt[i][1]*NZ*NT+xt[i][2]*NT+xt[i][3];
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
							for(thi =0;thi<NT;thi++) if( gg_get( sk->thsel[ll],i,ti ) > cmtprob[thi] ) ++xt[i][3];
							int iM = At*NP*NG*NS*NZ*NT + Pt[jt[i]]*NG*NS*NZ*NT + xt[i][0]*NS*NZ*NT+xt[i][1]*NZ*NT+xt[i][2]*NT+xt[i][3];
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
			printmat_int("thhist.csv",ht->xhist[ll][3],append );

			printmat("swprob_hist.csv", swprob_hist[ll],append );
			printmat("WUhist.csv", WU_hist[ll],append );
			printmat("WEhist.csv", WE_hist[ll],append );
			printmat_int("J2Jhist.csv", ht->J2Jhist[ll],append );
			append = 1;
		}
	}

	for(ll=0;ll<Npaths;ll++){
		gsl_matrix_free(swprob_hist[ll]);
		gsl_matrix_free(WE_hist[ll]);gsl_matrix_free(WU_hist[ll]);
	}
	free(swprob_hist);
	free(WE_hist);
	free(WU_hist);

}

int sum_stats(   struct cal_params * par, struct valfuns *vf, struct polfuns *pf, struct hists *ht, struct shocks *sk, struct stats *st  ){

	int ll, i,ti;

	int Nemp=0, Nunemp = 0, Nfnd =0, Nsep =0, NswU = 0, NswE=0, NJ2J=0, Nspell=0, NswSt=0;


#pragma parallel for private(ll,ti,i) reduction( +: Nemp ) reduction( +:Nunemp) reduction( +: Nfnd) reduction( +: Nsep) reduction( +: NswU) reduction( +: NswE) reduction( +: Nspell) reduction( +: NswSt)
	for(ll=0;ll<Npaths;ll++){

		int Nemp_ll = 0, Nunemp_ll = 0, Nfnd_ll =0, Nsep_ll = 0, NswU_ll = 0, NswE_ll=0, NJ2J_ll=0,Nspell_ll=0, NswSt_ll =0;
		int sw_spell=0;
		for(i=0;i<Nsim;i++){
			for(ti=0;ti<TT-1;ti++){
				if( ggi_get(ht->uhist[ll],i,ti) ==0 ){
					Nemp_ll += 1;
					if( ggi_get(ht->uhist[ll],i,ti+1) ==1 ){
						Nsep_ll +=1;
						sw_spell = ggi_get(ht->jhist[ll],i,ti);
					}
					if( ggi_get(ht->J2Jhist[ll],i,ti+1) ==1 ){
						NJ2J_ll +=1;
						if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti) ) NswE_ll +=1;
					}else{  // not EE
						if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti) ) NswSt_ll += 1;
					}
				}else{
					Nunemp_ll +=1;

					if( ggi_get(ht->uhist[ll],i,ti+1) ==0 ){
						Nfnd_ll +=1;
						Nspell_ll += 1;
						if( ggi_get(ht->jhist[ll],i,ti+1) != sw_spell ) NswU_ll += 1; //only count switches at the end of the spell.
					}
				}
			}
		}

		Nemp += Nemp_ll;
		Nunemp += Nunemp_ll;
		Nfnd += Nfnd_ll;
		Nsep += Nsep_ll;
		NswU += NswU_ll;
		NswE += NswE_ll;
		NJ2J += NJ2J_ll;
		Nspell += Nspell_ll;
		NswSt  += NswSt_ll;
	}

    st->J2Jprob  = (double) NJ2J / (double) Nemp ;
    st->findrate = (double) Nfnd / (double) Nunemp;
    st->seprate  = (double) Nsep / (double) Nemp ;
    st->swProb_EE = (double) NswE / (double) NJ2J;
    st->swProb_U = (double) NswU / (double) Nspell;
    st->unrate =  (double) Nunemp / (double) ( Nunemp + Nemp );


    double * w_stns = malloc(sizeof(double) *(Nemp-Nsep-NJ2J-NswSt) );
	double * w_stsw = malloc(sizeof(double) *(NswSt) );
	double * w_EEsw = malloc(sizeof(double) *(NswE) );
	double * w_EEns = malloc(sizeof(double) *(NJ2J - NswE) );
	double * w_EUsw = malloc(sizeof(double) *(NswU) );
	double * w_EUns = malloc(sizeof(double) *(Nspell - NswU) );
	double * w_UEsw = malloc(sizeof(double) *(NswU) );
	double * w_UEns = malloc(sizeof(double) *(Nspell - NswU) );


	int idx_EUns =0, idx_UEns =0, idx_EEns=0, idx_EUsw=0,idx_UEsw=0,idx_EEsw=0;
	int idx_stns =0, idx_stsw =0;
	for(ll=0;ll<Npaths;ll++){

		int sw_spell=0,uspell=0;
		double wlast =0., wnext=0.;
		for(i=0;i<Nsim;i++){
			for(ti=0;ti<TT-1;ti++){
				if( ti>3 && ti< TT-2  ){
					wlast = gg_get( ht->whist[ll] ,i,ti-1) + gg_get( ht->whist[ll] ,i,ti-2)+ gg_get( ht->whist[ll] ,i,ti-3);
					wnext = gg_get( ht->whist[ll] ,i,ti) + gg_get( ht->whist[ll] ,i,ti+1)+ gg_get( ht->whist[ll] ,i,ti+2);
				}else{
					wlast = 0.; wnext=0.;
				}

				double w_EU=0.;
				if( ggi_get(ht->uhist[ll],i,ti) ==0 ){
					if( ggi_get(ht->uhist[ll],i,ti+1) ==1 ){
						w_EU = wnext- wlast; //not yet sure if this will be a switch or not, so just store it for now.
						sw_spell = ggi_get(ht->jhist[ll],i,ti);
					}
					else{
						if( ggi_get(ht->J2Jhist[ll],i,ti+1) ==1 ){

							if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti) ){
	                            w_EEsw[idx_EEsw] = wnext - wlast;
	                            idx_EEsw ++;
							}else{
	                            w_EEns[idx_EEns] = wnext - wlast;
	                            idx_EEns ++;
							}
						}else{  // not EE
							if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti) ){
	                            w_stsw[idx_stsw] = wnext - wlast;
							    idx_stsw ++;
							}else{
	                            w_stns[idx_stns] = wnext - wlast;
	                            idx_stns ++;
							}
						}
					}
				}else{ //unemployed

					if( ggi_get(ht->uhist[ll],i,ti+1) ==0 ){
						if( ggi_get(ht->jhist[ll],i,ti+1) != sw_spell){
							w_UEsw[idx_UEsw] = wnext-wlast;
							w_EUsw[idx_EUsw] = w_EU;
							idx_UEsw ++;
							idx_EUsw ++;
						}else{
							w_UEns[idx_UEns] = wnext-wlast;
							w_EUns[idx_EUns] = w_EU;
							idx_UEns ++;
							idx_EUns ++;
						}
					}
//					if( ggi_get(ht->jhist[ll],i,ti+1) !=ggi_get(ht->jhist[ll],i,ti) && sw_spell == 0 ){
//						sw_spell = 1; //only count the switch once per unemployment spell
//					}
				}
			}
		}
	}

	qsort(w_stns, (size_t)idx_stns,sizeof(double), &comp_dble_asc);
	qsort(w_stsw, (size_t)idx_stsw,sizeof(double), &comp_dble_asc);

	qsort(w_EEns, (size_t)idx_EEns,sizeof(double), &comp_dble_asc);
	qsort(w_EEsw, (size_t)idx_EEsw,sizeof(double), &comp_dble_asc);

	qsort(w_EUns, (size_t)idx_EUns,sizeof(double), &comp_dble_asc);
	qsort(w_EUsw, (size_t)idx_EUsw,sizeof(double), &comp_dble_asc);

	qsort(w_UEns, (size_t)idx_UEns,sizeof(double), &comp_dble_asc);
	qsort(w_UEsw, (size_t)idx_UEsw,sizeof(double), &comp_dble_asc);

	w_qtls(w_stns,1,idx_stns,st->stns_qtls);
	w_qtls(w_stsw,1,idx_stsw,st->stsw_qtls);
	w_qtls(w_EEns,1,idx_EEns,st->EEns_qtls);
	w_qtls(w_EEsw,1,idx_EEsw,st->EEsw_qtls);
	w_qtls(w_EUns,1,idx_EUns,st->EUns_qtls);
	w_qtls(w_EUsw,1,idx_EUsw,st->EUsw_qtls);
	w_qtls(w_UEns,1,idx_UEns,st->UEns_qtls);
	w_qtls(w_UEsw,1,idx_UEsw,st->UEsw_qtls);




	free(w_stns);free(w_stsw);free(w_EEns);free(w_EEsw);free(w_EUns);free(w_EUsw);free(w_UEns);free(w_UEsw);
}

double param_dist(  struct cal_params *par , int Npar, double * err_vec , int Nerr){

	// Takes in a parameters vector and pumps out an error vector

	int i;

	struct valfuns vf;
	struct polfuns pf;
	struct shocks sk;
	struct hists ht;

	struct stats st;

	st.EEns_qtls = malloc(Nqtls*sizeof(double));st.EEsw_qtls = malloc(Nqtls*sizeof(double));
	st.EUns_qtls = malloc(Nqtls*sizeof(double));st.EUsw_qtls = malloc(Nqtls*sizeof(double));
	st.UEns_qtls = malloc(Nqtls*sizeof(double));st.UEsw_qtls = malloc(Nqtls*sizeof(double));
	st.stns_qtls = malloc(Nqtls*sizeof(double));st.stsw_qtls = malloc(Nqtls*sizeof(double));

	allocate_mats(&vf,&pf,&ht,&sk);
	allocate_pars(&par);

	init_pf(&pf,&par);
	init_vf(&vf,&par);


	double quad_dist =0;
	for(i=0;i<Nerr;i++)
		quad_dist += err_vec[i]*err_vec[i];


	return(quad_dist);

}


int draw_shocks(struct shocks * sk){
	int ji,i,ti,ll;

	int seed;

#pragma parallel for private(ll, ti,i,ji)
	for(ll=0;ll<Npaths;ll++){
		gsl_rng * rng0 = gsl_rng_alloc( gsl_rng_default );
		seed = 12281951 + 10*omp_get_thread_num();
		gsl_rng_set( rng0, seed );

		for(ti=0;ti<TTT;ti++){
			gsl_vector_set(sk->Asel[ll],ti, gsl_rng_uniform(rng0));
			for(i=0;i<Nsim;i++){
				gg_set( sk->zsel[ll]      ,i,ti, gsl_rng_uniform(rng0) );
				gg_set( sk->thsel[ll]     ,i,ti, gsl_rng_uniform(rng0) );
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
	par->tlev    = gsl_vector_calloc(NT) ;

	par->AloadP = gsl_vector_calloc(JJ);
	par->Ptrans = malloc(sizeof(gsl_matrix*)*JJ);
	for(ji=0;ji<JJ;ji++){
		par-> Ptrans[ji] = gsl_matrix_calloc(NP,NP);
	}
	par->Atrans = gsl_matrix_calloc(NA,NA);
	par->xGtrans = gsl_matrix_calloc(NG,NG);
	par->xStrans = gsl_matrix_calloc(NS,NS);
	par->ztrans = gsl_matrix_calloc(NZ,NZ);
	par->tprob = gsl_vector_calloc(NT);
	par->zprob = gsl_vector_calloc(NZ);
	par->jprob = gsl_vector_calloc(JJ);

}


void allocate_mats( struct valfuns * vf, struct polfuns * pf, struct hists * ht, struct shocks * sk ){

	int j;

	NN = NA*NP*NG*NS*NZ*NT;
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
	sk->thsel      = malloc(sizeof(gsl_matrix*)*Npaths);
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
		sk->thsel[j]      = gsl_matrix_calloc(Nsim,TTT);
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
		gsl_matrix_free(sk->thsel[j]);
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

	free(sk->thsel);
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
	gsl_vector_free(par->tlev);
	gsl_vector_free(par->AloadP);
	for(ji=0;ji<JJ;ji++){
		gsl_matrix_free(par-> Ptrans[ji]);
	}
	free(par->Ptrans);
	gsl_matrix_free(par->Atrans);
	gsl_matrix_free(par->xGtrans);
	gsl_matrix_free(par->xStrans);
	gsl_matrix_free(par->ztrans);
	gsl_vector_free(par->tprob);
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
