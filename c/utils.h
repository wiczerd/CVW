#ifndef _utils_h
#define _utils_h

// NOTE, to use MKL LAPACK define _MKL_USE
#ifdef _MKL_USE
#include "mkl_lapacke.h"
#endif
// is this syncing?
// is it syncing yet?
// now?

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>

void printarray(char* name, double* vec, int rows, int append);

int comp_dble_asc( const void * a, const void * b ){
	double x = *(double *)a;
	double y = *(double *)a;
	if( x<y ) return -1; //flip the sign for descending order
	if( x>y ) return  1;
	return 0;
}
int comp_dble_desc( const void * a, const void * b ){
	double x = *(double *)a;
	double y = *(double *)a;
	if( x<y ) return  1; //flip the sign for descending order
	if( x>y ) return -1;
	return 0;
}

int fact_int(int n){
    int i;
    int fact = n;
    for(i=1;i<n;i++) fact *= i;
    return( fact);
}

int set_matrix_block(gsl_matrix * dest,gsl_matrix *source,int left,int top){
	int i,j;
	int right = (int)source->size1 + left;
	int bottom = (int)source->size2 + top;
	for(i=left;i<right;i++){
		for(j=top;j<bottom;j++){
			double sourceij = gsl_matrix_get(source,i-left,j-top);
			gsl_matrix_set(dest,i,j,sourceij);
		}
	}
	return 0;
}
int outer(gsl_vector * a,gsl_vector* b, gsl_matrix* C){
	// computes a*b' = C
	int Na,Nb,r,c;
	Na = (int)a->size;
	Nb = (int)b->size;
	if (C->size1!=Na || C->size2!=Nb){
		printf("Matrix C is the wrong size, Na,Nb=(%d,%d), C=(%d,%d)\n",Na,Nb,(int)C->size1,(int)C->size2);
		return 1;
	}
	for(r=0;r<Na;r++){
		for(c=0;c<Nb;c++){
			double prod = gsl_vector_get(a,r)*gsl_vector_get(b,c);
			gsl_matrix_set(C,r,c,prod);
		}
	}
	return 0;
}

int find_eigens(gsl_matrix * subject, gsl_vector * eigenvals, gsl_matrix * eigenvecs){
	// note that subject is destroyed in here!
	int status=0;
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(subject->size1);
	status += gsl_eigen_symmv(subject, eigenvals, eigenvecs, w);
	gsl_eigen_symmv_free(w);
	//gsl_matrix_free(subject);
	return status;
}

double norm1mat(gsl_matrix *m){
	// sum of colmumn norms
	int nc = m->size2;
	int nr = m->size1;
	int c,r;
	double csum,mn = 0.0;
	for (c=0;c<nc;c++){
		csum = 0.0;
		for(r=0;r<nr;r++)
			csum+= fabs(gsl_matrix_get(m,r,c));
		mn = csum>mn ? csum : mn;
	}
	return mn;
}

int ergod_dist( gsl_matrix * Pi, gsl_vector * dist){
    int N,i;

    N = (int)Pi->size1;
    if(Pi->size2 != N){
        printf(" Transition matrix is not square ");
    }
    gsl_matrix * PiPi = gsl_matrix_calloc(N,N);
    gsl_matrix * PiPiPi = gsl_matrix_calloc(N,N);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,Pi,Pi,0.,PiPi);
    for(i=0;i<N*N*N*N*N*N*N*N;i++) {
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,Pi,PiPi,0.,PiPiPi);
        gsl_matrix_memcpy(PiPi,PiPiPi);
    }

    for(i=0;i<N;i++) gsl_vector_set(dist,i,gsl_matrix_get(PiPi,0,i));

    gsl_matrix_free(PiPi);gsl_matrix_free(PiPiPi);

    return 0;

}

int kron(int TransA,int TransB,gsl_matrix* A,gsl_matrix* B,int TransC,gsl_matrix* C){
//computes kronecker products of A \oprod B = C
	int rA,rB,cA,cB, ri1,ri2,ci1,ci2;
	rA = TransA==0 ? A->size1 : A->size2;
	cA = TransA==0 ? A->size2 : A->size1;
	rB = TransB==0 ? B->size1 : B->size2;
	cB = TransB==0 ? B->size2 : B->size1;

	if(rA*rB == C->size1 && cA*cB == C->size2){

	if (TransC == 0){
		for (ri1 = 0;ri1<rA;ri1++){
			for(ci1=0;ci1<cA;ci1++){

		for (ri2 = 0;ri2<rA;ri2++){
			for(ci2=0;ci2<cA;ci2++){
				gsl_matrix_set(C,rB*ri1 + ri2,cB*ci1 + ci2, gsl_matrix_get(A,ri1,ci1)*gsl_matrix_get(B,ri2,ci2) );
			}
		}

			}
		}
	}
	else{
		for (ri1 = 0;ri1<rA;ri1++){
			for(ci1=0;ci1<cA;ci1++){

		for (ri2 = 0;ri2<rA;ri2++){
			for(ci2=0;ci2<cA;ci2++){
				gsl_matrix_set(C,cB*ci1 + ci2, rB*ri1 + ri2, gsl_matrix_get(A,ri1,ci1)*gsl_matrix_get(B,ri2,ci2) );
			}
		}

			}
		}	
	}
	return 0;
	}
	else{
	return 1;
	}
}



double cov(gsl_vector* x1, gsl_vector * x2){
	double m1,m2,len,cov;
	int ri;
	if(x1->size != x2->size){
		printf("Cov vectors are the wrong size");
		return 0.0;
	}
	len = (double)x1->size;
	m1=0.0; m2=0.0;cov=0.0;
	for(ri =0;ri<x1->size;ri++){
		if(gsl_finite(x1->data[ri])==1 ) m1 += x1->data[ri];
		else len -= 1.;
	}
	m1 /= len;

	len = (double)x2->size;
	for(ri =0;ri<x2->size;ri++){
		if(gsl_finite(x2->data[ri])==1 ) m2 += x2->data[ri];
		else len -= 1.;
	}
	m2 /= len;

	len = (double)x1->size;
	for(ri =0;ri<x2->size;ri++)
		if(gsl_finite(x1->data[ri])==1 && gsl_finite(x2->data[ri])==1)
			cov +=(x2->data[ri] -m2)*(x1->data[ri] -m1);
		else
			len -= 1.;
	cov /= len;
	return cov;
}




int WOLS(const gsl_vector* y, const gsl_matrix* x,const gsl_matrix* w, gsl_vector * coef, gsl_vector * e){
	int status =0;
	int r = x->size1;
	int c = x->size2;
	gsl_set_error_handler_off();
	if (r<2){
		gsl_vector_set(coef,0,0.0);
		gsl_vector_set(coef,1,0.0);
		return 0;
	}

	gsl_matrix * xx 	= gsl_matrix_alloc(c,c);
	gsl_vector * xpwy	= gsl_vector_alloc(c);
	gsl_vector * p 	= gsl_vector_alloc (c);
	gsl_matrix * xpw		= gsl_matrix_alloc(c,r);

	gsl_blas_dgemm( CblasTrans, CblasNoTrans,1.0,x,w,0.0,xpw);
	gsl_blas_dgemm( CblasNoTrans, CblasNoTrans,1.0,xpw,x,0.0,xx);
	status = gsl_linalg_QR_decomp(xx, p);
	status+= gsl_blas_dgemv (CblasNoTrans, 1.0, xpw, y, 0.0, xpwy);
	status+= gsl_linalg_QR_solve(xx, p, xpwy,coef);
	gsl_blas_dgemv(CblasNoTrans,1.0,x,coef,0.0,e);

	gsl_vector_free(p);
	gsl_matrix_free(xx);
	gsl_vector_free(xpwy);
	gsl_matrix_free(xpw);

	return status;
}


int WOLS_LU(const gsl_vector* y, const gsl_matrix* x,const gsl_matrix* w, gsl_vector * coef, gsl_vector * e){
	int status =0,s;
	int r = x->size1;
	int c = x->size2;
	gsl_set_error_handler_off();
	if (r<2){
		gsl_vector_set(coef,0,0.0);
		gsl_vector_set(coef,1,0.0);
		return 0;
	}

	gsl_matrix * xx 	= gsl_matrix_alloc(c,c);
	gsl_vector * xpwy	= gsl_vector_alloc(c);
	gsl_permutation * p 	= gsl_permutation_alloc (c);
	gsl_matrix * xpw		= gsl_matrix_alloc(c,r);
	
	gsl_blas_dgemm( CblasTrans, CblasNoTrans,1.0,x,w,0.0,xpw);
	gsl_blas_dgemm( CblasNoTrans, CblasNoTrans,1.0,xpw,x,0.0,xx);
	status = gsl_linalg_LU_decomp(xx, p, &s);
	status+= gsl_blas_dgemv (CblasNoTrans, 1.0, xpw, y, 0.0, xpwy);
	status+= gsl_linalg_LU_solve(xx, p, xpwy,coef);
	gsl_blas_dgemv(CblasNoTrans,1.0,x,coef,0.0,e);	
	gsl_permutation_free(p);

	gsl_matrix_free(xx);
	gsl_vector_free(xpwy);
	gsl_matrix_free(xpw);

	return status;
}

int OLS(const gsl_vector* y, const gsl_matrix* x, gsl_vector * coef, gsl_vector * e){
	int status =0,s;
//	gsl_matrix * W = gsl_matrix_alloc(x->size1,x->size1);
//	gsl_matrix_set_identity(W);
//	status = WOLS_LU(y,x,W,coef,e);

//	gsl_matrix_free(W);
	int r = x->size1;
	int c = x->size2;

	if (r<2){
		gsl_vector_set(coef,0,0.0);
		gsl_vector_set(coef,1,0.0);
		return 0;
	}

	gsl_vector * xpy	= gsl_vector_alloc(c);
	gsl_matrix * xx 	= gsl_matrix_alloc(c,c);
	gsl_permutation * p 	= gsl_permutation_alloc (c);

	gsl_blas_dgemm( CblasTrans, CblasNoTrans,1.0,x,x,0.0,xx);
	status = gsl_linalg_LU_decomp(xx, p, &s);
	status+= gsl_blas_dgemv (CblasTrans, 1.0, x, y, 0.0, xpy);
	status+= gsl_linalg_LU_solve(xx, p, xpy,coef);
	gsl_blas_dgemv(CblasNoTrans,1.0,x,coef,0.0,e);
	gsl_permutation_free(p);

	gsl_matrix_free(xx);
	gsl_vector_free(xpy);

	return status;
}

int inv_wrap(gsl_matrix * invx, gsl_matrix * x){
	int r = x->size1;
	int c = x->size2;
	int status =0,s;

#ifdef _MKL_USE
	int* ipiv = (int*)malloc(sizeof(int)*r);
	for(s=0;s<r;s++) ipiv[s] =0; //initialize
	double * xdat;
	if(invx->tda!=c){// will have to recopy into a row-major array
		xdat = (double*)malloc(r*c*sizeof(double));
		int ri,ci;
		for(ri=0;ri<r;ri++){
			for(ci=0;ci<c;ci++)
				xdat[ri*c+ci]=x->data[ri*x->tda + ci];
		}
	}
	else{
		gsl_matrix_memcpy(invx,x);
		xdat = invx->data;
	}
	status += LAPACKE_dgetrf( LAPACK_ROW_MAJOR, r, c, xdat, c, ipiv );
	status += LAPACKE_dgetri( LAPACK_ROW_MAJOR, c, xdat, c, ipiv );
	if(invx->tda!=c){
		gsl_matrix_view invvw = gsl_matrix_view_array(xdat,r,c);
		gsl_matrix_memcpy(invx,&invvw.matrix);
		free(xdat);
	}
	free(ipiv);
#endif


#ifndef _MKL_USE
	gsl_permutation * p 	= gsl_permutation_alloc (c);
	status = gsl_linalg_LU_decomp(x, p, &s);
	status = gsl_linalg_LU_invert(x, p, invx);
	gsl_permutation_free(p);
#endif

	return status;
}

int pinv_sol(gsl_vector * x, const gsl_matrix * A, const gsl_vector * b){
	int status=0;
	int Nr,Nc,ri,ci;
	gsl_matrix * U, * V;
	gsl_vector * S, *work;
	
	Nr = A->size1;
	Nc = A->size2;
	
	U 	= gsl_matrix_alloc(Nr,Nc);
	V 	= gsl_matrix_alloc(Nc,Nc);
	S 	= gsl_vector_alloc(Nc);
	work 	= gsl_vector_alloc(Nc);

	gsl_matrix_memcpy(U,A);
	
	if(Nr>Nc){
		gsl_matrix * Xwork  = gsl_matrix_alloc(Nc,Nc);
		status += gsl_linalg_SV_decomp_mod(U, Xwork, V, S, work);
		gsl_matrix_free(Xwork);
	}
	else
		status += gsl_linalg_SV_decomp(U, V, S, work);	
	status += gsl_linalg_SV_solve(U,V,S,b,x);
	
	gsl_vector_free(S);
	gsl_matrix_free(V);
	gsl_matrix_free(U);
	gsl_vector_free(work);
	
	return status;
}

int sol_Axb_tri(gsl_matrix * A, gsl_vector * b, gsl_vector * x, gsl_permutation * p){
	// if tri==1 then assumes already LU decomp
	int s,status=0;
	int r = A->size1;
	int c = b->size;
	gsl_matrix* AA;

	if(r!=c)
		status = 1;
	int tri = 1;
	if(p==0){
		tri =0;
		p = gsl_permutation_alloc (c);
		AA = gsl_matrix_alloc(r,c);
		gsl_matrix_memcpy(AA,A);
		status = gsl_linalg_LU_decomp (AA, p, &s);
		status += gsl_linalg_LU_solve (AA, p, b, x);
	}
	else
		status = gsl_linalg_LU_solve (A, p, b, x);

	if(tri!=1){
		gsl_permutation_free (p);
		gsl_matrix_free(AA);
	}
	return status;
}
int sol_Axb(gsl_matrix * A, gsl_vector * b, gsl_vector * x){
	// this is just a wrapper that assumes that not already triangularized A
	int status = sol_Axb_tri(A, b, x, 0);
	return status;
}

int sol_AXM(gsl_matrix * A, gsl_matrix * M, gsl_matrix *X, void* tri_v){
	// if int tri==0 (null pointer) then assume not yet LU decomposed
	int status, r, c, rm,cm;
	r = A->size1;
	c = A->size2;
	rm= M->size1;
	cm= M->size2;



#ifndef _MKL_USE

	gsl_matrix * AA;
	gsl_permutation * p = (gsl_permutation*)tri_v;
	int s;
	AA = gsl_matrix_calloc(r,c);
	gsl_matrix_memcpy(AA,A);
	status = gsl_linalg_LU_decomp (AA, p, &s);

	int ci;
	for(ci=0;ci<cm;ci++){
		gsl_vector_view x = gsl_matrix_column(X,ci);
		gsl_vector_view b = gsl_matrix_column(M,ci);
		status = gsl_linalg_LU_solve (A, p, &b.vector, &x.vector);
	}
	gsl_matrix_free(AA);

#endif

#ifdef _MKL_USE
	int* tri = (int*)tri_v;
	int * ipiv;
	double *xdat,*ydat;
	if(A->tda!=c){// will have to recopy into a row-major array
		xdat = (double*)malloc(r*c*sizeof(double));
		int ri,ci;
		for(ri=0;ri<r;ri++){
			for(ci=0;ci<c;ci++)
				xdat[ri*c+ci]=gsl_matrix_get(A,ri,ci);
		}
	}
	else
		xdat = A->data;
	if(M->tda!=c){// will have to recopy into a row-major array
		ydat = (double*)malloc(r*c*sizeof(double));
		int ri,ci;
		for(ri=0;ri<rm;ri++){
			for(ci=0;ci<cm;ci++)
				ydat[ri*cm+ci]=gsl_matrix_get(M,ri,ci);
		}
	}
	else
		ydat = M->data;
	if(tri== 0){
		ipiv = (int*)malloc(sizeof(int)*c);
		status =LAPACKE_dgetrf(LAPACK_ROW_MAJOR, r,c, xdat, c,ipiv );
	}
	else
		ipiv = tri;

	status += LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', rm, c, xdat, c, ipiv, ydat, cm);
	if(A->tda!=c)
		free(xdat);
	if(M->tda!=c)
		free(ydat);

#endif


	return status;
}

int pca(const gsl_matrix * data, const int dimensions_we_want, const int demean, const int stdize,gsl_matrix * pc_space, gsl_matrix * projected){
	int ds = data->size2;
	int i,c,status=0;
	gsl_matrix * x;
	x = gsl_matrix_calloc(data->size1,data->size2);

	if (demean>0){
		gsl_matrix_memcpy(x,data);
		for(i=0;i<data->size1;i++){
			double rmean = 0.0;
			for (c=0;c<data->size2;c++)
				rmean += data->data[i*data->tda + c];
			rmean /= (double)data->size2;
			for (c=0;c<data->size2;c++)
				x->data[i*data->tda + c] -= rmean;
		}
	}
	else
		gsl_matrix_memcpy( x,data);
	if(stdize >0){
		for(i=0;i<data->size1;i++){
			double rvar =  0.;
			double rmean = 0.0;
			for (c=0;c<data->size2;c++)
				rmean += data->data[i*data->tda + c];
			rmean /= (double)data->size2;
			for (c=0;c<data->size2;c++)
				rvar += pow(data->data[i*data->tda + c] -rmean ,2);
			rvar /= (double)data->size2;
			for (c=0;c<data->size2;c++)
				x->data[i*data->tda + c] /= pow(rvar,0.5);
		}
	}

	gsl_matrix *cov_matrix = gsl_matrix_calloc(ds, ds);
	gsl_vector *eigenvals = gsl_vector_alloc(ds);
	gsl_matrix *eigenvecs = gsl_matrix_alloc(ds, ds);
	status += gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0, x, x, 0.0, cov_matrix);
	status += find_eigens(cov_matrix, eigenvals, eigenvecs);

	status += gsl_eigen_symmv_sort(eigenvals, eigenvecs,GSL_EIGEN_SORT_ABS_DESC);


	for (i=0;i<dimensions_we_want; i++){
		gsl_vector_view evec = gsl_matrix_column(eigenvecs, i);
		status += gsl_matrix_set_col(pc_space, i, &evec.vector);
	}
	status += gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0, x, pc_space, 0.0, projected);

	gsl_matrix_free(cov_matrix);
	gsl_vector_free(eigenvals);
	gsl_matrix_free(eigenvecs);
	gsl_matrix_free(x);
	return status;
}



int isfinmat(gsl_matrix * mat){
	// returns 0 if entire matrix is finite and o.w. the number of non-finite entries
	int r,c,Nr,Nc,status;
	status = 0;
	Nr = mat->size1;
	Nc = mat->size2;

	#pragma omp parallel for default(shared) private(r,c)
	for(r=0;r<Nr;r++){
		for(c=0;c<Nc;c++){
			if(gsl_finite(mat->data[r*mat->tda +c])==0){
				mat->data[r*mat->tda +c] = 0.0;
				status++;
			}
		}
	}
	return status;
}

int isfinvec(gsl_vector * vec){
	int r,Nr,status;
	status =0;
	Nr = (int)vec->size;
	//#pragma omp parallel for default(shared) private(r)
	for(r=0;r<Nr;r++){
		if(gsl_finite(vec->data[r])==0){
			vec->data[r] = 0.0;
			status++;
		}
	}
	return status;
}


int rouwenhorst(double ar1, double sig, gsl_matrix * pi_out, gsl_vector * grid_inout){
// ar1 is the auto correlation and
// sig is the conditinal standard deviation
	gsl_matrix * PP = gsl_matrix_calloc(2,2);
	gsl_matrix *PPm1= gsl_matrix_calloc(2,2);
	gsl_matrix * qPP= gsl_matrix_calloc(2,2);
	int i,ii, ri,ci;
	int N = (int)grid_inout->size;

	double q = (ar1+1.)/2.;double p =(ar1+1.)/2;
	double nu = pow((double)(N-1)/(1.-ar1*ar1),0.5)*sig;

	gsl_matrix_set(PP,0,0, p);gsl_matrix_set(PP,0,1, 1.-p);
	gsl_matrix_set(PP,1,1, q);gsl_matrix_set(PP,1,0, 1.-q);

	gsl_matrix_memcpy(  PPm1,PP);

	for(i=2;i<N;i++){
		gsl_matrix_free(PP);
		gsl_matrix_free(qPP);

		PP = gsl_matrix_calloc(i+1,i+1);
		qPP = gsl_matrix_calloc(i+1,i+1);

		set_matrix_block(PP,PPm1,0,0);
		gsl_matrix_scale(PP,p);
		set_matrix_block(qPP,PPm1,1,0);
		gsl_matrix_scale(qPP,1.-q);
		gsl_matrix_add(PP,qPP );
		gsl_matrix_set_zero(qPP);
		set_matrix_block(qPP,PPm1,0,1);
		gsl_matrix_scale(qPP,1.-p);
		gsl_matrix_add(PP,qPP);
		gsl_matrix_set_zero(qPP);
		set_matrix_block(qPP,PPm1,1,1);
		gsl_matrix_scale(qPP,q);
		gsl_matrix_add(PP,qPP);

		for(ri=1;ri<i;ri++){
			for(ci=0;ci<i;ci++) PP->data[ri*PP->tda+ci] *= 0.5;
		}

		gsl_matrix_free(PPm1);
		PPm1 = gsl_matrix_calloc(i+1,i+1);
		gsl_matrix_memcpy( PPm1,PP);
	}

	for(i=0;i<N;i++){
		gsl_vector_set(grid_inout,i, 2.*nu*((double)i)/((double)N-1.) -nu);
	}
	// make sure rowsums are 1
	for(ri=0;ri<N;ri++){
		double rsum=0;
		for(ci=0;ci<N;ci++) rsum+= PP->data[ri*PP->tda+ci];
		gsl_vector PPr = (gsl_matrix_row(PP,ri).vector);
		gsl_vector_scale( &PPr , 1./rsum);
	}

	gsl_matrix_memcpy(pi_out,PP);
	gsl_matrix_free(PP);gsl_matrix_free(PPm1);gsl_matrix_free(qPP);

}
int disc_Weibull(double * zprob, double * zlev, const int NZ, const double mean_z, const double scale_z, const double shape_z){
	int success = 0;
	int ii;
	double zub[NZ-1];
	for(ii=0;ii<(NZ-1);ii++)
		zub[ii] = gsl_cdf_weibull_Pinv( ((double)ii+1.)/(double)NZ , scale_z,shape_z);
	for(ii=1;ii<NZ;ii++)
		zlev[ii] =(zub[ii]+zub[ii-1])/2.;
	zlev[0]  = gsl_cdf_weibull_Pinv( 1./(double)NZ /2., scale_z,shape_z) ;
	zlev[NZ-1] = gsl_cdf_weibull_Pinv( ((double)NZ-1.)/(double)NZ  + 1./(double)NZ /2.,  scale_z,shape_z );
	for(ii=0;ii<NZ;ii++)  zprob[ii] = 1./(double)NZ;
	double meanz_tmp = scale_z*gsl_sf_gamma(1.+1./shape_z);
	for(ii=0;ii<NZ;ii++) zlev[ii] += mean_z - meanz_tmp;

	return success;
}

int disc_2emg(double * epsprob, double * epslev, const int NE, const double mean_eps, const double var_eps, const double lshape_eps, double rshape_eps ){
/*
Requires and odd number of points to have proper placement.

epsprob    : Output. Array of length NE with probabilities of cells
epslev     : Output. Center point of cells (unweighted by the probability). The exception are the top and bottom cells, which are essentially extrapolants from the rest given that long tails may have very far center points
NE         : Input. Number of points
mean_eps   : Input. Mean of normal
var_eps    : Input. Variance of normal
lshape_eps : Input. Left-side exponential parameter (larger numbers mean shorter tails)
rshape_eps : Input. Right-side exponential parameter (larger numbers mean shorter tails)

 Returns success = 0 if all good
*/

	int ii,ri;
	int success=1;
	int const Nevalpts=500;
	double epslev0[NE];
	double sd_eps = pow(var_eps,0.5);

	// heuristic for placement of endpoints: where there's very little mass even on the exponential:
	double epsUB = gsl_cdf_gaussian_Pinv( 1./3., sd_eps  ) +
	               gsl_cdf_exponential_Pinv( 0.995, 1./rshape_eps ) +mean_eps;
	double epsLB = -gsl_cdf_gaussian_Pinv( 1./3., sd_eps  ) -
	               gsl_cdf_exponential_Pinv( 0.995, 1./lshape_eps ) +mean_eps;

	double leps_frac = lshape_eps/(lshape_eps+rshape_eps);
	double reps_frac = rshape_eps/(lshape_eps+rshape_eps);

	double evalgrid[Nevalpts]; double epsPrGrid[Nevalpts];
	for(ii=0;ii<Nevalpts;ii++) evalgrid[ii] = (epsUB-epsLB)*(double)ii/(double)Nevalpts + epsLB;

	double sumprobs =0.;
	for(ii=1;ii < Nevalpts-1;ii++){
		double epsL = evalgrid[ii-1]/2 + evalgrid[ii]/2;
		double epsH = evalgrid[ii+1]/2 + evalgrid[ii]/2;
		double RcdfH = gsl_cdf_ugaussian_P( (epsH-mean_eps)/sd_eps  - rshape_eps*sd_eps );
		double LcdfH = gsl_cdf_ugaussian_P(-(epsH-mean_eps)/sd_eps  - lshape_eps*sd_eps );
		double RcdfL = gsl_cdf_ugaussian_P( (epsL-mean_eps)/sd_eps  - rshape_eps*sd_eps );
		double LcdfL = gsl_cdf_ugaussian_P(-(epsL-mean_eps)/sd_eps  - lshape_eps*sd_eps );
		double Pr_hr = gsl_cdf_ugaussian_P( (epsH-mean_eps)/sd_eps ) -
		               leps_frac*exp(- rshape_eps*(epsH-mean_eps) + var_eps/2 * rshape_eps * rshape_eps )*RcdfH+
		               reps_frac*exp(  lshape_eps*(epsH-mean_eps) + var_eps/2 * lshape_eps * lshape_eps )*LcdfH
		               -(	gsl_cdf_ugaussian_P(epsL/sd_eps ) -
		                     leps_frac*exp(- rshape_eps*(epsL-mean_eps) + var_eps/2 * rshape_eps * rshape_eps )*RcdfL+
		                     reps_frac*exp(  lshape_eps*(epsL-mean_eps) + var_eps/2 * lshape_eps * lshape_eps )*LcdfL
		               );
		epsPrGrid[ii]= Pr_hr;
		sumprobs += Pr_hr;
	}
	ii=0;
	double epsH = evalgrid[ii+1]/2 + evalgrid[ii]/2;
	double RcdfH = gsl_cdf_ugaussian_P( epsH/sd_eps  - rshape_eps*sd_eps );
	double LcdfH = gsl_cdf_ugaussian_P(-epsH/sd_eps  - lshape_eps*sd_eps );
	double Pr_hr = gsl_cdf_ugaussian_P(epsH/sd_eps ) -
	               leps_frac*exp(-rshape_eps*epsH + var_eps/2 * rshape_eps * rshape_eps )*RcdfH+
	               reps_frac*exp( lshape_eps*epsH + var_eps/2 * lshape_eps * lshape_eps )*LcdfH;
	epsPrGrid[ii]=Pr_hr;
	sumprobs += Pr_hr;
	ii=Nevalpts-1;
	double epsL  = evalgrid[ii-1]/2+ evalgrid[ii]/2;
	double RcdfL = gsl_cdf_ugaussian_P( epsL/sd_eps  - rshape_eps*sd_eps );
	double LcdfL = gsl_cdf_ugaussian_P(-epsL/sd_eps  - lshape_eps*sd_eps );
	Pr_hr = 1.-
	        (	gsl_cdf_ugaussian_P(epsL/sd_eps ) -
	              leps_frac*exp(- rshape_eps*epsL + var_eps/2 * rshape_eps * rshape_eps )*RcdfL+
	              reps_frac*exp(  lshape_eps*epsL + var_eps/2 * lshape_eps * lshape_eps )*LcdfL
	        );
	epsPrGrid[ii] = Pr_hr;
	sumprobs += Pr_hr;

	for(ii=0;ii<NE;ii++) epsprob[ii]  *= 1./sumprobs ;
	if(sumprobs>.999999 && sumprobs<1.000001) success=0;
	// pick points by inverting the distribution:

	double epsub[NE-1];
	double cumepsPrGrid[Nevalpts]; cumepsPrGrid[0] = epsPrGrid[0];
	for( ii=1;ii<Nevalpts;ii++) cumepsPrGrid[ii] = epsPrGrid[ii] + cumepsPrGrid[ii-1];
	// first pick the upper-bounds
	for( ii=0; ii<NE -1; ii++ ){
		int iL =0;
		while( (double)(ii+1)/(double)NE > cumepsPrGrid[iL+1] && iL<Nevalpts-1){
			iL++;
		}
		double wtL = (cumepsPrGrid[iL+1] - (double)(ii+1)/(double)NE)/(cumepsPrGrid[iL+1] - cumepsPrGrid[iL]);

		epsub[ii] = (wtL * evalgrid[iL] + (1.-wtL)*evalgrid[iL+1]);
		epsprob[ii] = 1./(double)NE;
	}
	// take midpoints:
	for(ii=1;ii<NE-1;ii++){
		epslev[ii] = (epsub[ii-1]+epsub[ii])/2;
	}
	// (epslev[0]+epslev[1])/2 = epsub[0]
	epslev[0]  = 2*epsub[0]-epslev[1];
	epsprob[0] = 1./(double)NE;
	// (epslev[NE-1] + epslev[NE-2])/2 = epsub[NE-1]
	epslev[NE-1] = 2*epsub[NE-2]-epslev[NE-2];
	epsprob[NE-1]= 1./(double)NE;

	// check for weirdness (non-monotone, etc and if so, go back to epslev0
	success =0;
	for(ii=0;ii<NE-1;ii++){
		if(epslev[ii]-epslev[ii+1]>0){

			success =100;
			break;
		}
	}

	return success;

}

void printmat(char* name, gsl_matrix* mat, int append){
	FILE* matfile;int yi,bi;
	int rows = (int)mat->size1;
	int cols = (int)mat->size2;
	if(append == 0){
		matfile = fopen(name, "w");
	}else{
		matfile = fopen(name, "a");
	}

	for(yi=0;yi<rows;yi++){
		for(bi=0;bi<cols-1;bi++){
			fprintf(matfile,"%f,",gsl_matrix_get(mat,yi,bi));
		}
		fprintf(matfile,"%f\n",gsl_matrix_get(mat, yi,cols-1));
	}
	if(append==0) printf("Printing matrix to, %s\n",name);
	fclose(matfile);
}


void printmat_int(char* name, gsl_matrix_int* mat, int append){
	FILE* matfile;int yi,bi;
	int rows = (int)mat->size1;
	int cols = (int)mat->size2;
	if(append == 0){
		matfile = fopen(name, "w");
	}else{
		matfile = fopen(name, "a");
	}
	for(yi=0;yi<rows;yi++){
		for(bi=0;bi<cols-1;bi++){
			fprintf(matfile,"%d,",gsl_matrix_int_get(mat,yi,bi));
		}
		fprintf(matfile,"%d\n",gsl_matrix_int_get(mat, yi,cols-1));
	}
	if(append==0) printf("Printing matrix to, %s\n",name);
	fclose(matfile);
}
void printvec(char* name, gsl_vector* vec, int append){
	FILE* matfile;int ri;
	int rows = (int)vec->size;
	if(append == 0){
		matfile = fopen(name, "w");
	}else{
		matfile = fopen(name, "a");
	}
	for(ri=0;ri<rows;ri++){
		fprintf(matfile,"%f\n",gsl_vector_get(vec,ri));
	}
	if(append==0) printf("Printing matrix to, %s\n",name);
	fclose(matfile);
}
void printarray(char* name, double* vec, int rows, int append){
	FILE* matfile;int ri;

	if(append == 0){
		matfile = fopen(name, "w");
	}else{
		matfile = fopen(name, "a");
	}
	for(ri=0;ri<rows;ri++){
		fprintf(matfile,"%f, ",vec[ri]);
	}
	fprintf(matfile,"\n");
	if(append==0) printf("Printing array to, %s\n",name);
	fclose(matfile);
}
void printarray_int(char* name, int* vec, int rows, int append){
    FILE* matfile;int ri;

    if(append == 0){
        matfile = fopen(name, "w");
    }else{
        matfile = fopen(name, "a");
    }
    for(ri=0;ri<rows;ri++){
        fprintf(matfile,"%d, ",vec[ri]);
    }
    fprintf(matfile,"\n");
    if(append==0) printf("Printing integer array to, %s\n",name);
    fclose(matfile);
}
void printvec_int(char* name, gsl_vector_int* vec, int append){
	FILE* matfile; int ri;
	int rows = (int)vec->size;
	if(append == 0){
		matfile = fopen(name, "w");
	}else{
		matfile = fopen(name, "a");
	}
	for(ri=0;ri<rows;ri++){
		fprintf(matfile,"%d\n",gsl_vector_int_get(vec,ri));
	}
	if(append==0) printf("Printing vector to, %s\n",name);
	fclose(matfile);
}

int readmat(char* name, gsl_matrix * mat){
	int status,Nr,Nc,r,c,rstatus,trashstatus;
	double trash;
	
	Nc = (int) mat->size2;
	Nr = (int) mat->size1;
	status =0;
	FILE * f;
	f = fopen(name,"r");
	if(f == NULL){
		status = 1;
		gsl_matrix_set_zero(mat);
		printf("Could not read from %s",name);
	}
	else{
		double dd;
		for(r=0;r<Nr;r++){
			for(c=0;c<Nc-1;c++){
				rstatus = fscanf(f,"%lf,",&dd);
				if(rstatus==0){
					rstatus = fscanf(f,"%lf\n",&dd);
				}
				if(rstatus==0){
					rstatus = fscanf(f,"%lf,\n",&dd);
				}
				if(rstatus==0){
					rstatus = fscanf(f,"%lf \n",&dd);
				}
				if(rstatus==0)
					printf("Error, did not read from %s\n",name);

				gsl_matrix_set(mat,r,c,(double) dd);
				
			}			
			rstatus = fscanf(f,"%lf\n",&dd);
			/* if this doesn't work, then keep reading --- This does not work
			if(rstatus==0){
				rstatus = fscanf(f,"%lf,",&dd);
				if(rstatus>0){// it wasn't the end of the line, so now read until the end of the line
					for(c=0;c<100;c++){
						trashstatus = fscanf(f,"%lf,",&trash);
						if(trashstatus ==0){
							trashstatus = fscanf(f,"%lf\n",&trash);
							break; //came to the end
						}
					}
				}
			}*/
			if(rstatus==0)
				printf("Error, did not read from %s\n",name);

			gsl_matrix_set(mat,r,Nc-1,dd);
		}
		

		fclose(f);
	}
	return status;
}

int readvec(char* name, gsl_vector * vec){
	int status,N,i,rstatus;

	N = (int) vec->size;

	status =0;
	FILE * f;
	f = fopen(name,"r");
	status = f == NULL ? 1 : 0;
	float dd;
	for(i=0;i<N;i++){
		rstatus = fscanf(f,"%f\n",&dd);
		if(rstatus==0){
			printf("Error, did not read from %s",name);
		}
		gsl_vector_set(vec,i,(double)dd);
	}
	fclose(f);
	return status;
}


void vec(gsl_matrix * mat, gsl_vector * vec_ret){
// stacks by column
// requires 2x memory, as this is not a copy in place transformation
	int r,c;
	for(c=0;c < mat->size2;c++){
		for(r=0; r< mat->size1; r++)
			gsl_vector_set(vec_ret,r+c*(mat->size1),
				gsl_matrix_get(mat,r,c)); 
	}

}

void vech(gsl_matrix * mat, gsl_vector * vech_ret ){
// stacks by column, only using upper triangle of mat
// requires 2x memory, as this is not a copy in place transformation

	int r,c;
	for(c=0;c < mat->size2;c++){
		for(r=0; r< c+1; r++)
			gsl_vector_set(vech_ret,r+c*(mat->size1),
				gsl_matrix_get(mat,r,c)); 
	}
}

int randn(gsl_matrix * mat,unsigned long int seed){
	// fills mat with i.i.d random normal draws
	int r,c,Nr,Nc,status;
	gsl_rng * rng;
	double draw;
	Nr = mat->size1;
	Nc = mat->size2;
	//gsl_rng_env_setup(); // Reads the environment variables GSL_RNG_TYPE and GSL_RNG_SEED

	rng  = gsl_rng_alloc(gsl_rng_mt19937);
	status = 0;
	gsl_rng_set (rng, seed);
	for(r=0;r<Nr;r++){
		for(c=0;c<Nc;c++){
			draw = gsl_ran_gaussian_ziggurat(rng,1.0);
			gsl_matrix_set(mat,r,c,draw);
		}
	}

	gsl_rng_free(rng);
	return status;
}

int randu(gsl_matrix * mat,unsigned long int seed){
	// fills mat with i.i.d random normal draws
	int r,c,Nr,Nc,status;
	gsl_rng * rng;
	double draw;
	Nr = mat->size1;
	Nc = mat->size2;
	//gsl_rng_env_setup(); // Reads the environment variables GSL_RNG_TYPE and GSL_RNG_SEED

	rng  = gsl_rng_alloc(gsl_rng_default);
	status = 0;
	gsl_rng_set (rng, seed);
	for(r=0;r<Nr;r++){
		for(c=0;c<Nc;c++){
			draw = gsl_rng_uniform(rng);
			gsl_matrix_set(mat,r,c,draw);
		}
	}

	gsl_rng_free(rng);
	return status;
}


double zero_brent ( double a, double b, double t, void * pin,
			  double f ( double x, void * params ) )

/******************************************************************************/
/*
  Purpose:

    ZERO seeks the root of a function F(X) in an interval [A,B].

  Discussion:

    The interval [A,B] must be a change of sign interval for F.
    That is, F(A) and F(B) must be of opposite signs.  Then
    assuming that F is continuous implies the existence of at least
    one value C between A and B for which F(C) = 0.

    The location of the zero is determined to within an accuracy
    of 6 * MACHEPS * abs ( C ) + 2 * T.

    Thanks to Thomas Secretin for pointing out a transcription error in the
    setting of the value of P, 11 February 2013.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2013

  Author:

    Original FORTRAN77 version by Richard Brent.
    C version by John Burkardt.

  Reference:

    Richard Brent,
    Algorithms for Minimization Without Derivatives,
    Dover, 2002,
    ISBN: 0-486-41998-3,
    LC: QA402.5.B74.

  Parameters:

    Input, double A, B, the endpoints of the change of sign interval.

    Input, double MACHEP, an estimate for the relative machine
    precision.

    Input, double T, a positive error tolerance.

    Input, double F ( double x ), a user-supplied function whose zero
    is being sought.

    Output, double ZERO, the estimated value of a zero of
    the function F.
*/
{
	double c;
	double d;
	double e;
	double fa;
	double fb;
	double fc;
	double m;
	double p;
	double q;
	double r;
	double s;
	double sa;
	double sb;
	double tol;
	double machep = 2.220446049250313E-016;
/*
  Make local copies of A and B.
*/
	sa = a;
	sb = b;
	fa = f ( sa ,pin );
	fb = f ( sb ,pin );

	c = sa;
	fc = fa;
	e = sb - sa;
	d = e;

	for ( ; ; )
	{
		if ( fabs ( fc ) < fabs ( fb ) )
		{
			sa = sb;
			sb = c;
			c = sa;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol = 2.0 * machep * fabs ( sb ) + t;
		m = 0.5 * ( c - sb );

		if ( fabs ( m ) <= tol || fb == 0.0 )
		{
			break;
		}

		if ( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
		{
			e = m;
			d = e;
		}
		else
		{
			s = fb / fa;

			if ( sa == c )
			{
				p = 2.0 * m * s;
				q = 1.0 - s;
			}
			else
			{
				q = fa / fc;
				r = fb / fc;
				p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
				q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
			}

			if ( 0.0 < p )
			{
				q = - q;
			}
			else
			{
				p = - p;
			}

			s = e;
			e = d;

			if ( 2.0 * p < 3.0 * m * q - fabs ( tol * q ) &&
				 p < fabs ( 0.5 * s * q ) )
			{
				d = p / q;
			}
			else
			{
				e = m;
				d = e;
			}
		}
		sa = sb;
		fa = fb;

		if ( tol < fabs ( d ) )
		{
			sb = sb + d;
		}
		else if ( 0.0 < m )
		{
			sb = sb + tol;
		}
		else
		{
			sb = sb - tol;
		}

		fb = f ( sb,pin );

		if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
		{
			c = sa;
			fc = fa;
			e = sb - sa;
			d = e;
		}
	}
	return sb;
}




#endif
