#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <gsl/gsl_vector.h>
#include "CDFt.h"

/* 
 * ObsRp, DataGp, DataGf, npas=1000, dev=2 
 void correct_cdft(
 	gsl_vector *obs_rp, gsl_vector *data_gp, gsl_vector *data_gf, int npas, double dev,
 	double *x, double *cdf_rp, double *cdf_gp, double *cdf_gf, double *cdf_rf, double * cdft
 );
*/
SEXP cdft(SEXP obs_rp, SEXP data_gp, SEXP data_gf, SEXP npas, SEXP dev){
	
	int length_obs_rp = length(obs_rp);
	int length_data_gp = length(data_gp);
	int length_data_gf = length(data_gf);
	int int_npas = INTEGER(npas)[0];
	double x[int_npas], cdf_rp[int_npas], cdf_gp[int_npas], cdf_gf[int_npas], cdf_rf[int_npas];
	double cdft[length_data_gf];

	gsl_vector 
		*gsl_obs_rp = gsl_vector_alloc(length_obs_rp),
		*gsl_data_gp = gsl_vector_alloc(length_data_gp),
		*gsl_data_gf = gsl_vector_alloc(length_data_gf);
  	
	{	
		double *r_obs_rp = REAL(obs_rp), *r_data_gp = REAL(data_gp), *r_data_gf = REAL(data_gf);
		for ( size_t index = 0; index < length_obs_rp; ++index ) {
			gsl_vector_set( gsl_obs_rp, index, r_obs_rp[index] );
		}
		for ( size_t index = 0; index < length_data_gp; ++index ) {
			gsl_vector_set( gsl_data_gp, index, r_data_gp[index] );
		}
		for ( size_t index = 0; index < length_data_gf; ++index ) {
			gsl_vector_set( gsl_data_gf, index, r_data_gf[index] );
		}
	}
 	correct_cdft(
 		gsl_obs_rp, gsl_data_gp, gsl_data_gf, int_npas, REAL(dev)[0],
 		x, cdf_rp, cdf_gp, cdf_gf, cdf_rf, cdft
 	);

	SEXP R_x = PROTECT(allocVector(REALSXP, int_npas));
	SEXP R_cdf_rp = PROTECT(allocVector(REALSXP, int_npas));
	SEXP R_cdf_gp = PROTECT(allocVector(REALSXP, int_npas));
	SEXP R_cdf_gf = PROTECT(allocVector(REALSXP, int_npas));
	SEXP R_cdf_rf = PROTECT(allocVector(REALSXP, int_npas));
	SEXP R_cdft = PROTECT(allocVector(REALSXP, length_data_gf));

	for(int i = 0; i < int_npas; i++){
		REAL(R_x)[i] = x[i];
		REAL(R_cdf_rp)[i] = cdf_rp[i];
		REAL(R_cdf_gp)[i] = cdf_gp[i];
		REAL(R_cdf_gf)[i] = cdf_gf[i];
		REAL(R_cdf_rf)[i] = cdf_rf[i];
	}
	for(int i = 0; i < length_data_gf; i++){
		REAL(R_cdft)[i] = cdft[i];
	}
  
	SEXP output = PROTECT(allocVector(VECSXP, 6));
	SET_VECTOR_ELT(output, 0, R_x);
	SET_VECTOR_ELT(output, 1, R_cdf_rp);
	SET_VECTOR_ELT(output, 2, R_cdf_gp);
	SET_VECTOR_ELT(output, 3, R_cdf_gf);
	SET_VECTOR_ELT(output, 4, R_cdf_rf);
	SET_VECTOR_ELT(output, 5, R_cdft);
  
	UNPROTECT(7);
	gsl_vector_free(gsl_obs_rp);
	gsl_vector_free(gsl_data_gp);
	gsl_vector_free(gsl_data_gf);

	return output;
}


