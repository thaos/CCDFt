#ifndef CDFT_INCLUDED
#define CDFT_UTILS_INCLUDED


#include <gsl/gsl_vector.h>

double gsl_stats_cdf_from_sorted_data( const double sorted_data[], size_t n, double x );

int count_unique( const double sorted_data[], size_t n );

void average_byx( const double sorted_x[], const double y[], size_t n,  double *xu, double *yu);

void correct_cdft(
	gsl_vector *obs_rp, gsl_vector *data_gp, gsl_vector *data_gf, int npas, double dev,
        double *x, double *cdf_rp, double *cdf_gp, double *cdf_gf, double *cdf_rf, double * cdft
);

#endif
