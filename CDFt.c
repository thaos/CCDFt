#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort_vector.h>

double gsl_stats_cdf_from_sorted_data( const double sorted_data[], size_t n, double x ){
	size_t cdf = 0;
	while ( cdf < n &&  sorted_data[cdf] <= x ) {
		cdf++;
	}
	return (double) cdf / (double) n;
}

int count_unique( const double sorted_data[], size_t n ){
	size_t 
		n_unique = n,
		icurr = 0,
	       	inext = icurr + 1;

	while ( inext < n ) {
		printf( "%i:%1i\n", (int) icurr, (int) inext );
		printf( "%3.1f:%3.1f\n", sorted_data[icurr], sorted_data[inext] );
		if( sorted_data[icurr] == sorted_data[inext] ){
			printf( "minus 1\n");
			n_unique--;
		} else {
			printf( "same\n");
			icurr = inext;
		}
		inext++;
	}
	return n_unique;
}

void average_byx( const double sorted_x[], const double y[], size_t n,  double *xu, double *yu){

	size_t 
		icurr = 0,
	       	inext = icurr + 1,
		ifill = 0;
	double 
		num = y[icurr],
		denom = 1;

	while ( inext < n ) {
		if( sorted_x[icurr] == sorted_x[inext] ){
			num += y[inext];
			denom++;
		} else {
			xu[ifill] = sorted_x[icurr];
			yu[ifill] = num / denom;
			ifill++;
			icurr = inext;
			num = y[icurr];
			denom = 1;

		}
		inext++;
	}
	xu[ifill] = sorted_x[icurr];
	yu[ifill] = num / denom;
}

	


void correct_cdft(
	gsl_vector *obs_rp, gsl_vector *data_gp, gsl_vector *data_gf, int npas, double dev,
        double *x, double *cdf_rp, double *cdf_gp, double *cdf_gf, double *cdf_rf, double * cdft
){

/* 
	ObsRp, data_gp, data_gf  gsl_vector
	npas int = 100
	dev float/double = 2
*/

/*	double gsl_stats_mean(const double data[], size_t stride, size_t n) 
	avec vect un gsl_vector
	double gsl_stats_mean(vect->data, vect->stride, vect->size)
	
	
	# create copy of DataGP, data_gf
	int gsl_vector_memcpy(gsl_vector *dest, const gsl_vector *src)
	int gsl_vector_add_constant(gsl_vector *a, const double x) 
*/

	double
		m_obs = gsl_stats_mean( obs_rp->data, obs_rp->stride, obs_rp->size ),
		m_gp = gsl_stats_mean( data_gp->data, data_gp->stride, data_gp->size ),
		m_gf = gsl_stats_mean( data_gf->data, data_gf->stride, data_gf->size );

	gsl_vector 
		*data_gp2 = gsl_vector_alloc( data_gp->size ),
		*data_gf2 = gsl_vector_alloc( data_gf->size ),
		*obs_rp_sorted = gsl_vector_alloc( obs_rp->size ),
		*data_gp2_sorted = gsl_vector_alloc( data_gp->size ),
		*data_gf2_sorted = gsl_vector_alloc( data_gf->size ),
	        *data_pooled = gsl_vector_alloc( data_gp->size + data_gf->size + obs_rp->size );

	gsl_vector_memcpy( data_gp2, data_gp );
	gsl_vector_memcpy( data_gf2, data_gf );
	gsl_vector_add_constant( data_gp2, (m_obs - m_gp) );
	gsl_vector_add_constant( data_gf2, (m_obs - m_gp) );
	gsl_vector_memcpy( obs_rp_sorted, obs_rp );
	gsl_vector_memcpy( data_gp2_sorted, data_gp2 );
	gsl_vector_memcpy( data_gf2_sorted, data_gf2 );
	gsl_sort_vector( obs_rp_sorted );
	gsl_sort_vector( data_gp2_sorted );
	gsl_sort_vector( data_gf2_sorted );

	/* perform the copy to portions of the output */
	{
		gsl_vector_view
			view_rp = gsl_vector_subvector( data_pooled, 0, obs_rp->size ),
			view_gp = gsl_vector_subvector( data_pooled, obs_rp->size, data_gp->size ),
			view_gf = gsl_vector_subvector( data_pooled, obs_rp->size + data_gp->size, data_gf->size );
		gsl_vector_memcpy( &view_rp.vector, obs_rp );
		gsl_vector_memcpy( &view_gp.vector, data_gp );
		gsl_vector_memcpy( &view_gf.vector, data_gf );
	}
	
  	double 
  		a = fabs( m_gf - m_gp),
		lb = gsl_vector_min( data_pooled ) - dev * a,
  		ub = gsl_vector_max( data_pooled ) + dev * a,
		by = (ub - lb) / (npas - 1),
		qntl;

	gsl_vector_free(data_pooled);
 	
  	for ( size_t i = 0; i < npas; i++ ) {
		x[i] = lb + i * by;
		cdf_gp[i] = gsl_stats_cdf_from_sorted_data( data_gp2_sorted->data, data_gp2_sorted->size, x[i] );
		cdf_gf[i] = gsl_stats_cdf_from_sorted_data( data_gf2_sorted->data, data_gf2_sorted->size, x[i] );
		cdf_rp[i] = gsl_stats_cdf_from_sorted_data( obs_rp_sorted->data, obs_rp_sorted->size, x[i] );
		cdf_rf[i] = gsl_stats_cdf_from_sorted_data(
				obs_rp_sorted->data,
			       	obs_rp_sorted->size,
			        gsl_stats_quantile_from_sorted_data( data_gp2_sorted->data, 1,  data_gp2_sorted->size, cdf_gf[i] )
		);
		printf( "%3.1f:%3.1f\n", x[i], cdf_rf[i] );
	}
	printf( "---------\n" );


	/* FRf=FRp with shift for x<min(data_gf) */
	double dgf2_min = gsl_vector_min( data_gf2 );

	if ( gsl_vector_min( obs_rp ) < dgf2_min ) {

		int i = 0;
		double qmin = gsl_stats_quantile_from_sorted_data( obs_rp_sorted->data, 1,  obs_rp->size, cdf_rf[0] );
		while(x[i] <= qmin) {
			i++;
		}
  
		int j = 0;
		while( x[j] < dgf2_min ) {
			j++;
		}

		while( j >= 0 && i >= 0 ) {
			cdf_rf[j] = cdf_rp[i];
			j--;
			i--;
		}

		while( j >= 0 ) {
			cdf_rf[j] = 0;
		}
	}

  	/* FRf=FRp with shift for x>max(data_gf) */
	if( cdf_rf[npas - 1] < 1 ){

		double qmax = gsl_stats_quantile_from_sorted_data( obs_rp_sorted->data, 1,  obs_rp->size, cdf_rf[npas] );
		int i = npas - 1;
		while( x[i] >= qmax){
			i--;
		}
		i++;

		int j = npas - 2;
		while(j >= 0 && cdf_rf[j] == cdf_rf[npas - 1]){
			j--;
		}

		if(j < 0){
			printf("In CDFt, dev must be higher");
			exit(EXIT_FAILURE);
		}

		int diff = npas - 1 - (i > j ? i : j);
		for ( size_t d = 0; d < diff; d++ ) {
			cdf_rf[j + d] = cdf_rp[i + d];
		}
		int k = j + diff;
		while ( k < npas ) {
			cdf_rf[k]=1;
			k++;
		}

	}
	printf( "---------\n" );
  	for ( size_t i = 0; i < npas; i++ ) {
		printf( "%3.1f:%3.1f\n", x[i], cdf_rf[i] );
	}
	printf( "---------\n" );
	size_t n_unique = count_unique(cdf_rf, npas);
	double 
		x_interp[n_unique],
		y_interp[n_unique];
	average_byx(cdf_rf, x, npas, x_interp, y_interp);
	printf( "%i\n", (int) n_unique);
	for ( size_t index = 0; index < n_unique; index++ ) {
		printf( "%3.1f:%3.1f\n", x_interp[index], y_interp[index] );
	}
	printf( "---------\n" );



	/* Quantile-matching based on the new large-scale CDF and downscaled local-scale CDF. */
	{
		/*
		for ( size_t i = 0; i < data_gf2->size; i++ ) {
			printf( "%3.1f\n", cdf_rf[i] );
		}
		*/
		gsl_interp_accel *acc = gsl_interp_accel_alloc();
		gsl_interp *interp = gsl_interp_alloc( gsl_interp_linear, n_unique);
		gsl_interp_init( interp, x_interp, y_interp, n_unique );

		for ( size_t i = 0; i < data_gf2->size; i++ ) {
			qntl = gsl_stats_cdf_from_sorted_data( data_gf2_sorted->data, data_gf2_sorted->size, data_gf2->data[i] );
			cdft[i] = gsl_interp_eval( interp, x_interp, y_interp, qntl, acc);
		}
		gsl_interp_free( interp );
		gsl_interp_accel_free( acc );
	}

	 
	gsl_vector_free( data_gp2 );
	gsl_vector_free( data_gf2 );
	gsl_vector_free( obs_rp_sorted );
	gsl_vector_free( data_gp2_sorted );
	gsl_vector_free( data_gf2_sorted );
}

int main () {
	int npas = 20;
	size_t length = 10;
	/* allocate all vectors */
	gsl_vector
		*vrp = gsl_vector_alloc( length ),
		*vgp = gsl_vector_alloc( length ),
		*vgf = gsl_vector_alloc( length );
	double x[npas], cdf_rp[npas], cdf_gp[npas], cdf_gf[npas], cdf_rf[npas];
	double cdft[length];

	printf( "%i\n", count_unique( vrp->data, length ) );
	/* fill input vectors with some test data */
	for ( size_t index = 0; index < length; ++index ) {
		gsl_vector_set( vrp, index, (double) index );
		gsl_vector_set( vgp, index, (double) index + 2 );
		gsl_vector_set( vgf, index, (double) index + 2 );
	}
	printf( "%i\n", count_unique( vrp->data, length ) );
	correct_cdft(
		vrp, vgp, vgf, npas, 2, 
		x, cdf_rp, cdf_gp, cdf_gf, cdf_rf, cdft	
	);

	/* display the result to see it is correct */
	for ( size_t index = 0; index < length; ++index ) {
		printf( "%3.1f\n", cdft[index] );
	}

	/* be nice and tidy: release resources after use */
	gsl_vector_free( vrp );
	gsl_vector_free( vgp );
	gsl_vector_free( vgf);
}

