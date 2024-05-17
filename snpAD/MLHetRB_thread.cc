
#include "MLHetRB_thread.h"
#include <numeric>
#include <thread>
#include <stdexcept>
#include <iomanip>
#include "log_add.h"


vector<string> genotypes_RBhet = { "AC", "AG", "AT", "CG", "CT", "GT" } ;

void MLHetRB_thread::write( ostream &os ) 
{
	double p_hom = 1.-accumulate( params.begin(), prev(params.end()), 0. ) ;
	
	os << setprecision(10) ;
	for ( size_t i = 0 ; i < genotypes_RBhet.size() ; ++i )
		os << "P(" << genotypes_RBhet[i] << ")\t" << params[i] << endl ;
	os << "P(AA)\t" << p_hom*data.basecomp[0] << endl ;
	os << "P(CC)\t" << p_hom*data.basecomp[1] << endl ;
	os << "P(GG)\t" << p_hom*data.basecomp[3] << endl ;
	os << "P(TT)\t" << p_hom*data.basecomp[2] << endl ;
	os << "P(ref)\t" << params[6] << endl ;
	
}

vector<double>* MLHetRB_thread::get_priors(  ) {
	double p_hom = 1.-accumulate( params.begin(), prev(params.end()), 0. ) ;
	
	// first homozygous then hets
	vector<double> *ret = new vector<double>(10,0.) ;
	for ( size_t i = 0 ; i < 4 ; i++ ) ret->at(i) = p_hom*data.basecomp[i] ;

	for ( size_t i = 0 ; i < genotypes_RBhet.size() ; ++i )
		ret->at(4+i) = params[i] ;

	return ret ;
}


void thread_worker_RBhet( vector<double> params, double p_hom, Data *d, size_t beg, size_t end, long double* return_val ) 
{
	
//	cerr << "Worker" << endl ;
	long double llik = 0.;
	double refbias = params[6] ;
	for( auto it = d->d.begin()+beg ; it != d->d.begin()+end ; ++it ) { // all sites
		long double llik_site = -numeric_limits<long double>::infinity() ;
		for ( size_t it_gts = 0 ; it_gts < genotypes_RBhet.size() ; ++it_gts ) { // all possible genotypes
			long double prior_gt = params[it_gts] ;
			double r_1 = 0.5 ; // GT1 
			double r_2 = 0.5 ; // GT2 
// XXX: this condition is not needed since r_1 + r_2 = 1 in all cases
//			if ( genotypes_RBhet[it_gts][0] != genotypes_RBhet[it_gts][1] ) {
				if ( genotypes_RBhet[it_gts][0] == d->poss[ it - d->d.begin() ].ref ) {
					r_1 += refbias ; r_2 -= refbias ;
				} else if ( genotypes_RBhet[it_gts][1] == d->poss[ it - d->d.begin() ].ref  ) {
					r_2 += refbias ; r_1 -= refbias ;
				}
//			}
			long double llik_gt = 0. ;
			for ( vector<Obs>::iterator it2 = it->begin() ; it2 != it->end() ; ++it2 ) { // all reads at this site
				if ( it2->strand ) 
					llik_gt += (1+it2->count) * log( r_1*d->e.probb( BASE_TO_NUM(genotypes_RBhet[it_gts][0]), it2->base, it2->pos )
							+ r_2*d->e.probb( BASE_TO_NUM(genotypes_RBhet[it_gts][1]), it2->base, it2->pos ) ) ;
				else
					llik_gt += (1+it2->count) * log( r_1*d->e.probb( COMP_NUM(BASE_TO_NUM(genotypes_RBhet[it_gts][0])), COMP_NUM(it2->base), it2->pos )
							+ r_2*d->e.probb( COMP_NUM(BASE_TO_NUM(genotypes_RBhet[it_gts][1])), COMP_NUM(it2->base), it2->pos ) ) ;
			}
			llik_site = log_add( llik_site, log(prior_gt) + llik_gt ) ;
		}
		for ( unsigned char homgt = 0 ; homgt < 4 ; ++homgt ) {
			long double prior_gt = d->basecomp[homgt] * p_hom ;
//			cerr << "Prior homozygous " << int(homgt) << " = " << prior_gt << endl ;
			long double llik_gt = 0. ;
			for ( vector<Obs>::iterator it2 = it->begin() ; it2 != it->end() ; ++it2 ) { // all reads at this site
				if ( it2->strand )
					llik_gt += (1+it2->count) * log( d->e.probb( homgt, it2->base, it2->pos ) ) ;
				else 
					llik_gt += (1+it2->count) * log( d->e.probb( COMP_NUM(homgt), COMP_NUM(it2->base), it2->pos ) ) ;
			}

			llik_site = log_add( llik_site, log(prior_gt) + llik_gt ) ;
		}
		llik += llik_site ;
	}
//	cerr << "Worker llik=" << llik << endl ;
	*return_val = llik ;
}

// Objective function
double obj_fun_RBhet_threaded( const vector<double> &params, vector<double> &grad, void *d_ ) 
{
        Data *d = reinterpret_cast<Data*>(d_) ;
 	
	// sum over probabilities of genotypes to get probability for homozygous
	double p_hom = 1.-accumulate( params.begin(), prev(params.end()), 0. ) ;
	if ( p_hom < 0. ) return -HUGE_VAL ;


	// calculate how to split workload based on the number of threads
	size_t total_data = d->d.size() ;
	size_t steps = total_data/d->threads ;
	vector<thread> threads ;
	vector<long double> return_vals(d->threads, 0.) ;

	// spawn threads
	for ( size_t i = 0 ; i < d->threads ; ++i ) {
		if ( i < d->threads -1 ) {
			threads.push_back( thread( thread_worker_RBhet, params, p_hom, d, i*steps, (i+1)*steps, &(return_vals[i]) ) ) ;
		} else {
			threads.push_back( thread( thread_worker_RBhet, params, p_hom, d, i*steps, total_data, &(return_vals[i]) ) ) ;
		}
//		cerr << "Thread created with id " << threads[i].get_id() << endl ;
	}

	// join threads
	for ( size_t i = 0 ; i < d->threads ; ++i ) {
		threads[i].join() ;
	}
	long double llik = accumulate( return_vals.begin(), return_vals.end(), 0. ) ;

	cerr << "logL( " ;
	cerr << "P(refbias)=" << params[6] << " " ;
	for( size_t i = 0 ; i < params.size()-1 ; ++i ) {
		cerr << "P(" << genotypes_RBhet[i] << ")=" << params[i] << " " ;
	}
	cerr << " P(AA)=" << p_hom*d->basecomp[0] ;
	cerr << " P(CC)=" << p_hom*d->basecomp[1] ;
	cerr << " P(GG)=" << p_hom*d->basecomp[3] ;
	cerr << " P(TT)=" << p_hom*d->basecomp[2] ;
	cerr << " | Data ) == " << setprecision(14) << llik << endl ;

        return llik ;
}

MLHetRB_thread::MLHetRB_thread( Data data_, double max_gtfreq, vector<double> *initial_params ) : data(data_)
{
	size_t dim = 6 + 1 ; // prob for each het genotype, then reference bias
//	nlopt::opt opt( nlopt::LN_NELDERMEAD, dim ) ; 
	nlopt::opt opt( nlopt::LN_BOBYQA, dim ) ; 

	// sanity check
	if ( data.d.size() != data.poss.size() ) {
		cerr << "Need reference information, alternative allele and alternative allele frequency for all sites in input file" << endl ;
		throw( invalid_argument( "Wrong inputfile format" ) ) ;
	}


	// objective function
	opt.set_max_objective( obj_fun_RBhet_threaded, reinterpret_cast<void*>(&data) ) ;

	// bounds
	vector<double> lb( dim, 0. ) ;
	vector<double> ub( dim, max_gtfreq ) ;
	// XXX: all of this stuff needs to be adjustable.
	ub[6] = 0.5 ; // refbias
	opt.set_lower_bounds( lb ) ;
	opt.set_upper_bounds( ub ) ;

	// set stopping criteria
//	vector<double> sc( dim, 1.0e-5 ) ; // precision on parameters: 1e-3 relative (arbitrary)
//	opt.set_xtol_abs( sc ) ;
	opt.set_xtol_rel( 1.0e-3 ) ;
//	opt.set_maxtime( 60 ) ; // do not run longer than a minute

	// set initial parameters to 0.
	params.resize( dim, 0. ) ;
	if ( initial_params ) {
		if ( initial_params->size() != dim ) {
			throw( logic_error( "wrong number of params" ) ) ;
		}
		params = *initial_params ;
	}

	// run optimization
	nlopt::result r = opt.optimize( params, ll ) ;

	cerr << "Return val: " << r << endl ;
}

