
#include <iostream>
#include <fstream>
#include <sstream>
#include <popt.h>
#include "Data.h"
#include "MLHet_thread.h"
#include "MLHetRB_thread.h"

using namespace std ;

void estimate_ci( LL &ll, vector<double> &priors,  vector<double> &ret_low,  vector<double> &ret_high, double delta, int multtest ) ;

int main( int argc, const char* argv[] ) {

	char *error_file = 0 ;
	int cpus = 40 ;
	char *priors_string ;
	double error_prior = 1e-4 ; // 1/10,000 probability of a base being misread (e.g. P(A->C) + P(A->G) + P(A->T) == 1e-4; other bases analogous)
	double precision = 0.01 ;
	int flag_multtesting = 0 ;

        struct poptOption optionsTable[] = {
                        { "error",                              'e', POPT_ARG_STRING, &error_file, 0, "error profile file; overrides prior_error", 0 },
                        { "multiple_testing",                              'M', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT, &flag_multtesting, 0, "Correct for multiple testing; recommended when summing values to get CI for heterozygosity", 0 },
                        { "genotype_frequencies",                             'f', POPT_ARG_STRING, &priors_string, 0, "best estimate for frequency of AC,AG,AT,CG,CT,GT (comma separated, in this order)", 0 },
                        { "cpus",                               'c', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT , &cpus, 0, "use N CPUs in parallel", "N" },
			{ "delta",      'D', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &precision, 0, "find parameters with likelihood \\approx best-1.92 +/- delta" },
                        { "flat_error",                              'E', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT | POPT_ARGFLAG_DOC_HIDDEN, &error_prior, 0, "prior probability for a base being misread; use -e for full error profile", 0 }, // doesn't make sense as an option, but I leave it in for testing.
                        POPT_AUTOHELP { NULL, 0, 0, NULL, 0 }
        };


        // parse options; do some sanity checking
        poptContext optCon ;
        optCon = poptGetContext( "snpAD", argc, argv, optionsTable, 0 ) ;
        poptSetOtherOptionHelp( optCon, "[OPTIONS]* <input-file>" ) ;
        int rc = poptGetNextOpt( optCon ) ;
        if (rc != -1) {
                poptPrintUsage( optCon, stderr, 0 ) ;
                exit( 1 ) ;
        }

        const char* fn = poptGetArg( optCon ) ;

        if ( !fn ) {
                cerr << "Error: need snpAD input-file as argument." << endl ;
                poptPrintUsage( optCon, stderr, 0 ) ;
                exit( 1 ) ;
        }

	// Error profile
	Error *e ;
	if ( error_file ) {
		e = new Error( error_file ) ; // profile
	} else {
		e = new Error( error_prior ) ; // flat error profile
	}

	// parse best estimate for genotype frequencies from cmdl
	vector<double> priors ;
	if ( ! priors_string ) {
		cerr << "Need priors for each genotype (-p option)." << endl ;
		exit( 1 ) ;
	} else {
		istringstream ss( priors_string ) ;
		string s ;
		while ( getline( ss, s, ',' ) ) {
			double x = atof( s.c_str() ) ;
			if ( x == 0.0  ) {
				cerr << "prior #" << priors.size()+1 << " is zero!" << endl ;
				exit( 1 ) ;
			} else {
				priors.push_back( x ) ;
			}
		}
	}
	if ( priors.size() != 6 ) {
		cerr << "Only the heterozygous genotype frequencies should be given for option -f!" << endl ;
		exit( 1 ) ;
	}

	Data d( fn, *e, *e, cpus ) ;
	LL ll( d ) ;

	// These are the return vectors. They are just initialized with /priors/ to give them the right size.
	vector<double> priors_low( priors ) ; 
	vector<double> priors_high( priors ) ;

	estimate_ci( ll, priors, priors_low, priors_high, precision, flag_multtesting ) ;

	const char* gts[6] = { "AC", "AG", "AT", "CG", "CT", "GT" } ;

	cout << "Param#\tEst\tLowCI\tHighCI" << endl ;
	for ( size_t i = 0 ; i < priors.size() ; ++i ) {
		cout << gts[i] << "\t" << priors[i] << "\t" << priors_low[i] << "\t" << priors_high[i] << endl ;
	}
	
}

// This is spectacularly ugly. It feels like I should have used a modified version of Brent's Method. 
void estimate_ci( LL &ll, vector<double> &priors,  vector<double> &ret_low,  vector<double> &ret_high, double delta, int multtest ) {  

	double pchisq_cutoff = 3.841459/2. ;
	if ( multtest ) 
		pchisq_cutoff = 6.960401/2. ;


	double best_ll = ll.run( &priors ) ;
	double target_ll = best_ll - pchisq_cutoff ; // we calculate the CI as those parameters that yield this target LL

	cerr << "Best ll: " << best_ll << endl ;
	cerr << "Target ll: " << target_ll << endl ;
	cerr << "priors:" ; for ( size_t i = 0 ; i < priors.size() ; ++i ) cerr << " " << priors[i] ; cerr << endl ;

	for ( size_t i = 0 ; i < priors.size() ; ++i ) {

		// Lower limit

		double downhill = priors[i] ;
		double downhill_ll = best_ll ;

		vector<double> test( priors ) ;
		test[i] *= 0.9 ; // start with 10% lower
		double uphill = test[i] ;
		double uphill_ll = ll.run( &test ) ;


		while ( abs(downhill_ll-target_ll) > delta && abs(uphill_ll-target_ll) > delta ) {
			// downhill is not right of target 
			if ( !(downhill_ll > target_ll) ) {
				// extend by 50% to right
				test = priors ;
				test[i] = downhill+abs(downhill-uphill)/2. ;
				// over the limit (shouldn't happen)
				if ( test[i] > priors[i] ) {
					test[i] = (priors[i]+downhill)/2. ;
				}
				downhill = test[i] ;
				downhill_ll = ll.run( &test ) ;
				continue ;
			} 
			// uphill is not left of target
			if ( !(uphill_ll < target_ll) ) {
				// extend by 50% to left
				test = priors ;
				test[i] = uphill-abs(downhill-uphill)/2. ;
				// over the limit (shouldn't happen)
				if ( test[i] < 1e-6 ) {
					test[i] = (uphill+1e-6)/2. ;
				}
				// can only happen over many iterations
				if ( test[i] < 1e-6+1e-10 ) {
					uphill = 0. ; 
					uphill_ll = target_ll ;
					continue ;
				}
				uphill = test[i] ;
				uphill_ll = ll.run( &test ) ;
				continue ;
			}
			// the target is in between these two. Extrapolate new point.
			test = priors ;
			test[i] = uphill+(target_ll-downhill_ll) * (downhill-uphill) / (uphill_ll-downhill_ll) ;
			double test_ll = ll.run( &test ) ;
			if ( test_ll < target_ll ) {
				uphill = test[i] ;
				uphill_ll = test_ll ;
			} else {
				downhill = test[i] ;
				downhill_ll = test_ll ;
			}
		}

		if ( abs(downhill_ll-target_ll) <= delta ) 
			ret_low[i] = downhill ;
		else if ( abs(uphill_ll-target_ll) <= delta ) 
			ret_low[i] = uphill ;

		// Upper limit
		downhill = priors[i] ;
		downhill_ll = best_ll ;

		test = priors ;
		test[i] *= 1.1 ; // start with 10% lower
		uphill = test[i] ;
		uphill_ll = ll.run( &test ) ;


		while ( abs(downhill_ll-target_ll) > delta && abs(uphill_ll-target_ll) > delta ) {
			// downhill is not left of target 
			if ( !(downhill_ll > target_ll) ) {
				// extend by 50% to left
				test = priors ;
				test[i] = downhill - abs(downhill-uphill)/2. ;
				// over the limit (shouldn't happen)
				if ( test[i] < priors[i] ) {
					test[i] = priors[i]+(priors[i]-downhill)/2. ;
				}
				downhill = test[i] ;
				downhill_ll = ll.run( &test ) ;
				continue ;
			} 
			// uphill is not right of target
			if ( !(uphill_ll < target_ll) ) {
				// extend by 50% to left
				test = priors ;
				test[i] = uphill + abs(downhill-uphill)/2. ;
				// over the limit (shouldn't happen)
				if ( test[i] > 0.2 ) { // XXX: this is fix and shouldn't be.
					test[i] = uphill+(.2-uphill)/2. ;
				}
				uphill = test[i] ;
				uphill_ll = ll.run( &test ) ;
				continue ;
			}
			// the target is in between these two. Extrapolate new point.
			test = priors ;
			test[i] = downhill+(target_ll-downhill_ll) * (uphill-downhill) / (uphill_ll-downhill_ll) ;
			double test_ll = ll.run( &test ) ;
			if ( test_ll < target_ll ) {
				uphill = test[i] ;
				uphill_ll = test_ll ;
			} else {
				downhill = test[i] ;
				downhill_ll = test_ll ;
			}
		}

		if ( abs(downhill_ll-target_ll) <= delta ) 
			ret_high[i] = downhill ;
		else if ( abs(uphill_ll-target_ll) <= delta ) 
			ret_high[i] = uphill ;

	}

}

