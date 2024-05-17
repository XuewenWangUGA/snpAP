
#include <iostream>
#include <fstream>
#include <popt.h>
#include "Data.h"
#include "MLHet_thread.h"
#include "MLHet_thread_tstv.h"
#include "MLHetRB_thread.h"

using namespace std ;

int main( int argc, const char* argv[] ) {

	char *error_file = 0 ;
	char *priors_out = 0 ;
	char *errors_out = 0 ;
	char *refbias_out = 0 ;

	// These 5 variables are for future development. They are not used.
	int flag_less_mem = 0 ;
	int flag_norepeat = 0 ;
	int flag_fast = 0 ;
	int verbose = 0 ;
	int flag_refbias = 0 ; 
	int flag_contam = 0 ; // XXX
	int full_error = 1 ;
	int no_full_error = 0 ;
	char *contam_error_file ; // XXX
	int flag_print_version = 0 ;
	
	double max_gtfreq = 0.001 ;

	int cpus = 40 ;
	
	double error_prior = 1e-4 ; // 1/10,000 probability of a base being misread (e.g. P(A->C) + P(A->G) + P(A->T) == 1e-4; other bases analogous)

        struct poptOption optionsTable[] = {
                        { "prior_error",                              'p', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &error_prior, 0, "prior probability for a base being misread", 0 },
                        { "error",                              'e', POPT_ARG_STRING, &error_file, 0, "error profile file; overrides prior_error", 0 },
                        { "cpus",                               'c', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT , &cpus, 0, "use N CPUs in parallel", "N" },
                        { "priors_out",                             'o', POPT_ARG_STRING, &priors_out, 0, "write best fitting priors to file (order: AA,CC,GG,TT,AC,AG,AT,CG,CT,GT)", 0 },
                        { "errors_out",                             'O', POPT_ARG_STRING, &errors_out, 0, "write best fitting error profile to file", 0 },
			{ "no_full_error",				'X', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT | POPT_ARGFLAG_DOC_HIDDEN , &no_full_error, 0, "Reestimate error profile by taking the most likely call instead of integrating over all possible genotypes", 0 },
			{ "lessmem",				'M', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT | POPT_ARGFLAG_DOC_HIDDEN , &flag_less_mem, 0, "Estimate priors with a less memory intensive procedure", 0 },
			{ "norepeat",				'R', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT , &flag_norepeat, 0, "Run one cycle of estimating priors and error profile (useful for when parameter error is given and closely matches true error profile)", 0 },
			{ "verbose",				'v', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT , &verbose, 0, "Print more intermediate messages on stdout", 0 },
			{ "fast",				'F', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT , &flag_fast, 0, "Estimate frequency of ts and tv heterozygotes instead of all 6 types of heterozygotes. Runs faster. Not implemented for ref-bias.", 0 },
			{ "refbias",				'B', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT , &flag_refbias, 0, "Estimate reference bias", 0 },
			{ "refbias_out",			0, POPT_ARG_STRING, &refbias_out, 0, "write best fitting reference bias estimate to file", 0 },
                        { "max_gt_freq",                              0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &max_gtfreq, 0, "Maximum genotype frequency for any individual heterozygous genotype (increase this parameter for SNP-capture data or species with heterozygosity >>1e-3)", 0 },
			{ "contam",				'C', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT | POPT_ARGFLAG_DOC_HIDDEN , &flag_contam, 0, "Estimate contamination. Requires --refbias. Incompatible with --lessmem. Option is under development, don't use!", 0 },
                        { "contam_error",                              'E', POPT_ARG_STRING | POPT_ARGFLAG_DOC_HIDDEN , &contam_error_file, 0, "Error profile for contamination (default: identical to error profile)", 0 },
			{ "version", 0, POPT_ARG_NONE, &flag_print_version, 0, "Print version", 0 },
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

	#include "../version.h"
	if ( flag_print_version ) {
		cout << "snpAD " << SNPAD_VERSION << endl ;
		exit( 0 ) ;
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
		e = new Error( error_file, error_prior ) ; // profile
	} else {
		e = new Error( error_prior ) ; // flat error profile
	}

	if ( no_full_error ) full_error = 0 ; // flip for historical reasons.

	// Best parameters:
	Error *best_error ;
	vector<double> *priors ; // priors are the genotype priors
	double best_rb ; //best reference bias, only when flag is set.

	if ( flag_refbias ) { 
		// Run first iteration
		Data d( fn, *e, *e, cpus ) ;
		MLHetRB_thread ml( d, max_gtfreq ) ;
		ml.write( cout ) ;
		cout << "ll: " << ml.ll << endl ;
		double prior_ll = ml.ll ;
		priors = ml.get_priors(  ) ;
		best_rb = ml.get_refbias(  ) ;
		if ( verbose ) cout << "Refbias: " << best_rb << endl ;

		// further iterations (if needed)
		while ( !flag_norepeat ) {
			best_error = d.get_error(  ) ;
			if ( full_error ) 
				d.reestimate_error_full( priors ) ;
			else 
				d.reestimate_error( priors ) ;
			if ( verbose ) d.error_print( cout ) ;
			MLHetRB_thread ml2( d, max_gtfreq ) ;
			ml2.write( cout ) ;
			cout << "ll: " << ml2.ll << endl ;
			if ( ml2.ll - prior_ll < 3. ) break ; 
			delete priors ;
			priors = ml2.get_priors(  ) ;
			prior_ll = ml2.ll ;
			best_rb = ml2.get_refbias(  ) ;
		}
		

	} else {  // No reference bias
	    	if ( flag_fast ) {
			// Run first iteration
			Data d( fn, *e, *e, cpus ) ;
			MLHet_thread_tstv ml( d, max_gtfreq ) ;
			ml.write( cout ) ;
			cout << "ll: " << ml.ll << endl ;
			double prior_ll = ml.ll ;
			priors = ml.get_priors(  ) ;

			// further iterations (if needed)
			while ( !flag_norepeat ) {
				best_error = d.get_error(  ) ;
				if ( full_error ) 
					d.reestimate_error_full( priors ) ;
				else 
					d.reestimate_error( priors ) ;
				if ( verbose ) d.error_print( cout ) ;
				MLHet_thread_tstv ml2( d, max_gtfreq ) ;
				ml2.write( cout ) ;
				cout << "ll: " << ml2.ll << endl ;
				if ( ml2.ll - prior_ll < 3. ) break ; 
				delete priors ;
				priors = ml2.get_priors(  ) ;
				prior_ll = ml2.ll ;
			}
		} else { 
			// Run first iteration
			Data d( fn, *e, *e, cpus ) ;
			MLHet_thread ml( d, max_gtfreq ) ;
			ml.write( cout ) ;
			cout << "ll: " << ml.ll << endl ;
			double prior_ll = ml.ll ;
			priors = ml.get_priors(  ) ;

			// further iterations (if needed)
			while ( !flag_norepeat ) {
				best_error = d.get_error(  ) ;
				if ( full_error ) 
					d.reestimate_error_full( priors ) ;
				else 
					d.reestimate_error( priors ) ;
				if ( verbose ) d.error_print( cout ) ;
				MLHet_thread ml2( d, max_gtfreq ) ;
				ml2.write( cout ) ;
				cout << "ll: " << ml2.ll << endl ;
				if ( ml2.ll - prior_ll < 3. ) break ; 
				delete priors ;
				priors = ml2.get_priors(  ) ;
				prior_ll = ml2.ll ;
			}
		}
	}

	// print priors
	if ( priors_out ) {
		// Note: we need to swap GG and TT which are in order TT, GG during estimation
		iter_swap( priors->begin()+2, priors->begin()+3 ) ;
		ofstream po( priors_out ) ;
		for ( size_t i = 0 ; i < priors->size( ) ; ++i ) {
			if ( i !=0 ) po << "," ;
			po << priors->at( i ) ;
		}
		po << endl ;
		po.close(  ) ;
	}
	
	// print error profile
	if ( errors_out ) {
		ofstream eo( errors_out ) ;
		best_error->print( eo ) ;
		eo.close(  ) ;
	}

	// print reference bias (if needed)
	if ( refbias_out && flag_refbias ) {
		ofstream rbo( refbias_out ) ;
		rbo << best_rb << endl ;
		rbo.close(  ) ;
	}

	cout << "Done." << endl ;

	exit( 0 ) ;

}

