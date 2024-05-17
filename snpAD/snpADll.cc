
#include <iostream>
#include <fstream>
#include <sstream>
#include <popt.h>
#include "Data.h"
#include "MLHet_thread.h"
#include "MLHetRB_thread.h"

using namespace std ;

int main( int argc, const char* argv[] ) {

	char *error_file = 0 ;
	int cpus = 40 ;
	double error_prior = 1e-4 ; // 1/10,000 probability of a base being misread (e.g. P(A->C) + P(A->G) + P(A->T) == 1e-4; other bases analogous)

        struct poptOption optionsTable[] = {
                        { "prior_error",                              'p', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &error_prior, 0, "prior probability for a base being misread", 0 },
                        { "error",                              'e', POPT_ARG_STRING, &error_file, 0, "error profile file; overrides prior_error", 0 },
                        { "cpus",                               'c', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT , &cpus, 0, "use N CPUs in parallel", "N" },
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
	e->print( cout ) ;

	Data d( fn, *e, *e, cpus ) ;
	LL ll( d ) ;
	cout << "Enter:" << endl ;
	string s ;
	while ( cin >> s ) {
		istringstream ss( s ) ;
		vector<double> params ;
		string x ;
		while ( getline( ss, x, ',' ) ) {
			double d = atof( x.c_str() ) ;
			params.push_back( d ) ;
			cerr << "Read: " << d << endl ;
		}
		cerr << "LL=" << ll.run( &params ) << endl ;
	}
		
}
