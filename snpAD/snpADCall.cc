
#include <iostream>
#include <fstream>
#include <popt.h>
#include <string>
#include <sstream>
#include "Error.h"
#include "VCFfac.h"
#include "VCFfac_GL.h"


int main( int argc, const char* argv[] ) 
{

	char *error_file = 0 ;
	char *priors_string = 0 ;
	char *sample_name = 0 ;
	
	double refbias = 0. ;

	int flag_gls = 0 ;
	int flag_pplog = 0 ;

	double set_zero_priors = 0. ;

	int flag_print_version = 0 ;

        struct poptOption optionsTable[] = {
                        { "error",                              'e', POPT_ARG_STRING, &error_file, 0, "error profile file", 0 },
                        { "priors",                             'p', POPT_ARG_STRING, &priors_string, 0, "priors for calling AA,CC,GG,TT,AC,AG,AT,CG,CT,GT (comma separated, in this order)", 0 },
			{ "refbias",				'B', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &refbias, 0, "Reference bias", 0 },
			{ "name", 				'N', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &sample_name, 0, "Name of sample", 0 },
			{ "gl",				'l', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT , &flag_gls, 0, "Print GLs", 0 },
			{ "pplog",				'l', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT , &flag_pplog, 0, "Print PP field as log", 0 },
			{ "set_zero_priors",	'Z', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &set_zero_priors, 0, "Set priors that are 0 to this value", 0 },
			{ "version", 0, POPT_ARG_NONE, &flag_print_version, 0, "Print version", 0 },
                        POPT_AUTOHELP { NULL, 0, 0, NULL, 0 }
        };


        // parse options; do some sanity checking
        poptContext optCon ;
        optCon = poptGetContext( "snpADCall", argc, argv, optionsTable, 0 ) ;
        poptSetOtherOptionHelp( optCon, "[OPTIONS]* <input-file>" ) ;
        int rc = poptGetNextOpt( optCon ) ;
        if (rc != -1) {
                poptPrintUsage( optCon, stderr, 0 ) ;
                exit( 1 ) ;
        }
        const char* fn = poptGetArg( optCon ) ;

	#include "../version.h"
	if ( flag_print_version ) {
		cout << "snpAD " << SNPAD_VERSION << endl ;
		exit( 0 ) ;
	}

        if ( !fn ) {
                cerr << "Error: need snpAD input-file as argument." << endl ;
                poptPrintUsage( optCon, stderr, 0 ) ;
                exit( 1 ) ;
        }

        if ( ! error_file ) {
                cerr << "Error: need error profile file as argument." << endl ;
                poptPrintUsage( optCon, stderr, 0 ) ;
                exit( 1 ) ;
        }

	// read Error model
	Error e( error_file ) ;

	// parse priors from cmdl
	vector<double> priors ;
	if ( ! priors_string ) {
		cerr << "Need priors for each genotype (-p option)." << endl ;
		exit( 1 ) ;
	} else {
		string priors_string_s( priors_string ) ;
		ifstream priorsfile( priors_string ) ;
		if ( priorsfile ) {
			getline( priorsfile, priors_string_s ) ;
			priorsfile.close() ;
		} 

		istringstream ss( priors_string_s ) ;
		string s ;
		while ( getline( ss, s, ',' ) ) {
			double x = atof( s.c_str() ) ;
			if ( x == 0.0 && set_zero_priors == 0.0 ) {
				cerr << "prior #" << priors.size()+1 << " is zero! Consider setting option --set_zero_priors to a resonable small value (e.g. 1/genome_size)." << endl ;
				exit( 1 ) ;
			} else {
		//		cerr << "read prior " << x << endl ;
				priors.push_back( x ) ;
			}
		}
		if ( set_zero_priors != 0. ) {
			for ( auto it = priors.begin() ; it != priors.end() ; ++it ) {
				if ( *it == 0. ) *it = set_zero_priors ;
			}
		}
// DEBUG:
//		for ( auto it = priors.begin() ; it != priors.end() ; ++it ) {
//			cerr << *it << endl ;
//		}
	
	}
	

	// read stdin line by line and print VCF
	if ( flag_gls ) { // GT likelihood format
		ifstream inf(fn) ;
		VCFfac_GL fac( inf, /*reffa,*/ e, priors, sample_name, refbias ) ;
		cout << fac.header() ;
		while ( fac.can_has_more_VCF() ) { // ^o.o^
			cout << fac.chr() << "\t" << fac.pos() << "\t" << "." << "\t" << fac.ref() 
				<< "\t" << fac.alt() << "\t" << fac.qual_alt() << "\t" << "." 
				<< "\t" << fac.info() << "\t" << fac.format_call() << '\n' ;
			fac.next() ;		
		}
	} else { // No GT likelihoods
//		cerr << "Here!" << endl ;
		ifstream inf(fn) ;
		VCFfac fac( inf, /*reffa,*/ e, priors, sample_name, refbias, flag_pplog ) ;
		cout << fac.header() ;
//		cerr << "Here2!" << endl ;
		while ( fac.can_has_more_VCF() ) { // ^o.o^
//			cerr << "Can have more" << endl; 
			cout << fac.chr() << "\t" << fac.pos() << "\t" << "." << "\t" << fac.ref() 
				<< "\t" << fac.alt() << "\t" << fac.qual_alt() << "\t" << "." 
				<< "\t" << fac.info() << "\t" << fac.format_call() << '\n' ;
			fac.next() ;		
		}
//		cerr << "Here3!" << endl ;
	}
	// done.

}
