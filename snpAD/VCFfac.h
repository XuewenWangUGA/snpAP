
#include <iostream>
#include <string>

#include "FastaSlurp.h"
#include "Error.h"

using namespace std ;

class VCFfac
{
	public:
		VCFfac( ifstream &in_, /* FastaSlurp &fa_, */ Error &e_, vector<double> &priors_, const char* sample_name_, double ref_bias_=0., int pplog_=0 ) : in(in_), /* fa(fa_), */ e(e_), 
									priors(priors_), sample_name(sample_name_), eof(0), PP( 10, 0 ), logPP( 10, 0. ), log_gt( 10, 0. ), 
									basesfwd( 4, 0 ), basesrev( 4, 0 ), ref_bias( ref_bias_ ), pplog( pplog_ )
			{ next() ; } ;
		void next() ; // parse the next line or set eof

		bool can_has_more_VCF() { return !eof ; }
		string chr() { return chr_ ; } 
		size_t pos() { return pos_ ; }
		char ref() { return ref_ ; }
		string alt() { return alt_ ; }
		size_t qual_alt() { return QUAL ; } 
		string header() ;

		string info() { return "." ; } 
		string format_call() ; 
		

	private:
		ifstream &in ;
//		FastaSlurp &fa ;
		Error &e ;
		vector<double> priors ; // priors for 10 genotypes
		const char* sample_name ;

		bool eof ;

		string chr_ ;
		size_t pos_ ;
		char ref_ ;
		string alt_ ;
		size_t QUAL ;
		string GT ;
		vector<size_t> PP ; // with prior and normalized by all GTs
		vector<double> logPP ; // with prior and normalized by all GTs
		size_t GQ ;
		size_t DP ;
		vector<long double> log_gt ;
		vector<size_t> basesfwd ;
		vector<size_t> basesrev ;

		double ref_bias ;
		int pplog ;
		
} ;

