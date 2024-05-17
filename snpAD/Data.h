
#ifndef DATA_H
#define DATA_H

#include <limits>
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include "Error.h" 
using namespace std ;

struct Obs
{
	Obs( char b, bool ori, unsigned char p ) : base(b), strand(ori), pos(p) {  } ;
	Obs( char b, bool ori, unsigned char p, unsigned char c ) : base(b), strand(ori), pos(p), count(c) {  } ;
	unsigned char pos : 7 ; // error profile for this base
	bool strand : 1 ; // 1 for pos, 0 for neg; base will be complemented before querying error when negative
	unsigned char base : 2 ; // code for base A = 0, C = 1 , G = 3, T = 4 (see BASE_TO_NUM() in Error.h)
	unsigned char count : 6 ; // how often did I see this pase at this position? (Note: 0 is interpreted as 1, so max is 64)
	// comparison operation (will not check count!)
	
} ;

struct CoordRefAltFreq
{
	string chr ;
	size_t pos ;
	char ref ;
	char alt ;
	double alt_f ;
} ;

class Data {
	public:
//		Data( const string &filename, Error &e_, unsigned int threads_ = 1, Error &e_c_ ) ;
		Data( const string &filename, Error e_, Error e_c_, unsigned int threads_ = 1, unsigned int maxsites=numeric_limits<unsigned int>::infinity() )  ;
		void reestimate_error( vector<double>* priors  ) ; // update e by comparing the data to a provisional consensus
		void reestimate_error_full( vector<double>* priors  ) ; // update e by comparing the data to a provisional consensus
		void error_print( ostream &os ) { e.print( os ) ; }
		Error* get_error(  ) { Error *ret = new Error( e ) ; return ret ; } 
		deque<vector<Obs> > d ;
		vector<CoordRefAltFreq> poss ;
		Error e ;
		Error e_c ;
		unsigned int threads ;
		vector<double> basecomp ;
		unsigned char max_ec ; // maximum error class we observed in the data
	private:
		void reestimate_basecomp( ) ; // update base composition with error profile (called as last step from reestimate_error)
} ;

#endif
