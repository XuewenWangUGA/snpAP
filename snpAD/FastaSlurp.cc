#include "FastaSlurp.h"
#include <stdexcept>

FastaSlurp::FastaSlurp( ifstream &in ) 
{
	string s ;
	getline( in, s ) ;
	if ( s[0] != '>' ) 
		throw invalid_argument( "Fasta file does not start with '>'" ) ;
	string chr = s.substr( 1, s.find_first_of(" \t")-1 ) ;
	while ( ! in.eof() ) 
	{
		getline( in, s ) ;
		if ( s[0] != '>' ) 
			chr_seq[chr] += s ;
		else 
			chr = s.substr( 1, s.find_first_of(" \t")-1 ) ;
	}
//	cerr << "read last '" << chr << "'" << endl ;
}

char FastaSlurp::fetch( string chr, size_t pos ) 
{
//	cerr << "Fetch " << chr << " " << pos << endl ;
	if ( chr_seq[chr].length() >= pos ) {
		return chr_seq.at(chr)[pos-1] ;
	} else {
		return 'N' ;
	}
}

