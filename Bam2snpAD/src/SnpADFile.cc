
#include "SnpADFile.h"
#include <sstream> 
#include <cstdlib>

SnpADFile::SnpADFile( istream &is_ ) : is(is_), eof_(0)
{
	string line ;
	getline( is, line ) ;
	istringstream iss( line, istringstream::in ) ;
	iss >> chr ;
	iss >> loc ;
	iss >> ref ;
	data.clear() ;
	string val ;
	while ( !iss.eof() ) { 
		iss >> val ; 
		data.push_back( val ) ;
	}
	if ( iss.fail() ) { 
		cerr << "File doesn't contain entries or irregular tab at end of line." << endl ;
		exit( 0 ) ;
	}

}

bool SnpADFile::next(  )
{
	string line ;
	getline( is, line ) ;
	istringstream iss( line, istringstream::in ) ;
	iss >> chr ;
	iss >> loc ;
	iss >> ref ;
	data.clear() ;
	string val ;
	while ( !iss.eof() ) { 
		iss >> val ; 
		data.push_back( val ) ;
	}
	if ( iss.fail() ) { eof_=1 ; return 0 ; }
	else return 1 ;
}
