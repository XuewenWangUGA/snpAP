
#include "Error.h"
#include <iostream>
#include <fstream>
#include <cmath>

char num2base( int x )
{
	switch (x) {
		case 0: return 'A' ;
		case 1: return 'C' ;
		case 2: return 'T' ;
		case 3: return 'G' ;
	}
}

double Error::prob( char from, char to, unsigned char pos )  
{ 
	return d[ BASE_TO_NUM(from) + BASE_TO_NUM(to)*4 + pos*16 ] ;  
}

double Error::probb( unsigned char from, unsigned char to, unsigned char pos )  
{ 
	return d[ from + to*4 + pos*16 ] ;  
}

Error::Error( double error_rate ) {
	default_error_rate = error_rate ;
	d.resize( 256*16 ) ;
	for ( size_t pos = 0 ; pos < 256 ; pos++ ) {
		for ( size_t from = 0 ; from < 4 ; from++ ) {
			for ( size_t to = 0 ; to < 4 ; to++ ) {
				d[ from + to*4 + pos*16 ] = from==to ? 1.-error_rate : error_rate/3. ;
			}
		}
	}
	
}

Error::Error( const string &filename, double error_rate )
{
	default_error_rate = error_rate ;
	ifstream in( filename ) ;

	// should be a #define somewhere. Matches the number of error-profiles that can be given for Data entries (see Data) * number of possible combinations of bases
	d.resize( 256*16 ) ;
	
	size_t pos ;
	char from, to ;
	double p ;
	while ( in >> pos >> from >> to >> p ) {
	//	cerr << "read " << pos << " " << from << " " << to << " " << p << endl ;
		d[ BASE_TO_NUM(from) + BASE_TO_NUM(to)*4 + pos*16 ] = p ;
	//	cerr << "extract: " << d[FTP_TO_FIELD( from, to, pos )] << endl ;
	}
}

void Error::print( ostream &os ) 
{
	for ( size_t pos = 0 ; pos < 256 ; ++pos ) {
		for( size_t from = 0 ; from < 4 ; ++from ) {
			for( size_t to = 0 ; to < 4 ; ++to ) { 
				if ( !std::isnan( d[ from + to*4 + pos*16 ] ) ) 
				os << pos << "\t" << num2base( from ) << "\t" << num2base( to ) << "\t" << d[ from + to*4 + pos*16 ] << endl ;
			}
		}
	}
}
