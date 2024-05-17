
#ifndef ERROR_H
#define ERROR_H

#include <iostream>
#include <string>
#include <vector>
using namespace std ;

// A = 0, C = 1, G = 3, T = 2
#define BASE_TO_NUM(base) (base >> 1 & 3)

#define BASE_A 0
#define BASE_C 1
#define BASE_G 3
#define BASE_T 2

// Toggle bit to complement
#define COMP_NUM(bnum) (bnum ^ 2)

char num2base( int x ) ;


/*
struct FromToPos
{
	FromToPos( char f, char t, unsigned char p ) : from(f), to(t), pos(p) {  } 
	char from ;
	char to ;
	unsigned char pos ;
	bool operator<( const FromToPos &x ) const {
		if ( this->pos != x.pos ) return this->pos < x.pos ;
		else if ( this->from != x.from ) return this->from < x.from ;
		else return this->to < x.to ;
	}
} ;
*/

class Error
{
	public:
		Error( const string &filename, double error_rate=1e-4 ) ;
		Error( double error_rate ) ;
		Error( vector<double> d_, double error_rate=1e-4 ) : d(d_), default_error_rate( error_rate ) {  } ;
		
		void print( ostream &os ) ;

		// prob(): from,to are given as actual bases 'A','C','G' or 'T':
		double prob( char from, char to, unsigned char pos ) ; 

		// probb(): from,to are already converted to number with BASE_TO_NUM():
		double probb( unsigned char from, unsigned char to, unsigned char pos ) ; 

		double default_error_rate ;

	private:
		vector<double> d ;
} ;

#endif

