
#ifndef FASTASLURP_H
#define FASTASLURP_H


#include <iostream>
#include <fstream>
#include <map>
using namespace std ;


class FastaSlurp
{
	public: 
		FastaSlurp( ifstream &in ) ;
		char fetch( string chr, size_t pos ) ;

	private:
		map<string,string> chr_seq ;
} ;

#endif
