
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdint.h>

using namespace std ;

class SnpADFile
{
	public:
		SnpADFile( istream &is_ ) ;
		bool next(  ) ;
		string get_chr(  ) const
			{ return chr ; } ;
		uint32_t get_loc(  ) const
			{ return loc ; } ;
		char get_ref(  ) const 
			{ return ref[0] ; } ;
		vector<string> &get_data(  ) 
			{ return data ; } ; 
		bool eof(  ) const
			{ return eof_ ; } ;
	private:
		istream &is ;
		bool eof_ ;
		string chr ;
		uint32_t loc ;
		string ref ;
		vector<string> data ;
} ;
