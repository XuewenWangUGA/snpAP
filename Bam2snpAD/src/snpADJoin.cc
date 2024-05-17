
#include "SnpADFile.h"
#include <limits.h>

int main( int argc, char *argv[] )
{
	vector<SnpADFile*> files ;
	for ( int i = 1 ; i < argc ; ++i ) {
		ifstream *in = new ifstream( argv[i] ) ;
		files.push_back( new SnpADFile( *in ) ) ;
	}

	// loop through lines in files until first file ends, then bail out
	while ( 1 ) {
		// identify the highest value in all files 
		// and check whether all have the same coordinate
		size_t lowest = INT_MAX ;
		size_t lowest_id = 0 ; // set to the entry that has lowest coordinate
		bool all_eof = 1 ;
		for ( size_t it = 0 ; it < files.size() ; ++it ) {
			if ( files[it]->eof() == false ) {
				all_eof = 0 ;
				if ( files[it]->get_loc() < lowest ) {
					lowest = files[it]->get_loc() ;
					lowest_id = it ;
				}
			}
		}
		if ( all_eof ) return 0 ; // all files ended

		// print all that have lowest coordinate, print 0 for rest
		// increase only those at lowest coordinate
	
		cout << files[lowest_id]->get_chr() << "\t" << files[lowest_id]->get_loc() << "\t" << files[lowest_id]->get_ref() ;
		for ( size_t it = 0 ; it < files.size() ; ++it ) {
			if ( files[it]->eof() == false && files[it]->get_loc() == lowest ) {
				for ( vector<string>::const_iterator jt = files[it]->get_data().begin() ; 
					jt < files[it]->get_data().end() ; ++jt ) {
					cout << "\t" << *jt ;
				}
				files[it]->next() ;
			}
		}
		cout << '\n' ;
	}
	cout.flush(  ) ;
}
