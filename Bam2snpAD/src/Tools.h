
#include <string>
#include <algorithm>

inline char rc( char c ) 
{
	switch (c) 
	{
		case 'A': return 'T' ;
		case 'C': return 'G' ;
		case 'G': return 'C' ;
		case 'T': return 'A' ;
		default: return c ;
	}
}

inline void revcomp( string &s ) 
{
	reverse( s.begin(), s.end() ) ;
	for ( size_t i = 0 ; i < s.length() ; ++i ) 
		s[i] = rc( s[i] ) ;
}


