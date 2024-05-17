
#include "log_add.h"
#include <cmath>
#include <algorithm>

long double log_add( long double x, long double y ) 
{
	if( x == -std::numeric_limits<long double>::infinity() ) return y ;
    	if( y == -std::numeric_limits<long double>::infinity() ) return x ;
	return std::max( x, y ) + log1p( exp( -fabs(x - y) ) ) ;
}
double log_add( double x, double y ) 
{
	if( x == -std::numeric_limits<double>::infinity() ) return y ;
    	if( y == -std::numeric_limits<double>::infinity() ) return x ;
	return std::max( x, y ) + log1p( exp( -fabs(x - y) ) ) ;
}

