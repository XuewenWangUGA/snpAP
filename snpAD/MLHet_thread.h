
#ifndef MLHETTHREAD_H
#define MLHETTHREAD_H

#include <iostream>
#include <vector>
#include <nlopt.hpp>
#include "Data.h"

class MLHet_thread
{
        public:
                MLHet_thread( Data data_, double max_gtfreq, vector<double> *initial_params = 0 ) ;
                double ll ; // log likelihood
                vector<double> params ;
		void write( ostream &os ) ; 
		vector<double>* get_priors(  ) ;
        private:
                Data data ;
        
} ;

class LL
{
	public:
		LL( Data data_ ) ;
		double run( vector<double> *params ) ;
	private:
		Data data ;
} ;

#endif

