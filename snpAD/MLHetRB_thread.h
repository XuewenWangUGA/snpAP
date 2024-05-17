
#ifndef MLHETRBTHREAD_H
#define MLHETRBTHREAD_H

#include <iostream>
#include <vector>
#include <nlopt.hpp>
#include "Data.h"

// Reference Bias 
class MLHetRB_thread
{
        public:
                MLHetRB_thread( Data data_, double max_gtfreq, vector<double> *initial_params = 0 ) ;
                double ll ; // log likelihood
                vector<double> params ;
		void write( ostream &os ) ;
		vector<double>* get_priors(  ) ;
		double get_refbias(  ) { return params[6] ; }
		
        private:
                Data data ;
        
} ;

#endif

