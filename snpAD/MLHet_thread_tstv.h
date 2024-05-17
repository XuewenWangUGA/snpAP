
#ifndef MLHETTHREADTSTV_H
#define MLHETTHREADTSTV_H

#include <iostream>
#include <vector>
#include <nlopt.hpp>
#include "Data.h"

class MLHet_thread_tstv
{
        public:
                MLHet_thread_tstv( Data data_, double max_gtfreq, vector<double> *initial_params = 0 ) ;
                double ll ; // log likelihood
                vector<double> params ;
		void write( ostream &os ) ; 
		vector<double>* get_priors(  ) ;
        private:
                Data data ;
        
} ;

#endif

