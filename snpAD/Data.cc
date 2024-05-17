
#include "Data.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <numeric>
#include <stdexcept>
#include "log_add.h"

// reverse complement
inline char rc( char c ) {
	switch (c) {
		case 'A': return 'T' ;
		case 'T': return 'A' ;
		case 'C': return 'G' ;
		case 'G': return 'C' ;
		default: return 'N' ;
	}
}


Data::Data( const string &filename, Error e_, Error e_c_, unsigned int threads_, unsigned int maxsites ) : e(e_), e_c(e_c_), threads(threads_), basecomp( 4, 0. ), max_ec( 0 ) 
{
	ifstream in( filename ) ;
	string l ;

	if ( !in ) {
		cerr << "Cannot read " << filename << endl ;
		exit( 1 ) ;
	}

	vector<size_t> ACTG_counts( 4, 0 ) ;

L:

	while (!in.eof()) {
		vector<Obs> tmp ;
		getline( in, l ) ;
		istringstream iss( l ) ;
		string c ; int pos ;
		iss >> c >> pos ;

		string rd ;
		unordered_map<string,size_t> m ;
		while ( iss >> rd ) {

			// file with ref or with ref, alt, freq
			if ( rd == "A" || rd == "C" || rd == "G" || rd == "T" ) {
				CoordRefAltFreq p ;
				p.chr = c ;
				p.pos = pos ;
				p.ref = rd[0] ;
				p.alt_f = 0. ; // default: no alternative freq.
				if ( iss >> rd ) {


					// with alt, freq
					if ( rd == "A" || rd == "C" || rd == "G" || rd == "T" ) {
						p.alt = rd[0] ;
						iss >> rd ;
						p.alt_f = atof( rd.c_str() ) ;
						iss >> rd ;
					} 
				} else {
					goto L ;
				}
				// only add if this site has sequences covering it.
				poss.push_back( p ) ;
//				cerr << "Saving to poss: " << p.chr << ":" << p.pos << " ref=" << p.ref << " alt=" << p.alt << " f=" << p.alt_f << endl ;
			} else if ( rd.length() == 1 ) { // anything else thats a single character is treated as "N"
				CoordRefAltFreq p ;
				p.chr = c ;
				p.pos = pos ;
				p.alt_f = 0. ; // default: no alternative freq.
				p.ref = rand()%4 ; // set to random base
				if ( ! (iss >> rd) ) goto L ;
				poss.push_back( p ) ;
			}
			m[rd]++ ;
		}
		
		vector<double> ACTGp( 4, 0. ) ;
		for ( auto it = m.begin() ; it != m.end() ; ++it ) {
//			cerr << "have " << it->first << " with count " << it->second << endl ;
			unsigned char ec = atoi(it->first.substr( 4, 10 ).c_str()) ; // error class
			if ( ec > max_ec ) max_ec = ec ;
			while ( it->second > 64 ) {
				tmp.push_back( Obs( BASE_TO_NUM(it->first[0]), it->first[2]=='+'?1:0, ec, 63 ) ) ;
				if ( it->first[2] == '+' ) {
					ACTGp[0] += log( 64*e.prob( 'A', it->first[0], ec ) ) ;
					ACTGp[1] += log( 64*e.prob( 'C', it->first[0], ec ) ) ;
					ACTGp[2] += log( 64*e.prob( 'T', it->first[0], ec ) ) ;
					ACTGp[3] += log( 64*e.prob( 'G', it->first[0], ec ) ) ;
				} else {
					ACTGp[0] += log( 64*e.prob( rc('A'), rc(it->first[0]), ec ) );
					ACTGp[1] += log( 64*e.prob( rc('C'), rc(it->first[0]), ec ) );
					ACTGp[2] += log( 64*e.prob( rc('T'), rc(it->first[0]), ec ) );
					ACTGp[3] += log( 64*e.prob( rc('G'), rc(it->first[0]), ec ) );
				}
//				cout << "Saved: " << (unsigned int)tmp.back().base << ":" << tmp.back().strand << ":" << (unsigned int)tmp.back().pos << " -- count=" << (unsigned int)tmp.back().count << endl ;
				it->second -= 64 ;
			}
			if ( it->second > 0 ) {
				if ( it->first[2] == '+' ) {
					ACTGp[0] += log( it->second*e.prob( 'A', it->first[0], ec ) );
					ACTGp[1] += log( it->second*e.prob( 'C', it->first[0], ec ) );
					ACTGp[2] += log( it->second*e.prob( 'T', it->first[0], ec ) );
					ACTGp[3] += log( it->second*e.prob( 'G', it->first[0], ec ) );
				} else {
					ACTGp[0] += log( it->second*e.prob( rc('A'), rc(it->first[0]), ec ) );
					ACTGp[1] += log( it->second*e.prob( rc('C'), rc(it->first[0]), ec ) );
					ACTGp[2] += log( it->second*e.prob( rc('T'), rc(it->first[0]), ec ) );
					ACTGp[3] += log( it->second*e.prob( rc('G'), rc(it->first[0]), ec ) );
				}
				tmp.push_back( Obs( BASE_TO_NUM(it->first[0]), it->first[2]=='+'?1:0, ec, it->second-1 ) ) ;
			}
//				cout << "Saved: " << (unsigned int)tmp.back().base << ":" << tmp.back().strand << ":" << (unsigned int)tmp.back().pos << " -- count=" << (unsigned int)tmp.back().count << endl ;
		}
		size_t largest_i = 0 ;
		for ( size_t i = 1 ; i < 4 ; ++i ) 
			if ( ACTGp[i] > ACTGp[largest_i] ) largest_i = i ;

		ACTG_counts[largest_i]++ ;
		
//		cerr << "saving " << c << ":" << pos << " with " << tmp.size() << " elements" << endl ;
		if ( ! tmp.empty() ) d.push_back( tmp ) ;
	}
	size_t sum = accumulate( ACTG_counts.begin(), ACTG_counts.end(), 0 ) ;
	basecomp[0] = static_cast<double>(ACTG_counts[0])/static_cast<double>(sum) ;
	basecomp[1] = static_cast<double>(ACTG_counts[1])/static_cast<double>(sum) ;
	basecomp[2] = static_cast<double>(ACTG_counts[2])/static_cast<double>(sum) ;
	basecomp[3] = static_cast<double>(ACTG_counts[3])/static_cast<double>(sum) ;

	cerr << "Basecomp: " << basecomp[0] << " " << basecomp[1] << " " << basecomp[2] << " " << basecomp[3] << " " << endl ;
	cerr << "total sites: " << d.size() << " total site information: " << poss.size() << endl ;
}

vector<string> Data_genotypes = { "AA", "CC", "TT", "GG", "AC", "AG", "AT", "CG", "CT", "GT" } ;

void Data::reestimate_error( vector<double>* priors ) 
{

	vector<size_t> e_counts( 256*16, 0 ) ;
	

	// for all sites
	for( auto it = d.begin() ; it != d.end() ; it++ ) {
		// call a genotype


		// log GT value
		vector<long double> log_gt ( Data_genotypes.size(), 0. ) ;

		// all reads at this site
		for ( auto it2 = it->begin() ; it2 != it->end() ; it2++ ) {
			for ( size_t i = 0 ; i < Data_genotypes.size() ; ++i ) {
				if ( it2->strand == 1 ) // plus
					log_gt[i] += (1+it2->count)*log( 0.5*e.probb( BASE_TO_NUM(Data_genotypes[i][0]), it2->base, it2->pos ) +
								0.5*e.probb( BASE_TO_NUM(Data_genotypes[i][1]), it2->base, it2->pos ) ) ;
				else // minus
					log_gt[i] += (1+it2->count)*log( 0.5*e.probb( COMP_NUM(BASE_TO_NUM(Data_genotypes[i][0])), COMP_NUM(it2->base), it2->pos ) +
								0.5*e.probb( COMP_NUM(BASE_TO_NUM(Data_genotypes[i][1])), COMP_NUM(it2->base), it2->pos ) ) ;
			}
		}
		// gt likelihood
		vector<long double> log_gt_prior( log_gt ) ;
		for ( size_t i = 0 ; i < Data_genotypes.size() ; ++i ) {
			log_gt_prior[i] += log( priors->at(i) ) ;
		}

		// what is the most likely genotype when using a prior?
		size_t largest_prior = 0 ;
		for( size_t i = 1 ; i < Data_genotypes.size() ; ++i ) {
			if ( log_gt_prior[largest_prior] < log_gt_prior[i] ) 
				largest_prior = i ;
		}
		
		// we only calculate error profiles from homozygous genotypes with at least 3x coverage
		if ( it->size() > 2 && Data_genotypes[largest_prior][0] == Data_genotypes[largest_prior][1] ) {
			for ( auto it2 = it->begin() ; it2 != it->end() ; it2++ ) {
				if ( it2->strand == 1 ) // plus
					e_counts[ BASE_TO_NUM(Data_genotypes[largest_prior][0]) + it2->base*4 + it2->pos*16 ] += it2->count+1 ;
				else // minus
					e_counts[ COMP_NUM(BASE_TO_NUM(Data_genotypes[largest_prior][0])) + COMP_NUM(it2->base)*4 + it2->pos*16 ] += it2->count+1;
			}
		}
	}

	vector<double> e_( 256*16, 0/.0 /*NaN*/ ) ;

	// initialize those up to max_ec properly
	for ( size_t pos = 0 ; pos <= max_ec ; ++pos ) {
		for ( size_t from = 0 ; from < 4 ; ++from ) {
			for ( size_t to = 0 ; to < 4 ; ++to ) {
				e_[ from + to*4 + pos*16 ] = from==to ? 1.-e.default_error_rate : e.default_error_rate/3. ; 
			}
		}
	}

	for ( size_t pos = 0 ; pos < 256 ; ++pos ) {
		for ( size_t from = 0 ; from < 4 ; ++from ) {
			// sum over all /to/ bases
			double sum = 0. ;
			for ( size_t to = 0 ; to < 4 ; ++to ) 
				sum += e_counts[ from + to*4 + pos*16 ] ;
			if ( sum == 0. ) continue ;
			size_t count0 = 0 ;
			for ( size_t to = 0 ; to < 4 ; ++to ) {
				if ( e_counts[ from + to*4 + pos*16 ] > 0. ) 
					e_[ from + to*4 + pos*16 ] = e_counts[ from + to*4 + pos*16 ] / sum ;
				else { // minimum probability is 1/10G 
					e_[ from + to*4 + pos*16 ] = 1e-10 ;
					count0++ ;
				}
			}
			// adjust the non-error probability for the adjustments to fix 0 counts
			if ( count0 > 0 ) e_[ from + from*4 + pos*16 ] -= count0 * 1e-10 ;

		}
	}

	e = Error( e_ ) ;

	e.print( cerr ) ;

	reestimate_basecomp(  ) ;
}

void Data::reestimate_error_full( vector<double>* priors ) 
{

	vector<long double> e_counts( 256*16, -numeric_limits<long double>::infinity() ) ;
	

	// for all sites
	for( auto it = d.begin() ; it != d.end() ; it++ ) {
		// call a genotype


		// log GT value
		vector<long double> log_gt ( Data_genotypes.size(), 0. ) ;

		// all reads at this site
		for ( auto it2 = it->begin() ; it2 != it->end() ; it2++ ) {
			for ( size_t i = 0 ; i < Data_genotypes.size() ; ++i ) {
				if ( it2->strand == 1 ) // plus
					log_gt[i] += (1+it2->count)*log( 0.5*e.probb( BASE_TO_NUM(Data_genotypes[i][0]), it2->base, it2->pos ) +
								0.5*e.probb( BASE_TO_NUM(Data_genotypes[i][1]), it2->base, it2->pos ) ) ;
				else // minus
					log_gt[i] += (1+it2->count)*log( 0.5*e.probb( COMP_NUM(BASE_TO_NUM(Data_genotypes[i][0])), COMP_NUM(it2->base), it2->pos ) +
								0.5*e.probb( COMP_NUM(BASE_TO_NUM(Data_genotypes[i][1])), COMP_NUM(it2->base), it2->pos ) ) ;
			}
		}
		// gt likelihood
		vector<long double> log_gt_prior( log_gt ) ;
		long double sum_log_gt = -numeric_limits<long double>::infinity() ;
		for ( size_t i = 0 ; i < Data_genotypes.size() ; ++i ) {
			log_gt_prior[i] += log( priors->at(i) ) ;
			sum_log_gt = log_add( sum_log_gt, log_gt_prior[i] ) ; // max(sum_log_gt, log_gt_prior[i]) + log1p( exp( -fabs( sum_log_gt - log_gt_prior[i] ) ) ) ;
		}

		// At least 3x cov
		if ( it->size() <= 2 ) continue ;

		// gt likelihood for each GT and prob of subst
		for ( size_t i = 0 ; i < Data_genotypes.size() ; ++i ) {
			long double log_gt_p = log_gt_prior[i] - sum_log_gt ;

			for ( auto it2 = it->begin() ; it2 != it->end() ; it2++ ) {
				// use the matching base of the GT or a random GT if neither matches
				char base = BASE_TO_NUM( Data_genotypes[i][0] ) ;
				if ( base != it2->base ) base = BASE_TO_NUM( Data_genotypes[i][1] ); 
				if ( base != it2->base ) base = BASE_TO_NUM( Data_genotypes[i][rand()%2] ) ; 
				if ( it2->strand == 1 ) { // plus
					e_counts[ base + it2->base*4 + it2->pos*16 ] = log_add( e_counts[ base + it2->base*4 + it2->pos*16 ], log_gt_p+log(static_cast<long double>(it2->count+1)) ) ;
				} else { // minus
					e_counts[ COMP_NUM(base) + COMP_NUM(it2->base)*4 + it2->pos*16 ] = log_add( e_counts[ COMP_NUM(base) + COMP_NUM(it2->base)*4 + it2->pos*16 ], log_gt_p+log(static_cast<long double>(it2->count+1)) ) ;
				}
			}
			
		}

	}

	vector<double> e_( 256*16, 0/.0 /*NaN*/ ) ;

	// initialize those up to max_ec properly
	for ( size_t pos = 0 ; pos <= max_ec ; ++pos ) {
		for ( size_t from = 0 ; from < 4 ; ++from ) {
			for ( size_t to = 0 ; to < 4 ; ++to ) {
				e_[ from + to*4 + pos*16 ] = from==to ? 1.-e.default_error_rate : e.default_error_rate/3. ; 
			}
		}
	}


	for ( size_t pos = 0 ; pos < 256 ; ++pos ) {
		for ( size_t from = 0 ; from < 4 ; ++from ) {
			// sum over all /to/ bases
			long double sum = -numeric_limits<long double>::infinity() ;
			for ( size_t to = 0 ; to < 4 ; ++to ) 
				sum = log_add( sum, e_counts[ from + to*4 + pos*16 ] );
			if ( sum == -numeric_limits<long double>::infinity() ) continue ;
			for ( size_t to = 0 ; to < 4 ; ++to ) {
					e_[ from + to*4 + pos*16 ] = exp( e_counts[ from + to*4 + pos*16 ] - sum ) ;
			}
		}
	}

	e = Error( e_ ) ;

	e.print( cerr ) ;

	reestimate_basecomp(  ) ;
}

void Data::reestimate_basecomp(  ) 
{
	// for all sites
	vector<size_t> ACTG_counts( 4, 0 ) ;
	for( auto it = d.begin() ; it != d.end() ; it++ ) {
		vector<double> ACTGp( 4, 0. ) ;
		// all reads
		for ( auto it2 = it->begin() ; it2 != it->end() ; it2++ ) {
			if ( it2->strand == 1 ) { // plus
				ACTGp[BASE_A] += log( (it2->count+1)*e.probb( BASE_A, it2->base, it2->pos ) ) ;
				ACTGp[BASE_C] += log( (it2->count+1)*e.probb( BASE_C, it2->base, it2->pos ) ) ;
				ACTGp[BASE_T] += log( (it2->count+1)*e.probb( BASE_T, it2->base, it2->pos ) ) ;
				ACTGp[BASE_G] += log( (it2->count+1)*e.probb( BASE_G, it2->base, it2->pos ) ) ;
			} else { // minus
				ACTGp[BASE_A] += log( (it2->count+1)*e.probb( COMP_NUM(BASE_A), COMP_NUM(it2->base), it2->pos ) ) ;
				ACTGp[BASE_C] += log( (it2->count+1)*e.probb( COMP_NUM(BASE_C), COMP_NUM(it2->base), it2->pos ) ) ;
				ACTGp[BASE_T] += log( (it2->count+1)*e.probb( COMP_NUM(BASE_T), COMP_NUM(it2->base), it2->pos ) ) ;
				ACTGp[BASE_G] += log( (it2->count+1)*e.probb( COMP_NUM(BASE_G), COMP_NUM(it2->base), it2->pos ) ) ;
			}

		}
		size_t largest_i = 0 ;
		for ( size_t i = 1 ; i < 4 ; ++i ) 
			if ( ACTGp[i] > ACTGp[largest_i] ) largest_i = i ;


		ACTG_counts[largest_i]++ ;
	}

	size_t sum = accumulate( ACTG_counts.begin(), ACTG_counts.end(), 0 ) ;
	basecomp[0] = static_cast<double>(ACTG_counts[0])/static_cast<double>(sum) ;
	basecomp[1] = static_cast<double>(ACTG_counts[1])/static_cast<double>(sum) ;
	basecomp[2] = static_cast<double>(ACTG_counts[2])/static_cast<double>(sum) ;
	basecomp[3] = static_cast<double>(ACTG_counts[3])/static_cast<double>(sum) ;
	cerr << "Re-estimated basecomp: " << basecomp[0] << " " << basecomp[1] << " " << basecomp[2] << " " << basecomp[3] << " " << endl ;
}


