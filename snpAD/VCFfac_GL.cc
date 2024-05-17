
#include "VCFfac_GL.h"
#include <sstream>
#include <cmath>
#include <limits>
#include "log_add.h"

vector<string> genotypes_GL = { "AA", "CC", "GG", "TT", "AC", "AG", "AT", "CG", "CT", "GT" } ;

size_t gt2pos( char a, char b ) 
{
	if ( a == 'A' ) {
		if ( b == 'A' ) {
			return 0 ;
		} else if ( b == 'C' ) {
			return 4 ;
		} else if ( b == 'G' ) {
			return 5 ;
		} else if ( b == 'T' ) {
			return 6 ;
		}
	} else if ( a == 'C' ) {
		if ( b == 'A' ) {
			return 4 ; // AC
		} else if ( b == 'C' ) {
			return 1 ;
		} else if ( b == 'G' ) {
			return 7 ;
		} else if ( b == 'T' ) {
			return 8 ;
		}

	} else if ( a == 'G' ) {
		if ( b == 'A' ) {
			return 5 ; // AG
		} else if ( b == 'C' ) {
			return 7 ; // CG
		} else if ( b == 'G' ) {
			return 2 ;
		} else if ( b == 'T' ) {
			return 9 ;
		}
	} else if ( a == 'T' ) {
		if ( b == 'A' ) {
			return 6 ; // AT
		} else if ( b == 'C' ) {
			return 8 ; // CT
		} else if ( b == 'G' ) {
			return 9 ; // GT
		} else if ( b == 'T' ) {
			return 3 ;
		}
	}
	return 12 ;
}

void VCFfac_GL::next(  ) 
{

   L:
	// Bail out if end of VCF
	string l ;
	getline( in, l ) ;
//	cerr << "Line: " ;
//	cerr << l << endl ;
	if ( l == "" ) {
		eof = 1 ;
		return ;
	}
	if ( in.eof() ) {
		eof=1 ;
		return ;
	}

	// read chr, pos, reference
	istringstream iss( l ) ;
	iss >> chr_ >> pos_ ;

	if ( iss.eof() ) {
		eof=1 ;
		return ;
	}

	string rd ;
	iss >> rd ;
	ref_ = rd[0] ;
	if ( iss.eof() ) goto L ;
	iss >> rd ;



//	if ( rd == "A" || rd == "C" || rd == "G" || rd == "T" ) {
//		ref_ = rd[0] ;
//		iss >> rd ;
//		if ( rd == "A" || rd == "C" || rd == "G" || rd == "T" ) {
//			iss >> rd ; iss >> rd ; // ignore alt + freq
//		}
//	} else {
//		ref_ = 'N' ;
//		iss >> rd ;
//		ref_ = toupper(fa.fetch( chr_, pos_ )) ;
//	}

	// reset counts 
	fill( basesfwd.begin(), basesfwd.end(), 0 ) ;
	fill( basesrev.begin(), basesrev.end(), 0 ) ;
	DP = 0 ;

	// log GT value
	vector<long double> log_gt ( genotypes_GL.size(), 0. ) ;

	// read bases in seqs individually
	do {
		char base = rd[0] ;
		char strand = rd[2] ;
		size_t errprofnum = atoi( rd.substr( 4,10 ).c_str() ) ;		
		
		DP++ ;

		// count bases
		if ( strand == '+' ) {
			switch( base ) {
				case 'A': basesfwd[0]++ ; break ;
				case 'C': basesfwd[1]++ ; break ;
				case 'G': basesfwd[2]++ ; break ;
				case 'T': basesfwd[3]++ ; break ;
			}
		} else { // rev
			switch( base ) {
				case 'A': basesrev[0]++ ; break ;
				case 'C': basesrev[1]++ ; break ;
				case 'G': basesrev[2]++ ; break ;
				case 'T': basesrev[3]++ ; break ;
			}
		}


		// assign the number of the base
		unsigned char basenum = BASE_TO_NUM(base) ;
		
		// for each read, add log(p)
		if ( ref_bias == 0. ) {
			for ( size_t i = 0 ; i < genotypes_GL.size() ; ++i ) {
				if ( strand == '+' ) 
					log_gt[i] += log( 0.5*e.probb( BASE_TO_NUM(genotypes_GL[i][0]), basenum, errprofnum ) +
								0.5*e.probb( BASE_TO_NUM(genotypes_GL[i][1]), basenum, errprofnum ) ) ;
				else 
					log_gt[i] += log( 0.5*e.probb( COMP_NUM(BASE_TO_NUM(genotypes_GL[i][0])), COMP_NUM(basenum), errprofnum ) +
								0.5*e.probb( COMP_NUM(BASE_TO_NUM(genotypes_GL[i][1])), COMP_NUM(basenum), errprofnum ) ) ;
			}
		} else { // add ref bias
			for ( size_t i = 0 ; i < genotypes_GL.size() ; ++i ) {
				double r_1 = 0.5 ;
				double r_2 = 0.5 ;
				if ( ref_ == genotypes_GL[i][0] ) {
					r_1 += ref_bias ; r_2 -= ref_bias ;
				} else if ( ref_ == genotypes_GL[i][1] ) {
					r_2 += ref_bias ; r_1 -= ref_bias ;
				}
				if ( strand == '+' ) 
					log_gt[i] += log( r_1*e.probb( BASE_TO_NUM(genotypes_GL[i][0]), basenum, errprofnum ) +
								r_2*e.probb( BASE_TO_NUM(genotypes_GL[i][1]), basenum, errprofnum ) ) ;
				else 
					log_gt[i] += log( r_1*e.probb( COMP_NUM(BASE_TO_NUM(genotypes_GL[i][0])), COMP_NUM(basenum), errprofnum ) +
								r_2*e.probb( COMP_NUM(BASE_TO_NUM(genotypes_GL[i][1])), COMP_NUM(basenum), errprofnum ) ) ;
			}
		}

	} while ( iss >> rd ) ;

	if ( DP==0 ) goto L ;

	// what is the most likely genotype?
	size_t largest = 0 ;
	for( size_t i = 1 ; i < genotypes_GL.size() ; ++i ) {
		if ( log_gt[largest] < log_gt[i] ) 
			largest = i ;
	}
	
	vector<long double> log_gt_prior( log_gt ) ;
//	long double sum = 0. ;
	for ( size_t i = 1 ; i < genotypes_GL.size() ; ++i ) {
		log_gt_prior[i] += log( priors[i] ) ;
	}

	// what is the most likely genotype when using a prior?
	size_t largest_prior = 0 ;
	for( size_t i = 1 ; i < genotypes_GL.size() ; ++i ) {
		if ( log_gt_prior[largest_prior] < log_gt_prior[i] ) 
			largest_prior = i ;
	}
	
	// calculate QUAL as the phred scaled probability of all other genotypes_GL than the one called
	// sum over all
	long double sum_log_gt = -numeric_limits<long double>::infinity() ;
	for ( size_t i = 0 ; i < log_gt_prior.size() ; ++i ) 
		sum_log_gt = log_add( sum_log_gt, log_gt_prior[i] ) ; 
	
	long double qual = 0 ;
	for ( size_t i = 0 ; i < log_gt_prior.size() ; ++i ) 
		if ( i != largest_prior ) qual += exp( log_gt_prior[i] - sum_log_gt ) ;
	
	QUAL = round( -10. * log10( qual ) ) ;
	

//	vector<long double> log_gt_prior_norm( log_gt_prior ) ;
//	for ( size_t i = 1 ; i < genotypes_GL.size() ; ++i ) 
//		log_gt_prior[i] += log( priors[i] ) ;

	// PP field and GQ field:
	// Posterior probabilities on a phred 
	size_t min = numeric_limits<size_t>::max() ;
	size_t min2 = numeric_limits<size_t>::max()-1 ;
	for ( size_t i = 0 ; i < genotypes_GL.size() ; ++i ) {
		PP[i] = round( -10. * (log_gt_prior[i]/log(10.)) ) ;
		if ( PP[i] < min ) {
			min2 = min ;
			min = PP[i] ;
		} else if ( PP[i] < min2 ) 
			min2 = PP[i] ;
	}
	GQ = min2-min ; 
	for ( size_t i = 1 ; i < genotypes_GL.size() ; ++i ) 
		PP[i] -= min ;
	
#ifdef DEBUG2
	cerr << chr_ << "\t" << pos_ << " Basesfwd:" << basesfwd[0] << "," << basesfwd[1] << "," << basesfwd[2] << "," << basesfwd[3] 
			<< " Basesrev:" << basesrev[0] << "," << basesrev[1] << "," << basesrev[2] << "," << basesrev[3] 
			<< " Genotype logl:" ;
	for( size_t i = 1 ; i < genotypes_GL.size() ; ++i ) {
		cerr << " LL(" << genotypes_GL[i] << "|D)=" << log_gt[i] ;
	}
	cerr << "::" ;
	for( size_t i = 1 ; i < genotypes_GL.size() ; ++i ) {
		cerr << " LL(" << genotypes_GL[i] << "|D)=" << log_gt_prior[i] ;
	}
	cerr << " -> Genotype:" << genotypes_GL[largest] ;
	cerr << " Genotype_prior:" << genotypes_GL[largest_prior] << endl ;
#endif

	vector<char> order_gtbases ; 
	order_gtbases.push_back( ref_ ) ;

	const char* bases = "ACGT" ;
	// GT and ALT:
	// check whether this is a homozygous or heterozygous genotype
	if ( genotypes_GL[largest_prior][0] != genotypes_GL[largest_prior][1] ) { // het
		// one of the two alleles matches reference
		if ( genotypes_GL[largest_prior][0] == ref_ ) {
			alt_ = genotypes_GL[largest_prior][1] ;
			order_gtbases.push_back( genotypes_GL[largest_prior][1] ) ;
			GT = "0/1" ;
			for ( int i = 0 ; i < 4 ; ++i ) {
				if ( bases[i] != ref_ && bases[i] != genotypes_GL[largest_prior][1] ) { 
					alt_ += ',' ; 
					alt_ += bases[i] ;
					order_gtbases.push_back( bases[i] ) ;
				}
			}
		} else if ( genotypes_GL[largest_prior][1] == ref_ ) {
			alt_ = genotypes_GL[largest_prior][0] ;
			order_gtbases.push_back( genotypes_GL[largest_prior][0] ) ;
			GT = "0/1" ;
			for ( int i = 0 ; i < 4 ; ++i ) {
				if ( bases[i] != ref_ && bases[i] != genotypes_GL[largest_prior][0] ) { 
					alt_ += ',' ; 
					alt_ += bases[i] ;
					order_gtbases.push_back( bases[i] ) ;
				}
			}
		} else { // neither matches reference
			ostringstream os ;
			os << genotypes_GL[largest_prior][0] ;
			os << "," << genotypes_GL[largest_prior][1] ;
			alt_ = os.str() ;
			GT = "1/2" ;
			order_gtbases.push_back( genotypes_GL[largest_prior][0] ) ;
			order_gtbases.push_back( genotypes_GL[largest_prior][1] ) ;
			for ( int i = 0 ; i < 4 ; ++i ) {
				if ( bases[i] != ref_ && bases[i] != genotypes_GL[largest_prior][0] && bases[i] != genotypes_GL[largest_prior][1] ) { 
					alt_ += ',' ; 
					alt_ += bases[i] ; 
					order_gtbases.push_back( bases[i] ) ;
				}
			}
		}
	} else { // homozygous
		if ( genotypes_GL[largest_prior][0] == ref_ ) { // equal to reference
			GT = "0/0" ;
//			alt_ = "." ;
			bool first = 1 ;
			for ( int i = 0 ; i < 4 ; ++i ) {
				if ( bases[i] == ref_ ) continue ;
				if ( first ) { alt_ = bases[i] ; first = 0 ; }
				else { 
					alt_ += ',' ;
					alt_ += bases[i] ;
				} 
				order_gtbases.push_back( bases[i] ) ;
			}
		} else { // different
			GT = "1/1" ;
			alt_ = genotypes_GL[largest_prior][0] ;
			order_gtbases.push_back( genotypes_GL[largest_prior][0] ) ;
			for ( int i = 0 ; i < 4 ; ++i ) {
				if ( bases[i] != ref_ && bases[i] != genotypes_GL[largest_prior][0] ) { 
					alt_ += "," ; 
					alt_ += bases[i] ;  
					order_gtbases.push_back( bases[i] ) ;
				}
			}
		}
	}
//	cerr << "GT Bases are: " << order_gtbases.size() << endl ;
//	for ( size_t i = 0 ; i < order_gtbases.size() ; ++i ) cerr << order_gtbases[i] << endl ;
	for ( size_t i = 0 ; i < 10 ; ++i ) GL[i] = 1. ;
	if ( ref_ == 'A' || ref_ == 'C' || ref_ == 'G' || ref_ == 'T' ) {
		size_t GLpos = 0 ;
		for ( size_t i2 = 0 ; i2 < 4 ; ++i2 ) {
			for ( size_t i1 = 0 ; i1 <= i2 ; ++i1 ) {
				char gt[2] ;
	//			cerr << "Checking: " << i1 << "," << i2 << endl ;
				gt[0] = order_gtbases[i1] ;
				gt[1] = order_gtbases[i2] ;
	//			cerr << "Bases: " << gt[0] << "," << gt[1] << endl ;
	//			cerr << "Pos: " << gt2pos( gt[0], gt[1] ) << endl ;
	//			cerr << "LogGT: " << log_gt[gt2pos( gt[0], gt[1] )] << endl ;
				GL[GLpos] = log_gt[gt2pos( gt[0], gt[1] )]/log(10.) ;
				GLpos++ ;
			}
		}
	}

}

string VCFfac_GL::format_call(  )
{
	ostringstream ss ;
	ss << "GT:DP:A:C:G:T:PP:GQ:GL\t" ;
	ss << GT << ":" ;
	ss << DP << ":" ;
	for ( size_t i = 0 ; i < 4 ; ++i ) 
		ss << basesfwd[i] << ',' << basesrev[i] << ':' ; 
	ss << PP[0] ;
	for ( size_t i = 1 ; i < PP.size() ; ++i ) 
		ss << ',' << PP[i] ;
	ss << ':' << GQ ; 
	ss << fixed ;
	ss.precision(2) ;
	ss << ':' << GL[0] ;
	for ( size_t i = 1 ; i < GL.size() ; ++i ) 
		ss << ',' << GL[i] ;

	return ss.str() ;
}

string VCFfac_GL::header(  )
{
	string s( 
	"##fileformat=VCFv4.1\n"
	"##FORMAT=<ID=A,Number=2,Type=Integer,Description=\"Number of A bases on forward and reverse strand\">\n"
	"##FORMAT=<ID=C,Number=2,Type=Integer,Description=\"Number of C bases on forward and reverse strand\">\n"
	"##FORMAT=<ID=G,Number=2,Type=Integer,Description=\"Number of G bases on forward and reverse strand\">\n"
	"##FORMAT=<ID=T,Number=2,Type=Integer,Description=\"Number of T bases on forward and reverse strand\">\n"
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
	"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
	"##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality (difference between lowest and second lowest PP value)\">\n"
	"##FORMAT=<ID=PP,Number=10,Type=Integer,Description=\"Normalized, Phred-scaled posteror probability for genotypes AA, CC, GG, TT, AC, AG, AT, CG, CT, GT, in this order\">\n"
	"##FORMAT=<ID=GL,Number=10,Type=Float,Description=\"Genotype Likelihoods (log10 in order AA,AB,BB,AC,BC,CC,AD,BD,CD,DD for ref:alt = A:B,C,D)\">\n"
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" ) ;
	if ( sample_name ) 
		s += sample_name ;
	else 
		s += "AncientSample" ;
	s += "\n" ;
	return s ;
}
