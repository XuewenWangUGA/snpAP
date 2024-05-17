
#include "VCFfac.h"
#include <sstream>
#include <cmath>
#include <limits>
#include "log_add.h"

vector<string> genotypes = { "AA", "CC", "GG", "TT", "AC", "AG", "AT", "CG", "CT", "GT" } ;

void VCFfac::next(  ) 
{

   L:
	// Bail out if end of VCF
	string l ;
	getline( in, l ) ;
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
	vector<long double> log_gt ( genotypes.size(), 0. ) ;

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
			for ( size_t i = 0 ; i < genotypes.size() ; ++i ) {
				if ( strand == '+' ) 
					log_gt[i] += log( 0.5*e.probb( BASE_TO_NUM(genotypes[i][0]), basenum, errprofnum ) +
								0.5*e.probb( BASE_TO_NUM(genotypes[i][1]), basenum, errprofnum ) ) ;
				else 
					log_gt[i] += log( 0.5*e.probb( COMP_NUM(BASE_TO_NUM(genotypes[i][0])), COMP_NUM(basenum), errprofnum ) +
								0.5*e.probb( COMP_NUM(BASE_TO_NUM(genotypes[i][1])), COMP_NUM(basenum), errprofnum ) ) ;
			}
		} else { // add ref bias
			for ( size_t i = 0 ; i < genotypes.size() ; ++i ) {
				double r_1 = 0.5 ;
				double r_2 = 0.5 ;
				if ( ref_ == genotypes[i][0] ) {
					r_1 += ref_bias ; r_2 -= ref_bias ;
				} else if ( ref_ == genotypes[i][1] ) {
					r_2 += ref_bias ; r_1 -= ref_bias ;
				}
				if ( strand == '+' ) 
					log_gt[i] += log( r_1*e.probb( BASE_TO_NUM(genotypes[i][0]), basenum, errprofnum ) +
								r_2*e.probb( BASE_TO_NUM(genotypes[i][1]), basenum, errprofnum ) ) ;
				else 
					log_gt[i] += log( r_1*e.probb( COMP_NUM(BASE_TO_NUM(genotypes[i][0])), COMP_NUM(basenum), errprofnum ) +
								r_2*e.probb( COMP_NUM(BASE_TO_NUM(genotypes[i][1])), COMP_NUM(basenum), errprofnum ) ) ;
			}
		}

	} while ( iss >> rd ) ;

	if ( DP==0 ) goto L ;

	// what is the most likely genotype?
	size_t largest = 0 ;
	for( size_t i = 1 ; i < genotypes.size() ; ++i ) {
		if ( log_gt[largest] < log_gt[i] ) 
			largest = i ;
	}

	vector<long double> log_gt_prior( log_gt ) ;
	for ( size_t i = 0 ; i < genotypes.size() ; ++i ) {
		log_gt_prior[i] += log( priors[i] ) ;
	}

	// what is the most likely genotype when using a prior?
	size_t largest_prior = 0 ;
	for( size_t i = 0 ; i < genotypes.size() ; ++i ) {
		if ( log_gt_prior[largest_prior] < log_gt_prior[i] ) 
			largest_prior = i ;
	}
	
	// calculate QUAL as the phred scaled probability of all other genotypes than the one called
	// sum over all
	long double sum_log_gt = -numeric_limits<long double>::infinity() ;
	for ( size_t i = 0 ; i < log_gt_prior.size() ; ++i ) 
		sum_log_gt = log_add( sum_log_gt, log_gt_prior[i] ) ; 
	
	long double qual = 0 ;
	for ( size_t i = 0 ; i < log_gt_prior.size() ; ++i ) 
		if ( i != largest_prior ) qual += exp( log_gt_prior[i] - sum_log_gt ) ;
	
	QUAL = round( -10. * log10( qual ) ) ;

	size_t min = numeric_limits<size_t>::max() ;
	size_t min2 = numeric_limits<size_t>::max()-1 ;
	for ( size_t i = 0 ; i < genotypes.size() ; ++i ) {
		logPP[i] = log_gt_prior[i] - sum_log_gt ;
		PP[i] = round( -10. * logPP[i]/log(10.) ) ;
		if ( PP[i] < min ) {
			min2 = min ;
			min = PP[i] ;
		} else if ( PP[i] < min2 ) 
			min2 = PP[i] ;
	}
// cerr << "min=" << min << " min2=" << min2 << " min2-min=" << min2-min << endl ;
	GQ = min2-min ;
	

//	vector<long double> log_gt_prior_norm( log_gt_prior ) ;
//	for ( size_t i = 1 ; i < genotypes.size() ; ++i ) 
//		log_gt_prior[i] += log( priors[i] ) ;

	// PP field and GQ field:
	// Posterior probabilities on a phred 
/*	size_t min = numeric_limits<size_t>::max() ;
	size_t min2 = numeric_limits<size_t>::max()-1 ;
	for ( size_t i = 0 ; i < genotypes.size() ; ++i ) {
		PP[i] = round( -10. * (log_gt_prior[i]/log(10.)) ) ;
		if ( PP[i] < min ) {
			min2 = min ;
			min = PP[i] ;
		} else if ( PP[i] < min2 ) 
			min2 = PP[i] ;
	}
	GQ = min2-min ; 
	for ( size_t i = 1 ; i < genotypes.size() ; ++i ) 
		PP[i] -= min ;
*/
	
#ifdef DEBUG2
	cerr << chr_ << "\t" << pos_ << " Basesfwd:" << basesfwd[0] << "," << basesfwd[1] << "," << basesfwd[2] << "," << basesfwd[3] 
			<< " Basesrev:" << basesrev[0] << "," << basesrev[1] << "," << basesrev[2] << "," << basesrev[3] 
			<< " Genotype logl:" ;
	for( size_t i = 0 ; i < genotypes.size() ; ++i ) {
		cerr << " LL(" << genotypes[i] << "|D)=" << log_gt[i] ;
	}
	cerr << "::" ;
	for( size_t i = 0 ; i < genotypes.size() ; ++i ) {
		cerr << " LL(" << genotypes[i] << "|D)=" << log_gt_prior[i] ;
	}
	cerr << endl ;
	for( size_t i = 0 ; i < genotypes.size() ; ++i ) {
		cerr << " log_P(" << genotypes[i] << "|D)=" << logPP[i] << endl ;
	}
	cerr << " -> Genotype:" << genotypes[largest] ;
	cerr << " Genotype_prior:" << genotypes[largest_prior] << endl ;
#endif

	// GT and ALT:
	// check whether this is a homozygous or heterozygous genotype
	if ( genotypes[largest_prior][0] != genotypes[largest_prior][1] ) { // het
		// one of the two alleles matches reference
		if ( genotypes[largest_prior][0] == ref_ ) {
			alt_ = genotypes[largest_prior][1] ;
			GT = "0/1" ;
		} else if ( genotypes[largest_prior][1] == ref_ ) {
			alt_ = genotypes[largest_prior][0] ;
			GT = "0/1" ;
		} else { // neither matches reference
			ostringstream os ;
			os << genotypes[largest_prior][0] ;
			os << "," << genotypes[largest_prior][1] ;
			alt_ = os.str() ;
			GT = "1/2" ;
		}
	} else { // homozygous
		if ( genotypes[largest_prior][0] == ref_ ) { // equal to reference
			GT = "0/0" ;
			alt_ = "." ;
		} else { // different
			GT = "1/1" ;
			alt_ = genotypes[largest_prior][0] ;
		}
	}
}

string VCFfac::format_call(  )
{
	ostringstream ss ;
	ss << "GT:DP:A:C:G:T:PP:GQ\t" ;
	ss << GT << ":" ;
	ss << DP << ":" ;
	for ( size_t i = 0 ; i < 4 ; ++i ) 
		ss << basesfwd[i] << ',' << basesrev[i] << ':' ; 
	if ( pplog ) {
		ss << logPP[0] ;
		for ( size_t i = 1 ; i < logPP.size() ; ++i ) 
			ss << ',' << logPP[i] ;
	} else {
		ss << PP[0] ;
		for ( size_t i = 1 ; i < PP.size() ; ++i ) 
			ss << ',' << PP[i] ;
	}
	ss << ':' << GQ ; 
	return ss.str() ;
}

string VCFfac::header(  )
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
	"##FORMAT=<ID=PP,Number=10,Type=Integer,Description=\"Phred-scaled posteror probability for genotypes AA, CC, GG, TT, AC, AG, AT, CG, CT, GT, in this order\">\n"
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" ) ;
	if ( sample_name ) 
		s += sample_name ;
	else 
		s += "AncientSample" ;
	s += "\n" ;
	return s ;
}
