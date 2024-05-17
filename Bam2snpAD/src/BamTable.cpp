#include <popt.h>
#include <limits.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <list>
#include <unordered_set>
#include <sstream>
using namespace std;

//#define DEBUG2;
//#define DEBUG;


#include "api/BamReader.h"
#include "api/BamWriter.h"
using namespace BamTools;

#include "BamAliIt.h"
#include "FastaSlurp.h"
// #include "SeqEditCounts.h"

// forward declarations

// parses region string of the format chr:start-end (or parts thereof) and sets the BamReader to that region
void parse_and_set_region( string region, BamReader &reader, int &left_bound, int &right_bound, string &refname ) ;

// prints lines for a specific region 
void process_region( BamReader &reader, unordered_set<string> &readgroups, string &refname, int left_bound, int right_bound, FastaSlurp &fasta, int offset_profilenum, 
		int mapqual_cutoff, int basequal_cutoff, int basequal_cutoff_max, 
		int ends_size, int merge_ends, int nopaired ) ;

// prints lines for a specific regions (faster version for arrays where we walk through one full chr for all regions on it)
void process_region_vec( BamReader &reader, unordered_set<string> &readgroups, string &refname, vector<pair<size_t,size_t>> &regions, FastaSlurp &fasta, int offset_profilenum, 
		int mapqual_cutoff, int basequal_cutoff, int basequal_cutoff_max, 
		int ends_size, int merge_ends, int nopaired ) ;
//
int main(int argc, const char* argv[]) {

	/////
	// parse options

	int mapqual_cutoff = 25 ; // mapping quality cutoff
	int basequal_cutoff = 30 ; // base quality min
	int basequal_cutoff_max = 100 ; // base quality max
	int ends_size = 15 ;
//	int snip_left = 0;
//	int snip_right = 0;

	char* index_fn = 0; // .bai file (if given)
	char* region = 0; // region given as "chr:start-end", "chr", or "chr:start-" 
	char* s_regionfile = 0 ;
	int regionfilebed = 0 ;
	char* fasta_file = 0 ;
	char* s_rgs = 0 ;
	int offset_profilenum = 0 ;

	int array_flag = 0 ;

	int nopaired = 0 ;
	int merge_ends = 0 ;

	int flag_print_version = 0 ;

	// popt parameters
	struct poptOption optionsTable[] = {
			{ "bam_index",				'i', POPT_ARG_STRING, &index_fn, 0, "BAM index FILE", "FILE" },
			{ "region",				'r', POPT_ARG_STRING, &region, 0,"restrict to CHR:START-END. Alternatively you may restrict to a specific chr by giving CHR or the remainder of a chromosome with CHR:START-", "CHR:START-END" },
			{ "region_file",		'F', POPT_ARG_STRING, &s_regionfile, 0, "file with regions to be parsed (format: chr, start, end  in tab separated format, one per line, 1-based closed interval)" },
			{ "region_file_is_bedformat",   'b', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT, &regionfilebed, 0, "The region file is in bed format (i.e. coordinates are 0-based open interval not 1-based closed interval)" }, 
                        { "region_file_is_array",       'A', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT, &array_flag, 0, "The regionfile contains lots of small regions that are sorted by start position within each chromosome. Setting this flag usually speeds up the process for these files." },
			{ "readgroups",			'R', POPT_ARG_STRING, &s_rgs, 0, "read groups, comma separated", 0 },
			{ "fasta", 				'f', POPT_ARG_STRING, &fasta_file, 0, "Read consensus or reference sequence from this fasta file", "FASTA" },
			{ "map_qual",				'Q', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &mapqual_cutoff, 0,"mapping quality >= MQ", "MQ" },
			{ "min_qual",			'q', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &basequal_cutoff, 0, "base quality >= Q", "Q" },
			{ "max_qual",			0, POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &basequal_cutoff_max, 0, "base quality <= Q", "Q" },
			{ "size",			's', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &ends_size, 0, "label N bases from either end of the sequence", "N" },
			{ "offset",			'o', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &offset_profilenum, 0, "Add N to the position-profile number", "N" },
//			{ "by_strand",				'S', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT, &by_strand, 0, "give separate columns for + (first four) and - strand (last four) molecules; uses only properly-paired and single reads", 0 },
//			{ "output_bed",				'B', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT, &bedformat, 0, "Print output in bed format", 0 },
			{ "nopaired",				'P', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT, &nopaired, 0, "Ignore paired reads", 0 },
			{ "merge_ends",				0, POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT, &merge_ends, 0, "Do not separate 3' and 5' ends for position-profiles (same distance to end == same profile)", 0 },
			// not urgently needed:
//			{ "snip_left",				0  , POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &snip_left, 0, "ignore n bases at the 5' end of the read", "n" },
//			{ "snip_right",				0  , POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,&snip_right, 0, "ignore n bases at the 3' end of the read", "n" }, 
			{ "version", 0, POPT_ARG_NONE, &flag_print_version, 0, "Print Version", 0 },
			POPT_AUTOHELP { NULL, 0, 0, NULL, 0 }
	};


	// parse options; do some sanity checking
	poptContext optCon ;
	optCon = poptGetContext( "Bam2snpAD", argc, argv, optionsTable, 0 ) ;
	poptSetOtherOptionHelp( optCon, "[OPTIONS]* <bamfile>" ) ;
	int rc = poptGetNextOpt( optCon ) ;
	if (rc != -1) {
		poptPrintUsage( optCon, stderr, 0 ) ;
		exit( 1 ) ;
	}

	#include "../../version.h" 
	if ( flag_print_version ) {
		cout << "snpAD " << SNPAD_VERSION << endl ;
		exit( 0 ) ;
	}

	const char* fn = poptGetArg( optCon ) ;
	if ( !fn ) {
		cerr << "Error: need bam file as argument." << endl ;
		poptPrintUsage( optCon, stderr, 0 ) ;
		exit( 1 ) ;
	}

	if ( !region && !s_regionfile) {
		cerr << "Error: Need either region (-r option) or region file (-f option)" << endl ;
		poptPrintUsage( optCon, stderr, 0 ) ;
		exit( 1 ) ;
	} 

	if ( ! fasta_file ) {
		cerr << "Error: need fasta file as argument." << endl ; 
		poptPrintUsage( optCon, stderr, 0 ) ;
		exit( 1 ) ;
	}

	// populate readgroup unordered set
	unordered_set<string> readgroups ;
	if ( s_rgs ) {
		istringstream iss( s_rgs ) ;
		string s ;
		while ( getline( iss, s, ',' ) ) {
			readgroups.insert( s ) ;
		}
	}

	// ios_base::sync_with_stdio( false ) ;


	////// 
	// Open Fasta file and read into mem
	ifstream fa_stream(fasta_file) ;
 	FastaSlurp fasta( fa_stream ) ;
//	cerr << "fetching chr1, 1000000: " << fasta.fetch( "1", 1000000 ) << endl ;


	///////
	// open bam file and seek

	BamReader reader ;
	reader.Open( string(fn) ) ;

	if (index_fn) {
		if ( !reader.OpenIndex( string(index_fn) ) ) {
			cerr << "Error: Not a valid index: " << index_fn << endl ;
			exit( 1 ) ;
		}
	} else {
		if ( !reader.LocateIndex() ) {
			cerr << "Error: Cannot find index. Either run `samtools index' on the bam file or supply the correct index name with -i" << endl ;
			exit( 1 ) ;
		}
	}

	// refname_id is the chromosome id
	int left_bound = 0, right_bound = 0 ; // start and end 1-based interval
	string refname ;

	ios_base::sync_with_stdio( false ) ;

	if ( region ) {
		parse_and_set_region( region, reader, left_bound, right_bound, refname ) ;
		process_region( reader, readgroups, refname, left_bound, right_bound, fasta, offset_profilenum, 
				mapqual_cutoff, basequal_cutoff, basequal_cutoff_max, 
				ends_size, merge_ends, nopaired ) ;
	} else if ( s_regionfile && array_flag ) {
	                // Read the regions
                ifstream in( s_regionfile ) ;
                if ( !in.is_open() ) {
                        cerr << "Error: cannot open region file '" << s_regionfile << "'." << endl ;
                        exit( 1 ) ;
                }
                vector<vector<pair<size_t,size_t> > > coords ; // coordinates, by chr, then vector of (sorted) coordinates
                vector<string> chrs ; // chromosomes in the order they occured. 
                vector<int> chrs_ids ; // chromosomes in the order they occured. 
                vector<size_t> chrs_sizes ;
                string line ;
                int last = -1 ;
                while ( getline( in, line ) ) {
                        istringstream ss(line) ;
                        string refname ; ss >> refname ;
                        int refname_id = reader.GetReferenceID( refname ) ;
                        if ( refname_id == -1 ) {
                                cerr << "Error: reference `" << refname << "' in region file not found." << endl;
                                exit( 1 );
                        }
                        int left_bound ; if ( !( ss >> left_bound ) ) {
                                cerr << "Cannot parse line '" << line << "' from region file." << endl ;
                                exit( 1 ) ;
                        }
                        if ( regionfilebed ) left_bound++ ;
                        int right_bound ; if ( !( ss >> right_bound ) ) {
                                cerr << "Cannot parse line '" << line << "' from region file." << endl ;
                                exit( 1 ) ;
                        }
                        if ( left_bound > right_bound ) {
                                cerr << "Start coordinate larger than end at line '" << line << "'." << endl ;
                                exit( 1 ) ;
                        }
                        if ( chrs.empty() || chrs.back() != refname ) {
                                chrs.push_back( refname ) ;
                                chrs_ids.push_back( refname_id ) ;
                                chrs_sizes.push_back( reader.GetReferenceData()[refname_id].RefLength ) ;
                                coords.push_back( vector<pair<size_t,size_t>>() ) ;
                                last = -1 ;
                        }
                        if ( left_bound < last ) {
                                cerr << "Regionfile not sorted. Either omit -A or sort file." << endl ;
                                exit( 1 ) ;
                        }
                        last = left_bound ;
                        coords.back().push_back( pair<size_t,size_t>(left_bound,right_bound) ) ;
                }
                for ( size_t i = 0 ; i < chrs.size() ; ++i ) {
                        // set chr
                        if ( !reader.SetRegion( chrs_ids[i], 0, chrs_ids[i], chrs_sizes[i] ) ) {
                                cerr << "Error: Cannot set bounds to " << reader.GetReferenceData()[chrs_ids[i]].RefName << ":" << left_bound << "-" << right_bound << endl;
                                exit( 1 ) ;
                        }
			process_region_vec( reader, readgroups, chrs[i], coords[i], fasta, offset_profilenum, 
					mapqual_cutoff, basequal_cutoff, basequal_cutoff_max, 
					ends_size, merge_ends, nopaired ) ;
                }

	} else if ( s_regionfile && !array_flag ) {
		ifstream in( s_regionfile ) ;
		if ( !in.is_open() ) {
			cerr << "Error: cannot open region file '" << s_regionfile << "'." << endl ;
			exit( 1 ) ;
		}
		string line ;
		while ( getline( in, line ) ) {
			istringstream ss(line) ;
			string refname ; ss >> refname ;
			int refname_id = reader.GetReferenceID( refname ) ;
			if ( refname_id == -1 ) {
				cerr << "Error: reference `" << refname << "' in region file not found." << endl;
				exit( 1 );
			}
			int left_bound ; if ( !( ss >> left_bound ) ) { 
				cerr << "Cannot parse line '" << line << "' from region file." << endl ;
				exit( 1 ) ;
			}
			if ( regionfilebed ) left_bound++ ;
			int right_bound ; if ( !( ss >> right_bound ) ) { 
				cerr << "Cannot parse line '" << line << "' from region file." << endl ;
				exit( 1 ) ;
			}
			if ( left_bound > right_bound ) {
				cerr << "Start coordinate larger than end at line '" << line << "'." << endl ;
				exit( 1 ) ;
			}
			if ( !reader.SetRegion( refname_id, left_bound-1, refname_id, right_bound ) ) {
				cerr << "Error: Cannot set bounds to " << reader.GetReferenceData()[refname_id].RefName << ":" << left_bound << "-" << right_bound << endl;
				exit( 1 ) ;
			}
			process_region( reader, readgroups, refname, left_bound, right_bound, fasta, offset_profilenum, 
					mapqual_cutoff, basequal_cutoff, basequal_cutoff_max, 
					ends_size, merge_ends, nopaired ) ;
		}
	}

}

void process_region( BamReader &reader, unordered_set<string> &readgroups, string &refname, int left_bound, int right_bound, FastaSlurp &fasta, int offset_profilenum, 
		int mapqual_cutoff, int basequal_cutoff, int basequal_cutoff_max, 
		int ends_size, int merge_ends, int nopaired ) 
{


	///// List of Alignments to process
	list<BamAlignmentIterator> ali_list ;
	///// current Alignment under consideration
	BamAlignment ali ;

	// get the first valid alignment ...
	bool validalis=false ;
	while ( reader.GetNextAlignmentCore( ali ) ) {
		// Need to check for correct readgroup?
		if ( !readgroups.empty() ) { // need to check read group
			ali.BuildCharData(  ) ;
			string s ;
			ali.GetTag( "RG", s ) ;
			if ( readgroups.count( s ) == 0 ) continue ; // no good, try next
		}
		if ( !ali.IsFailedQC() && !ali.IsDuplicate() && ali.IsMapped() && ali.MapQuality >= mapqual_cutoff && (!nopaired || !ali.IsPaired()) ) {
			validalis=true ;
			break ;
		}
	}
	// Not a single valid alignment. Just bail.
	if ( !validalis ) {
		return ; 
	}
	int pos = ali.Position ;
	if ( readgroups.empty() ) ali.BuildCharData(  ) ;
	// ... and add to list of alignments to check
	ali_list.push_back( BamAlignmentIterator(ali) ) ;


	//////////
	// read through file until we hit the first aligment that starts outside of the interval or we have no more reads
	//
	while ( reader.GetNextAlignmentCore( ali ) ) {
		if ( ali.IsFailedQC() ) continue ;
		if ( !ali.IsMapped() ) continue ;
		if ( ali.IsDuplicate() ) continue ;
		if ( ali.MapQuality < mapqual_cutoff ) continue ; 
		if ( nopaired && ali.IsPaired() ) continue ;
		ali.BuildCharData(  ) ; 
		if ( !readgroups.empty() ) { // need to check read group
			string s ;
			ali.GetTag( "RG", s ) ;
			if ( readgroups.count( s ) == 0 ) continue ;
		}

		// alignments at the same position are added to the list (for processing this round)
		if ( pos == ali.Position ) {
			ali_list.push_back( BamAlignmentIterator(ali) ) ;
			continue ;
		}
		// alignments with start outside of our desired range do not need to be processed (see below for rest of printing...)
		if ( pos > right_bound ) break ; 
		
//		cerr << " [ " << pos << " - " << ali.Position-1 << " ]" <<  endl ;
		
		// go through and process all alignments up to ali.Position
		// (since the file is sorted, we know that we have all overlapping alignments up to that point)
		for ( int i = pos ; i < ali.Position ; i++ ) {
			if ( i+1 < left_bound ) continue ; // +1 because BamTools is 0-based for the start coordinate (meh.)
			if ( i+1 > right_bound ) break ;


//			if ( i % 100000 == 0 ) cerr << i << '\r' ;

			bool header = 0 ; // "header" for the line (chr + pos) not yet printed

			// go through ali_list, count 
			list<BamAlignmentIterator>::iterator it = ali_list.begin() ;
			while ( it != ali_list.end() ) {
				it->go( i ) ;
				
				if ( it->gimme_qual() >= basequal_cutoff && it->gimme_qual() <= basequal_cutoff_max ) {
					char c = toupper( it->gimme_base(  ) );
					if ( c == 'A' || c == 'C' || c == 'G' || c == 'T' ) {
						if ( ! header ) { 
							header = 1 ;
							char ref = toupper( fasta.fetch( refname, i+1 ) ) ;
							cout << refname << "\t" << i+1 << "\t" << ref ;
						}
						int dist_5 = it->gimme_dist5( ) ;
						int dist_3 = it->gimme_dist3( ) ;
						if ( ! merge_ends ) { // default
							if ( it->orientation == 1 ) {
								if ( dist_5 >= 0 && dist_5 < ends_size ) 
									cout << "\t" << c << ":+:" << dist_5 + offset_profilenum ; 
								else if ( dist_3 >= 0 && dist_3 < ends_size ) 
									cout << "\t" << c << ":+:" << ends_size+dist_3 + offset_profilenum ;
								else 
									cout << "\t" << c << ":+:" << ends_size*2 + offset_profilenum ;
							} else if ( it->orientation == -1 ) {
								if ( dist_5 >= 0 && dist_5 < ends_size ) 
									cout << "\t" << c << ":-:" << dist_5 + offset_profilenum ;
								else if ( dist_3 >= 0 && dist_3 < ends_size ) 
									cout << "\t" << c << ":-:" << ends_size+dist_3 + offset_profilenum ;
								else 
									cout << "\t" << c << ":-:" << ends_size*2 + offset_profilenum ;
							}
						} else { // merge end
							if ( it->orientation == 1 ) {
								if ( dist_5 >= 0 && dist_5 < ends_size ) 
									cout << "\t" << c << ":+:" << dist_5 + offset_profilenum ; 
								else if ( dist_3 >= 0 && dist_3 < ends_size ) 
									cout << "\t" << c << ":+:" << dist_3 + offset_profilenum ;
								else 
									cout << "\t" << c << ":+:" << ends_size + offset_profilenum ;
							} else if ( it->orientation == -1 ) {
								if ( dist_5 >= 0 && dist_5 < ends_size ) 
									cout << "\t" << c << ":-:" << dist_5 + offset_profilenum ;
								else if ( dist_3 >= 0 && dist_3 < ends_size ) 
									cout << "\t" << c << ":-:" << dist_3 + offset_profilenum ;
								else 
									cout << "\t" << c << ":-:" << ends_size + offset_profilenum ;
							}
						}
					}
				}

				if ( it->end() )  
					it = ali_list.erase( it ) ;
				else 
					++it ;
			}
			
			if ( header ) {
				cout << endl ;
			}

			if ( ali_list.empty() ) break ;

		}
		// add alignment that ended at pos+1 for next round.
		ali_list.push_back( BamAlignmentIterator(ali) ) ;
		pos = ali.Position ;
	}


	// keep going until our list is empty or end is reached 
	while ( !ali_list.empty() && pos < right_bound ) {

		if ( pos+1 < left_bound ) {
			pos++ ;
			continue ;
		}


		// go through ali_list, count 
		list<BamAlignmentIterator>::iterator it = ali_list.begin() ;
		bool header = 0 ; // "header" for the line (chr + pos) not yet printed
		while ( it != ali_list.end() ) {

			it->go( pos ) ;
			
			
			if ( it->gimme_qual() >= basequal_cutoff && it->gimme_qual() <= basequal_cutoff_max ) {
				char c = toupper( it->gimme_base(  ) );
				if ( c == 'A' || c == 'C' || c == 'G' || c == 'T' ) {
					if ( ! header ) { 
						header = 1 ;
						char ref = toupper( fasta.fetch( refname, pos+1 ) ) ;
						cout << refname << "\t" << pos+1 << "\t" << ref ;
					}
					int dist_5 = it->gimme_dist5( ) ;
					int dist_3 = it->gimme_dist3( ) ;
					if ( ! merge_ends ) { // default
						if ( it->orientation == 1 ) {
							if ( dist_5 >= 0 && dist_5 < ends_size ) 
								cout << "\t" << c << ":+:" << dist_5 + offset_profilenum ;
							else if ( dist_3 >= 0 && dist_3 < ends_size ) 
								cout << "\t" << c << ":+:" << ends_size+dist_3 + offset_profilenum ;
							else 
								cout << "\t" << c << ":+:" << ends_size*2 + offset_profilenum ;
						} else if ( it->orientation == -1 ) {
							if ( dist_5 >= 0 && dist_5 < ends_size ) 
								cout << "\t" << c << ":-:" << dist_5 + offset_profilenum ;
							else if ( dist_3 >= 0 && dist_3 < ends_size ) 
								cout << "\t" << c << ":-:" << ends_size+dist_3 + offset_profilenum ;
							else 
								cout << "\t" << c << ":-:" << ends_size*2 + offset_profilenum ;
						}
					} else {
						if ( it->orientation == 1 ) {
							if ( dist_5 >= 0 && dist_5 < ends_size ) 
								cout << "\t" << c << ":+:" << dist_5 + offset_profilenum ;
							else if ( dist_3 >= 0 && dist_3 < ends_size ) 
								cout << "\t" << c << ":+:" << dist_3 + offset_profilenum ;
							else 
								cout << "\t" << c << ":+:" << ends_size + offset_profilenum ;
						} else if ( it->orientation == -1 ) {
							if ( dist_5 >= 0 && dist_5 < ends_size ) 
								cout << "\t" << c << ":-:" << dist_5 + offset_profilenum ;
							else if ( dist_3 >= 0 && dist_3 < ends_size ) 
								cout << "\t" << c << ":-:" << dist_3 + offset_profilenum ;
							else 
								cout << "\t" << c << ":-:" << ends_size + offset_profilenum ;
						}
					}
				}
			}
			if ( it->end() )  
				it = ali_list.erase( it ) ;
			else 
				++it ;
		}
		if ( header ) {
			cout << endl ;
		}
		pos++ ;		
	}
}

void process_region_vec( BamReader &reader, unordered_set<string> &readgroups, string &refname, vector<pair<size_t,size_t>> &regions, FastaSlurp &fasta, int offset_profilenum, 
		int mapqual_cutoff, int basequal_cutoff, int basequal_cutoff_max, 
		int ends_size, int merge_ends, int nopaired ) 
{


	///// List of Alignments to process
	list<BamAlignmentIterator> ali_list ;
	///// current Alignment under consideration
	BamAlignment ali ;


	auto region_it = regions.begin() ;

	// get the first valid alignment ...
	bool validalis=false ;
	while ( reader.GetNextAlignmentCore( ali ) ) {
		// Need to check for correct readgroup?
		if ( !readgroups.empty() ) { // need to check read group
			ali.BuildCharData(  ) ;
			string s ;
			ali.GetTag( "RG", s ) ;
			if ( readgroups.count( s ) == 0 ) continue ; // no good, try next
		}
		if ( !ali.IsFailedQC() && !ali.IsDuplicate() && ali.IsMapped() && ali.MapQuality >= mapqual_cutoff && (!nopaired || !ali.IsPaired()) ) {
			validalis=true ;
			break ;
		}
	}
	// Not a single valid alignment. Just bail.
	if ( !validalis || region_it == regions.end() ) {
		return ; 
	}
	size_t pos = ali.Position ;
	if ( readgroups.empty() ) ali.BuildCharData(  ) ;
	// ... and add to list of alignments to check
	ali_list.push_back( BamAlignmentIterator(ali) ) ;


	//////////
	// read through file until we hit the first aligment that starts outside of the interval or we have no more reads
	//
	while ( reader.GetNextAlignmentCore( ali ) ) {
		if ( ali.IsFailedQC() ) continue ;
		if ( !ali.IsMapped() ) continue ;
		if ( ali.IsDuplicate() ) continue ;
		if ( ali.MapQuality < mapqual_cutoff ) continue ; 
		if ( nopaired && ali.IsPaired() ) continue ;
		ali.BuildCharData(  ) ; 
		if ( !readgroups.empty() ) { // need to check read group
			string s ;
			ali.GetTag( "RG", s ) ;
			if ( readgroups.count( s ) == 0 ) continue ;
		}

		// alignments at the same position are added to the list (for processing this round)
		if ( pos == static_cast<size_t>(ali.Position) ) {
			ali_list.push_back( BamAlignmentIterator(ali) ) ;
			continue ;
		}
		// alignments with start outside of our desired range do not need to be processed (see below for rest of printing...)
		if ( pos > regions.back().second ) break ; 
		
//		cerr << " [ " << pos << " - " << ali.Position-1 << " ]" <<  endl ;
		
		// go through and process all alignments up to ali.Position
		// (since the file is sorted, we know that we have all overlapping alignments up to that point)
		for ( size_t i = pos ; i < static_cast<size_t>(ali.Position) ; i++ ) {
//			if ( i+1 < left_bound ) continue ; // +1 because BamTools is 0-based for the start coordinate (meh.)
//			if ( i+1 > right_bound ) break ;

			while( region_it != regions.end() && i+1 > region_it->second ) region_it++ ;

			if ( region_it == regions.end() ) return ;

//			if ( i % 100000 == 0 ) cerr << i << '\r' ;


			if ( i+1 >= region_it->first && i+1 <= region_it->second ) {

				bool header = 0 ; // "header" for the line (chr + pos) not yet printed

				// go through ali_list, count 
				list<BamAlignmentIterator>::iterator it = ali_list.begin() ;
				while ( it != ali_list.end() ) {
					it->go( i ) ;
					
					if ( it->gimme_qual() >= basequal_cutoff && it->gimme_qual() <= basequal_cutoff_max ) {
						char c = toupper( it->gimme_base(  ) );
						if ( c == 'A' || c == 'C' || c == 'G' || c == 'T' ) {
							if ( ! header ) { 
								header = 1 ;
								char ref = toupper( fasta.fetch( refname, i+1 ) ) ;
								cout << refname << "\t" << i+1 << "\t" << ref ;
							}
							int dist_5 = it->gimme_dist5( ) ;
							int dist_3 = it->gimme_dist3( ) ;
							if ( ! merge_ends ) { // default
								if ( it->orientation == 1 ) {
									if ( dist_5 >= 0 && dist_5 < ends_size ) 
										cout << "\t" << c << ":+:" << dist_5 + offset_profilenum ; 
									else if ( dist_3 >= 0 && dist_3 < ends_size ) 
										cout << "\t" << c << ":+:" << ends_size+dist_3 + offset_profilenum ;
									else 
										cout << "\t" << c << ":+:" << ends_size*2 + offset_profilenum ;
								} else if ( it->orientation == -1 ) {
									if ( dist_5 >= 0 && dist_5 < ends_size ) 
										cout << "\t" << c << ":-:" << dist_5 + offset_profilenum ;
									else if ( dist_3 >= 0 && dist_3 < ends_size ) 
										cout << "\t" << c << ":-:" << ends_size+dist_3 + offset_profilenum ;
									else 
										cout << "\t" << c << ":-:" << ends_size*2 + offset_profilenum ;
								}
							} else { // merge end
								if ( it->orientation == 1 ) {
									if ( dist_5 >= 0 && dist_5 < ends_size ) 
										cout << "\t" << c << ":+:" << dist_5 + offset_profilenum ; 
									else if ( dist_3 >= 0 && dist_3 < ends_size ) 
										cout << "\t" << c << ":+:" << dist_3 + offset_profilenum ;
									else 
										cout << "\t" << c << ":+:" << ends_size + offset_profilenum ;
								} else if ( it->orientation == -1 ) {
									if ( dist_5 >= 0 && dist_5 < ends_size ) 
										cout << "\t" << c << ":-:" << dist_5 + offset_profilenum ;
									else if ( dist_3 >= 0 && dist_3 < ends_size ) 
										cout << "\t" << c << ":-:" << dist_3 + offset_profilenum ;
									else 
										cout << "\t" << c << ":-:" << ends_size + offset_profilenum ;
								}
							}
						}
					}

					if ( it->end() )  
						it = ali_list.erase( it ) ;
					else 
						++it ;
				}
				
				if ( header ) {
					cout << endl ;
				}

				if ( ali_list.empty() ) break ;
			}

		}
		// add alignment that ended at pos+1 for next round.
		ali_list.push_back( BamAlignmentIterator(ali) ) ;
		pos = ali.Position ;
	}


	// keep going until our list is empty or end is reached 
	while ( !ali_list.empty() ) {

		while ( region_it != regions.end() && pos+1 > region_it->second ) region_it++ ;
		if ( region_it == regions.end() ) return ;
/*		if ( pos+1 < left_bound ) {
			pos++ ;
			continue ;
		}*/


		if ( pos+1 >= region_it->first && pos+1 <= region_it->second ) {

			// go through ali_list, count 
			list<BamAlignmentIterator>::iterator it = ali_list.begin() ;
			bool header = 0 ; // "header" for the line (chr + pos) not yet printed
			while ( it != ali_list.end() ) {

				it->go( pos ) ;
				
				
				if ( it->gimme_qual() >= basequal_cutoff && it->gimme_qual() <= basequal_cutoff_max ) {
					char c = toupper( it->gimme_base(  ) );
					if ( c == 'A' || c == 'C' || c == 'G' || c == 'T' ) {
						if ( ! header ) { 
							header = 1 ;
							char ref = toupper( fasta.fetch( refname, pos+1 ) ) ;
							cout << refname << "\t" << pos+1 << "\t" << ref ;
						}
						int dist_5 = it->gimme_dist5( ) ;
						int dist_3 = it->gimme_dist3( ) ;
						if ( ! merge_ends ) { // default
							if ( it->orientation == 1 ) {
								if ( dist_5 >= 0 && dist_5 < ends_size ) 
									cout << "\t" << c << ":+:" << dist_5 + offset_profilenum ;
								else if ( dist_3 >= 0 && dist_3 < ends_size ) 
									cout << "\t" << c << ":+:" << ends_size+dist_3 + offset_profilenum ;
								else 
									cout << "\t" << c << ":+:" << ends_size*2 + offset_profilenum ;
							} else if ( it->orientation == -1 ) {
								if ( dist_5 >= 0 && dist_5 < ends_size ) 
									cout << "\t" << c << ":-:" << dist_5 + offset_profilenum ;
								else if ( dist_3 >= 0 && dist_3 < ends_size ) 
									cout << "\t" << c << ":-:" << ends_size+dist_3 + offset_profilenum ;
								else 
									cout << "\t" << c << ":-:" << ends_size*2 + offset_profilenum ;
							}
						} else {
							if ( it->orientation == 1 ) {
								if ( dist_5 >= 0 && dist_5 < ends_size ) 
									cout << "\t" << c << ":+:" << dist_5 + offset_profilenum ;
								else if ( dist_3 >= 0 && dist_3 < ends_size ) 
									cout << "\t" << c << ":+:" << dist_3 + offset_profilenum ;
								else 
									cout << "\t" << c << ":+:" << ends_size + offset_profilenum ;
							} else if ( it->orientation == -1 ) {
								if ( dist_5 >= 0 && dist_5 < ends_size ) 
									cout << "\t" << c << ":-:" << dist_5 + offset_profilenum ;
								else if ( dist_3 >= 0 && dist_3 < ends_size ) 
									cout << "\t" << c << ":-:" << dist_3 + offset_profilenum ;
								else 
									cout << "\t" << c << ":-:" << ends_size + offset_profilenum ;
							}
						}
					}
				}
				if ( it->end() )  
					it = ali_list.erase( it ) ;
				else 
					++it ;
			}
			if ( header ) {
				cout << endl ;
			}
		}
		pos++ ;		
	}
}

void parse_and_set_region( string region, BamReader &reader, int &left_bound, int &right_bound, string &refname )
{

	// refname_id is the chromosome id
	int refname_id ;

	string sregion( region ) ;
	size_t colpos = sregion.find( ':' ) ;
	// not found, must be just chr then.
	if ( colpos == std::string::npos ) {
		refname_id = reader.GetReferenceID( sregion ) ;
		if ( refname_id == -1 ) {
			cerr << "Error: reference `" << sregion << "' in region argument not found." << endl;
			exit( 1 );
		}
		left_bound = 1 ;
		right_bound = reader.GetReferenceData()[refname_id].RefLength ;
		refname = sregion ;
	} else { // there are coordinates
		refname_id = reader.GetReferenceID( sregion.substr(0, colpos) ) ;
		refname = sregion.substr( 0,colpos ) ;
		if ( refname_id == -1 ) {
			cerr << "Error: reference `" << sregion.substr(0,colpos) << "' in region argument not found." << endl;
			exit( 1 );
		}
		size_t dashpos = sregion.find( '-' ) ;
		if ( dashpos == std::string::npos ) {
			cerr << "Error: Malformed region argument: colon after chr, but no coordinate given. Use `chr:start-end', `chr' or `chr:start-'." << endl ;
			exit( 1 ) ;
		}
		left_bound = std::stoi( sregion.substr(colpos+1) ) ;
		if ( dashpos+1 < sregion.length() ) 
			right_bound = std::stoi( sregion.substr(dashpos+1) ) ;
		if ( left_bound <= 0 ) {
			cerr << "Error: Malformed region argument: Start coordinate seems wrong. Use `chr:start-end', `chr' or `chr:start-'." << endl ;
			exit( 1 ) ;
		}
		if (!right_bound) {
			right_bound = reader.GetReferenceData()[refname_id].RefLength;
		}
		if ( right_bound < left_bound ) {
			cerr << "Error: Malformed region argument: End coordinate smaller than start. Use `chr:start-end', `chr' or `chr:start-'." << endl ;
			exit( 1 ) ;
		}
	}

	if ( !reader.SetRegion( refname_id, left_bound-1, refname_id, right_bound ) ) {
		cerr << "Error: Cannot set bounds to " << reader.GetReferenceData()[refname_id].RefName << ":" << left_bound << "-" << right_bound << endl;
		exit( 1 ) ;
	}


}


