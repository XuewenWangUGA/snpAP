
#ifndef BAMALIIT_H
#define BAMALIIT_H

#include <vector>
using namespace std;

#include "api/BamReader.h"
#include "api/BamWriter.h"
using namespace BamTools;


struct s_indel 
{
	bool is_insertion ; // if the read has an insertion compared to the reference
	bool is_reverse ; // strand orientation
	bool eof ; // is set if this is the last indel
	unsigned char size ; // size of insertion or deletion
	// distances are the coordinates of the 5' and 3' most base next to the insertion/deletion
	int dist5 ; // distance to 5' end of sequence 
	int dist3 ; // distance to 3' end of sequence   
	string seq ; // sequence in case of an insertion
	size_t pos ; // position (5' most adjacent basepair in the genome)
} ;


// iterates over an alignment to a reference genome.
class BamAlignmentIterator
{
	public:
		BamAlignmentIterator( BamAlignment& in ) ;
		void go( int genome_position ) ; // go to genomic position (must be > prev. pos)
		void skip_genome( int count ) ; // skip over count genome bases 
		char gimme_base(  ) ; // returns base aligning to current position
		char gimme_qual(  ) ; // returns qual for base aligning to current position
		int gimme_dist5(  ) ; // returns distance to 5' end in fragment (or -1 if n.a.)
		int gimme_dist3(  ) ; // returns distance to 3' end in fragment (or -1 if n.a.)
		bool end()
		{ return endofali ; } // we stepped past the last character

		bool next_indel( s_indel& indel ) ;

	private:
		char gimme( char ) ;
		BamAlignment bAli ;
		// where we are on the genome:
		size_t position_genome ;
		// where we are in unaligned sequence & quality scores:
		size_t position_seq ;
		//  where we are in the cigar line:
		size_t position_cigar ;
		size_t cigar_char_count ;
		// we are at the end of the sequence
		bool endofali ;
		
		// small nonsense that makes things easier
		bool ignore5, ignore3 ; // ignore 5' resp. 3' end distance if we are in a paired end mate.
		bool had_indel ; // did we see an indel?
		bool parse_indels ; // is set when the Cigar line was parsed once and we are fetching indels

	public:
		bool has_indel() const { 
			if ( endofali ) return had_indel ; 
			else throw logic_error( "should never be asked for indel before iterating once through entire Alignment." ) ; 
		} 
		// orientation of the molecule (not necessarily equal with the reads orientation in the case of paired-end)
		char orientation ; // { 1, -1, 0 } for { +, -, undef }
		string name ;


} ;

#endif
