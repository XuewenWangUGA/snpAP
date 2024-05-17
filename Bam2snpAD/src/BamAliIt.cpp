
#include "api/BamReader.h"
using namespace BamTools;

#include "BamAliIt.h"
#include "Tools.h"
#include <iostream>
using namespace std;

BamAlignmentIterator::BamAlignmentIterator( BamAlignment &in ) : bAli(in), position_seq(0), position_cigar(0), cigar_char_count(0), endofali(false), ignore5(false), ignore3(false), had_indel(0), orientation(0)
{
	position_genome = bAli.Position ;
	parse_indels = false ;
	name = bAli.Name ;

	// jump over initial padding and insertions
	while (  bAli.CigarData[position_cigar].Type == 'I' || bAli.CigarData[position_cigar].Type == 'S' || 
		bAli.CigarData[position_cigar].Type == 'P' ) 
	{
		if ( bAli.CigarData[position_cigar].Type == 'I' ) {
			had_indel = true ;
			position_seq += bAli.CigarData[position_cigar].Length ;
		} else if ( bAli.CigarData[position_cigar].Type == 'S' ) {
			position_seq += bAli.CigarData[position_cigar].Length ;
		}
		position_cigar++ ;
	}
	// This is hypothetical: It can only happen if a CigarLine consists completely of P, S and I.
	if ( position_cigar == bAli.CigarData.size() ) {
		endofali = true ;
	}


	// Check orientation of the molecule that was the basis for this read
	// default = 0, i.e. undefined strand 
	if ( bAli.IsPaired() && bAli.IsProperPair() ) {
		if ( bAli.IsFirstMate() ) { // First mate is always in correct orientation
			ignore3 = true ; // we ignore the 3' distance if asked (see gimme_dist3)
			if ( bAli.IsReverseStrand() ) 
				orientation = -1 ;
			else 
				orientation = 1 ;
		} else if ( bAli.IsSecondMate() ) { // Its opposite day for second mate
			ignore5 = true ; // we ignore the 5' distance if asked (see gimme_dist3)
			if ( bAli.IsReverseStrand() ) 
				orientation = 1 ;
			else 
				orientation = -1 ;
		}
	} else if ( !bAli.IsPaired( ) ) { // not paired
		if ( bAli.IsReverseStrand() ) 
			orientation = -1 ;
		else 
			orientation = 1 ;
	 }
	
}

void BamAlignmentIterator::go( int newpos ) 
{

	this->skip_genome( newpos - position_genome ) ;
}

void BamAlignmentIterator::skip_genome( int count )
{
	if ( count < 0 ) {
		throw range_error( "Negative distance to next tested genome position." ) ;
	}
	// move along sequence and cigar string until we skipped over /count/ genome bases
	while ( count > 0 && position_cigar < bAli.CigarData.size() ) {

		// move ahead in positions on the sequence and along the genome depending on the cigar op
		switch( bAli.CigarData[position_cigar].Type ) {
			case 'M':
			case '=':
			case 'X':
				position_seq++ ;
				position_genome++ ;
				count-- ;
				break ;
			case 'S':
				position_seq++ ;
				break ;
			case 'I':
				had_indel = true ;
				position_seq++ ;
				break ;
			case 'D':
				had_indel = true ;
				position_genome++ ;
				count-- ;
				break ;
			case 'N':
			case 'H':    
				position_genome++ ;
				count-- ;
				break ;
		}

		// go to next cigar-aligned character
		cigar_char_count++ ;
		if ( cigar_char_count >= bAli.CigarData[position_cigar].Length ) {
			cigar_char_count = 0 ;
			position_cigar++ ;
		}
	}

	// We could now be in an insertion or a padded stretch: jump over these
	while ( position_cigar < bAli.CigarData.size() && 
		( bAli.CigarData[position_cigar].Type == 'P' || bAli.CigarData[position_cigar].Type == 'I' ) )
	{
		if ( bAli.CigarData[position_cigar].Type == 'I' ) {
			had_indel = true ;
			if ( cigar_char_count != 0 ) throw logic_error( "non-start position in indel stretch in go()" ) ;
			position_seq += bAli.CigarData[position_cigar].Length ;
		}
		position_cigar++ ;
	}
	
	// Are we at the end?
	if ( position_cigar == bAli.CigarData.size() )
		endofali = true ;
	else if ( bAli.CigarData[position_cigar].Type == 'S' || bAli.CigarData[position_cigar].Type == 'H' ) {
		position_cigar = bAli.CigarData.size() ;
		endofali = true ; // if we see one of these clip-types, we must have hit the end 
	}

}

char BamAlignmentIterator::gimme( char what ) 
{
    	if ( endofali ) return 0 ;
	switch ( bAli.CigarData[position_cigar].Type ) {
		case 'H':
			return 0 ; // not exactly a gap, so we give this back...
			break ;
		case 'S':
		case 'P': 
		case 'I':
#ifdef DEBUG
			cerr << "Pos cigar: " << position_cigar << endl ;
			cerr << "Pos seq: " << position_seq << endl ;
			for ( int i = 0 ; i < bAli.CigarData.size() ; i++ ) {
				cerr << bAli.CigarData[i].Type << bAli.CigarData[i].Length << " " ;
			}
			cerr << endl ;
			cerr << bAli.CigarData[position_cigar].Type << endl ;
#endif
			throw logic_error( "Hit P, S or I in Cigar Line which shouldn't happen (ever!)" ) ;
			break ;
		case 'D':
			had_indel = true ; // XXX: probably not needed; 
			return '-' ;
			break ;
		case 'N':
			return '-' ;
			break ;
		case 'M':
		case 'X':
		case '=':
			if ( what == 's' ) // sequence
				return bAli.QueryBases[position_seq] ;
			else if ( what == 'q' ) // quality
				return bAli.Qualities[position_seq]-33 ;
			break ;
		
	}
#ifdef DEBUG
	cerr << "Pos cigar: " << position_cigar << endl ;
	cerr << "Pos seq: " << position_seq << endl ;
	for ( int i = 0 ; i < bAli.CigarData.size() ; i++ ) {
		cerr << bAli.CigarData[i].Type << bAli.CigarData[i].Length << " " ;
	}
	cerr << endl ;
	cerr << bAli.CigarData[position_cigar].Type << endl ;
#endif
	throw logic_error( "Hit an unknown character in a Cigar Line" ) ;
	
}

char BamAlignmentIterator::gimme_base(  ) 
{
	return gimme( 's' ) ;
}

char BamAlignmentIterator::gimme_qual(  ) 
{

	return gimme( 'q' ) ;

}

int BamAlignmentIterator::gimme_dist5(  ) 
{
	if ( ignore5 ) return -1 ;
	else {
		if ( orientation == 1 ) return position_seq ;
		else if ( orientation == -1 ) return bAli.Length - 1 - position_seq ;
		else return -1 ;
	}
}

int BamAlignmentIterator::gimme_dist3(  )
{
	if ( ignore3 ) return -1 ;
	else {
		if ( orientation == -1 ) return position_seq ;
		else if ( orientation == 1 ) return bAli.Length -1 - position_seq ;
		else return -1 ;
	}
}


bool BamAlignmentIterator::next_indel( s_indel& indel ) 
{
	// no need to bother with reads of unknown orientation
	if ( orientation == 0 ) 
		return false ;

	// if at the end, we reset 
	if ( parse_indels == false ) {
		if ( position_cigar != bAli.CigarData.size() ) 
			throw logic_error( "parsing indels in read that hasn't reached end." ) ;
		parse_indels = true ;
		position_cigar = 0 ;
		position_seq = 0 ;
		position_genome = bAli.Position ;
	}	

	while ( bAli.CigarData[position_cigar].Type != 'I' && 
		bAli.CigarData[position_cigar].Type != 'D' &&
		position_cigar != bAli.CigarData.size() ) {
		// move positions according to CigarData entry
		switch (bAli.CigarData[position_cigar].Type) 
		{
			case 'M':
			case '=':
			case 'X':
				position_seq += bAli.CigarData[position_cigar].Length ;
				position_genome += bAli.CigarData[position_cigar].Length ;
				break ;
			case 'S':
				position_seq += bAli.CigarData[position_cigar].Length ;
				break ;
			case 'I':
				// cant happen
				position_seq++ ;
				break ;
			case 'D':
				// cant happen
				position_genome++ ;
				break ;
			case 'N':
			case 'H':    
				position_genome += bAli.CigarData[position_cigar].Length ;
				break ;
			// 'P' doesn't do anything to the iterators
		}
		position_cigar++ ;
	}

	if ( position_cigar == bAli.CigarData.size() ) 
		return false ;
	else {
		if ( bAli.CigarData[position_cigar].Type == 'I' ) {
			indel.is_insertion = true ;
			indel.size = bAli.CigarData[position_cigar].Length ;
			indel.seq = bAli.QueryBases.substr( position_seq, indel.size ) ;
			indel.pos = position_genome ;
			if ( orientation == 1 ) {
				indel.dist5 = position_seq ;
				indel.dist3 = bAli.Length - position_seq - indel.size ;
				indel.is_reverse = false ;
			} else if ( orientation == -1 ) {
				indel.dist3 = bAli.Length - position_seq - indel.size ;
				indel.dist5 = position_seq ;
				revcomp( indel.seq ) ;
				indel.is_reverse = true ;
			} //else  throw logic_error( "Reads with unknown orientation should not be queried for insertions" ) ;
			
			// move iterators
			position_seq += bAli.CigarData[position_cigar].Length ;
			position_cigar++ ;
		} else { // must be 'D'
			indel.is_insertion = false ;
			indel.size = bAli.CigarData[position_cigar].Length ;
			indel.pos = position_genome ;
			// indel.seq = "" ; // no need to set "". Will be overwritten with next insertion
			if ( orientation == 1 ) {
				indel.dist5 = position_seq ;
				indel.dist3 = bAli.Length - position_seq ;
				indel.is_reverse = false ;
			} else if ( orientation == -1 ) {
				indel.dist5 = bAli.Length - position_seq ;
				indel.dist3 = position_seq ;
				indel.is_reverse = true ;
			} else if ( orientation == -1 ) {
			} //else  throw logic_error( "Reads with unknown orientation should not be queried for deletions" ) ;
			// move iterators
			position_genome += bAli.CigarData[position_cigar].Length ;
			position_cigar++ ;
		}
	} 
	return true ;
}

