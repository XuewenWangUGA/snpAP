#!/usr/bin/perl

use strict ;
use Getopt::Long ;

my $help ;
my $winlength = 1_000_000 ;
my $jump = 100_000 ;
my $inputf="" ;
my $errorf="" ;

GetOptions( 'help|?' => \$help, 'window_length|l=i' => \$winlength, 'window_jump|j=i' => \$jump, 'input|i=s' => \$inputf, 'error|e=s' => \$errorf ) ;

if ( $help ) { print_help() ; exit(0) }
if ( $inputf eq "" ) { print_help() ; print "Please specify inputfile with -i\n" ; exit( 1 ) ; }
if ( $errorf eq "" ) { print_help() ; print "Please specify error profile with -e\n" ; exit( 1 ) ; }

my $IN ;
if ( $inputf =~ /\.gz$/ ) {
	open $IN, "gunzip -c $inputf|" ;
} else {
	open $IN, "<$inputf" ;
}

my @buf ;
my $l = <$IN> ;
my ( @r ) = split /\t/, $l ;
push @buf, [ @r ] ;
my $pos = $r[1] ;

while (<$IN>) {
	my ( @r ) = split /\t/ ;
	# we have not filled this window
	if ( $r[1] < $pos+$winlength ) {
		push @buf, [ @r ] ;
		next ;
	}
	while ( $r[1] >= $pos+$winlength ) {
		# we have filled this window, process
		my $het = snpAD( \@buf ) ;
		print join "\t", ( $pos, scalar @buf, $het ) ; print "\n" ;

		# shorten array, adjust pos
		if ( @buf == 0 ) { $pos += $jump ; next ; }
		my $n = 0 ;
		while ( $n < @buf && $buf[$n]->[1] < $pos+$jump ) { $n++ ; }
		@buf = splice( @buf, $n, $#buf ) ; # shorten
		$pos += $jump ;
	}
	push @buf, [ @r ] ;
}
# run the leftovers
my $het = snpAD( \@buf ) ;
print join "\t", ( $pos, scalar @buf, $het ) ; print "\n" ;

# Run snpAD
sub snpAD
{
	my ( $buf ) = @_ ;
	return -1 if @$buf == 0 ;

	my $X = $$ ;

	open SNPAD, "|/home/pruefer/svn/snpAD-0.3/snpAD/snpAD -e $errorf -R -F -o $X.est /dev/stdin 2>$X.err >$X.out" ;
	foreach ( my $i = 0 ; $i < @$buf ; $i++ ) {
		print SNPAD join "\t", @{$buf->[$i]} ;
	}
	close SNPAD ;

	open IN, "<$X.est" ;
	my $EST = <IN> ;
#	print STDERR $EST ;
	chomp $EST ;
	return $EST ;
	my @est = split /,/, $EST ;
	return ( $est[4]+$est[5]+$est[6]+$est[7]+$est[8]+$est[9] ) ;

#	open IN, "<$X.err" ;
#	while (<IN>) {
#		print ;
#	}
#	close IN ;

}


sub print_help
{
	print "Usage: $0 -l [window_length] -j [window_jump] -e [error profile] -i [inputfile]\n" ;
	print " Note: snpAD must be in \$PATH\n" ;
}

