#!/usr/bin/perl

# "Algorithm R" (Jeffrey Vitter 1985)

use strict ;

my @res ;
my $n = $ARGV[0] ;
my $i = 0 ;

while (<STDIN>) {
        $i++ ;
        if ( @res < $n ) {
                push @res, $_ ;
        } else {
                my $r = rand( $i ) ;
                $res[$r] = $_ if ( $r < $n ) ;
        }
}

foreach my $x (@res) { print $x ; }

