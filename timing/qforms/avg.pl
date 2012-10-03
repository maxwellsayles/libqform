#!/usr/bin/perl -w

use strict;
use vars;
use warnings;

# ranges are inclusive
my @ranges = ([16, 59], [60, 118], [119, 140]);

sub sumFileWithRanges {
    my ($filename, $ranges) = @_;
	
    my @sums = map { 0 } 1..scalar(@$ranges);
	
    open FILE, "<$filename";
	
    while (<FILE>) {
	if ($_ =~ m/(\d+),? (\d+)/) {
	    my $i = 0;
	    while ($i < scalar(@$ranges)) {
		if ($1 >= $ranges->[$i]->[0] and $1 <= $ranges->[$i]->[1]) {
		    $sums[$i] += $2;
		}
		$i ++;
	    } 
	}
    }	
    
    close FILE;
    
    return \@sums;
}

sub composeAndSquare {
    open COMPOSE, "<compose.dat";
    open SQUARE, "<square.dat";
    open CUBE, "<cube.dat";

    open COMPOSE_SQUARE, ">compose_square.dat";

    my $composeLine = <COMPOSE>;
    my $squareLine = <SQUARE>;
    my $cubeLine = <CUBE>;
    
    while ($composeLine && $squareLine && $cubeLine) {
	
	last if $composeLine !~ m/(\d+),? (\d+)/;
	my $nbits = $1;
	my $composeTime = $2;
	
	last if $squareLine !~ m/\d+,? (\d+)/;
	my $squareTime = $1;

	last if $cubeLine !~ m/\d+,? (\d+)/;
	my $cubeTime = $1;

	my $sum = $composeTime + $squareTime;

	print COMPOSE_SQUARE "$nbits, $sum\n";

	$composeLine = <COMPOSE>;
	$squareLine = <SQUARE>;
	$cubeLine = <CUBE>;
    }


    close COMPOSE;
    close SQUARE;
    close CUBE;
}

eval {
    my $compose_sums = sumFileWithRanges 'compose.dat', \@ranges;
    my $square_sums = sumFileWithRanges 'square.dat', \@ranges;
    my $cube_sums = sumFileWithRanges 'cube.dat', \@ranges;
    
    my $i = 0;
    while ($i < scalar(@ranges)) {
	print "$ranges[$i]->[0]-$ranges[$i]->[1] compose: $compose_sums->[$i]\n";
	print "$ranges[$i]->[0]-$ranges[$i]->[1] square: $square_sums->[$i]\n";
	print "$ranges[$i]->[0]-$ranges[$i]->[1] cube: $cube_sums->[$i]\n";
	print "\n";
	$i ++;
    }

    $i = 0;
    while ($i < scalar(@ranges)) {
	my $b = $square_sums->[$i];
	my $a = $compose_sums->[$i] / $b;
	my $c = $cube_sums->[$i] / $b;
	print "$ranges[$i]->[0]-$ranges[$i]->[1] compose: $a\n";
	print "$ranges[$i]->[0]-$ranges[$i]->[1] square: 1\n";
	print "$ranges[$i]->[0]-$ranges[$i]->[1] cube: $c\n";
	print "\n";
	$i ++;
    }


    composeAndSquare;
};
if ($@) {
    print STDERR "$@\n";
    exit 1;
};
exit 0;
