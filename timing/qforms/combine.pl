#!/usr/bin/perl -w

use strict;
use vars;
use warnings;

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

=pod
This reads the compose and square files and generates a compose-square
file that is the sum of the two.  It represents the time of the
alternative to using the cube routine.
=cut
sub combine {
    my $ext = shift;
    open COMPOSE, "<compose-${ext}.dat";
    open SQUARE, "<square-${ext}.dat";
    open COMPOSE_SQUARE, ">compose_square-${ext}.dat";

    my $composeLine = <COMPOSE>;
    my $squareLine = <SQUARE>;
    while ($composeLine && $squareLine) {
	last if $composeLine !~ m/(\d+),? (\d+\.\d+)/;
	my $nbits = $1;
	my $composeTime = $2;

	last if $squareLine !~ m/\d+,? (\d+\.\d+)/;
	my $squareTime = $1;

	my $sum = $composeTime + $squareTime;
	print COMPOSE_SQUARE "$nbits, $sum\n";

	$composeLine = <COMPOSE>;
	$squareLine = <SQUARE>;
    }

    close COMPOSE;
    close SQUARE;
}

eval {
    combine "64";
    combine "128";
    combine "mpz";
    combine "mpir";
    combine "pari";
};
if ($@) {
    print STDERR "$@\n";
    exit 1;
};
exit 0;
