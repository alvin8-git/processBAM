#!/usr/bin/perl -w
# =============================================================================
# writeTMSPDepthtoXLS.pl - Write TMSP depth data to Excel
# =============================================================================
# Usage: writeTMSPDepthtoXLS.pl <Depth.TMSP.txt>

use strict;
use Excel::Writer::XLSX;
use File::Basename;

my ($file)=@ARGV;
my $sample = basename("$file",".txt");

open(COVERAGE, "$sample.txt" ) or die "$sample.txt: $!";

my $workbook  = Excel::Writer::XLSX->new( "CoverageTMSP.xlsx" );
my $coverage = $workbook->add_worksheet("$sample");

# Row and column are zero indexed
my $a = 0;
while ( <COVERAGE> ) {
    chomp;

    # Split on single tab
    my @fields = split( '\t', $_ );

    my $b = 0;
    for my $c ( @fields ) {
        $coverage->write( $a, $b, $c );
        $b++;
    }
    $a++;
}

$workbook->close();
