#!/usr/bin/perl -w
# =============================================================================
# writeCEBNXDepthtoXLS.pl - Write CEBNX depth data to Excel
# =============================================================================
# Usage: writeCEBNXDepthtoXLS.pl <Depth.CEBNX.txt>

use strict;
use Excel::Writer::XLSX;
use File::Basename;

my ($file)=@ARGV;
my $sample = basename("$file",".txt");

open(COVERAGE, "$sample.txt" ) or die "$sample.txt: $!";

my $workbook  = Excel::Writer::XLSX->new( "CoverageCEBNX.xlsx" );
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
