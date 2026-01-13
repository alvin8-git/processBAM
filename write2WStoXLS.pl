#!/usr/bin/perl -w
# =============================================================================
# write2WStoXLS.pl - Write two worksheets to Excel
# =============================================================================
# Usage: write2WStoXLS.pl <XLSname> <WS1 data> <WS2 data>

use strict;
use Excel::Writer::XLSX;
use File::Basename;

my ($XLS,$file1,$file2)=@ARGV;
my $data1 = basename("$file1",".txt");
my $data2 = basename("$file2",".txt");

open(WS1, "$data1.txt" ) or die "$data1.txt: $!";
open(WS2, "$data2.txt" ) or die "$data2.txt: $!";

my $workbook  = Excel::Writer::XLSX->new( "$XLS.xlsx" );
my $worksheet1 = $workbook->add_worksheet("$data1");
my $worksheet2 = $workbook->add_worksheet("$data2");

# Row and column are zero indexed
my $a = 0;
while ( <WS1> ) {
    chomp;

    # Split on single tab
    my @fields = split( '\t', $_ );

    my $b = 0;
    for my $c ( @fields ) {
        $worksheet1->write( $a, $b, $c );
        $b++;
    }
    $a++;
}

my $d = 0;
while ( <WS2> ) {
    chomp;

    # Split on single tab
    my @fields = split( '\t', $_ );

    my $e = 0;
    for my $f ( @fields ) {
        $worksheet2->write( $d, $e, $f );
        $e++;
    }
    $d++;
}

$workbook->close();
