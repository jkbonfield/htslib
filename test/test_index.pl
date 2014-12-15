#! /usr/bin/env perl
#
#    Copyright (C) 2014 Genome Research Ltd.
#
#    Author: Rob Davies <rmd+git@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

use strict;
use warnings;
use Getopt::Long;

my $passed = 0;
my $failed = 0;
my $framework = $ENV{TEST_FRAMEWORK} ? "$ENV{TEST_FRAMEWORK} " : '';
my $tidy = 0;
my $fail_early = 0;

GetOptions("tidy!" => \$tidy, "fail_early!" => \$fail_early);

my $small_ref = 'index_small_ref';
my $big_ref   = 'index_big_ref';

if (0 == make_bam($small_ref)) {
    if (!$fail_early || !$failed) { test_index($small_ref, 'b'); }
    if (!$fail_early || !$failed) { test_index($small_ref, 'c'); }
}

if ((!$fail_early || !$failed) && 0 == make_bam($big_ref)) {
    if (!$fail_early || !$failed) { test_index($big_ref,   'c'); }
}

print "\nNumber of tests:\n";
printf "\ttotal  .. %d\n", $passed + $failed;
printf "\tpassed .. %d\n", $passed;
printf "\tfailed .. %d\n", $failed;
print "\n";

exit($failed ? 1 : 0);

sub test {
    my ($cmd) = @_;

    print "  $cmd\n";
    if (system("$cmd || exit 1") != 0) {
        print "FAIL\n";
        $failed++;
	return -1;
    } else {
        $passed++;
	return 0;
    }
}

sub make_bam {
    my ($name) = @_;

    print "\n === Making BAM file $name.tmp.bam ===\n\n";

    return test("${framework}./test_view -S -b '$name.sam' > '$name.tmp.bam'");
}

sub test_index {
    my ($name, $type) = @_;

    printf("\n === Test %s index on $name.tmp.bam ===\n\n",
	   $type eq 'b' ? 'bai' : 'csi');

    my $res = 0;

    # Extract the locations of each read (NB the test sam file must
    # be constructed so none overlap for the tests to work.)
    my @headers;
    my @lines;
    my @regions;
    open(my $sam, '<', "$name.sam") || die "Couldn't open $name.sam : $!\n";
    while (<$sam>) {
	if (/^@/) {
	    push(@headers, $_);
	    next;
	}
	my @F = split;
	push(@lines, $_);
	push(@regions, [$F[2], $F[3], $F[3] + length($F[9])]);
    }
    close($sam) || die "Error reading $name.sam : $!\n";

    unlink("$name.tmp.bam.bai", "$name.tmp.bam.csi");
    if (0 != test("${framework}./test_index -$type '$name.tmp.bam'")) {
	return -1;
    }
    
    # Try to extract each sequence
    for (my $seq = 0; $seq < @regions; $seq++) {
	if (0 != test_extract($name, $seq, $seq, \@headers, \@regions, \@lines)){
	    if ($fail_early) { return -1; }
	    $res = -1;
	}
    }

    return $res;
}

sub test_extract {
    my ($name, $start, $end, $headers, $regions, $lines) = @_;

    my $reg = "$regions->[$start]->[0]:$regions->[$start]->[1]-$regions->[$end]->[2]";
    my $bam = "$name.tmp.bam";
    my $expected = "$bam.expected.$reg";
    my $out = "$bam.got.$reg";

    open(my $exp, '>', $expected) || die "Couldn't open $expected : $!\n";
    print $exp @$headers, @$lines[$start..$end];
    close($exp) || die "Error writing to $exp : $!\n";

    
    if (0 != test("${framework}./test_view '$bam' $reg > $out")) { return -1; }
    if (0 != test("./compare_sam.pl '$expected' '$out'")) { return -1; }
    if ($tidy) {
	unlink($expected, $out);
    }

    return 0;
}
