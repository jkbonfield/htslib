#    Copyright (C) 2026 Genome Research Ltd.
#
#    Author: James Bonfield <jkb@sanger.ac.uk>
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

# First field:
#   INIT = initialisation, not counted in testing
#   P = expected to pass (zero return; expected output matches, if present)
#   N = expected to return non-zero
#   F = expected to fail
#
# Second field (P/N/F only):
#   Filename of expected output.  If '.', output is not checked
#
# Rest:
#   Command to execute.  $pileup is replaced with the path to the pileup test
# program

# We test BGZF2 using a small block size so we can test threading and
# indexing properly without needing large files


# ---- bzst_cli
# Check we can decompress with zstd as well as bzst, region queries,
# and multiple threads.
P empty $bzst -@4 -c -b 10000 ../ce#1000.sam > _tmp.zstd
P ../ce#1000.sam zstd -dc _tmp.zstd
P ../ce#1000.sam $bzst -dc -@4 _tmp.zstd

INIT head -c200000 ../ce#1000.sam | tail -c100000 > ce#1000.segment
P ce#1000.segment  $bzst -dc -@4 -r 100000-199999 _tmp.zstd

# ---- SAM

# FIXME: -x is to force test_view to write an index.  (It should be always
# on for zstd, but we haven't added that yet)
P empty $test_view -@4 -zzz -o block_size=10000 -l1 -p _sam.bzst -x /dev/null ../ce#1000.sam

# Test file byte regions, same as bzst
P ce#1000.segment  $bzst -dc -@4 -r 100000-199999 _sam.bzst

# Test genomic regions
INIT $test_view -bz -x _sam.gz.bai -p _sam.gz ../ce#1000.sam
INIT $test_view -p _sam.140-150 _sam.gz CHROMOSOME_I:140-150
P _sam.140-150 $test_view -@4 _sam.bzst CHROMOSOME_I:140-150

# ---- BAM
# As per SAM, but using bam binary data instead
# FIXME: -o block_size=100 fails region query.  Index issue?
P empty $test_view -@4 -bzz -o block_size=10000 -l1 -p _bam.bzst -x /dev/null ../ce#1000.sam
P _sam.140-150 $test_view -@4 _bam.bzst CHROMOSOME_I:140-150

# ---- VCF
P empty $test_view -@4 -zz -o block_size=1000 -l1 -p _vcf.bzst -x /dev/null ../index.vcf
P ../index.vcf $test_view -@4 _vcf.bzst
INIT $test_view -z -p _vcf.gz -x _vcf.gz.csi ../index.vcf
INIT $test_view -p _vcf.part _vcf.gz 2:5000000-5000029
P _vcf.part $test_view     _vcf.bzst 2:5000000-5000029
P _vcf.part $test_view -@4 _vcf.bzst 2:5000000-5000029

# ---- BCF
P empty $test_view -@4 -bzz -o block_size=1000 -l1 -p _bcf.bzst -x /dev/null ../index.vcf
P ../index.vcf $test_view -@4 _bcf.bzst
P _vcf.part $test_view     _bcf.bzst 2:5000000-5000029
P _vcf.part $test_view -@4 _bcf.bzst 2:5000000-5000029
