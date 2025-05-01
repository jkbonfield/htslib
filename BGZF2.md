Newest todos

- Make the magic number block also have a predefined magic number for
  SAM, BAM, VCF, BCF, BED, FASTQ, etc as well as the generic detection
  via string.

- Use a uniform BGZF2 frame type where the first byte is sub-type.
  0=file format magic (first X number of bytes for detection)
  1=per-frame meta-data (generic frame info, like cram's container)
  2=genomic index
  3=per-file meta-data (a generic idxstats)

  Per frame meta-data, use name space like SAM tags.
  Official:
  - Number of mapped sequences
  - Number of unmapped sequences
  - Total number of bases
  - Local seq checksum?
  - Local qual checksum?
  
  Unofficial:
  - Whatever we like (lowercase key?)
  
  Per file meta-data.
  Official:
  - Everything in IDX stats
  - Samtools checksum outputs if BAM

  Unmapped data? Name sorted data?

  Genomic index meta-data.
  - Per frame (eg number of records, mapped, unmapped, bases, etc)
  - Per chromosome (as above)
  - Per file
  Perhaps we have key=value here too, with fixed keys that always
  occur in that order and non-fixed keys which are explicit.
  Eg:
  Fixed_per_frame=number_of_mapped, number_of_unmapped
  A frame with 990 mapped, 10 unmapped can then have:
      <990> <10> <read1=490> <read2=510>
  So we see standard fixed items not needing the key duplicated, and
  some that have an explicit key=value text.
 
  
  Unofficial:
  - Whatever we wish

- Genomic index frequency.  We should stipulate at least one gindex
  entry per data frame and at least one gindex per chromosome, but we
  can have more if we need them.  Whatever frequency we wish
  basically.

- If seekable index or genomic index don't fit within a single zstd
  frame then we need to be able to concatenate them together to form a
  single index.  How likely is this to happen?  Unknown, but we can
  write it in the spec even if we don't support it in code yet.

- Add checksum to meta-data frames too.

=============================================================================


TODO
----

- Index.
  - Also need to consider multi-ref frames, for when we get
    many small contigs.  So multiple index entries per frame.

- Break from pzstd and use a generic next frame meta-data block.
  This will allow us to encode arbitrary key-value pairs.
  Like CRAM header (from cram_dump):

  Container pos 9495 size 607898
      Ref id:            0
      Ref pos:           10001 + 594302
      Rec counter:       0
      No. recs:          10000
      No. bases          1000000

  But we'd permit more than one ref:start-end range per meta-data block.
  This means we can switch chromosome in a middle of a data frame *if
  we wish* (but it may be ideal not to).

  We should also consider things like no. (un)mapped reads / bases for
  crude depth plots, and possibly linear indices.

- Add more key/value pair fields in the header.
  That includes a mandatory file format version number!



BGZF2
-----

The BGZF2 format is based on ZSTD, using a block based format much
like BGZF(1) does for Deflate.

ZSTD has data frames and skippable frames.  We use skippable frames
for meta-data and a self contained index mapping uncompressed offsets
to compressed offsets.

BGZF2 can be viewed as a union of the pzstd and seekable_format
"contrib" formats as found at the zstd Github page.
    https://github.com/facebook/zstd/tree/dev/contrib


File Layout
-----------

[bgzf2 header]
[bgzf2 block meta data]
[zstd data frame]
[bgzf2 block meta data]
[zstd data frame]
...
[bgzf2 file meta data]
[bgzf2 genomic index]
[seekable_format index]

Extensive use of zstd skippable frames are used.  These are ZSTD
compliant frames which carry no compressed payload, so can be skipped
over by any ZSTD decoder.  Instead they carry meta-data useful for
file type detection, rapid parallel decoding, and internal indices.

In order for standard zstd tools to be able to skip past these custom
meta-data frames, the skippable frames have a specific data layout of
a 4 byte magic number (0x184D2A50 to 0x184D2A5F), a 4 byte frame
length N, and then N bytes of frame specific payload.

We use the skippable magic number 0x184D2A5B as a BGZF2 frame.
We further subdivide this into BGZF2 frame types as follows:

- 4 bytes BGZF2 frame magic number: 0x184D2A5B
- 4 bytes of remaining frame length "N"
- 1 byte of BGZF2 sub-type
- 1 byte of sub-type specific version number
- N-2 bytes of meta-data

This combination of skippable frames means the following tools can
decompress a bgzf2 file:

- "zstd".
  This provides single threaded streaming decode.

- "seekable_decompression" from Zstd's contrib/seekable_format directory.
  This provides random access via byte ranges in the uncompressed data.

- "bgzip2" and htslib.
  These provide all of the above, but also offers random access by
  genomic region or a range of record numbers.


Bgzf2 header frame
------------------

[TODO: Maybe we want to make it more feature
rich, and to add specific BGZF2 wrapper version numbering in there
too, so we can extend it with additional key/value meta-data (one of
which is file type).]

The data contents are an uncompressed copy of the first data bytes,
used for file format detection.  There is no fixed limit on this
length, but it is recommended to be not significantly more than is
required to accurately determine file type and some basic versioning.

A header frame is formatted as:

   4: 0x184D2A5B (bgzf2 magic number)
   4: N (remaining size of skippable frame)
   1: 0 (BGZF2 Header type)
   1: 0 (BGZF2 header format version)
   4: "BGZ2" identifier
N-10: "BAM\x01????@HD.VN:1.4.SO:coordinate\n"
   4: XXHash-64

For VCF we need to get more header so we can decode
"##fileformat=VCF4.2".  The meta-data length is therefore not fixed.

Tools should be capable of working with reduced meta-data here, for
example having just "BAM\x01".

The checksum at the end uses the Zstd XXHash-64 used elsewhere in zstd
and is the hash of the entire frame.


BGZF2 block meta data
---------------------

This holds meta-data about the next zstd compressed frame.

The zstd data frames do not have length information, so this is an
important part of the block meta-data (inspired by pzstd) and is
necessary to permit parallel processing.

Additionally we support arbitrary key/value pairs.  This could include
coverage data for alignments, such as a chromosome range, checksums,
or arbitrary user defined data.

The block meta-data is formatted as:

   4: 0x184D2A5B (bgzf2 magic number)
   4: N (remaining size of skippable frame)
   1: 1 (BGZF2 block meta-data type)
   1: 0 (BGZF2 header format version)
   4: size of next zstd data frame
N-10: meta-data matching ((key=value)(;key=value)*)?
   4: XXHash-64

TODO: consider the alternatives on meta-data encoding.
- Binary too, BAM style?
  ID length of 2 is a bit restrictive perhaps.  Maybe 4 bytes instead?
- Textual, VCF INFO style?
- Nested, so maybe JSON?

Also learn from SAM and have specific name spaces, so uppercase
starting letter for key is reserved, lowercase is user controlled.



Zstd data frame
---------------

This is a standard data frame holding compressed data to be decoded
and returned by the zstd library.  It starts with the magic number
0xFD2FB528.

See RFC8478 for further details.


BGZF2 file meta data
---------------------

This holds aggregated meta-data about the entire file.
It has similar to the block meta-data but lacks the size of the next
zstd data frame.

The file meta-data is formatted as:

   4: 0x184D2A5B (bgzf2 magic number)
   4: N (remaining size of skippable frame)
   1: 2 (BGZF2 file meta-data type)
   1: 0 (BGZF2 header format version)
N-10: meta-data matching ((key=value)(;key=value)*)?
   4: XXHash-64


Indexing
--------

BGZF version 1 indices can be BAI, CSI and TBI.  They are also
external files.  This myriad of index formats leads to a variety of
problems:

- The naming differs between tools, sometimes it is foo.bam.bai and
  sometimes foo.bai.

- Multiple files causes problems when downloading from object stores
  where the filenames are content hashes.

- If multiple indices exist (.bai and .csi) there are undocumented
  precedence rules.

- There is potential for catastrophy when the bam file is rewritten
  without the index being recreated.  Indices do not have time stamps
  or shared secrets with their BAM counterparts, making it hard to
  detect this problem.

- These indexes are only for aligned data, so there is no index
  capability on unmapped records, name-sorted or unsorted data.  Some
  of these are covered by yet another index format: GZI.

- There is insufficient flexibility in what is indexed, leading to
  custom indices such as PBI.

So in BGZF2 the index is embedded within the file.

[TODO: It is also partially distributed.  Useful? Probably not]
[TODO: Consider and index on the index to speed up searches on tiny
regions.]

The BGZF1 indices combined mapping of genomic index to compressed
offsets together.  We split this in two, with a genomic index to map
chromosome ranges to an uncompressed offsets, and an uncompressed to
compressed offset which is covered by the seekable index format.

BGZF2 Genomic index frame
-------------------------

This is a genomic coordinate index used for random access by
chromosome and position sorted data, or by record number if unsorted.

[TODO: The index itself is self-indexing, meaning it is possible to
hop around within it via precomputed offsets.  In practice however it is
likely the most performant use is to load the entire index into
memory, especially if compressed.]

The genomic index is formatted as:

   4: 0x184D2A5B (bgzf2 magic number)
   4: N (remaining size of skippable frame)
   1: 3 (BGZF2 file meta-data type)
   1: 0 (BGZF2 header format version)
   1: encoding flags
   ?: genomic index data
   4: N-8 (to permit reverse reading)
   4: "BG2I" magic number
   4: XXHash-64

Encoding flags:

Bit 0:    1 if "genomic index data" is zstd compressed, 0 if uncompressed
Bits 1-7: 0 (reserved)

[ FIXME: we should have encoding flag per reference, ie NR times), so
we can do random access on our index permits its own compressed block.
Or is it worth compressing the index at all given larger block sizes?
Possibly not. ]

[ Should we use variable sized integer encoding instead?  It's a bit
more faffing, but not too much and it makes things totally data size
agnostic. However it makes decoding also trickier, especially random
access to within an index.  An sparse index on the index solves that
however. ]

The genomic index section is formatted as:

[ Per NR references]
1: reference flag
   Bit 0: 1 if long ref (>= 4GB)
   Bit 1-7: 0 (reserved)
4 or 8: length of reference (see ref flag)
?: Nul-terminated name of reference.
4: byte offset into NI index table for start of this reference
4: number of bytes in NI index table for this reference
4: size of meta-data for this reference MD
MD: meta-data

[ per NI index entry]
1: flag
   Bit 0: 1 if multiple references are within this data frame
   Bit 1: 1 if all records in data frame are unmapped (also set if
          (global.flag & 1) == 0)
   Bit 2: 1 if multiple references are present and data is mapped
   Bit 3-7: 0 (reserved)
4: NIR: Number of references (if NI.flag bit 2 set),
   otherwise 1 if NI.flag bit 1 is clear,
   otherwise 0.
[ Per NIR ]
4: Reference number
4 or 8: reference start (size from NR.flag bit 0)
4 or 8: reference span (end-start+1)
4: Number of records aligned
4: Number of records unaligned


[ Footer ]
4: distance back to start of zstd frame
4: 0x8F92EABB: BGZF2 index magic number

[ TODO: maybe a self-describing index format so we can add extra
meta-data columns. Eg number of QC failures, or number of secondary
alignments. ]

[ FIXME: do we want meta-data here too?  Is it just a copy of the
per-frame meta-data, or a user-controlled subset of it?  It can be
useful to aggregate it together in one place perhaps, permitting
additional querying and index abilities.

Global and per-ref meta-data is key\0value\0 pairs.  It replaces the
magic bin numbers used in BAI and CSI.  We can view it aggregated
together for the entire file, or on a per-reference basis. ]

Tips on usage:

- The index can be loaded into a nested containment list.  This is a
  list of items where all curr.start >= last.start and curr.end >=
  last.end. Ie we never have items that are smaller than their previous
  one and entirely contained within it.  If that happens, the item
  itself becomes a new nested containment list, in a recursive manner.

- We can then do binary searching to find overlaps.  Search on start
  to find first item overlapping the range.  We can then linearly step
  right until we get to the first item beyond the range.

  - For items containing sub-lists (containments), recurse.

Question: is this still inefficient in a meaningful manner when
compared to an R-tree binning system.  Imagine a mix of very
long-reads and short-reads.  Querying a region may return multiple
frames.

For example a file of 8 frames numbered 0 to 7 with long and short
reads shown mapped against a chromosome (left to right), with their
frame numbers listed, along with a region to query:

                1         2         3         4         5         6
Pos:   123456789012345678901234567890123456789012345678901234567890123456
                                  |-----| query
long:              11111111111111111111  44444444444444444444
        000000000000000000000   3333333333333333333  55555555555555555555
                                  |-----|
short: 000 000 111 111 222 222 333 333 444 444 555 555 666 666 777 777
        000 000 111 111 222 222 333 333 444 444 555 555 666 666 777 777
          000 000 111 111 222 222 333 333 444 444 555 555 666 666 777 777
                                  |-----|

Our index would be
frame 0: chr:1-22  
frame 1: chr:9-32
frame 2: chr:17-26 (region contained within frame 1 region)
frame 3: chr:25-44
frame 4: chr:33-54
frame 5: chr:41-66
frame 6: chr:49-58 (region contained within frame 5 region)
frame 7: chr:57-66 (region contained within frame 5 region)

As a nested containment list we would have frames 0, 1, 3, 4, 5 with
frames 1 and 5 containing sub-lists.  The benefit of enforced nesting
of chr:start-end is as start increases so does end, which in turn
makes binary searching trivial and yields the same complexity as a
binary R-tree.

Our query for chr:28-34 hits long reads in frame 1 and 3 and short
reads in frame 3 and 4.

We could either explicitly seek to frame 1, 3, and 4, or we could
simply seek to frame 1 and read/skip as apprpriate (skipping frame 2,
based on frame meta-data) until we are beyond the end of our query
region.

TODO: We may wish to add per-frame linear indices, like the BAI linear
index.  Consider data within a single frame:

                                           |---| query
      --------------------------------
                       -----------------------------------------
      ...   ...   ...   ...   ...   ...   ...   ...   ...   ...  
        ...   ...   ...   ...   ...   ...   ...   ...   ...   ...  
          ...   ...   ...   ...   ...   ...   ...   ...   ...   ...  

                       ^                  ^
                       data offsets

This returns 1 long read and several short reads.  As we've stored
data sorted by the left end, this means our region has two offsets
within the frame that need to be decoded which are separated by data
we will discard.

If we implement this, it should be in the per-frame meta-data block,
so it's distributed and we only pay the cost of loading the per-frame
index for frames that we decode.  A sparse map of poosition to offset
within the uncompressed frame would suffice.

TODO: Similar to above, we can also have an index per frame that
records chromosomes listed as a way of handling multiple small
references or the junction from migrating from the end of one
chromosome to the start of the next.

[Both of these above issues requires our iterator to cache a copy of
the region(s), as an index query alone won't have sufficient data to
perform this decode optimisartion]


Seekable index frame
--------------------

We have one seekable-format index frame as dodcumented in

https://github.com/facebook/zstd/blob/dev/contrib/seekable_format/zstd_seekable_compression_format.md

This is always the last frame in the file.

The "seekable" index frame consists of the standard skippable frame
header (using the magic number 0x184D2A5E and a length holding the
remaining size of the skippable frame), a series of index entries
(compressed size, uncompressed size and an optional CRC) and a footer
with the number of index entries, a flag byte, and ending with a
magic number (0x8F92EAB1).

The seekable index trailing magic number is also used as an EOF detector.
