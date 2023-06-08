TODO
----

- Consider a small meta-data skippable frame at the start.
  - This can explicitly label it as BGZF2 and not just arbitrary pzstd
    output.
  - It can help inform the decoder of the content type, by replicating
    the first e.g. 64 bytes of the data stream verbatim.

- Index.
  - This needs to go after the last zstd compressed block and before
    the seekable index.
  - Track genomic chr:start-end range (eg as in CRAM)
    - Ref 1-N or 0 for unmapped?  Or 0-N and -ve?
    - Also need to consider multi-ref frames, for when we get
      many small contigs.  So multiple index entries per frame.
  - Track number of sequences, maybe also bases.
    - Permits crude index based depth analysis.
    - Also permits easy data segmentation, even on unsorted or
      unaligned data.
  - R-tree?  Is it still helpful for a mix of wildly different seq lengths?
    Probably. (CRAM misses this)
    - Long reads indicate which frames have data spanning a query region
    - Short reads may need skipping within a frame.

                                           |---|
      --------------------------------
                       -----------------------------------------
      ...   ...   ...   ...   ...   ...   ...   ...   ...   ...  
        ...   ...   ...   ...   ...   ...   ...   ...   ...   ...  
          ...   ...   ...   ...   ...   ...   ...   ...   ...   ...  

                       ^                  ^

      Optimally, two key points to seek to within this frame.
    - Possibly multiple overlapping indices covering different length
      data (via RG -> LB?).
      - Tracks start locations for long reads
      - May be "good enough".


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

[bgzf2 header skippable frame]
[pzstd skippable frame]
[zstd data frame]
[pzstd skippable frame]
[zstd data frame]
...
[bgzf2 index skippable frame]
[seekable_format skippable frame (index)]

Extensive use of zstd skippable frames are used.  These are ZSTD
compliant frames which carry no compressed payload, so can be skipped
over by any ZSTD decoder.  Instead they carry meta-data useful for
file type detection, rapid parallel decoding, and internal indices.

In order for standard zstd tools to be able to skip past these custom
meta-data frames, the skippable frames have a specific data layout:

- 4 bytes of magic number: 0x184D2A50 to 0x184D2A5F
- 4 bytes of remaining frame length "N"
- N bytes of meta-data


This combination of skippable frames means the following tools can
decompress a bgzf2 file:
bgzf2 data:

- "zstd".
  This provides single threaded streaming decode.

- "pzstd", from Zstd's contrib/pzstd directory.
  This provides a parallel decompression capability.

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

[TODO: this breaks pzstd detection.  Maybe we should make the initial
pzstd frame bigger and merge in with this?  Does that work?  It looks
likely not as it only reads 12 bytes in contrib/pzstd/SkippableFrame.cpp]

The header frame has magic number 0x184D2A5B.
The data contents are an uncompressed copy of the first data bytes,
used for file format detection.  There is no fixed limit on this
length, but it is recommended to be not significantly more than is
required to accurately determine file type and some basic versioning.

An example header frame:

 4: 0x184D2A5B (bgzf2 magic number)
 4: 23 (length of meta-data in header)
 4: "BGZ2" identifier
19: "BAM\x01????@HD.VN:1.4.SO:coordinate\n"

Tools should be capable of working with reduced meta-data here, for
example having just "BAM\x01".


PZstd skippable frames
----------------------

The pzstd skippable frames hold the size of the next compressed zstd
data frame (as this is not part of the zstd format).  They are 12
bytes long, consisting of 3 little-endian values.

- 0x184D2A50: the magic number for pzstd's skippable frame
- 4: the length of the remaining skippable frame data
- comp_sz: the compressed size of the next data frame

We have one pzstd skippable frame preceeding each zstd compressed data
frame.

TODO: can we augment this with additional data, such as the region
being used (eg see the CRAM container struct).  This would provide a
way to gather information about the upcoming frame as we're streaming,
without needing the index.  This in turn can lead to efficient partial
decode where we skip over data that's going to be filtered out as it
doesn't match a region filter.  (An index and random access is
preferable, but not always feasible.)

If not supported by pzstd, it may still be useful to do and just break
from the pzstd standard.  Should use a generic key-value pair
mechanism, with BGZF2 indicating specific keys in use.  That makes
this format a generic one for any data type.

Examples: count records (mapped and unmapped), skip data outside of a
region, distributed processing to turn one BGZF2 to multi sub-BGZF2
without any decompression and recompression.


BGZF2 index frame
-----------------

BGZF indices can be BAI, CSI and TBI.  They are also external files.
This myriad of index formats leads to a variety of problems:

- The naming differs between tools, sometimes it is foo.bam.bai and
  sometimes foo.bai.

- If multiple indices exist (.bai and .csi) there are undocumented
  precedence rules.

- Multiple files causes problems when downloading from object stores
  where the filenames are content hashes.

- There is potential for catastrophy when the bam file is rewritten
  without the index being recreated.  Indices do not have time stamps
  or shared secrets with their BAM counterparts, making it hard to
  detect this problem.

- These indexes are only for aligned data, so there is no index
  capability on unmapped records, name-sorted or unsorted data.  Some
  of these are covered by yet another index format: GZI.

So in BGZF2 the index is embedded within the file. It is also
partially distributed. [TODO]

The main global index is at the end of the file.  The purpose of our
index is two fold:

- To map genomic coordinates to file offsets, for random access.

- For parallel processing by splitting data into discrete work units,
  regardless of whether it has been aligned and sorted.

  This latter use case is already covered by the seekable index. So is
  ignored here.

For genomic sorted data, we are given a region and turn it into one or
more zstd frames that hold data covering that region.  Minimally, we
must do a zstd decompression of either an entire frame, or from the
start of a frame until we have finished the region query.  Hence the
global index simply needs to map region to frame offset, or region to
uncompressed offset (and use the seekable index to convert).

(If we're using uncompressed offsets, our resolution of genomic region
queries do not need to match the resolution of zstd data frames.)

Prior to each compressed data frame we have an uncompressed meta-data
frame, currently pzstd and holding only the next frame size.  This
will be extended to also hold arbitrary meta-data such as chromosome
and range.  This permits a streaming mode where we skip data.  We
still have the I/O requirement, but we do not need to decompress data
if it doesn't match our desired range.

Furthermore, we can also provide additional genomic indexing within
the uncompressed frame, so for example a multiple-reference frame may
list the offset within the uncompressed frame for each new chromosome,
or perhaps every 16kb into that chromosome.  This is like the BAI
linear index, but instead it can be stream inline with the data.  This
distributed nature makes the indexing capabilities more efficient when
seeking is unavailable.

Hence in this section, we can focus on a simplistic index capability.
How to map a genomic range to the start of a zstd frame.

[TODO]

This is a genomic coordinate index used for random access by
chromosome and position sorted data, or by record number if unsorted.
The index itself is self-indexing, meaning it is possible to hop
around within it via precomputed offsets.  In practice however it is
likely the most performant use is to load the entire index into
memory, especially if compressed.

This is another skippable frame with the same magic numebr as the
bgzf2 header.

 4: 0x184D2A5B (bgzf2 magic number)
 4: N+1 (length of meta-data in header)
 1: index flags
 N: An array of index entries
 8: Footer (another magic number and size, to permit reverse reading)

Index flags use the following bit-field

Bit 0:    1 if the remaining bytes are zstd compressed, 0 if uncompressed
Bits 1-7: 0 (reserved)

The remaining N bytes may be zstd compressed.  Once decompressed, the
format of the index is as follows.

[ Should we use variable sized integer encoding instead?  It's a bit
more faffing, but not too much and it makes things totally data size
agnostic. However it makes decoding also trickier, especially random
access to within an index.]

1: global flag
   Bit 0: 1 if aligned, 0 otherwise
   Bit 1: 1 if genomic sorted 
   Bit 2-7: 0 (reserved)
4: number of references NR
4: number of index entries NI (matching number of compressed data frames)
4: size of meta-data MDG
MDG: meta-data

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

Global and per-ref meta-data is key\0value\0 pairs.  It replaces the
magic bin numbers used in BAI and CSI.  We can view it aggregated
together for the entire file, or on a per-reference basis.

TODO: add a controlled dictionary here.  Eg nrec_mapped,
nrec_unmapped, nbase_mapped, nbase_unmapped, etc.

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

The seekable index trailing magic number is used as an EOF detector.
