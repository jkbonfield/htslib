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

[bgzip2 header skippable frame]
[pzstd skippable frame]
[zstd data frame]
[pzstd skippable frame]
[zstd data frame]
...
[bgzip2 index skippable frame]
[seekable_format skippable frame (index)]

The first zstd data frame should be end within the first 4KB of the
file, so that htslib hpeek() function can provide enough data to
decompress the frame and determine the internal format.

[TODO: consider an additional meta-data skippable frame that
duplicates this, avoiding this peek issue?]

The pzstd skippable frames hold the size of the next compressed zstd
data frame (as this is not part of the zstd format).  They are 12
bytes long, consisting of 3 little-endian values.

- 0x184D2A50: the magic number for pzstd's skippable frame
- 4: the length of the remaining skippable frame data
- comp_sz: the compressed size of the next data frame.

The "seekable" index frame consists of the standard skippable frame
header (using the magic number 0x184D2A5E and a length holding the
remaining size of the skippable frame), a series of index entries
(compressed size, uncompressed size and an optional CRC) and a footer
with the number of index entries, a flag byte, and ending with a
magic number (0x8F92EAB1).

The Seekable format skippable frame is documented in
https://github.com/facebook/zstd/blob/dev/contrib/seekable_format/zstd_seekable_compression_format.md

The seekable index trailing magic number is used as an EOF detector.
