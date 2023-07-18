# SAM API

## HTSLib, SAM APIs and samtools

HTSLib is a C library implementation used to access and process the genome sequence data. HTSLib implements multiple API interfaces, HTS API, VCF API and SAM API.
HTS API provides a framework for use by other APIs and applications, implements bgzf compression, htscodecs and provides CRAM format support.
VCF APIs work with variant data in VCF and BCF format.

SAM API works with sequence data of different formats, SAM / BAM / CRAM / FASTA / FASTQ, and provides methods to do operations on the data. It uses methods from HTS API.

‘samtools’ is the utility used to read and modify sequence data. It uses SAM APIs from HTSLib to work on the sequence data.

## About this document

There are a number of demonstration utilities and their source code in ‘samples’ directory of HTSLib and this document gives the description of them and the usage of SAM API of HTSLib.
The samples are for demonstration purposes only and proper error handling is required for actual usage. This document is based on HTSLib version 1.17.

## The sample apps

Flags - This application shows the basic read of alignment files and flag access. This utility shows the alignment count against different flags in it.

Split - This application shows the basic read and write of alignment data. This utility saves the read1 and read2 as separate files in a given directory.

Read_header - This application shows the read of header line fields.

Update_header - This application shows the update of header line fields, where update is allowed.

Read_bam - This application shows how the alignment data fields are accessed.

Read_aux - This application shows how the auxiliary data are read from bam data.

Flags_tpool - This application is similar to the ‘flags’ and uses HTSLib’s thread pool facility.

Flags_field - This application is similar to the ‘flags’ and better, as it reads only the selected fields it needs, with CRAM file types.

Idx_on_write / Idx_by_read - These applications show how to create index while creating sequence data and on existing files.

Read_with_idx / read_multi_index - These applications show how to read data using index file and numerical position / positions.

Read_reg / read_multireg - These applications show how to read data using index file and interested region / regions through textual description.

Pileup / mpileup  - These applications show how to have the sequence data in transposed fashion from a usual SAM style representation and lets access each read and position for different processing.


## Building the sample apps

The samples expect the HTSLib is installed, libraries and header file path are part of the PATH environment variable. If not, these paths need to be explicitly passed during the build time.

Gcc and compatible compilers can be used to build the samples.

Along with the htslib and gcc, other libraries and/or headers required to build are, math, pthread, curl, lzma, z and bz2 libraries.

The sample alignment data file present is in SAM format and it need to be sorted and saved as compressed file, to be used with the samples.

## Usage of SAM APIs
### Sequence data file access for read

The sequence data file for read may be opened using the sam_open method. 
It opens the file and returns samFile (htsFile) pointer on success or NULL on failure. 
The input can be path to a file in disk, network, cloud or “-” designating the standard input.

SAM, BAM and CRAM file formats are supported and the input file format is detected from the file content.

Once done with the file, it needs to be closed with sam_close.

Many times, header details would be required and can be read using sam_hdr_read api.
It returns sam_hdr_t pointer or NULL. The returned header needs to be destroyed using sam_hdr_destroy when no longer required.

The sequence data may be compressed or uncompressed on disk and on memory it is read and kept as uncompressed BAM format. It can be read from a file using sam_read1 api. samFile pointer, header and bam storage are to be passed as argument and it returns 0 on success, -1 on end of file and < -1 in case of errors.

The bam storage has to be initialised using bam_init1 api before the call and can be reused for successive reads. Once done, it needs to be destroyed using bam_destroy1.
The member field named core - bam1_core_t - in bam storage, bam1_t, has the sequence data in an easily accessible way. Using the fields and macros, data can easily be read from it.

    #include <htslib/sam.h>
    
    int main(int argc, char *argv[])
    {
    
       //initialise
       bamdata = bam_init1();
    
    
       //open input files
       infile = sam_open(inname, "r");
    
    
       //read header
       in_samhdr = sam_hdr_read(infile);
    
    
       //read data, check flags and update count
       while (0 <= (c = sam_read1(infile, in_samhdr, bamdata))) {
           //check flag and update count
           if (bamdata->core.flag & BAM_FQCFAIL) {
               flgcnt.qcfailed++;
           }
           else {
               if (bamdata->core.flag & BAM_FPAIRED)
                   flgcnt.paired++;
               .
           }
       }
       if (-1 == c) {
           //EOF
           printf("File %s have %ld alignments\n \
           
           inname, cnt, flgcnt.paired, flgcnt.proppaired, flgcnt.unmapped, flgcnt.mateunmapped, flgcnt.reveresed, flgcnt.matereversed, flgcnt.read1, flgcnt.read2,
           flgcnt.secondary, flgcnt.supplementary, flgcnt.duplicate, flgcnt.qcfailed);
    
    
           ret = EXIT_SUCCESS;
       }
       //else not fully read
    end:
       //clean up
       if (in_samhdr)
           sam_hdr_destroy(in_samhdr);
       if (infile)
           sam_close(infile);
       if (bamdata)
           bam_destroy1(bamdata);
       return ret;
    }

Refer: flags_demo.c

This shows the count of alignments against different flags.

    ./flags -i /tmp/sample.sam.gz

### Sequence data file access for write

File access for write is similar to read with a few additional steps.

sam_open_mode method sets the output file type and compression based on the filename. This method expects a buffer to append type and compression flags. Usually a buffer with standard file open flag is used, the buffer past the flag is passed to the method to ensure existing flags and updates from this method are present in the same buffer without being overwritten. This method will add more flags indicating file type and compression based on name.

sam_open_format method may also be used to open the file for output as more information on the output file can be specified using this.

The header data can be written using the sam_hdr_write api. When the header data is copied to another variable and has different lifetime, it is good to increase the reference count of the header using sam_hdr_incr_ref and sam_hdr_dstroy called as many times as required.

The alignment data can be written using the sam_write1 api. It takes a samFile pointer, header pointer and the alignment data. The header data is required to set the reference name in the alignment. It returns -ve value on error.

    int main(int argc, char *argv[])
    {   
       outfile1 = sam_open(file1, "w");`
       outfile2 = sam_open(file2, "w");
       
       //read header, required to resolve the target names to proper ids
       in_samhdr = sam_hdr_read(infile);
       
       //write header
       if (-1 == sam_hdr_write(outfile1, in_samhdr) || -1 == sam_hdr_write(outfile2, in_samhdr)) {
       
       while (0 <= (c = sam_read1(infile, in_samhdr, bamdata))) {
           cnt++;
           if (bamdata->core.flag & BAM_FREAD1) {
               if (0 > sam_write1(outfile1, in_samhdr, bamdata)) {
                   fprintf(stderr, "Failed to write output data\n");
       
       if (-1 == c) {	
           //EOF
           fprintf(stderr, "Done");
       
       if (outfile1)
           sam_close(outfile1);
   
Refer: split.c

This creates 1.sam and 2.sam in /tmp/ containing read1 and read2 respectively.

    ./split -i /tmp/sample.sam.gz -o /tmp/

### Header data read

The header gives the version, reference details, read group, change history and comments. These data are stored inside the sam_hdr_t. Each of these entries, except comments, have their unique identifier and it is required to access different fields of them.
The api sam_hdr_count_lines gives the count of the specified type of header line.
The value of a unique identifier to a specific type of header line can be retrieved with sam_hdr_line_name api.
The api sam_hdr_find_tag_id and sam_hdr_find_tag_pos can get the field data from a header line. The full header line can be retrieved using sam_hdr_find_line_pos.

    int main(int argc, char *argv[])
    {  
       infile = sam_open(inname, "r");
       
       //read header
       in_samhdr = sam_hdr_read(infile);
       
           //get count of given header type
           linecnt = sam_hdr_count_lines(in_samhdr, *(header + c));
           if (linecnt) {
               //get unique identifier field's value for given type
               id = sam_hdr_line_name(in_samhdr, *(header + c), pos);
               if (tag) {
                   //get value of given tag for given unique identifier
                   sam_hdr_find_tag_id(in_samhdr, *(header + c), *(idf + c), id, tag, &data);
               }
               else {
                   //get line at given position
                   sam_hdr_find_line_pos(in_samhdr, *(header + c), pos, &data);
               }
        
       }

Refer: read_header.c

This will show the VN tag’s value from HD header.

    ./read_header -i /tmp/sample.sam.gz -h HD -n 0 -t VN

Shows the 2nd SQ line’s LN field value.

    ./read_header -i /tmp/sample.sam.gz -h SQ -n 1 -t LN

### Header data update

The unique identifier for the line needs to be found to update a field in line, though not all types in the header may be modifiable.
The api sam_hdr_update_line takes the unique identifier for the header line type, its value, the field which needs to be modified and the value with which to modify it, followed by a NULL.
Eg. To change LN field from 2000 to 2250 in SQ line with unique identifier SN as “chr1”, sam_hdr_update_line( header, “SQ”, “SN”, “chr1”, “LN”, “2250”, NULL).
To change PP field from ABC to DEF in PG line with ID APP.10,
sam_hdr_update_line( header, “PG”, “ID”, “APP.10”, “PP”, “DEF”, NULL).

    int main(int argc, char *argv[])
    {  
       //open / input files
       infile = sam_open(inname, "r");
       
       outfile = sam_open(outname, "w");
       
       //read header
       in_samhdr = sam_hdr_read(infile);
       
       //count of lines of given type  
       linecnt = sam_hdr_count_lines(in_samhdr, header);
       if (pos >= linecnt) {
       
       //get unique identifier value for the required line
       id = sam_hdr_line_name(in_samhdr, header, pos);
       //update with new data
       if (0 > sam_hdr_update_line(in_samhdr, header, idf, id, tag, val, NULL)) {
       
       if (0 > sam_hdr_write(outfile, in_samhdr)) {
   
Refer: update_header.c

Saves new sam file with 2nd SQ line having length as 38.

    ./update_header -i /tmp/sample.sam.gz -o /tmp/sam.sam -h SQ -n 1 -t LN -v 38


### Alignment data read

The alignment / sequence data contains many fields. Mainly the read/query name, flags indicating the properties of the read, reference sequence name, position in reference to which it matches, quality of the read, CIGAR string indicating the match status, position of mate / reverse strand, name of reference sequence to which mate matches, the insert length, base sequence, quality value of each base and auxiliary fields.

Header data would be required to retrieve the reference names as alignment contains the position of the reference in the header.

A few of the data are directly visible in bam1_t and the rest are hidden inside data member of bam1_t and can easily be retrieved using macros.
bam_get_qname gives the name of the read, sam_hdr_tid2name gives the reference name. bam_get_cigar retrieves the cigar operation array, which can be decoded using bam_cigar_oplen to get count of bases to which that operation applicable and bam_cigar_opchr to get the cigar operation.
bam_seqi retrieves the base data at a given position in alignment and it can be converted to character by indexing the seq_nt16_str array.

    int main(int argc, char *argv[])
    {    
       while (0 <= (ret_r = sam_read1(infile, in_samhdr, bamdata)))	
       {   
           //QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAG:TYPE:VALUE]
           fprintf(stdout, "%s", bam_get_qname(bamdata));                                  //get the query name using the macro
           fprintf(stdout, "\t%d", bamdata->core.flag);                                    //flag is available in core structure
           fprintf(stdout, "\t%s", sam_hdr_tid2name(in_samhdr, bamdata->core.tid));        //retrieves the target name using the value in bam and by referring the header
           fprintf(stdout, "\t%d", bamdata->core.pos + 1);                                 //internally position is 0 based and on text output / SAM it is 1 based
           fprintf(stdout, "\t%d", bamdata->core.qual);                                    //map quality value
           cigar = bam_get_cigar(bamdata);                                                 //retrieves the cigar data
           fprintf(stdout, "\t");
           for (i = 0; i < bamdata->core.n_cigar; ++i) {                                   //no. of cigar data entries
               fprintf(stdout, "%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));  //the macros gives the count of operation and the symbol of operation for given cigar entry
           }
           fprintf(stdout, "\t%s", bamdata->core.mtid == bamdata->core.tid ? "=" : sam_hdr_tid2name(in_samhdr, bamdata->core.mtid));   // = if rnext and pnext are same and the target name otherwise
           fprintf(stdout, "\t%d", bamdata->core.mpos+1);
           fprintf(stdout, "\t%d", bamdata->core.isize);
           data = bam_get_seq(bamdata);                                                    //get the sequence data
           if (bamdata->core.l_qseq != bam_cigar2qlen(bamdata->core.n_cigar, cigar)) {     //checks the length with CIGAR and query
    	
           for (i = 0; i < bamdata->core.l_qseq ; ++i) {                                   //sequence length
               fprintf(stdout, "%c", seq_nt16_str[bam_seqi(data, i)]);                     //retrieves the base from (internal compressed) sequence data
           }
           fprintf(stdout, "\t");
           for (int i = 0; i < bamdata->core.l_qseq ; ++i) {
               fprintf(stdout, "%c", bam_get_qual(bamdata)[i]+33);                         //retrives the quality value
       
Refer: read_bam.c

Shows the 2nd alignment data fields on screen.

    ./read_bam -i /tmp/sample.sam.gz -n 1


### Aux data read

Auxiliary data gives more information about the alignment. There can be a number of such data and can be accessed by iterating one by one through then once the alignment is read as bam1_t. The auxiliary data are stored along with the variable length data in the data field of bam1_t. There are macros defined to retrieve information about auxiliary data from the data field of bam1_t.

The start of aux data is retrieved using macro bam_aux_first and successive ones using bam_aux_next. Macro bam_aux_tag gives the tag of the aux field and bam_aux_type gives the information about type of the aux field.

Bam_aux2i, bam_aux2f,bam_aux2Z macros retrieve the aux data’s value as integer, float and string respectively. The integer value may be of different precision / size and the bam_aux_type character indicates how to use the value. The string/hex data are NULL terminated.

For array data, bam_aux_type will return ‘B’ and bam_auxB_len gives the length of the array. bam_aux_type with the next byte will give the type of data in the array. bam_auxB2i, bam_auxB2f will give integer and float data from a given position of the array.

    int printauxdata(FILE *fp, char type, int32_t idx, const uint8_t *data)
    {  
       switch(type) {
       case 'A':
           fprintf(fp, "%c", bam_aux2A(data));                                                 //byte data
           break;
       case 'c':
           fprintf(fp, "%d", (int8_t)(idx > -1 ? bam_auxB2i(data, idx) : bam_aux2i(data)));    //signed 1 byte data; bam_auxB2i - from array or bam_aux2i - non array data
        
       case 'f':
       case 'd':
           fprintf(fp, "%g", (float)(idx > -1 ? bam_auxB2f(data, idx) : bam_aux2f(data)));     //floating point data, 4 bytes
           break;
       case 'H':
       case 'Z':
           fprintf(fp, "%s", bam_aux2Z(data));                                                 //array of char or hex data
           break;
       case 'B':                                                                               //array of char/int/float
           auxBcnt = bam_auxB_len(data);                                                       //length of array
           auxBType = bam_aux_type(data + 1);                                                  //type of element in array
           fprintf(fp, "%c", auxBType);
           for (i = 0; i < auxBcnt; ++i) {                                                     //iterate the array
               fprintf(fp, ",");
               printauxdata(fp, auxBType, i, data);                                            //calling recurssively  with indexto reuse a few lines
    
    }
    int main(int argc, char *argv[])
    {  
       while (0 <= (ret_r = sam_read1(infile, in_samhdr, bamdata)))
       {
       
           errno = 0, i = 0;
           data = NULL;
           data = bam_aux_first(bamdata);                                              //get the first aux data
           for (i = 0; data; ++i) {
               fprintf(stdout, "%.2s:%c:", bam_aux_tag(data), NULL != strchr("cCsSiI", bam_aux_type(data)) ? 'i' : bam_aux_type(data));  //macros gets the tag and type of aux data
               printauxdata(stdout, bam_aux_type(data), -1, data);                     //dump the data
               fprintf(stdout, "\t");
               data = bam_aux_next(bamdata, data);                                     //get the next aux data
           }
           if (ENOENT != errno) {

Refer: read_aux.c

Shows the aux tags from 2nd alignment line from the file.

    ./read_aux -i ../../samtools/test/mpileup/mpileup.1.bam -n 1


### Threadpool to read / write

The HTSLib api provides a thread pool for better performance. The number of threads that need to be used can be configured and the thread pool can be associated with required samFile handles that the operation on those files will be carried out using the threadpool. The thread pool can be associated with different file handles at the same time.
The hts_tpool_init api initialises the thread pool, it takes the number of threads in the pool.
The hts_set_opt api takes the samfile pointer and thread pool address and associates the file with thread pool. Further operations in the file will be made through the thread pool.
Once done with the pool, it needs to be destroyed using the hts_tpool_destroy method.
Utilities like ‘time’ will show the change in performance with use of thread pool.

    int main(int argc, char *argv[])
    { 
       htsThreadPool tpool = {NULL, 0};
       
       //initialise bam data holder
       bamdata = bam_init1();
       
       //open input files
       infile = sam_open(inname, "r");
       
       //set 4 threads for thread pool
       if (!(tpool.pool = hts_tpool_init(4)))
       {
           fprintf(stderr, "\nFailed to setup thread pool\n");
           goto end;
       }
    
    
       //map threadpool with infile
       hts_set_opt(infile, HTS_OPT_THREAD_POOL, &tpool);
    
    
       //read header
       in_samhdr = sam_hdr_read(infile);
       
       //check flags and update count
       while (0 <= (c = sam_read1(infile, in_samhdr, bamdata))) {
           cnt++;
           if (bamdata->core.flag & BAM_FQCFAIL) {
               flgcnt.qcfailed++;
           }
           
       if (tpool.pool)
           hts_tpool_destroy(tpool.pool);
   
Refer: flags_tpool.c

### Read selected fields

At times the whole alignment data may not be of interest and it would be better to read required fields alone from the alignment data. CRAM file format supports such specific field data read and HTSLib provides an option to use this. This can improve the data read operation.
The hts_set_opt method does the selection of specified fields. There are flags indicating specific fields, , like SAM_FLAG, SAM_SEQ, SAM_QNAME, in alignment data and a combination of flags for the required fields shall be passed with CRAM_OPT_REQUIRED_FIELDS to this api.

    int main(int argc, char *argv[])
    {  
       bamdata = bam_init1();
       
       //open input files
       infile = sam_open(inname, "r");
       
       //select required field alone, this is useful for CRAM alone
       if (0 > hts_set_opt(infile, CRAM_OPT_REQUIRED_FIELDS,SAM_FLAG)) {
           fprintf (stderr, "Failed to set htsoption\n");
           goto end;
       }
       //read header
       in_samhdr = sam_hdr_read(infile);
       
       //check flags and update count
       while (0 <= (c = sam_read1(infile, in_samhdr, bamdata))) {
           cnt++;
           if (bamdata->core.flag & BAM_FQCFAIL) {
               flgcnt.qcfailed++;
           }

   
Refer: flags_htsopt_field.c

### Create an index

Indexes help to read data faster without iterating sequentially through the file. Indexes contain the position information about alignments and that they can be read easily. They are usually used with iterators.

Indexing of plain SAM files is not supported. Compressed SAM, BAM and CRAM files can be indexed. CRAM files are indexed as .crai and the other two can be indexed as .bai or .csi files. Each of these types have different internal representations of the index information. Bai uses a fixed configuration values where as csi has them dynamically updated based on the alignment data.

Indexes can be created either with the alignment data save or explicitly by reading existing alignment data file.

To create index along with alignment write, the sam_idx_init api need to be invoked before the start of alignment data. This api takes the output samFile pointer, header pointer, minimum shift and index file path. For BAI index, the min shift has to be 0. 

At the end of write, sam_idx_save api need to be invoked to save the index.

    int main(int argc, char *argv[])
    {  
       if (sam_hdr_write(outfile, in_samhdr)) {
       
       // initialize indexing, before start of write
       if (sam_idx_init(outfile, in_samhdr, size, fileidx)) {
           fprintf (stderr, "idx initialization failed\n");
           goto exit;
       }
       //read and write alignments
       while (0 <= (c = sam_read1(infile, in_samhdr, bamdata))) {
           if (0 > sam_write1(outfile, in_samhdr, bamdata)) {
               fprintf(stderr, "Failed to create gz file\n");
               goto exit;
           }
           cnt++;
       }
       if (-1 == c) {
           //EOF, save index
           if (sam_idx_save(outfile)) {
               fprintf(stderr, "Could not save index\n");
               goto exit;
           }
           ret = EXIT_SUCCESS;
   
Refer:index_write.c

Creates mpileup.1.bam and mpileup.1.bam.bai in /tmp/.

    ./idx_on_write -i ../../samtools/test/mpileup/mpileup.1.bam -o /tmp

To create index explicitly on an existing alignment data file, the sam_index_build api or its alikes can be used. sam_index_build takes the alignment file path, min shift for the index and creates the index file in same path. The output name will be based on the alignment file format and min shift passed.
The sam_index_build2 api takes the index file path as well and gives more control than the previous one.
The sam_index_build3 api provides an option to configure the number of threads in index creation.

    int main(int argc, char *argv[])
    {  
       if (!outdir) {
           //have the index in same location, with similarly named
           if (sam_index_build(inname, size)) {
               fprintf (stderr, "idx creation failed\n");
               goto exit; 
           }
       } else {
           //have the index in given location
        
           // build index by explicitly naming index file
           if (sam_index_build2(inname, fileidx, size)) {
               fprintf (stderr, "idx creation failed\n");
               goto exit;
           }
   
Refer:index_create.c

Creates /tmp/mp.bai.

    ./idx_by_read -i /tmp/mpileup.1.bam -n /tmp/mp.bai

### Read with iterators

Index file helps to read required data without sequentially accessing the file and are required to use iterators. The interested reference, start and end position etc. are required to read data with iterators. With index and these information, an iterator is created and relevant alignments can be accessed by iterating through it.

The api sam_index_load and it alikes does the index loading. It takes input samFile pointer and file path. It loads the index file based on the input file name, from the same path and with implicit index file extension - cram file with .crai and others with .bai. The sam_index_load2 api accepts explicit path to index file, which allows loading it from a different location and explicit extensions. The sam_index_load3 api supports download/save of the index locally from a remote location. These apis returns NULL on failure and index pointer on success.

The sam_iter_queryi or sam_iter_querys apis may be used to create an iterator and sam_itr_next api does the alignment data retrieval. Along with retrieval of current data, it advances the iterator to next relevant data. The sam_iter_queryi takes the interested positions as indexes and sam_iter_querys takes the interested position as a string.

With sam_iter_queryi, the reference id can be the 0 based index of reference data, -2 for unmapped alignments, -3 to start read from beginning of file, -4 to continue from current position, -5 to return nothing. Based on the reference id given, alignment covering the given start and end positions will be read with sam_iter_next api.

    int main(int argc, char *argv[])
    {  
       hts_idx_t *idx = NULL;
       hts_itr_t *iter = NULL;
       
       //load index file
       if (idxfile) {
           //index file provided
           idx = sam_index_load2(infile, inname, idxfile);
       }
       else {
           //assume it to be present in same location
           idx = sam_index_load(infile, inname);
       }   
       
       //create iterator
       if (!(iter = sam_itr_queryi(idx, tid, start, end))) {
           fprintf(stderr, "Failed to get iterator\n");
           goto exit;
       }
       //read using iterator
       while ((c = sam_itr_next(infile, iter, bamdata) >= 0)) {
       
       if (iter)
           sam_itr_destroy(iter);
       if (idx)
           hts_idx_destroy(idx);
   
Refer: index_read.c

With sample.sam, reference -2 will show alignments with name UNMAP2 and UNMAP3
With reference -3, it shows all alignments

    ./read_with_idx  /tmp/sample.sam.gz -t -2 
    ./read_with_idx  /tmp/sample.sam.gz -t -3

With reference 0, start 1 and end 4 it shows nothing and with start 1 and end 5 it shows alignment with name ITR1.

    ./read_with_idx  /tmp/sample.sam.gz -t 0 -s 1 -e 5

With reference 1, start 30 and end 100, it shows alignment with name ITR2M which refers the 2nd reference data 

    ./read_with_idx  /tmp/sample.sam.gz -t 1 -s 30 -e 100


With sam_iter_querys, the reference sequence is identified with the name and interested positions can be described with start and end separated by ‘-’ as string. When sequence is identified as ‘.’, it begins from the start of file and when it is ‘*’, unmapped alignments are read. Reference with <name>[:], <name>:S, <name>:S-E, <name>:-E retrieves all data, all data covering position S onwards, all data covering position S to E, all data covering upto position E of reference with ID <name> respectively on read using sam_iter_next.

       //load index file
       if (idxfile) {
           //index file provided
           idx = sam_index_load2(infile, inname, idxfile);
       }
       else {
           //assume it to be present in same location
           idx = sam_index_load(infile, inname);
       }   
       
       //create iterator	
       if (!(iter = sam_itr_querys(idx, in_samhdr, region))) {
           fprintf(stderr, "Failed to get iterator\n");
           goto exit;
       }
       //read using iterator
       while ((c = sam_itr_next(infile, iter, bamdata) >= 0)) {
      
       if (iter)
           sam_itr_destroy(iter);
       if (idx)
           hts_idx_destroy(idx);
 


Refer:index_reg_read.c

With sample.sam, region as \* will show alignments with name UNMAP2 and UNMAP3

    ./read_reg -i /tmp/sample.sam.gz -r \*

With region as \., it shows all alignments

    ./read_reg -i /tmp/sample.sam.gz -r \.

With region as T1:1-4, start 1 and end 4 it shows nothing and with T1:1-5 it shows alignment with name ITR1.

    ./read_reg -i /tmp/sample.sam.gz -r T1:1-5

With region as T2:30-100, it shows alignment with name ITR2M which refers the reference data T2.

    ./read_reg -i /tmp/sample.sam.gz -r T2:30-100

The iterator and index need to be destroyed once done with it, using sam_itr_destroy and hts_idx_destroy apis.	 

### Read multiple regions

Multiple interested regions can be specified for read using apis sam_iter_regions and sam_itr_regarray.
sam_itr_regions takes index file path, header, count of regions and an array of region description. The region description includes the reference indicator string, count of regions interested in this reference and an array of start end positions. sam_itr_multinext is used to get data through this method. 

    int main(int argc, char *argv[])
    { 
       hts_idx_t *idx = NULL;
       hts_itr_t *iter = NULL;
       hts_reglist_t *reglist = NULL, *tmp = NULL;
    
       //load index file
       //assume it to be present in same location
       idx = sam_index_load(infile, inname);
    
       if (!(iter = sam_itr_regions(idx, in_samhdr, reglist, regcnt))) {
     
       //read using iterator
       while ((c = sam_itr_multi_next(infile, iter, bamdata) >= 0)) {
    
       if (iter)
           sam_itr_destroy(iter);
       if (idx)
           hts_idx_destroy(idx);



Refer:index_multi_read.c

With compressed sample.sam and 2 regions from reference T1 (1-5 and 30-100) and 1 region from T2 (30-100), alignments with ITR1, A2, ITR1M and ITR2M will be shown.

    ./read_multi_index -i /tmp/sample.sam.gz -N 2 -d T1,2,1,5,30,100 -d T2,1,30,100


    int main(int argc, char *argv[])
    { 
       hts_idx_t *idx = NULL;
       hts_itr_t *iter = NULL;
       char **regions = NULL;
       
       //load index file
       //assume it to be present in same location
       idx = sam_index_load(infile, inname);
       
       if (!(iter = sam_itr_regarray(idx, in_samhdr, regions, regcnt))) {
           fprintf(stderr, "Failed to get iterator\n");
           goto exit;
       }
    
       //read using iterator
       while ((c = sam_itr_multi_next(infile, iter, bamdata) >= 0)) {
    
       if (iter)
           sam_itr_destroy(iter);
       if (regions)
           free(regions);
    
      if (idx)
           hts_idx_destroy(idx);

Refer:index_multireg_read.c

With compressed sample.sam and 2 regions from reference T1 (30 to 32) and 1 region from T2 (30 onwards), alignments with name A1, B1, A2 and ITR2M would be shown.

    ./read_multireg -i /tmp/sample.sam.gz -N2 -d T1:30-32,T2:34

The index and iterators are to be destroyed using the sam_itr_destroy and hts_idx_destroy. The hts_reglist_t* array(reglist) passed is destroyed by the library on iterator destroy. The regions array need to be destroyed by the user itself.

### Pileup and MPileup

Pileup shows the transposed view of the SAM alignment data, i.e. it shows the the reference positions and bases which cover that position through different reads. MPileup facilitates the piling up of multiple sam files against each other and same reference at the same time.

Mpileup has replaced the pileup. The input expects the input to be sorted by position.

Pileup needs to be initialized with bam_pileup_init method which takes pointer to a method, which will be called by pileup to read data from required files, and pointer to data which might be required for this read method to do the read operation. It returns a pointer to the pileup iterator.

User can specify methods which need to be invoked during the load and unload of an alignment, like constructor and destructor of objects. Bam_plp_constructor and bam_plp_destructor methods does the setup of these methods in the pileup iterator. During invocation of these methods, the pointer to data passed in the initialization is passed as well. If user want to do any custom status handling or actions during load or unload, in these methods it can be done. Alignment specific data can be created and stored in constructor and the same will be accessible during pileup status return. The same will be accessible during destructor as well where any deallocation can be made.

User is expected to invoke bam_plp_auto api to get the pileup status. It returns the pileup status or NULL on end. During this all alignments are read one by one, using the method given in initialization for data read, until one for a new reference is found or all alignment covering a position is read. On such condition, the pileup status is returned and the same continuous on next bam_plp_auto call.
The pileup status returned is an array for all positions for which the processing is completed. Along with the result, the reference index, position in reference data and number of alignments which covers this position are passed. User can iterate the result array and get bases from each alignment which covers the given reference position. The alignment specific custom data which were created in constructor function will also be available in the result.
 
The bam_plp_auto api invokes the data read method to load an alignment and the constructor method is invoked during the load. Once the end of alignment is passed, it is removed from the processing and destructor method is invoked, that user could do deallocations and custom actions as in load during this time. The custom data passed during the initialization is passed to the constructor and destructor methods during invocation. 

One the forward and reverse strands are identified, the better of the quality is identified and used. Both reads are required for this and hence reads are cached until its mate is read. The maximum number of reads that can be cached is controlled by bam_plp_set_maxcnt. Reads covering a position are cached and as soon as mate is found, quality is adjusted and is removed from cache. Reads above the cache limit are discarded.

Once done, the pileup iterator to be discarded by sam_plp_destroy api.

    int main(int argc, char *argv[])
    { 
       bam_plp_t plpiter = NULL;
       bam_pileup1_t *plp = NULL;
       
       if (!(plpiter = bam_plp_init(readdata, &conf))) {
           fprintf(stderr, "Failed to initialize pileup data\n");
           goto end;
       }
    
    
       bam_plp_constructor(plpiter, plpconstructor);
       bam_plp_destructor(plpiter, plpdestructor);
    
    
       //set depth / max count if given or != 0
       if (cnt) {
           bam_plp_set_maxcnt(plpiter, cnt);
           printf("Set max cnt as %d\n", cnt);
       }
       while ((plp = bam_plp_auto(plpiter, &tid, &pos, &n))) {
           fprintf(stdout, "%d\t%d\tN\t", tid+1, pos+1);
           for (j = 0; j < n; ++j) {
               fprintf(stdout, "%s%c%s", plp[j].is_head? "^" : "", seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)], plp[j].is_tail? "$" : "" );
     
         if (plpiter) {
           bam_plp_destroy(plpiter);

Refer:pileup.c

This sample shows log on the constructor, destructor method invocation and shows pileup like output.
The read method may use a simple read or it could be an advanced read using indexes, iterators and region specifications based on the need.
The constructor method may create any custom data and store it in the pointer passed to it. The same need to be released by use on destructor method.

Pileup needs to be initialized with bam_mpileup_init method which takes pointer to a method, which will be called by pileup to read data from required files, and an array of pointer to data which might be required for this read method to do the read operation. It returns a pointer to the pileup iterator.

User can specify methods which need to be invoked during the load and unload of an alignment, like constructor and destructor of objects. bam_mplp_constructor and bam_mplp_destructor methods does the setup of these methods in the pileup iterator. During invocation of these methods, the pointer to data passed in the initialization is passed as well. If user want to do any custom status handling or actions during load or unload, in these methods it can be done. Alignment specific data can be created and stored in the custom data pointer and the same will be accessible during pileup status return. The same will be accessible during destructor as well where any deallocation can be made.

User is expected to invoke bam_mplp_auto api to get the pileup status. It returns the pileup status. During this all alignments are read one by one, using the method given in initialization for data read, until one for a new reference is found or all alignment covering a position is read. On such condition, the pileup status is returned and the same continuous on next bam_mplp_auto call.

The pileup status is returned through a parameter in the method itself, is an array for all inputs, each containing array for positions on which the processing is completed. Along with the result, the reference index, position in reference data and number of alignments which covers this position are passed. User can iterate the result array and get bases from each alignment which covers the given reference position. The alignment specific custom data which were created in constructor function will also be available in the result.
 
Once the forward and reverse strands are identified, the better of the quality is identified and used. Both reads are required for this and hence reads are cached until its mate is read. The maximum number of reads that can be cached is controlled by bam_mplp_set_maxcnt. Reads covering a position are cached and as soon as mate is found, quality is adjusted and is removed from cache. Reads above the cache limit are discarded.

Once done, the pileup iterator to be discarded by sam_mplp_destroy api.

    int main(int argc, char *argv[])
    { 
       bam_mplp_t plpiter = NULL;
       bam_pileup1_t **plp = NULL;
        
       if (!(plpiter = bam_mplp_init(ipcnt, readdata, confarray))) {
       
       bam_mplp_constructor(plpiter, plpconstructor);
       bam_mplp_destructor(plpiter, plpdestructor);
        //set depth / max count if given or != 0
       if (cnt) {
           //bam_plp_set_maxcnt(plpiter, cnt);
           bam_mplp_set_maxcnt(plpiter, cnt);
           printf("Set max cnt as %d\n", cnt);
       }
    
       while ((bam_mplp_auto(plpiter, &tid, &pos, depth, plp)) > 0) {
           fprintf(fp, "%d\t%d\tN\t", tid+1, pos+1);
           for (k = 0; k < ipcnt; ++k) {
               for (j = 0; j < depth[k]; ++j) {
                   fprintf(fp, "%s%c%s", plp[k][j].is_head? "^" : "", seq_nt16_str[bam_seqi(bam_get_seq(plp[k][j].b), plp[k][j].qpos)], plp[k][j].is_tail? "$" : "" );
    
       if (plpiter) {
           bam_mplp_destroy(plpiter);
       }

Refer:mpileup.c

This sample takes multiple sam files and shows the pileup of data side by side.

    ./mpileup -i /tmp/mp.bam -i /tmp/mp.sams


## More Information

### CRAM reference files
The cram reference data is required for the read of sequence data in CRAM format. The sequence data file may have it as embedded or as a reference to the actual file. When it is a reference, it is downloaded locally, in the cache directory for later usage. It will be stored in a directory structure based on the MD5 checksum in the cache directory.

Each chromosome in a reference file gets saved as a separate file with md5sum as its path and name. The initial 4 numerals make the directory name and rest as the file name (<cache path>/<1st 2 of md5sum>/<2nd 2 of md5sum>/<rest of md5sum>). 

The download would be attempted from standard location, EBI ENA (https://www.ebi.ac.uk/ena).

### Bam1_t
This structure holds the sequence data in BAM format.There are fixed and variable size fields, basic and extended information on sequence data. Variable size data and extended information are kept together in a buffer, named data in bam1_t. Fields in the member named core, bam1_core_t, and a few macros together support the storage and handling of the whole sequence data. 

core has a link to reference as a 0 based index in field tid. The mate / reverse strand’s link to reference is given by mtid. 
Field pos and mpos gives the position in reference to which the sequence and its mate / reverse strand match.
Field flag gives the properties of the given alignment. It shows the alignment’s orientation, mate status, read order etc.
Field qual gives the quality of the alignment read.
l_qname gives the length of the name of the alignment / read, l_extranul gives the extra space used internally in the data field.
l_qseq gives the length of the alignment / read in the data field.
n_cigar gives the number of CIGAR operations for the given alignment.
isize gives the insert size of the read / alignment.
The bases in sequence data are stored by compressing 2 bases together in a byte.
When the reverse flag is set, the base data is reversed and complemented from the actual read (i.e. if the read is ACTG, with a reverse flag, it will be stored as CAGT).

Macros bam_get_qname, bam_get_seq, bam_get_qual, bam_get_aux, bam_get_l_aux, bam_seqi etc access the data field and retrieve the required data. The aux macros support the retrieval of auxiliary data from the data field.
 
### Sam_hdr_t
This structure holds the header information. This holds the number of targets / SQ lines in the file, each one's length, name and reference count to this structure. It also has this information in an internal data structure for easier access of each field of this data.
When this data is shared or assigned to another variable of a different scope or purpose, the reference count needs to be incremented to ensure that it is valid till the end of the variable’s scope. sam_hdr_incr_ref and it needs to be destroyed as many times with sam_hdr_destroy api.

### Index
Indexes need the data to be sorted by position. 
They can be of different types with extension .bai, .csi or .tbi for compressed SAM/BAM files and .crai for CRAM files.
The index name can be passed along with the alignment file itself by appending a specific character sequence. The apis can detect this sequence and extract the index path. ##idx## is the sequence which separates the file path and index path.

### Data files
The data files can be a local file, a network file, a file accessible through the web or in cloud storage like google and amazon.
The data files can be represented with URIs like file://, file://localhost/.., ,ftp://.., gs+http[s].., s3+http[s]://

