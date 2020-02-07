# SNP-Pipeline

A pipeline written using bash scripting to implement SNP (Single Nucleotide Polymorphism) Calling. The pipeline maps genomic reads to a reference genome and calls SNPs from the mapping.

Tools used:
- bwa (for performing the alignment)
- samtools (processing the alignment file and variant calling)
- GATK v3.7.0 (Alignment improvement)

It follows the following steps:

- alignment of FASTQ reads to a reference genome to create an alignment file (using bwa)
- processing the alignment file: file format conversion, sorting, alignment improvement
- variant calling

Command-line arguments to execute the bash script:
- -a: Input reads file -pair1
- -b: Input reads file -pair2
- -r: Reference genome file
- -e: Select to perform read re-alignment
- -o: Output VCF file name
- -f: Mills file location
- -z: Output VCF file will be gunzipped (.gz)
- -v: Verbose mode
- -i: Select to index the output BAM file (uses samtools index)
- -h: prints usage information and exits

**Input files:** Input reads file - pair1, Input reads file - pair2, Reference genome file, Mills file

**Execution:** ./snp_pipeline.bash -a <Input reads file -pair1> -b <Input reads file -pair2> -r <Reference genome file> -o <Output VCF file name> -f <Mills file location>

Other command-line arguments mentioned above can also be added depending on whether read re-alignment, indexing the BAM file, etc need to be performed.

**Output:** VCF file  
