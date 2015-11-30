# blobtools-light pipeline
Light version of the blobtools package, with additional pipeline scripts.

The core blobtools-light package allows the visualisation of (draft) genome assemblies using TAGC (Taxon-Annotated Gc-Coverage) plots (<a href="http://www.ncbi.nlm.nih.gov/pubmed/24348509">Kumar et al. 2012</a>).

In addition, the blobtools-light pipeline scripts perform the steps required to create input for the core package, as well as filtering of read data based on blobtools-light output. 

## Data Requirements

The core blobtools-light scripts require these input files:
- Assembly file* (FASTA)
- Coverage information (one or more):
  - COV (Tab-separated: "header\tcoverage")
  - BAM/SAM**
  - CAS format** (CLC mapper)
- BLAST result file

Most of these files can be created using blobtools-light pipeline scripts (see next section).

--------------------------------------------------------------------------------
###### Notes

\* If assembly was generated through Spades, Abyss or Velvet, coverage can be parsed from contig headers (e.g. -spades SPADESASSEMBLY); if coverage from assembly is not needed -exclude_assembly_cov can be specified.
\*\* The first time BAM/SAM/CAS file(s) gets parsed by blobtools-light, a COV file will be created for each.

--------------------------------------------------------------------------------

If running BLAST directly, the suggested command is:
```
blastn -task megablast -query assembly.fa -db blast_db/nt \
  -outfmt '6 qseqid staxids bitscore evalue std sscinames sskingdoms stitle' \
  -max_target_seqs 25 -culling_limit 2 -num_threads 16 -evalue 1e-25 \
  -out ASSEMBLY.vs.nt.25cul1.1e25.megablast.out
```
An essential part of this command is the output format parameter, in which the first four fields must be as shown (i.e. "-outfmt '6 qseqid staxids bitscore evalue'").

## Data Preparation

Input for blobtools-light can be created with the following pipeline commands.

Creating a blob contig assembly:

```
python assemble_blobs.py -1 reads_1.fastq -2 reads_2.fastq -a assembly.fa
```

Mapping sequencing reads to a blob contig assembly:

```
python map_blobs.py -1 reads_1.fastq -2 reads_2.fastq -a assembly.fa -bam mapping.bam
```

Screening of blob contigs against NCBI BLAST 'nt' database:

```
python blast_blobs.py -a assembly.fa -blastdb blast_db/ -blast blast.tsv
```

The previous three steps can be combined into a single call to the blobtools-light prep script:

```
python prep_blobs.py -1 reads_1.fastq -2 reads_2.fastq -a assembly.fa \
  -bam mapping.bam -blastdb blast_db/ -blast blast.tsv
```

## Making Blobs

The core blobtools-light pipeline is as follows.

Making blobplot file:
```
python makeblobs.py -a assembly.fa -bam mapping.bam -blast blast.tsv -o test
```

The resulting blob table file ('test.blobplot.txt') can be input into another script to plot the results.
```
python plotblobs.py test.blobplot.txt -o test
```

This will produce output as in the example below:

--------------------------------------------------------------------------------

![Example](example.blobplot.png?raw=true "Example Blobplot")

--------------------------------------------------------------------------------

## Filtering Blobs

Blob contigs (and their corresponding reads) can be filtered using the command:

```
python filter_blobs.py -target_taxa 6231 -foreign_taxa 1760 -blastdb blast_db/ \
   -1 reads_1.fastq -2 reads_2.fastq -a assembly.fa -blob_table test.blobplot.txt \
   -filter_plot filter_plot -email john.smith@example.com -bam mapping.bam \
   -f1 filtered_1.fastq -f2 filtered_2.fastq
```

This script generates filter information and uses that information to filter 
read files. The target and foreign taxa are specified as NCBI Taxonomy IDs, each 
of which uniquely identifies the given taxon. Target and foreign Taxonomy IDs 
should be chosen based on the information in the blob table and blob plot files 
that are output by the core blobtools-light scripts.

## Software Requirements
The blobtools-light pipeline package requires Python 2.7+ with the following packages:
- Biopython 1.65
- ete2 2.3.7
- matplotlib 1.4
- NumPy 1.9.1+
- PySAM 0.8.3
- scikit-learn 0.15.2

All the Python packages listed can be installed using pip.

One or more blobtools-light pipeline scripts require the following:
- BLAST 2.2.29+
- Bowtie2 2.2.5
- SGA 0.10.13
- samtools 1.2

