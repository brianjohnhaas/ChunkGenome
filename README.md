# Chunking large Axo genome into smaller chunks

Code breaks large chromosomes (and corresponding annotations) into smaller chunks based on:

- a target chunk size (ie. 500M)
- positions of annotated intergenic regions and sequencing gaps (N-stretches).

A sequencing gap within an intergenic region that is nearest the target chunk size is chosen as a breakpoint.

A chunked genome.fa and annotation.gtf file is generated.

>Note, any contigs with lengths below the target chunk size are intouched but carried into the final outputs.


# Running

## Defining the genome chunks
    
First, run define_chunks.py to generate the list of target breakpoints based on the intergenic regions and identified sequencing gaps.

```
    usage: define_chunks.py [-h] --genome_fa GENOME_FA --out OUT --max_chunk_size MAX_CHUNK_SIZE --gene_annot_gtf GENE_ANNOT_GTF --gene_spans GENE_SPANS
                        --N_regions N_REGIONS

define genome chunks

options:
  -h, --help            show this help message and exit
  --genome_fa GENOME_FA
                        genome in fasta format (default: None)
  --out OUT, -o OUT     chunks output filename (default: None)
  --max_chunk_size MAX_CHUNK_SIZE
                        maximum chunk size (default: None)
  --gene_annot_gtf GENE_ANNOT_GTF
                        annotation gtf file (default: None)
  --gene_spans GENE_SPANS
                        gene spans (default: None)
  --N_regions N_REGIONS
                        N regions (default: None)
```

An example:

```
    define_chunks.py --genome_fa my.genome.fa \
                     --max_chunk_size 2000000 \
                     --gene_annot_gtf annotation.gtf \
                     --gene_spans annotation.gene_spans \
                     --N_regions  genome_Nregions.tsv \
                     -o chunks.tsv
```    

>see the example/ and *.sh commands for more info.
    
    
See below for info on gene_spans and Nregions files, utilities, and formatting.

    
## Creating the chunked genome based on defined chunks
    
Second, given the above defined chunks, run create_chunks.py to produce the chunked genome and annotation file.

```
    usage: create_chunks.py [-h] --genome_fa GENOME_FA --out_prefix OUT_PREFIX --gene_annot_gtf GENE_ANNOT_GTF --chunks CHUNKS [--gtf_only]

chunk_genome

options:
  -h, --help            show this help message and exit
  --genome_fa GENOME_FA
                        genome in fasta format (default: None)
  --out_prefix OUT_PREFIX, -o OUT_PREFIX
                        prefix for output files (default: None)
  --gene_annot_gtf GENE_ANNOT_GTF
                        annotation gtf file (default: None)
  --chunks CHUNKS       chunks tsv file (default: None)
  --gtf_only            only write chunked annotation, not the genome (default: False)

```

An example:

```

       ../create_chunks.py --genome_fa my.genome.fa \
                           --chunks chunks.tsv \
                           -o chunked \
                           --gene_annot_gtf annotation.gtf
    
```
    
>see the example/ and *.sh commands for more info.


# Helper utilities

## Creating a gene spans file:

The following will generate a gene spans file based on the annotation.gtf file:

     ./util/gtf_to_gene_spans.pl annotation.gtf > annotation.gene_spans

and the format will look like so (tab-delimited):

```
ENSG00000151689.11  chr2   190343470  190371665  +  INPP1                         protein_coding
ENSG00000211978.2   chr14  106851123  106851417  -  IGHV5-78                      IG_V_pseudogene
ENSG00000228737.2   chr5   141618414  141626481  +  AC008781.7                    antisense
ENSG00000250205.1   chr4   155381042  155381239  -  YWHAEP4                       processed_pseudogene
ENSG00000152939.13  chr5   69415112   69444330   +  MARVELD2                      protein_coding
ENSG00000278135.1   chr4   25250545   25250639   +  U2                            snRNA
ENSG00000267667.1   chr17  61130704   61135964   -  RP11-136H19.1                 antisense
ENSG00000251680.4   chr5   129500361  129905917  -  CTC-575N7.1                   antisense
ENSG00000151135.8   chr12  106955719  106978778  +  TMEM263                       protein_coding
ENSG00000225752.1   chr6   143094034  143099645  -  RP1-45I4.2                    antisense
...
```    

## Defning the N-regions

Sequencing gaps in genome assemblies are usually represented by inserting a series of N characters into the genome at the positions of the sequencing gaps.  We can find these by searching out strings of N's within the genome fasta file.  The following script will do this and report the sequencing regions represented by N characters:

```
    util/getN_seq_ranges.pl genome.fasta > N_regions.tsv
```

and the output will be formatted like so (tab-delimited):

```
chr1  1          10000
chr1  207667     257666
chr1  297969     347968
chr1  535989     585988
chr1  2702782    2746290
chr1  12954385   13004384
chr1  16799164   16849163
chr1  29552234   29553835
chr1  121757469  121757469
chr1  121757509  121757509
...
```
