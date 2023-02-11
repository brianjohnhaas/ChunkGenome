# Chunking large Axo genome into smaller chunks

Code breaks large chromosomes (and corresponding annotations) into smaller chunks based on:

- a target chunk size (ie. 500M)
- positions of annotated intergenic regions and sequencing gaps (N-stretches).

A sequencing gap within an intergenic region that is nearest the target chunk size is chosen as a breakpoint.

A chunked genome.fa and annotation.gtf file is generated.

>Note, any contigs with lengths below the target chunk size are intouched but carried into the final outputs.


# Running

First, run define_chunks.py to generate the list of target breakpoints based on the intergenic regions and identified sequencing gaps.

Second, given the above defined chunks, run create_chunks.py to produce the chunked genome and annotation file.


>see the example/ and *.sh commands for more info.



