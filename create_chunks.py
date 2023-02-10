#!/usr/bin/env python3


import sys, os, re
import pyranges as pr
import pandas as pd
import logging
import argparse
import math
import subprocess

logging.basicConfig(level=logging.INFO,
                                        format='%(asctime)s : %(levelname)s : %(message)s',
                                        datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

def main():

    parser = argparse.ArgumentParser(description="chunk_genome",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--genome_fa", type=str, required=True,
                        help="genome in fasta format")

    parser.add_argument("--out_prefix", "-o", type=str, help="prefix for output files", required=True)
    parser.add_argument("--gene_annot_gtf", type=str, help="annotation gtf file", required=True)
    parser.add_argument("--chunks", type=str, help="chunks tsv file", required=True)
    
    args = parser.parse_args()

    genome_fa = args.genome_fa
    out_prefix = args.out_prefix
    chunks_file = args.chunks

    
    logger.info("-parsing chromosome lengths")
    genome_fai_file = genome_fa + ".fai"
    genome_contig_info = pd.read_csv(genome_fai_file, sep="\t", names=["Chromosome", "Length", "Offset", "Linebases", "Linewidth"])
    
    chromosomes = genome_contig_info['Chromosome'].values.tolist()
    chr_lengths = dict()
    for x in genome_contig_info.itertuples():
        chr_lengths[ x.Chromosome]  = x.Length


    logger.info("-reading in chunks")
    chunks = pd.read_csv(chunks_file, sep="\t")
    chromosomes_need_chunking = chunks['Chromosome'].unique().tolist()
    logger.info("-chromosomes needing chunking: {}".format(chromosomes_need_chunking))
    chromosomes_no_chunking_needed = [chrom for chrom in chromosomes if chrom not in chromosomes_need_chunking]
    

    logger.info("-parsing genome annotations in gtf")
    annotation_gtf = pd.read_csv("AmexT_v47-AmexG_v6.0-DD.gtf", sep="\t", names=[
                        "Chromosome", "Source", "ev_type", "Start", "End", "Score", "Strand", "dot", "info"])

    # beginning chunking
    logger.info("Starting to chunk data.\n")

    chunked_genome_filename = out_prefix + ".chunked.genome.fa"
    genome_ofh = open(chunked_genome_filename, "wt")

    
    

    chunked_annotation_df = None
    for chromosome_to_chunk in chromosomes_need_chunking:
        logger.info("-processing chromosome: {}".format(chromosome_to_chunk))
        chr_gtf = annotation_gtf[ annotation_gtf.Chromosome == chromosome_to_chunk ]

        breakpoints = chunks[chunks.Chromosome == chromosome_to_chunk ]['brkpt'].values
        breakpoints = sorted(breakpoints)
        logger.info(breakpoints)
        breakpoints.insert(0, 1)
        breakpoints.append(chr_lengths[chromosome_to_chunk]+1)

    
        for i in range(len(breakpoints)-1):
            lend_brkpt = breakpoints[i]
            rend_brkpt = breakpoints[i+1]
            rend_brkpt -= 1
            print(lend_brkpt, rend_brkpt)
            
            chrom_partition_name = chromosome_to_chunk + f"^c{i}^o{lend_brkpt}"
            
            chr_gtf_partition = chr_gtf[ (chr_gtf.Start >= lend_brkpt) & (chr_gtf.End <= rend_brkpt) ].copy()
            print(chr_gtf_partition.shape)
            
            # make adjustment to the gtf partition
            chr_gtf_partition['Chromosome'] = chrom_partition_name
            if lend_brkpt != 1:
                chr_gtf_partition['Start'] = chr_gtf_partition.Start - lend_brkpt + 1  # lend_brkpt becomes position 1
                chr_gtf_partition['End'] = chr_gtf_partition.End - lend_brkpt + 1

            
            chunked_annotation_df = pd.concat([chunked_annotation_df, chr_gtf_partition])
                
            # get genome chunk

            genome_fa_chunk = "\n".join(subprocess.check_output(f"samtools faidx {genome_fa} {chromosome_to_chunk}:{lend_brkpt}-{rend_brkpt}", shell=True).decode().split("\n")[1:])
            genome_fa_chunk = f">{chrom_partition_name}\n" + genome_fa_chunk
            genome_ofh.write(genome_fa_chunk)


    ## get rest of them, no chunking required.
    logger.info("-Done working on chunks.  Now adding in the rest that didn't need chunking....")

    for chrom in chromosomes_no_chunking_needed:
        # get chrom seq
        chrom_fa = subprocess.check_output(f"samtools faidx {genome_fa} {chrom}", shell=True).decode()
        genome_ofh.write(chrom_fa)


    # get gtf annotations
    rest_gtf = annotation_gtf[ annotation_gtf.Chromosome.isin(chromosomes_no_chunking_needed) ]
    chunked_annotation_df = pd.concat([chunked_annotation_df, rest_gtf])
    chunked_genome_annotation_filename = out_prefix + ".chunked.gtf"
    chunked_annotation_df.to_csv(chunked_genome_annotation_filename, sep="\t", header=False, index=False)

    logger.info("Done. See outputs: {}.chunked.*".format(out_prefix))


    sys.exit(0)

if __name__=='__main__':
    main()
                                    
