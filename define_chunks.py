#!/usr/bin/env python3


import sys, os, re
import pyranges as pr
import pandas as pd
import logging
import argparse
import math

logging.basicConfig(level=logging.INFO,
                                        format='%(asctime)s : %(levelname)s : %(message)s',
                                        datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

def main():

    parser = argparse.ArgumentParser(description="chunk_genome",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--genome_fa", type=str, required=True,
                        help="genome in fasta format")

    parser.add_argument("--out", "-o", dest="out", type=str, help="chunks output filename", required=True)
    parser.add_argument("--max_chunk_size", type=int, help="maximum chunk size", required=True)
    parser.add_argument("--gene_annot_gtf", type=str, help="annotation gtf file", required=True)
    parser.add_argument("--gene_spans", type=str, help="gene spans", required=True)
    parser.add_argument("--N_regions", type=str, help="N regions", required=True)
    
    
    args = parser.parse_args()

    genome_fa = args.genome_fa
    output_filename = args.out
    max_chunk_size = args.max_chunk_size
    gene_annot_gtf = args.gene_annot_gtf
    gene_spans_file = args.gene_spans
    N_regions_file = args.N_regions
    


    logger.info("-parsing chromosome lengths")
    genome_fai_file = genome_fa + ".fai"
    genome_contig_info = pd.read_csv(genome_fai_file, sep="\t", names=["Chromosome", "Length", "Offset", "Linebases", "Linewidth"])
    

    chromosomes_require_chunking = genome_contig_info[ genome_contig_info['Length'] > max_chunk_size]['Chromosome'].values.tolist()
    logger.info("-chromosomes requiring chunking include: {}".format(chromosomes_require_chunking))
    
    gene_spans_info = pd.read_csv(gene_spans_file, sep="\t", names=["gene_id", "Chromosome", "Start", "End", "Strand", "gene_sym", "gene_type"]) 
    gene_spans_info = gene_spans_info.drop(["Strand", "gene_sym", "gene_type"], axis=1)
    
    logger.info("-parsing N regions")
    N_regions_info = pd.read_csv(N_regions_file, sep="\t", names=["Chromosome", "Start", "End"])
    N_regions_info['Length'] = N_regions_info['End'] - N_regions_info['Start'] + 1
    N_regions_info = N_regions_info[N_regions_info['Length'] > 10 ]
    
    pr_N_regions = pr.PyRanges(N_regions_info)

    # merge gene spans
    logger.info("-merging gene span intervals")
    pr_gene_spans = pr.PyRanges(gene_spans_info)
    pr_gene_spans = pr_gene_spans.merge()
    df_gene_spans_merged = pr_gene_spans.as_df()

    # want intergenic regions:
    logger.info("-defining intergenic regions")
    df_intergenic = df_gene_spans_merged.copy()
    df_intergenic['intergenic_lend'] = df_intergenic['End'] + 1
    df_intergenic['intergenic_rend'] = df_intergenic.groupby('Chromosome').Start.shift(-1) -1

    df_intergenic = df_intergenic.drop(['Start','End'], axis=1)
    df_intergenic = df_intergenic.rename(columns={'intergenic_lend' : 'Start', 'intergenic_rend' : 'End'})
    df_intergenic = df_intergenic[ ~ df_intergenic['End'].isna() ]


    # find intergenic regions that encompass N regions.
    logger.info("-finding N regions within intergenic regions")
    pr_intergenic = pr.PyRanges(df_intergenic)
    pr_intergenic_N = pr_N_regions.join(pr_intergenic)
    pr_intergenic_N = pr_intergenic_N.drop(['Start_b', 'End_b'])

    df_chunks_defined = None
    
    genome_contigs_to_chunk_info = genome_contig_info[ genome_contig_info['Chromosome'].isin(chromosomes_require_chunking)]
    for contig_row in genome_contigs_to_chunk_info.itertuples():
        chromosome = contig_row.Chromosome

        logger.info("-processing {}".format(chromosome))

        chr_len = contig_row.Length
        # estimate number of chunks:
        num_chunks = math.ceil(chr_len/max_chunk_size)
        chunk_size = round(chr_len/num_chunks)
        chunk_brkpts = [ chunk_size * i for i in range(1, num_chunks) ]

        

        pr_brkpts = pr.from_dict({ 'Chromosome' : chromosome,
                                   'Start' : chunk_brkpts,
                                   'End' : chunk_brkpts })


        pr_brkpts_nearest = pr_brkpts.nearest(pr_intergenic_N)

        df_brkpts_nearest = pr_brkpts_nearest.as_df().copy()
        df_brkpts_nearest['brkpt'] = df_brkpts_nearest.apply(lambda x: int( (x['Start_b'] + x['End_b'])/2), axis=1)
        
        df_brkpts_nearest['chunksize_offset_frac'] = df_brkpts_nearest['Distance'] / max_chunk_size
        
        df_chunks_defined = pd.concat([df_chunks_defined, df_brkpts_nearest])


    df_chunks_defined[ ['Chromosome', 'brkpt', 'chunksize_offset_frac'] ].to_csv(output_filename, sep="\t", index=False)
    

if __name__=='__main__':
    main()
                                    
