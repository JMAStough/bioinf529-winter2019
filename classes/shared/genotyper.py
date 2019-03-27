import bamnostic as bs
from collections import Counter
from scipy.stats import binom_test
from data_readers import get_fasta

def get_pileup(alignments, region_start = None, region_end = None):
    ''' Function build a read pileup list
        this is implemented as a list of Counter() from region_start to region_end
        with our small genome, it is reasonable to cover the entire genome
        but for larger genomes a smaller window is required.
    
    Args:
        alignments (str): bam file of alignments
        region_start (int): start position to build pileup
        region_end (int): end position to build pileup
    
    Returns:
        genome (list of Counter()): a list from region_start to region_end of
            Counters of allele frequencies
        
    Example:
        >>> get_pileup('sample_reads.bam') #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE 
        [Counter({'G': 5}), Counter({'G': 8}), Counter({'T': 12}), ...] 
    '''
    bam = bs.AlignmentFile(alignments, 'rb')

    if(region_start == None):
        region_start = 0
        
    if(region_end == None):
        region_end = bam.header['SQ'][0]['LN']
        
    genome = [Counter() for _ in range(region_start, region_end)]
    
    for read in bam:
        for i, base in enumerate(read.seq, read.pos):
            if i >= region_start and i <= region_end:
                genome[i-region_start].update(base)
  
    return genome

def binomial_test(major, minor):
    ''' Function to perform binomial test
        We will consider a Pvalue threshold of 0.10: 
        SNPs for which the P value of the binomial test < 0.10 failed the heterozygosity test.

    Args:
        major (int): count of most frequent allele
        minor (int): count of second most frequent allele
    
    Returns:
        is_above_threshold (bool): true if passes heterozygosity test, otherwise false
        
    Example:
        >>> binomial_test(8, 4)
        True
    '''
    p_val = binom_test(major, major+minor, 1/2)
    if p_val >= .1:
        return True

    return False

def find_variants(reference, alignments):
    ''' Function to find variants given sequencing alignments and a reference
        Identify variants that are heterozygous using heterozygosity test
        Identify variants that are homozygous alternative allele
        
        Note: Variants are reported as 1-based coordinates
        
    Args:
        reference (str): genome reference fasta file
        alignments (str): sequence alignments
    
    Returns:
        variant_list (list of tuples): list of variants as (position, allele1, allele2)
        
    Example:
        >>> find_variants(reference = 'data/sample_genome.fna', alignments = 'sample_reads.bam') #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE 
        [(240, 'A', 'G'), (354, 'G', 'A'), (803, 'C', 'A'), ...]
    '''
    
    variant_list = []
    
    # Get our genome to track variants
    genome = None
    for name, seq in get_fasta(reference):
        genome = seq

    # Build pileup
    genome_pileup = get_pileup(alignments, 0, len(genome))
    
    # Iterate through genome to identify variants
    for i, counter in enumerate(genome_pileup, 1): # start at 1 since variants are reported as 1-based
        ref_allele = genome[i-1] # Expected allele is our reference allele:
        top_alleles = counter.most_common(2) # Get two top alleles
        
        if len(counter) > 1: #test for heterozygous alleles
            if binomial_test(top_alleles[0][1], top_alleles[1][1]):
                # This is a SNP:
                variant_list.append((i, top_alleles[0][0], top_alleles[1][0]))
        elif len(counter) == 0: #no reads here
            next
        else: # test for homozygous alternate site
            if counter[ref_allele] == 0:
                variant_list.append((i, top_alleles[0][0], top_alleles[0][0]))

    return variant_list