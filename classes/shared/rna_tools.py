import bamnostic as bs
from data_readers import get_gff

def calculate_RPKM(gene_reads, gene_length, sequence_depth):
    ''' Calculate RPKM value 
        RPKM = 𝑅𝐺/(𝐿𝐺/1000∗𝐷/1000000)
        Where, 𝑅𝐺 are the number of reads 𝑅 mapping to gene 𝐺, 𝐿𝐺 is then length of gene 𝐺, and 𝐷 is the total read count.
        
    Args:
        gene_reads (int): number of reads overlapping a gene
        gene_length (int): size of the gene
        sequence_depth (int): total number of reads in our sample
    
    Returns:
        rpkm (int): rpkm value
        
    Example:
        >>> calculate_RPKM(1953, 1301, 6167) #doctest: +ELLIPSIS
        243417.05...
    '''
    #return gene_reads / (gene_length / 1000 * sequence_depth/123)
    return gene_reads / (gene_length/1000 * sequence_depth/1000000)

def get_read_count(rna_alignments, seqid, start_position, end_position):
    ''' Calculates the read depth of a bam in a region
    
    Args:
        rna_alignments (str): bam file of aligned reads (requires a bai file to have been created)
        seqid (str): name of the contig for region (ie. chromosome or 'sample')
        start_position (int): start position of region
        end_position (int): end position of region
        
    Returns:
        read_count (int): total number of reads overlapping region
        
    Example:
        >>> get_read_count('sample_RNA_reads.bam', 'sample', 336, 1637)
        1953
    '''
    bam = bs.AlignmentFile(rna_alignments, 'rb')
    return sum(1 for read in enumerate(bam.fetch(seqid, start_position, end_position)))

def get_transcript_levels(rna_alignments, gff_file):
    ''' Gets RPKM values for all genes in a gff_file
    
    Args:
        rna_alignments (str): bam file of aligned reads
        gff_file (str): gff formatted gene annotation file
        
    Returns:
        transcript_expression (list of tuples): list of tuples containing contig name, 
            start position of gene, end position of gene, and RPKM for that gene
            
    Example:
        >>> get_transcript_levels('sample_RNA_reads.bam', 'data/sample_genomic.gff')  #doctest: +ELLIPSIS
        [('sample', 336, 1637, 243417.05...
    '''
    
    transcript_expression = []
    
    bam = bs.AlignmentFile(rna_alignments, 'rb')
    read_depth = sum(1 for read in bam)

    for gff_entry in get_gff(gff_file):
        if gff_entry.type == 'CDS':
            gene_reads = get_read_count(rna_alignments, gff_entry.seqid, gff_entry.start, gff_entry.end)
            rpkm = calculate_RPKM(gene_reads, gff_entry.end-gff_entry.start, read_depth)
            transcript_expression.append((gff_entry.seqid, gff_entry.start, gff_entry.end, rpkm))
        
    return transcript_expression