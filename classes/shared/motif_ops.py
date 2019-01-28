import numpy as np
from seq_ops import reverse_complement

def build_pfm(sequences, length):
    """Function to build a PFM using entries from a fasta file

    Args:
        sequences (list): list of sequence strings
        length (int): size of pfm we are building

    Yields:
        pfm (numpy array): dimensions are 4xlength
    """
    
    # Initialize an empty numpy array
    pfm = np.zeros((4,length), dtype=int)

    # Add base-wise counts to the numpy array to build PFM matrix
    for seq in sequences:

        # Counter will track our position along the sequence (j)
        counter = 0

        # For each sequence we count which bases we see at each position (j)
        for char in list(seq):

            # An if/switch statement would work here, but this is an ASCII trick
            #  to convert A, C, G, T to integers 0, 1, 3, 2 respectively based
            #  on their binary representations
            pfm[((ord(char) >> 1 & 3), counter)] += 1

            # Increment counter along the sequence
            counter+=1

    # Because our trick above results in G and T in the wrong positions in the matrix
    #  we need to swap T and G rows for consistency
    pfm[3,:], pfm[2,:] = pfm[2,:], pfm[3,:].copy()

    # return PFM
    return pfm

def build_pwm(pfm):
    """Function to build a PWM from a PFM

    Args:
        pfm (numpy array): dimensions are 4xlength

    Yields:
        pwm (numpy array): dimensions are 4xlength
    """

    # Initialize an empty numpy array
    pwm = np.zeros((4,len(pfm[0])), dtype=float)
    
    # Calculate the sums of each column (Xij)
    sums = np.sum(pfm, axis = 0)
    
    # For simplicity, pseudocounts and background are set to .25
    p = .25
    bg = .25
    
    # For each position in the matrix, apply the PWM formula
    for i in range(4):
            for j in range(len(pwm[0])):
                    pwm[i,j] = np.log2((pfm[i,j] + p) / (sums[j]+p*4)) - np.log2(bg)

    # Return pwm
    return pwm

def score_kmer(seq, pwm):
    """Function to score a kmer with a pwm
        kmer length is expected to be the same as pwm length

    Args:
        seq(str): kmer to score
        pwm (numpy array): pwm for scoring

    Yields:
        score (float): PWM score for kmer
    """
    
    # Initialize score to 0
    score = 0
    
    if len(seq) != len(pwm[0]):
        raise ValueError('K-mer and PWM are different lengths!')
    
    # Translator for DNA to numeric indices
    dna_to_index = str.maketrans('ACGT', '0123')
    
    # Iterate across kmer and sum log likelihoods
    for j, val in enumerate(list(seq.translate(dna_to_index)), 0):
        score += pwm[int(val), j]

    # Return score
    return score

def score_nmer(seq, pwm):
    """Function to score a nmer with a pwm
        This will scan sequence and score all 
        subsequences of length k with a pwm
        and return the maximum score

    Args:
        seq(str): nmer to score
        pwm (numpy array): pwm for scoring

    Yields:
        score (float): PWM score for nmer
        position (int): 0-based index of the best match location
    """
    
    # Initialize score to -100
    max_score = -100
    max_index = 'None'
    max_strand = 'None'
    
    # Get PWM length
    pwm_len = len(pwm[0])
    
    # Make sure that our full sequence is at least the PWM length
    if len(seq) < pwm_len:
        raise ValueError('N-mer shorter than PWM!')
    
    # Iterate through the length of the sequence and test to see if the score is above max
    # Note: this is greedy and only the first max hit will be recorded
    for i in range(0, len(seq)-pwm_len+1):
        
        # Score the subsequence starting at i
        if score_kmer(seq[i:pwm_len+i], pwm) > max_score:
            max_score = score_kmer(seq[i:pwm_len+i], pwm)
            max_index = i
            max_strand = 0
            
        # Also test the reverse complement
        if score_kmer(reverse_complement(seq[i:pwm_len+i]), pwm) > max_score:
            max_score = score_kmer(reverse_complement(seq[i:pwm_len+i]), pwm)
            max_index = i
            max_strand = 1
    
    # Return maximum score and index
    return max_score, max_index, max_strand

def pfm_ic(pfm):
    ic = 0
    
    # Calculate the sums of each column (Xij)
    sums = np.sum(pfm, axis = 0)
    
    # For simplicity, pseudocounts are set to .25
    p = .25
    
    # For each position in the matrix, apply the IC formula
    for j in range(len(pfm[0])):
        ic += 2
        for i in range(4):
            fij = (pfm[i,j] + p) / (sums[j] + p*4)
            ic += fij * np.log2(fij)

    # Return ic
    return ic
