#!/usr/bin/env python
"""
Author: Lily Sakellaridi
Partial implementation of the SSAHA algorithm (doi:10.1101/gr.194201).
"""

# implement your functions here
def build_hash_table(dataset, k):
    """
    Return a hash table of kmers and their positions in the dataset.
    
    dataset: a list of DNA sequences. Each DNA sequence is a string.
    k: positive integer; the length of a kmer.
    
    The hash table only contains the kmers that occur in the dataset,
    not all possible kmers.
    The hash table is implemented as a dictionary. The keys
    are the kmers. The values are lists that contain all the positions
    where a kmer occurs.
    Each position is a list of two integers, i and j.
    i is the 1-based index of the sequence.
    j is the 1-based index of the kmer's position within the sequence.
    """
    
    hash_table = {}
    
    for i, seq in enumerate(dataset):
        for j in range(0, len(seq) - k + 1, k):
            kmer = seq[j: j + k]
            if kmer not in hash_table:
                hash_table[kmer] = [[i + 1, j + 1]] # 1-based indexing
            else:
                hash_table[kmer].append([i + 1, j + 1])
                
    return hash_table
    
def create_master_list(query, hash_table, k):
    """
    Return hit list of a query, sorted by index and shift.
    
    query: the query sequence; str
    hash_table: dict where keys: kmers and values: hits
    k: the length of each kmer; positive int
    
    The hit list is a list of lists. Each sublist contains
    three integers: index, shift, and offset.
    index: 1-based index of the sequence in the database
    (integer i in the output of build_hash_table)
    offset: 1-based position of the kmer within the sequence
    (integer j in the output of build_hash_table)
    The shift is calculated by substracting the current base, t,
    from the offset.
    The current base t is 0-indexed.
    
    The list is sorted first by index, then by shift.
    The sorted hit list is called the master list.
    """
    hit_list = []
    
    for t in range(len(query) - k + 1):
        kmer = query[t: t + k]
        if kmer in hash_table:
            # position list is a list of two-element sublists,
            # where each sublist is a position.
            position_list = hash_table[kmer]
            
            # Calculate the shift for each position
            for position in position_list:
                shift = position[1] - t
                hit = [position[0], shift, position[1]]
                hit_list.append(hit)
                
    hit_list.sort(key = lambda hit: (hit[0], hit[1]))
                
    return hit_list
    
def create_master_dict(master_list):
    """
    Return the master list as a nested dictionary.
    
    master_list: list of lists
    
    The master_dict is a dictionary with two levels of nesting.
    There is an outer dictionary with many inner dictionaries.
    
    The keys of the outer dictionary are the unique indices:
    so 1, 2 and 3 for the seqs dataset.
    
    The values of the outer dictionary are lists of inner dictionaries.
    
    For every inner dictionary, the keys are unique shifts,
    and the values are lists of offsets.
    
    For example, let the list of lists representation be:
    [[1, 3, 5], [1, 3, 7], [1, 4, 6], [1, 4, 9], [2, 4, 5], [2, 4, 9]]
    Then the nested dict representation is:
    {1: {3: [5, 7], 4: [6, 9]}, 2: {4: [5, 9]}}
    """
    
    outer_dict = {}
    
    for sublist in master_list:
        index = sublist[0]
        
        if index not in outer_dict:
            outer_dict[index] = {}
            shift = sublist[1]
            offset = sublist[2]
            
            if shift not in outer_dict[index]:
                outer_dict[index][shift] = [offset]
            else:
                outer_dict[index][shift].append(offset)
                
        else:
            inner_dict = outer_dict[index]
            shift = sublist[1]
            offset = sublist[2]
            
            if shift not in inner_dict:
                inner_dict[shift] = [offset]
            else:
                inner_dict[shift].append(offset)
            
                
    return outer_dict
    
                
            
                
            
    
def longest_match(master_dict, dataset, k):
    """
    Return the longest match between the query and the database.
    
    master_dict: nested dictionary representation of hits master list.
    dataset: list of lists; the set of sequences to be searched
    k: length of kmer; positive integer
    
    There are two outputs: longest and match_length
    
    longest: the longest match; a list of three integers.
    The first integer is the index of the sequence where the 
    longest match is found (1-indexing)
    The second integer is the first position of the match
    The third integer is the last position of the match 
    (non inclusive, so if a match is from 2 to 15, the 
    real last position is 14.)
    
    match_length: integer, the length of the match.
    """
    longest = []
    match_length = 0
    
    for i, seq in enumerate(dataset):
        inner_dict = master_dict[i + 1] # matches have been stored with 1-indexing
        for shift in inner_dict:
            offsets = inner_dict[shift]
            current_match = [i + 1, offsets[0], offsets[-1] + k]
            current_length = offsets[-1] + k - offsets[0]
            if current_length > match_length:
                longest = current_match
                match_length = current_length
                
    return longest, match_length
    
def reverse_complement(seq):
    """
    Return the reverse complement of a sequence.
    
    seq: dna sequence; str
    
    To implement SSAHA matches in the reverse direction,
    we take the reverse complement of the query and then
    follow the same procedure.
    """
    base_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    revseq = seq[::-1]
    revcomp = []
    
    for base in revseq:
        revcomp.append(base_dict[base])
        
    revcomp = ''.join(revcomp)
        
    return revcomp
    
def fasta_parser(fasta_file):
    """
    Return a dictionary where keys: headers and values: sequences.
    
    fasta_file: a fasta file to be parsed
    """
    with open(fasta_file) as f:
        records = f.read().split('>')
        
        headers = []
        seqs = []
        
        for record in records[1:]:
            lines = record.split('\n')
            headers.append('>' + lines[0])
            seqs.append("".join(lines[1:]))
            
    seq_dict = dict(zip(headers, seqs))
    
    return seq_dict
    
def print_max_alignment(query, dataset, longest):
    """
    Print the longest exact alignment.
    
    query: the query sequence; string
    dataset: database of sequences; list of strings
    
    The -1's are to compensate for 1-indexing
    used in the hash table, master list
    and master dict.
    """
    
    # the index of the sequence where the longest match is found
    db_idx = longest[0] - 1
    db_seq = str(dataset[db_idx])
    
    start = longest[1] - 1
    end = longest[2] - 1
    
    match = db_seq[start: end]
            
    print(' ' * start + match)
    print(db_seq)
    
    
    
def print_ssaha_report(query, db_headers, query_headers, i, match_length, direction):
    """
    Prints exact matches in the form of a SSAHA report.
    
    query: the query sequence; string
    dataset: the database to be searched; list of strings
    direction: 'F' -> forward, 'C' -> reverse complement; string
    """
    
    db_idx = longest[0]
    db_header = str(db_headers[db_idx])
    query_header = str(query_headers[i])
    
    start = str(longest[1])
    end = str(longest[2])
    
    if direction == 'C':
        end, start = start, end
    
    print(db_header)
    print("=" * 51)
    print("Matches for query " + str(i) + "( " + str(len(query)) + " bases):\n")
    print("=" * 51)
    print("\t".join(["Score", "Q_Name", "S_Name", "Q_Start", "Q_End", "S_Start", "S_End", "Direction", "#Bases", "identity"]))
    print("\t".join(['ALIGNMENT::50', query_header, db_header, start, end, start, end, direction, str(match_length), '100.0']))
    
    
    
