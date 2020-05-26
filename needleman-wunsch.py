#!/usr/bin/env python3

"""
Author: Lily Sakellaridi

Description: this is a script to implement the Needleman-Wunsch algorithm for
global sequence alignment. 

This implementation first initializes the alignment and traceback matrix.
Both matrices are filled. Then, the sequences are aligned based on the full matrices,
the alignment score is returned and the percent identity is calculated.
Finally the alignment and its statistics are printed.

The gap penalties are treated as negative integers because it is more convenient
to code them that way. In the report, they are treated as positive integers because that
is the proper definition. So a -1 in the code is a +1 in the report.

Referenced figures are from Chapter 5 of
'Understanding Bioinformatics' by Marketa Zvelebil & Jeremy O.Baum.
"""
#import statements here

from sys import exit


# functions between here and __main__
blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""

def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int,parts[1:])))
    return order, blosum_matrix

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()

def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]

# write your own functions below here

def zerosMat(seq1, seq2):
    """
    Return a matrix of zeros.
    
    seq1: str; the first sequence to be aligned
    seq2: str; the second sequence to be aligned.
    
    This is a helper function for 'initializeAlignmentMat' and
    'initializeTracebackMat'.
    The number of rows is equal to the length of seq2 plus one.
    The number of columns is equal to the length of seq1 plus one.
    The 'plus one' is because an extra zero is needed at the beginning of the matrix.
    """
    return [[0 for i in range(len(seq2) + 1)] for j in range(len(seq1) + 1)]
    
def initializeAlignmentMat(seq1, seq2, egp):
    """
    Return the initial alignment matrix.
    
    seq1: str; the first sequence to be aligned
    seq2: str; the second sequence to be aligned
    zerosMat: matrix of all zeros, output of zerosMat() function
    egp: negative int; end gap penalty
    
    The initial alignment matrix is all zeros, except the first row and first column.
    The values first row and first column are equal to the index of the current position
    multiplied by the end gap penalty.
    """
    
    init_align_mat = zerosMat(seq1, seq2)
    for i in range(len(seq1) + 1):
        init_align_mat[i][0] = egp * i
    for j in range(len(seq2) + 1):
        init_align_mat[0][j] = egp * j
    
    
    return init_align_mat
    
def initializeTracebackMat(seq1, seq2):
    """
    Return the initial traceback matrix.
    
    seq1: str; the first sequence to be aligned
    seq2: str; the second sequence to be aligned
    zeros_mat: matrix of zeros
    
    The traceback matrix is filled with pointers that 
    correspond to the possible origins. 'h' stands for
    'horizontal', 'v' for 'vertical' and 'd' for
    'diagonal'. Initially, we fill the first row with 'h', 
    and the first column with 'v'.
    """
    init_traceback_mat = zerosMat(seq1, seq2)
    
    for i in range(1, len(seq1) + 1):
        init_traceback_mat[i][0] = "v"
    for j in range(1, len(seq2) + 1):
        init_traceback_mat[0][j] = "h"
    return init_traceback_mat
    
def setOrigin(ver, hor, diag, max_score):
    """
    Set the origin of an alignment step.
    
    ver: the score of the vertical direction
    hor: the score of the horizontal direction
    diag: the score of the diagonal direction
    max_score: the maximum score
    
    This is a helper function for 'fillMatrices'.
    Pointers are: 'd' for diagonal, 'v' for vertical,
    and 'h' for horizontal. These are the same pointers
    used in the initialization of the traceback matrix.
    'Vertical' is the same as 'top'. 'Horizontal' is the same as 'left'.
    
    In case of ties, the diagonal origin is preferred,
    then the vertical, and finally the horizontal.
    This decision is completely arbitrary. I could not find
    a computational reason to favor one specific direction.
    Splitting of ties should be done based on biological reasons
    on a case by case basis.
    """
    if diag == max_score:
        return "d"
    elif ver == max_score:
        return "v"
    elif hor == max_score:
        return "h"
    else:
        print("None of your origins correspond to the maximum score.")
        exit(0)
    
def fillMatrices(seq1, seq2, lgp, egp):
    """
    Return the full alignment and traceback matrices.
    
    seq1: str; the first sequence to be aligned
    seq2: str; the second sequence to be aligned
    lgp: negative int; the linear gap penalty
    egp: negative int; the end gap penalty
    
    
    """
    alignment_mat = initializeAlignmentMat(seq1, seq2, egp)
    traceback_mat = initializeTracebackMat(seq1, seq2)
    
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            # calculate top, left and diagonal scores
            diag = alignment_mat[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1])
            
            # For top and left origin, decide whether to
            # apply end or linear gap penalty
            # If we're at the end of the sequence, use end penalty
            # Else, linear.
            if j == len(seq2):
                ver = alignment_mat[i - 1][j] + egp
            else:
                ver = alignment_mat[i - 1][j] + lgp
                
            if i == len(seq1):
                hor = alignment_mat[i][j - 1] + egp
            else:
                hor = alignment_mat[i][j - 1] + lgp
                
            max_score = max(ver, hor, diag)
            origin = setOrigin(ver, hor, diag, max_score)
            
            alignment_mat[i][j] = max_score
            traceback_mat[i][j] = origin
            
    return alignment_mat, traceback_mat
    
def align(seq1, seq2, lgp, egp):
    """
    Return the alignment of two sequences.
    
    seq1: str; the first sequence to be aligned
    seq2: str; the second sequence to be aligned
    lgp: negative int; linear gap penalty
    egp: negative int; end gap penalty
    
    The output alignment is a list of four elements. 
    The first element is aligned seq1, with gaps as '-'.
    The second element is aligned seq2, with gaps as '-'.
    The third element is a string of symbols: '|' means match,
    '.' means mismatch, and ' ' means gap.
    The fourth element is the score of the alignment.
    """
    
    alignment_mat, traceback_mat = fillMatrices(seq1, seq2, lgp, egp)
    
    i = len(seq1)
    j = len(seq2)
    
    aligned_seq1 = []
    aligned_seq2 = []
    
    # The score of the alignment is in the bottom right cell of
    # the filled alignment matrix.
    alignment_score = alignment_mat[i][j]
    
    
    while i > 0 or j > 0:
        origin = traceback_mat[i][j]
        
        if origin == 'd':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif origin == 'v':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif origin == 'h':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        else:
            print("Your traceback matrix contains elements other than 'v', 'h' and 'd'.")
            exit(0)
    
    # The residues/gaps are appended from the end of the alignment to the beginning,
    # so the sequences must be reversed to get the actual alignment.
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]
    
    
    symbols = []
    count_identity = 0
    
    for pos in range(len(aligned_seq1)):
        # This conditional assumes that gaps are never aligned with gaps.
        # We can make this assumption because, when the alignment matrix is filled,
        # gaps aligning with gaps is not an option.
        if aligned_seq1[pos] == aligned_seq2[pos]:
            symbols.append('|')
            count_identity += 1
        elif aligned_seq1[pos] == '-' or aligned_seq2[pos] == '-':
            symbols.append(' ')
        else:
            symbols.append('.')
            
    alignment_length = len(symbols)
    percent_identity = count_identity / alignment_length * 100
    symbols = "".join(symbols)
    alignment = ["".join(aligned_seq1), "".join(aligned_seq2), symbols, alignment_score, percent_identity]
    
    return alignment
    
    
def printAlignment(seq1, seq2, lgp, egp):
    """
    Print the alignment of two sequences.
    
    seq1: the first sequence to be aligned.
    seq2: the second sequence to be aligned.
    lgp: negative int; the linear gap penalty.
    egp: negative int; the end gap penalty.
    
    """
    alignment = align(seq1, seq2, lgp, egp)
    
    print("This is the optimal global alignment under the current settings: " + "\n")
    print(alignment[0] + "\n")
    print(alignment[1] + "\n")
    print(alignment[2] + "\n")
    print("The alignment score is: " + (str(alignment[3]) + "\n"))
    print("The percent identity of the alignment is: " + str(round(alignment[4], 2)) + "\n")
    
    
if __name__ == "__main__":
    
    

    seq1 = "THISLINE"
    seq2 = "ISALIGNED"
    # seq3: GPA1_ARATH
    seq3 = "MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRLDYPRLTKDIAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA"
    # seq4: GPA1_ORYSI
    seq4 = "MGSSCSRSHSLSEAETTKNAKSADIDRRILQETKAEQHIHKLLLLGAGESGKSTIFKQIKLLFQTGFDEAELRSYTSVIHANVYQTIKILYEGAKELSQVESDSSKYVISPDNQEIGEKLSDIDGRLDYPLLNKELVLDVKRLWQDPAIQETYLRGSILQLPDCAQYFMENLDRLAEAGYVPTKEDVLYARVRTNGVVQIQFSPVGENKRGGEVYRLYDVGGQRNERRKWIHLFEGVNAVIFCAAISEYDQMLFEDETKNRMMETKELFDWVLKQRCFEKTSFILFLNKFDIFEKKIQKVPLSVCEWFKDYQPIAPGKQEVEHAYEFVKKKFEELYFQSSKPDRVDRVFKIYRTTALDQKLVKKTFKLIDESMRRSREGT"
    
    # Reproduce alignments in fig.5.9, fig.5.11 and fig.5.12.
    
    print("Trying to reproduce the alignment in fig. 5.9.")
    printAlignment(seq1, seq2, -8, -8)
    print("\n")
    
    print("Trying to reproduce the alignment in fig. 5.11.")
    printAlignment(seq1, seq2, -4, -4)
    print("\n")
    
    print("Trying to reproduce the alignment in fig. 5.12.")
    printAlignment(seq1, seq2, -8, 0)
    print("\n")
    
    # Find alignments with linear gap penalties from 1-20, and no separate end gap penalty.
    
    print("Finding the alignments with linear gap penalties ranging from 1 to 20, and no separate end gap penalty.\n")
    print("\n")
    for i in range(1, 21):
        print("This is the alignment with linear gap penalty equal to " + str(i) + "\n")
        printAlignment(seq1, seq2, -i, -i)
        print("\n")
    
    # Align seq3 and seq4 with a linear gap penalty of 5 and an end gap penalty of 1.
    
    print("Aligning seq3 and seq4 with a linear gap penalty of 5 and an end gap penalty of 1.\n")
    printAlignment(seq3, seq4, -5, -1)
    print("\n")
    
    print("Aligning seq3 and seq4 with a linear gap penalty of 5 and an end gap penalty of 10.\n")
    printAlignment(seq3, seq4, -5, -10)
    print("\n")
    
    
    
    
