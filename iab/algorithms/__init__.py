#!/usr/bin/env python

#-----------------------------------------------------------------------------
# This work is licensed under the Creative Commons
# Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
# copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/.
#-----------------------------------------------------------------------------
from __future__ import division
from random import choice
from skbio.core.sequence import BiologicalSequence

blosum50 = {'A': {'A': 5, 'C': -1, 'D': -2, 'E': -1, 'F': -3, 'G': 0, 'H': -2, 'I': -1, 'K': -1, 'L': -2, 'M': -1, 'N': -1, 'P': -1, 'Q': -1, 'R': -2, 'S': 1, 'T': 0, 'V': 0, 'W': -3, 'Y': -2},
'C': {'A': -1, 'C': 13, 'D': -4, 'E': -3, 'F': -2, 'G': -3, 'H': -3, 'I': -2, 'K': -3, 'L': -2, 'M': -2, 'N': -2, 'P': -4, 'Q': -3, 'R': -4, 'S': -1, 'T': -1, 'V': -1, 'W': -5, 'Y': -3},
'D': {'A': -2, 'C': -4, 'D': 8, 'E': 2, 'F': -5, 'G': -1, 'H': -1, 'I': -4, 'K': -1, 'L': -4, 'M': -4, 'N': 2, 'P': -1, 'Q': 0, 'R': -2, 'S': 0, 'T': -1, 'V': -4, 'W': -5, 'Y': -3},
'E': {'A': -1, 'C': -3, 'D': 2, 'E': 6, 'F': -3, 'G': -3, 'H': 0, 'I': -4, 'K': 1, 'L': -3, 'M': -2, 'N': 0, 'P': -1, 'Q': 2, 'R': 0, 'S': -1, 'T': -1, 'V': -3, 'W': -3, 'Y': -2},
'F': {'A': -3, 'C': -2, 'D': -5, 'E': -3, 'F': 8, 'G': -4, 'H': -1, 'I': 0, 'K': -4, 'L': 1, 'M': 0, 'N': -4, 'P': -4, 'Q': -4, 'R': -3, 'S': -3, 'T': -2, 'V': -1, 'W': 1, 'Y': 4},
'G': {'A': 0, 'C': -3, 'D': -1, 'E': -3, 'F': -4, 'G': 8, 'H': -2, 'I': -4, 'K': -2, 'L': -4, 'M': -3, 'N': 0, 'P': -2, 'Q': -2, 'R': -3, 'S': 0, 'T': -2, 'V': -4, 'W': -3, 'Y': -3},
'H': {'A': -2, 'C': -3, 'D': -1, 'E': 0, 'F': -1, 'G': -2, 'H': 10, 'I': -4, 'K': 0, 'L': -3, 'M': -1, 'N': 1, 'P': -2, 'Q': 1, 'R': 0, 'S': -1, 'T': -2, 'V': -4, 'W': -3, 'Y': 2},
'I': {'A': -1, 'C': -2, 'D': -4, 'E': -4, 'F': 0, 'G': -4, 'H': -4, 'I': 5, 'K': -3, 'L': 2, 'M': 2, 'N': -3, 'P': -3, 'Q': -3, 'R': -4, 'S': -3, 'T': -1, 'V': 4, 'W': -3, 'Y': -1},
'K': {'A': -1, 'C': -3, 'D': -1, 'E': 1, 'F': -4, 'G': -2, 'H': 0, 'I': -3, 'K': 6, 'L': -3, 'M': -2, 'N': 0, 'P': -1, 'Q': 2, 'R': 3, 'S': 0, 'T': -1, 'V': -3, 'W': -3, 'Y': -2},
'L': {'A': -2, 'C': -2, 'D': -4, 'E': -3, 'F': 1, 'G': -4, 'H': -3, 'I': 2, 'K': -3, 'L': 5, 'M': 3, 'N': -4, 'P': -4, 'Q': -2, 'R': -3, 'S': -3, 'T': -1, 'V': 1, 'W': -2, 'Y': -1},
'M': {'A': -1, 'C': -2, 'D': -4, 'E': -2, 'F': 0, 'G': -3, 'H': -1, 'I': 2, 'K': -2, 'L': 3, 'M': 7, 'N': -2, 'P': -3, 'Q': 0, 'R': -2, 'S': -2, 'T': -1, 'V': 1, 'W': -1, 'Y': 0},
'N': {'A': -1, 'C': -2, 'D': 2, 'E': 0, 'F': -4, 'G': 0, 'H': 1, 'I': -3, 'K': 0, 'L': -4, 'M': -2, 'N': 7, 'P': -2, 'Q': 0, 'R': -1, 'S': 1, 'T': 0, 'V': -3, 'W': -4, 'Y': -2},
'P': {'A': -1, 'C': -4, 'D': -1, 'E': -1, 'F': -4, 'G': -2, 'H': -2, 'I': -3, 'K': -1, 'L': -4, 'M': -3, 'N': -2, 'P': 10, 'Q': -1, 'R': -3, 'S': -1, 'T': -1, 'V': -3, 'W': -4, 'Y': -3},
'Q': {'A': -1, 'C': -3, 'D': 0, 'E': 2, 'F': -4, 'G': -2, 'H': 1, 'I': -3, 'K': 2, 'L': -2, 'M': 0, 'N': 0, 'P': -1, 'Q': 7, 'R': 1, 'S': 0, 'T': -1, 'V': -3, 'W': -1, 'Y': -1},
'R': {'A': -2, 'C': -4, 'D': -2, 'E': 0, 'F': -3, 'G': -3, 'H': 0, 'I': -4, 'K': 3, 'L': -3, 'M': -2, 'N': -1, 'P': -3, 'Q': 1, 'R': 7, 'S': -1, 'T': -1, 'V': -3, 'W': -3, 'Y': -1},
'S': {'A': 1, 'C': -1, 'D': 0, 'E': -1, 'F': -3, 'G': 0, 'H': -1, 'I': -3, 'K': 0, 'L': -3, 'M': -2, 'N': 1, 'P': -1, 'Q': 0, 'R': -1, 'S': 5, 'T': 2, 'V': -2, 'W': -4, 'Y': -2},
'T': {'A': 0, 'C': -1, 'D': -1, 'E': -1, 'F': -2, 'G': -2, 'H': -2, 'I': -1, 'K': -1, 'L': -1, 'M': -1, 'N': 0, 'P': -1, 'Q': -1, 'R': -1, 'S': 2, 'T': 5, 'V': 0, 'W': -3, 'Y': -2},
'V': {'A': 0, 'C': -1, 'D': -4, 'E': -3, 'F': -1, 'G': -4, 'H': -4, 'I': 4, 'K': -3, 'L': 1, 'M': 1, 'N': -3, 'P': -3, 'Q': -3, 'R': -3, 'S': -2, 'T': 0, 'V': 5, 'W': -3, 'Y': -1},
'W': {'A': -3, 'C': -5, 'D': -5, 'E': -3, 'F': 1, 'G': -3, 'H': -3, 'I': -3, 'K': -3, 'L': -2, 'M': -1, 'N': -4, 'P': -4, 'Q': -1, 'R': -3, 'S': -4, 'T': -3, 'V': -3, 'W': 15, 'Y': 2},
'Y': {'A': -2, 'C': -3, 'D': -3, 'E': -2, 'F': 4, 'G': -3, 'H': 2, 'I': -1, 'K': -2, 'L': -1, 'M': 0, 'N': -2, 'P': -3, 'Q': -1, 'R': -1, 'S': -2, 'T': -2, 'V': -1, 'W': 2, 'Y': 8}}

nt_substitution_matrix = {'A': {'A':  1, 'C': -2, 'G': -2, 'T': -2, 'N': 0},
                          'C': {'A': -2, 'C':  1, 'G': -2, 'T': -2, 'N': 0},
                          'G': {'A': -2, 'C': -2, 'G':  1, 'T': -2, 'N': 0},
                          'T': {'A': -2, 'C': -2, 'G': -2, 'T':  1, 'N': 0},
                          'N': {'A':  0, 'C':  0, 'G':  0, 'T':  0, 'N': 0 }}


def hamming_distance(s1, s2):
    s1 = BiologicalSequence(s1)
    s2 = BiologicalSequence(s2)
    return s1.distance(s2)


def format_matrix(row_headers, col_headers, data, hide_zeros=False, cell_width=3):
    result = []
    cell_format = "%" + str(cell_width) + "s"
    line_format = cell_format * (len(row_headers) + 1)
    
    # print a header row 
    result.append(line_format % tuple([' '] + list(row_headers)))
    
    # print the data rows
    for b2, row in zip(col_headers,data):
        if hide_zeros:
            display_row = []
            for v in row:
                if v == 0:
                    display_row.append('')
                else:
                    display_row.append(v)
        else:
            display_row = row
        result.append(line_format % tuple([b2] + display_row))
    
    return '\n'.join(result)

def format_dynamic_programming_matrix(seq1, seq2, matrix, cell_width=4):
    """ define a function for formatting dynamic programming matrices
    """
    lines = []

    cell_format = "%" + str(cell_width) + "s"
    line_format = cell_format * (len(seq1) + 2)
    # print seq1 (start the line with two empty strings)
    lines.append(line_format % tuple([' ',' '] + map(str,list(seq1))))

    # iterate over the rows and print each (starting with the
    # corresponding base in sequence2)
    for row, base in zip(matrix,' ' + seq2):
        lines.append(line_format % tuple([base] + map(str,row)))
    
    return '\n'.join(lines)


def generate_score_matrix(seq1,seq2,substitution_matrix):
    # Initialize a matrix to use for storing the scores
    score_matrix = []
    # Iterate over the amino acids in sequence two (which will correspond
    # to the vertical sequence in the matrix)
    for aa2 in seq2:
        # Initialize the current row of the matrix
        current_row = []
        # Iterate over the amino acids in sequence one (which will
        # correspond to the horizontal sequence in the matrix)
        for aa1 in seq1:
            # score as 1 if the bases are equal and 0 if they're not
            current_row.append(substitution_matrix[aa1][aa2])
        # append the current row to the matrix
        score_matrix.append(current_row)
    return score_matrix

##
# Needleman-Wunsch alignment
##

def generate_nw_and_traceback_matrices(seq1,seq2,gap_penalty,substitution_matrix):

    # Initialize a matrix to use for scoring the alignment and for tracing
    # back the best alignment
    nw_matrix = [[-gap_penalty * i for i in range(0,len(seq1)+1)]]
    traceback_matrix = [[None] + ['-' for i in range(0,len(seq1))]]
    # Iterate over the amino acids in sequence two (which will correspond
    # to the vertical sequence in the matrix)
    # Note that i corresponds to column numbers, as in the 'Biological Sequence
    # Analysis' example
    for i,aa2 in zip(range(1,len(seq2)+1),seq2):
        # Initialize the current row of the matrix
        current_row = [i * -gap_penalty]
        current_traceback_matrix_row = ['|']
        # Iterate over the amino acids in sequence one (which will
        # correspond to the horizontal sequence in the matrix)
        # Note that j corresponds to row numbers, as in the 'Biological Sequence
        # Analysis' example from class
        for j,aa1 in zip(range(1,len(seq1)+1),seq1):
            substitution_score = substitution_matrix[aa1][aa2]
            diag_score = (nw_matrix[i-1][j-1] + substitution_score,'\\')
            up_score = (nw_matrix[i-1][j] - gap_penalty,'|')
            left_score = (current_row[-1] - gap_penalty,'-')
            best_score = max(diag_score,up_score,left_score)
            current_row.append(best_score[0])
            current_traceback_matrix_row.append(best_score[1])
        # append the current row to the matrix
        nw_matrix.append(current_row)
        traceback_matrix.append(current_traceback_matrix_row)
    return nw_matrix, traceback_matrix


def nw_traceback(traceback_matrix,nw_matrix,seq1,seq2,gap_character='-'):

    aligned_seq1 = []
    aligned_seq2 = []

    current_row = len(traceback_matrix) - 1
    current_col = len(traceback_matrix[0]) - 1

    best_score = nw_matrix[current_row][current_col]

    while True:
        current_value = traceback_matrix[current_row][current_col]

        if current_value == '\\':
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
            current_col -= 1
        elif current_value == '|':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
        elif current_value == '-':
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append('-')
            current_col -= 1
        elif current_value == None:
            break
        else:
            raise ValueError, "Invalid value in traceback matrix: %s" % current_value

    return ''.join(aligned_seq1[::-1]), ''.join(aligned_seq2[::-1]), best_score


def nw_align(seq1, seq2, gap_penalty, substitution_matrix):
    """ Perform Needleman-Wunsch alignment of seq1 and seq2
    """
    nw_matrix, traceback_matrix = generate_nw_and_traceback_matrices(
                                    seq1, seq2, gap_penalty, substitution_matrix)
    aligned_seq1, aligned_seq2, score = nw_traceback(traceback_matrix,nw_matrix,seq1,seq2)
    return aligned_seq1, aligned_seq2, score

def nw_align_nt(seq1, seq2, gap_penalty=8, substitution_matrix=nt_substitution_matrix):
    """Globally align two nucleotide seqs (Needleman-Wunsch w single gap scoring)
    
       Parameters
       ----------
       sequence1 : string
           The first unaligned sequence
       sequence2 : string
           The second unaligned sequence
       gap_penalty : int, float, optional
           penalty for inserting a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar), optional
           lookup for substitution scores (these values are added to the 
           previous best alignment score)
        
       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment
         
    """
    return nw_align(seq1, seq2, gap_penalty, substitution_matrix)


##
# Smith-Waterman alignment
##


def generate_sw_and_traceback_matrices(seq1,seq2,gap_penalty,substitution_matrix):
    # Initialize a matrix to use for scoring the alignment and for tracing
    # back the best alignment
    sw_matrix = [[0 for i in range(0,len(seq1)+1)]]
    traceback_matrix = [[None] + [None for i in range(0,len(seq1))]]
    # Iterate over the amino acids in sequence two (which will correspond 
    # to the vertical sequence in the matrix)
    # Note that i corresponds to column numbers, as in the 'Biological Sequence 
    # Analysis' example from class
    for i,aa2 in zip(range(1,len(seq2)+1),seq2):
        # Initialize the current row of the matrix
        current_row = [0]
        current_traceback_matrix_row = [None]
        # Iterate over the amino acids in sequence one (which will 
        # correspond to the horizontal sequence in the matrix)
        # Note that j corresponds to row numbers, as in the 'Biological Sequence 
        # Analysis' example from class
        new_alignment_score = (0,None)
        for j,aa1 in zip(range(1,len(seq1)+1),seq1):
            substitution_score = substitution_matrix[aa1][aa2]
            diag_score = (sw_matrix[i-1][j-1] + substitution_score,'\\')
            up_score = (sw_matrix[i-1][j] - gap_penalty,'|')
            left_score = (current_row[-1] - gap_penalty,'-')
            best_score = max(diag_score,up_score,left_score,new_alignment_score)
            current_row.append(best_score[0])
            current_traceback_matrix_row.append(best_score[1])
        # append the current row to the matrix
        sw_matrix.append(current_row)
        traceback_matrix.append(current_traceback_matrix_row)
    return sw_matrix, traceback_matrix


def sw_traceback(traceback_matrix,sw_matrix,seq1,seq2,gap_character='-'):
    
    aligned_seq1 = []
    aligned_seq2 = []
    
    current_row = None 
    current_col = None
    best_score = 0
    for i in range(len(sw_matrix[0])):
        for j in range(len(sw_matrix)):
            current_score = sw_matrix[j][i]
            if current_score > best_score:
                best_score = current_score
                current_row = j
                current_col = i
    
    while True:
        current_value = traceback_matrix[current_row][current_col]
        
        if current_value == '\\':
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
            current_col -= 1
        elif current_value == '|':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
        elif current_value == '-':
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append('-')
            current_col -= 1
        elif current_value == None:
            break
        else:
            raise ValueError, "Invalid value in traceback matrix: %s" % current_value
        
    return ''.join(aligned_seq1[::-1]), ''.join(aligned_seq2[::-1]), best_score, current_col, current_row

def sw_multiple_traceback(traceback_matrix,sw_matrix,seq1,seq2,gap_character='-'):
    
    aligned_seq1 = []
    aligned_seq2 = []
    
    current_row = None 
    current_col = None
    best_scores = []
    for i in range(len(sw_matrix[0])):
        for j in range(len(sw_matrix)):
            try:
                new_alignment = traceback_matrix[j+1][i+1] == None
            except IndexError:
                new_alignment = True
            if new_alignment:
                current_score = sw_matrix[j][i]
                current_row = j
                current_col = i
                best_scores.append((current_score, current_row, current_col))
    best_scores.sort()
    best_scores.reverse()
    
    results = []
    for current_score, current_row, current_col in best_scores[:10]:
        while True:
            current_value = traceback_matrix[current_row][current_col]
        
            if current_value == '\\':
                aligned_seq1.append(seq1[current_col-1])
                aligned_seq2.append(seq2[current_row-1])
                current_row -= 1
                current_col -= 1
            elif current_value == '|':
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[current_row-1])
                current_row -= 1
            elif current_value == '-':
                aligned_seq1.append(seq1[current_col-1])
                aligned_seq2.append('-')
                current_col -= 1
            elif current_value == None:
                break
            else:
                raise ValueError, "Invalid value in traceback matrix: %s" % current_value
        results.append((''.join(aligned_seq1[::-1]), ''.join(aligned_seq2[::-1]), current_score, current_col, current_row))
        
    return results


def sw_align(sequence1, sequence2, gap_penalty, substitution_matrix):
    sw_matrix, traceback_matrix = generate_sw_and_traceback_matrices(sequence1,
                                                                 sequence2,
                                                                 gap_penalty,
                                                                 substitution_matrix)

    return sw_traceback(traceback_matrix,sw_matrix,sequence1,sequence2)


def generate_sw_and_traceback_matrices_affine_gap(seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix):
    # Initialize a matrix to use for scoring the alignment and for tracing
    # back the best alignment
    sw_matrix = [[0 for i in range(0,len(seq1)+1)]]
    traceback_matrix = [[None] + [None for i in range(0,len(seq1))]]
    # Iterate over the amino acids in sequence two (which will correspond 
    # to the vertical sequence in the matrix)
    # Note that i corresponds to column numbers, as in the 'Biological Sequence 
    # Analysis' example from class
    for i,aa2 in zip(range(1,len(seq2)+1),seq2):
        # Initialize the current row of the matrix
        current_row = [0]
        current_traceback_matrix_row = [None]
        # Iterate over the amino acids in sequence one (which will 
        # correspond to the horizontal sequence in the matrix)
        # Note that j corresponds to row numbers, as in the 'Biological Sequence 
        # Analysis' example from class
        new_alignment_score = (0,None)
        for j,aa1 in zip(range(1,len(seq1)+1),seq1):
            substitution_score = substitution_matrix[aa1][aa2]
            diag_score = (sw_matrix[i-1][j-1] + substitution_score,'\\')
            if traceback_matrix[i-1][j] == '|':
                # gap extend, because the cell above was also a gap
                up_score = (sw_matrix[i-1][j] - gap_extend_penalty,'|')
            else:
                # gap open, because the cell above was not a gap
                up_score = (sw_matrix[i-1][j] - gap_open_penalty,'|')
            if current_traceback_matrix_row[-1] == '-':
                # gap extend, because the cell to the left was also a gap
                left_score = (current_row[-1] - gap_extend_penalty,'-')
            else:
                # gap open, because the cell to the left was not a gap
                left_score = (current_row[-1] - gap_open_penalty,'-')
            best_score = max(diag_score,up_score,left_score,new_alignment_score)
            current_row.append(best_score[0])
            current_traceback_matrix_row.append(best_score[1])
        # append the current row to the matrix
        sw_matrix.append(current_row)
        traceback_matrix.append(current_traceback_matrix_row)
    return sw_matrix, traceback_matrix

def sw_align_affine_gap(sequence1, sequence2, gap_open_penalty, gap_extend_penalty, substitution_matrix):
    sw_matrix, traceback_matrix = generate_sw_and_traceback_matrices_affine_gap(sequence1,
                                                                 sequence2,
                                                                 gap_open_penalty,
                                                                 gap_extend_penalty,
                                                                 substitution_matrix)

    return sw_traceback(traceback_matrix,sw_matrix,sequence1,sequence2)

def sw_align_affine_gap_nt(sequence1, sequence2, gap_open_penalty=5,
                        gap_extend_penalty=2, substitution_matrix=None):
    """Locally align two nucleotide seqs (Smith-Waterman w affine gap scoring)
    
       Parameters
       ----------
       sequence1 : string
           The first unaligned sequence
       sequence2 : string
           The second unaligned sequence
       gap_open_penalty : int, float, optional
           penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive)
       gap_extend_penalty : int, float, optional
           penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar), optional
           lookup for substitution scores (these values are added to the 
           previous best alignment score); default is nt_substitution_matrix
        
       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment
       int
          The start position of the alignment in sequence 1
       int
          The start position of the alignment in sequence 2
    
       Examples
       --------
       >>> from iab.algorithms import sw_align_affine_gap_nt
       >>> s1 = "GCGTGCCTAAGGTATGCAAG"
       >>> s2 = "ACGTGCCTAGGTACGCAAG"
       >>> a1, a2, score, a1_start, a2_start = sw_align_affine_gap_nt(s1, s2)
       >>> print a1
       CGTGCCTAAGGTATGCAAG
       >>> print a2
       CGTGCCT-AGGTACGCAAG
       
    """
    if substitution_matrix is None:
        substitution_matrix = nt_substitution_matrix
        
    return sw_align_affine_gap(sequence1, sequence2, gap_open_penalty,
                               gap_extend_penalty, substitution_matrix)

def sw_multiple_align_affine_gap_nt(sequence1, sequence2, gap_open_penalty=5,
                        gap_extend_penalty=2,
                        substitution_matrix=nt_substitution_matrix):
    sw_matrix, traceback_matrix = generate_sw_and_traceback_matrices_affine_gap(sequence1,
                                                                 sequence2,
                                                                 gap_open_penalty,
                                                                 gap_extend_penalty,
                                                                 substitution_matrix)

    return sw_multiple_traceback(traceback_matrix,sw_matrix,sequence1,sequence2)

def sw_align_affine_gap_pr(sequence1, sequence2, gap_open_penalty=11,
                        gap_extend_penalty=1,
                        substitution_matrix=blosum50):
    """Locally align two protein seqs (Smith-Waterman w affine gap scoring)
    
       Parameters
       ----------
       sequence1 : string
           The first unaligned sequence
       sequence2 : string
           The second unaligned sequence
       gap_open_penalty : int, float, optional
           penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive)
       gap_extend_penalty : int, float, optional
           penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar), optional
           lookup for substitution scores (these values are added to the 
           previous best alignment score)
        
       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment
       int
          The start position of the alignment in sequence 1
       int
          The start position of the alignment in sequence 2
         
    """
    return sw_align_affine_gap(sequence1, sequence2, gap_open_penalty,
                               gap_extend_penalty, substitution_matrix)

def sw_align_nt(seq1, seq2, gap_penalty=8, substitution_matrix=nt_substitution_matrix):
    """Locally align two nucleotide seqs (Smith-Waterman w single gap scoring)
    
       Parameters
       ----------
       sequence1 : string
           The first unaligned sequence
       sequence2 : 
           The second unaligned sequence
       gap_penalty : int, float, optional
           penalty for inserting a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar), optional
           lookup for substitution scores (these values are added to the 
           previous best alignment score)
        
       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment
       int
          The start position of the alignment in sequence 1
       int
          The start position of the alignment in sequence 2
         
    """
    return sw_align(seq1, seq2, gap_penalty, substitution_matrix)

###
# Convenience wrapper
###


def align(sequence1, sequence2, gap_penalty, substitution_matrix, local):
    if local:
        return sw_align(sequence1, sequence2, gap_penalty, substitution_matrix)
    else:
        return nw_align(sequence1, sequence2, gap_penalty, substitution_matrix)


from scipy.cluster.hierarchy import average, dendrogram, to_tree

from bipy.core.distance import SymmetricDistanceMatrix

from iab.algorithms import nt_substitution_matrix, hamming_distance

def get_k_words(s, k, overlapping=True):
    result = []
    len_s = len(s)
    if overlapping:
        step = 1
    else:
        step = k
    for i in range(0,len_s,step):
        if i+k > len_s:
            # if there are no more k-mers left
            break
        else:
            result.append(s[i:i+k])
    return result

def fraction_unique_words(words1, words2):
    words1_set = set(words1)
    words2_set = set(words2)
    all_words = words1_set | words2_set
    shared_words = words1_set & words2_set
    number_unique = len(all_words) - len(shared_words)
    result = number_unique / len(all_words)
    return result

def kmer_distance(seq1, seq2, k):
    seq1_k_words = get_k_words(seq1, k)
    seq2_k_words = get_k_words(seq2, k)
    return fraction_unique_words(seq1_k_words, seq2_k_words)

def three_mer_distance(seq1, seq2):
    return kmer_distance(seq1, seq2, k=3)

def guide_tree_from_query_sequences(query_sequences, 
                                    distance_fn=three_mer_distance,
                                    display_tree = False):
    guide_dm = []
    seq_ids = []
    for seq_id1, seq1 in query_sequences:
        seq_ids.append(seq_id1)
        row = []
        for seq_id2, seq2 in query_sequences:
            row.append(kmer_distance(seq1, seq2, k=3))
        guide_dm.append(row)
    
    guide_dm = SymmetricDistanceMatrix(guide_dm, seq_ids)
    guide_lm = average(guide_dm.condensed_form())
    guide_tree = to_tree(guide_lm)
    if display_tree:
        guide_d = dendrogram(guide_lm, labels=guide_dm.ids, orientation='right', 
               link_color_func=lambda x: 'black')
    return guide_tree

def msa_generate_nw_and_traceback_matrices(aln1,aln2,gap_open_penalty,
                                           gap_extend_penalty,substitution_matrix):
    gap_open_penalty = float(gap_open_penalty)
    gap_extend_penalty = float(gap_extend_penalty)
    # Initialize a matrix to use for scoring the alignment and for tracing
    # back the best alignment
    nw_matrix = [[0, -1 * gap_open_penalty]]
    for i in range(2,len(aln1[0])+1):
        nw_matrix[0].append(nw_matrix[0][-1] - gap_extend_penalty)
    traceback_matrix = [[None] + ['-' for i in range(0,len(aln1[0]))]]
    # Iterate over the amino acids in sequence two (which will correspond
    # to the vertical sequence in the matrix)
    # Note that i corresponds to column numbers, as in the 'Biological Sequence
    # Analysis' example
    for i in range(1,len(aln2[0])+1):
        # Initialize the current row of the matrix
        if i == 1:
            current_row = [nw_matrix[i-1][0] - gap_open_penalty]
        else:
            current_row = [nw_matrix[i-1][0] - gap_extend_penalty]
        current_traceback_matrix_row = ['|']
        # Iterate over the amino acids in sequence one (which will
        # correspond to the horizontal sequence in the matrix)
        # Note that j corresponds to row numbers, as in the 'Biological Sequence
        # Analysis' example from class
        for j in range(1,len(aln1[0])+1):
            # computing the subsitution score is different when aligning alignments
            substitution_score = 0
            aln2_aas = [seq[i-1] for seq in aln2]
            aln1_aas = [seq[j-1] for seq in aln1]
            for aa2 in aln2_aas:
                for aa1 in aln1_aas:
                    if aa1 == "-" or aa2 == "-":
                        substitution_score += 0
                    else:
                        substitution_score += substitution_matrix[aa1][aa2]
            substitution_score /= (len(aln1) * len(aln2))

            # everything else is the same as for pairwise nw
            diag_score = (nw_matrix[i-1][j-1] + substitution_score,'\\')

            # affine gaps
            # up_score = (nw_matrix[i-1][j] - gap_penalty,'|')
            if traceback_matrix[i-1][j] == '|':
                # gap extend, because the cell above was also a gap
                up_score = (nw_matrix[i-1][j] - gap_extend_penalty,'|')
            else:
                # gap open, because the cell above was not a gap
                up_score = (nw_matrix[i-1][j] - gap_open_penalty,'|')
                
            #left_score = (current_row[-1] - gap_penalty,'-')
            if current_traceback_matrix_row[-1] == '-':
                # gap extend, because the cell to the left was also a gap
                left_score = (current_row[-1] - gap_extend_penalty,'-')
            else:
                # gap open, because the cell to the left was not a gap
                left_score = (current_row[-1] - gap_open_penalty,'-')
            
            best_score = max(diag_score,up_score,left_score)
            current_row.append(best_score[0])
            current_traceback_matrix_row.append(best_score[1])
        # append the current row to the matrix
        nw_matrix.append(current_row)
        traceback_matrix.append(current_traceback_matrix_row)
    return nw_matrix, traceback_matrix

def msa_nw_traceback(traceback_matrix,nw_matrix,aln1,aln2,gap_character='-'):

    # initialize the result alignments
    len_aln1 = len(aln1)
    aligned_seqs1 = []
    for e in range(len_aln1):
        aligned_seqs1.append([])
        
    len_aln2 = len(aln2)
    aligned_seqs2 = []
    for e in range(len_aln2):
        aligned_seqs2.append([])

    current_row = len(traceback_matrix) - 1
    current_col = len(traceback_matrix[0]) - 1

    best_score = nw_matrix[current_row][current_col]

    while True:
        current_value = traceback_matrix[current_row][current_col]

        if current_value == '\\':
            for i in range(len_aln1):
                aligned_seqs1[i].append(aln1[i][current_col-1])
            for i in range(len_aln2):
                aligned_seqs2[i].append(aln2[i][current_row-1])
            current_row -= 1
            current_col -= 1
        elif current_value == '|':
            for i in range(len_aln1):
                aligned_seqs1[i].append('-')
            for i in range(len_aln2):
                aligned_seqs2[i].append(aln2[i][current_row-1])
            current_row -= 1
        elif current_value == '-':
            for i in range(len_aln1):
                aligned_seqs1[i].append(aln1[i][current_col-1])
            for i in range(len_aln2):
                aligned_seqs2[i].append('-')
            current_col -= 1
        elif current_value == None:
            break
        else:
            raise ValueError, "Invalid value in traceback matrix: %s" % current_value
    
    for i in range(len_aln1):
        aligned_seqs1[i] = ''.join(aligned_seqs1[i][::-1])
    for i in range(len_aln2):
        aligned_seqs2[i] = ''.join(aligned_seqs2[i][::-1])
    
    return aligned_seqs1, aligned_seqs2, best_score

def msa_nw_align(aln1, aln2, gap_open_penalty=8, gap_extend_penalty=1, substitution_matrix=nt_substitution_matrix):
    """ Perform Needleman-Wunsch alignment of seq1 and seq2
    """
    nw_matrix, traceback_matrix = msa_generate_nw_and_traceback_matrices(
                                    aln1, aln2, gap_open_penalty, gap_extend_penalty, substitution_matrix)
    aligned_seq1, aligned_seq2, score = msa_nw_traceback(traceback_matrix,nw_matrix,aln1,aln2)
    return aligned_seq1, aligned_seq2, score

def progressive_msa(query_sequences, guide_tree=None, gap_open_penalty=8, gap_extend_penalty=1, 
                    substitution_matrix=nt_substitution_matrix):
    if guide_tree == None:
        # create a guide tree if one was not provided
        guide_tree = guide_tree_from_query_sequences(query_sequences,display_tree=False)
        
    left = guide_tree.get_left()
    if left.is_leaf():
        left_id = left.get_id()
        left_aln = [query_sequences[left_id]]
    else:
        left_aln, _ = progressive_msa(query_sequences, left,
                                   gap_open_penalty, gap_extend_penalty, substitution_matrix)
    
    right = guide_tree.get_right()
    if right.is_leaf():
        right_id = right.get_id()
        right_aln = [query_sequences[right_id]]
    else:
        right_aln, _ = progressive_msa(query_sequences, right,
                                    gap_open_penalty, gap_extend_penalty, substitution_matrix)
    
    aln1_ids = [s[0] for s in left_aln]
    aln1_seqs = [s[1] for s in left_aln]
    aln2_ids = [s[0] for s in right_aln]
    aln2_seqs = [s[1] for s in right_aln]
    aln1, aln2, score = msa_nw_align([s[1] for s in left_aln], 
                                 [s[1] for s in right_aln], 
                                 gap_open_penalty, gap_extend_penalty, substitution_matrix)
    msa = zip(aln1_ids, aln1) + zip(aln2_ids, aln2)
    msa.sort()
    return msa, guide_tree

def compute_aligned_sequence_distances(seqs, distance_fn=hamming_distance):
    dm = []
    ids = []
    for id1, seq1 in seqs:
        ids.append(id1)
        row = []
        for id2, seq2 in seqs:
            row.append(hamming_distance(seq1, seq2))
        dm.append(row)
    return SymmetricDistanceMatrix(dm, ids)

def progressive_msa_and_tree(query_sequences, gap_open_penalty=8, gap_extend_penalty=1, 
                             substitution_matrix=nt_substitution_matrix,
                             msa_distance_fn=compute_aligned_sequence_distances,
                             guide_tree=None, display_aln=False, display_tree=False):
    msa, guide_tree = progressive_msa(query_sequences, guide_tree, 
                                      gap_open_penalty, gap_extend_penalty, substitution_matrix)
    if display_aln:
        print "Multiple sequence alignment:\n"
        for seq_id, seq in msa:
            print seq, "(%s)" % seq_id
    
    dm = msa_distance_fn(msa)
    lm = average(dm.condensed_form())
    tree = to_tree(lm)
    if display_tree:
        print "\nOutput tree:"
        d = dendrogram(lm, labels=dm.ids, orientation='right', 
                   link_color_func=lambda x: 'black', leaf_font_size=24)
    return msa, tree

def iterative_msa_and_tree(query_sequences, num_iterations, 
                           gap_open_penalty=8, gap_extend_penalty=1, 
                           substitution_matrix=nt_substitution_matrix, 
                           guide_tree_fn=guide_tree_from_query_sequences,
                           msa_distance_fn=compute_aligned_sequence_distances,
                           guide_tree=None, display_aln=False, display_tree=False):
    if num_iterations > 5:
        raise ValueError("A maximum of five iterations is allowed. You requested %d." % num_iterations)
    previous_iter_tree = None
    for i in range(num_iterations):
        if i == (num_iterations - 1):
            # only display the last iteration
            display_iter = True
        else:
            display_iter = False
        previous_iter_msa, previous_iter_tree = progressive_msa_and_tree(query_sequences, 
            gap_open_penalty=gap_open_penalty, gap_extend_penalty=gap_extend_penalty, 
            substitution_matrix=substitution_matrix, msa_distance_fn=msa_distance_fn,
            guide_tree=previous_iter_tree, display_aln=display_aln and display_iter,
            display_tree=display_tree and display_iter)
                                                                    
    return previous_iter_msa, previous_iter_tree
