#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2014, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from __future__ import division
from random import choice
from bipy.core.sequence import BiologicalSequence

blosum50 = {'A': {'A': 5.0, 'C': -1.0, 'D': -2.0, 'E': -1.0, 'F': -3.0, 'G': 0.0, 'H': -2.0, 'I': -1.0, 'K': -1.0, 'L': -2.0, 'M': -1.0, 'N': -1.0, 'P': -1.0, 'Q': -1.0, 'R': -2.0, 'S': 1.0, 'T': 0.0, 'V': 0.0, 'W': -3.0, 'Y': -2.0},
'C': {'A': -1.0, 'C': 13.0, 'D': -4.0, 'E': -3.0, 'F': -2.0, 'G': -3.0, 'H': -3.0, 'I': -2.0, 'K': -3.0, 'L': -2.0, 'M': -2.0, 'N': -2.0, 'P': -4.0, 'Q': -3.0, 'R': -4.0, 'S': -1.0, 'T': -1.0, 'V': -1.0, 'W': -5.0, 'Y': -3.0},
'D': {'A': -2.0, 'C': -4.0, 'D': 8.0, 'E': 2.0, 'F': -5.0, 'G': -1.0, 'H': -1.0, 'I': -4.0, 'K': -1.0, 'L': -4.0, 'M': -4.0, 'N': 2.0, 'P': -1.0, 'Q': 0.0, 'R': -2.0, 'S': 0.0, 'T': -1.0, 'V': -4.0, 'W': -5.0, 'Y': -3.0},
'E': {'A': -1.0, 'C': -3.0, 'D': 2.0, 'E': 6.0, 'F': -3.0, 'G': -3.0, 'H': 0.0, 'I': -4.0, 'K': 1.0, 'L': -3.0, 'M': -2.0, 'N': 0.0, 'P': -1.0, 'Q': 2.0, 'R': 0.0, 'S': -1.0, 'T': -1.0, 'V': -3.0, 'W': -3.0, 'Y': -2.0},
'F': {'A': -3.0, 'C': -2.0, 'D': -5.0, 'E': -3.0, 'F': 8.0, 'G': -4.0, 'H': -1.0, 'I': 0.0, 'K': -4.0, 'L': 1.0, 'M': 0.0, 'N': -4.0, 'P': -4.0, 'Q': -4.0, 'R': -3.0, 'S': -3.0, 'T': -2.0, 'V': -1.0, 'W': 1.0, 'Y': 4.0},
'G': {'A': 0.0, 'C': -3.0, 'D': -1.0, 'E': -3.0, 'F': -4.0, 'G': 8.0, 'H': -2.0, 'I': -4.0, 'K': -2.0, 'L': -4.0, 'M': -3.0, 'N': 0.0, 'P': -2.0, 'Q': -2.0, 'R': -3.0, 'S': 0.0, 'T': -2.0, 'V': -4.0, 'W': -3.0, 'Y': -3.0},
'H': {'A': -2.0, 'C': -3.0, 'D': -1.0, 'E': 0.0, 'F': -1.0, 'G': -2.0, 'H': 10.0, 'I': -4.0, 'K': 0.0, 'L': -3.0, 'M': -1.0, 'N': 1.0, 'P': -2.0, 'Q': 1.0, 'R': 0.0, 'S': -1.0, 'T': -2.0, 'V': -4.0, 'W': -3.0, 'Y': 2.0},
'I': {'A': -1.0, 'C': -2.0, 'D': -4.0, 'E': -4.0, 'F': 0.0, 'G': -4.0, 'H': -4.0, 'I': 5.0, 'K': -3.0, 'L': 2.0, 'M': 2.0, 'N': -3.0, 'P': -3.0, 'Q': -3.0, 'R': -4.0, 'S': -3.0, 'T': -1.0, 'V': 4.0, 'W': -3.0, 'Y': -1.0},
'K': {'A': -1.0, 'C': -3.0, 'D': -1.0, 'E': 1.0, 'F': -4.0, 'G': -2.0, 'H': 0.0, 'I': -3.0, 'K': 6.0, 'L': -3.0, 'M': -2.0, 'N': 0.0, 'P': -1.0, 'Q': 2.0, 'R': 3.0, 'S': 0.0, 'T': -1.0, 'V': -3.0, 'W': -3.0, 'Y': -2.0},
'L': {'A': -2.0, 'C': -2.0, 'D': -4.0, 'E': -3.0, 'F': 1.0, 'G': -4.0, 'H': -3.0, 'I': 2.0, 'K': -3.0, 'L': 5.0, 'M': 3.0, 'N': -4.0, 'P': -4.0, 'Q': -2.0, 'R': -3.0, 'S': -3.0, 'T': -1.0, 'V': 1.0, 'W': -2.0, 'Y': -1.0},
'M': {'A': -1.0, 'C': -2.0, 'D': -4.0, 'E': -2.0, 'F': 0.0, 'G': -3.0, 'H': -1.0, 'I': 2.0, 'K': -2.0, 'L': 3.0, 'M': 7.0, 'N': -2.0, 'P': -3.0, 'Q': 0.0, 'R': -2.0, 'S': -2.0, 'T': -1.0, 'V': 1.0, 'W': -1.0, 'Y': 0.0},
'N': {'A': -1.0, 'C': -2.0, 'D': 2.0, 'E': 0.0, 'F': -4.0, 'G': 0.0, 'H': 1.0, 'I': -3.0, 'K': 0.0, 'L': -4.0, 'M': -2.0, 'N': 7.0, 'P': -2.0, 'Q': 0.0, 'R': -1.0, 'S': 1.0, 'T': 0.0, 'V': -3.0, 'W': -4.0, 'Y': -2.0},
'P': {'A': -1.0, 'C': -4.0, 'D': -1.0, 'E': -1.0, 'F': -4.0, 'G': -2.0, 'H': -2.0, 'I': -3.0, 'K': -1.0, 'L': -4.0, 'M': -3.0, 'N': -2.0, 'P': 10.0, 'Q': -1.0, 'R': -3.0, 'S': -1.0, 'T': -1.0, 'V': -3.0, 'W': -4.0, 'Y': -3.0},
'Q': {'A': -1.0, 'C': -3.0, 'D': 0.0, 'E': 2.0, 'F': -4.0, 'G': -2.0, 'H': 1.0, 'I': -3.0, 'K': 2.0, 'L': -2.0, 'M': 0.0, 'N': 0.0, 'P': -1.0, 'Q': 7.0, 'R': 1.0, 'S': 0.0, 'T': -1.0, 'V': -3.0, 'W': -1.0, 'Y': -1.0},
'R': {'A': -2.0, 'C': -4.0, 'D': -2.0, 'E': 0.0, 'F': -3.0, 'G': -3.0, 'H': 0.0, 'I': -4.0, 'K': 3.0, 'L': -3.0, 'M': -2.0, 'N': -1.0, 'P': -3.0, 'Q': 1.0, 'R': 7.0, 'S': -1.0, 'T': -1.0, 'V': -3.0, 'W': -3.0, 'Y': -1.0},
'S': {'A': 1.0, 'C': -1.0, 'D': 0.0, 'E': -1.0, 'F': -3.0, 'G': 0.0, 'H': -1.0, 'I': -3.0, 'K': 0.0, 'L': -3.0, 'M': -2.0, 'N': 1.0, 'P': -1.0, 'Q': 0.0, 'R': -1.0, 'S': 5.0, 'T': 2.0, 'V': -2.0, 'W': -4.0, 'Y': -2.0},
'T': {'A': 0.0, 'C': -1.0, 'D': -1.0, 'E': -1.0, 'F': -2.0, 'G': -2.0, 'H': -2.0, 'I': -1.0, 'K': -1.0, 'L': -1.0, 'M': -1.0, 'N': 0.0, 'P': -1.0, 'Q': -1.0, 'R': -1.0, 'S': 2.0, 'T': 5.0, 'V': 0.0, 'W': -3.0, 'Y': -2.0},
'V': {'A': 0.0, 'C': -1.0, 'D': -4.0, 'E': -3.0, 'F': -1.0, 'G': -4.0, 'H': -4.0, 'I': 4.0, 'K': -3.0, 'L': 1.0, 'M': 1.0, 'N': -3.0, 'P': -3.0, 'Q': -3.0, 'R': -3.0, 'S': -2.0, 'T': 0.0, 'V': 5.0, 'W': -3.0, 'Y': -1.0},
'W': {'A': -3.0, 'C': -5.0, 'D': -5.0, 'E': -3.0, 'F': 1.0, 'G': -3.0, 'H': -3.0, 'I': -3.0, 'K': -3.0, 'L': -2.0, 'M': -1.0, 'N': -4.0, 'P': -4.0, 'Q': -1.0, 'R': -3.0, 'S': -4.0, 'T': -3.0, 'V': -3.0, 'W': 15.0, 'Y': 2.0},
'Y': {'A': -2.0, 'C': -3.0, 'D': -3.0, 'E': -2.0, 'F': 4.0, 'G': -3.0, 'H': 2.0, 'I': -1.0, 'K': -2.0, 'L': -1.0, 'M': 0.0, 'N': -2.0, 'P': -3.0, 'Q': -1.0, 'R': -1.0, 'S': -2.0, 'T': -2.0, 'V': -1.0, 'W': 2.0, 'Y': 8.0}}

nt_substitution_matrix = {'A': {'A':  1, 'C': -2, 'G': -2, 'T': -2, 'N': 0},
                          'C': {'A': -2, 'C':  1, 'G': -2, 'T': -2, 'N': 0},
                          'G': {'A': -2, 'C': -2, 'G':  1, 'T': -2, 'N': 0},
                          'T': {'A': -2, 'C': -2, 'G': -2, 'T':  1, 'N': 0},
                          'N': {'A':  0, 'C':  0, 'G':  0, 'T':  0, 'N': 0 }}


def hamming_distance(s1, s2):
    s1 = BiologicalSequence(s1)
    s2 = BiologicalSequence(s2)
    return s1.distance(s2)


def format_matrix(row_headers, col_headers, data, hide_zeros=False):
    result = []
    line_format = "%6s" * (len(row_headers) + 1)
    
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


def format_dynamic_programming_matrix(seq1,seq2,matrix):
    """ define a function for formatting dynamic programming matrices
    """
    lines = []

    line_format = "%6s" * (len(seq1) + 2)
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
        
    return ''.join(aligned_seq1[::-1]), ''.join(aligned_seq2[::-1]), best_score


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


def sw_align_affine_gap_nt(sequence1, sequence2, gap_open_penalty=5,
                        gap_extend_penalty=2,
                        substitution_matrix=nt_substitution_matrix):
    """
        sequence1: the first unaligned sequence
        sequence2: the second unaligned sequence
        gap_open_penalty: penalty for opening a gap (this is substracted from 
         previous best alignment score, so is typically positive)
        gap_extend_penalty: penalty for opening a gap (this is substracted from 
         previous best alignment score, so is typically positive)
        substitution_matrix: 2D-dict-like object for looking up substitution 
         scores
    """
    sw_matrix, traceback_matrix = generate_sw_and_traceback_matrices_affine_gap(sequence1,
                                                                 sequence2,
                                                                 gap_open_penalty,
                                                                 gap_extend_penalty,
                                                                 substitution_matrix)

    return sw_traceback(traceback_matrix,sw_matrix,sequence1,sequence2)

def sw_align_affine_gap_pr(sequence1, sequence2, gap_open_penalty=11,
                        gap_extend_penalty=1,
                        substitution_matrix=blosum50):
    """
        sequence1: the first unaligned sequence
        sequence2: the second unaligned sequence
        gap_open_penalty: penalty for opening a gap (this is substracted from 
         previous best alignment score, so is typically positive)
        gap_extend_penalty: penalty for opening a gap (this is substracted from 
         previous best alignment score, so is typically positive)
        substitution_matrix: 2D-dict-like object for looking up substitution 
         scores
    """
    sw_matrix, traceback_matrix = generate_sw_and_traceback_matrices_affine_gap(sequence1,
                                                                 sequence2,
                                                                 gap_open_penalty,
                                                                 gap_extend_penalty,
                                                                 substitution_matrix)

    return sw_traceback(traceback_matrix,sw_matrix,sequence1,sequence2)

###
# Convenience wrapper
###


def align(sequence1, sequence2, gap_penalty, substitution_matrix, local):
    if local:
        return sw_align(sequence1, sequence2, gap_penalty, substitution_matrix)
    else:
        return nw_align(sequence1, sequence2, gap_penalty, substitution_matrix)


