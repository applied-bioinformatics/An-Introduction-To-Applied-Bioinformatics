#!/usr/bin/env python

# -----------------------------------------------------------------------------
# This work is licensed under the Creative Commons
# Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
# copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/.
# -----------------------------------------------------------------------------
from __future__ import division
from random import choice, random, shuffle

import numpy as np
from scipy.cluster.hierarchy import average, dendrogram, to_tree
from skbio.sequence import BiologicalSequence
from skbio.stats.distance import DistanceMatrix
from skbio.alignment import local_pairwise_align_ssw, Alignment
from skbio import TreeNode

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

traceback_decoding = {1: '\\', 2:'|', 3: '-', -1: 'E', 0: '*'}
###
# pairwise alignment notebook
###

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

def format_dynamic_programming_matrix(seq1, seq2, matrix, cell_width=6):
    """ define a function for formatting dynamic programming matrices
    """
    lines = []

    if isinstance(seq1, Alignment):
        seq1 = str(seq1[0])
    if isinstance(seq2, Alignment):
        seq2 = str(seq2[0])
    cell_format = "%" + str(cell_width) + "s"
    line_format = cell_format * (len(seq1) + 2)
    # print seq1 (start the line with two empty strings)
    lines.append(line_format % tuple([' ',' '] + [str(s) for s in list(seq1)]))

    # iterate over the rows and print each (starting with the
    # corresponding base in sequence2)
    for row, base in zip(matrix,' ' + seq2):
        row_list = [base]
        for s in row:
            if isinstance(s, np.float):
                s = str(s)
            else:
                s = s.decode('ascii')
            row_list.append(s)
        line = line_format % tuple(row_list)
        lines.append(line)

    return '\n'.join(lines)

def format_traceback_matrix(seq1, seq2, matrix, cell_width=6):
    if isinstance(seq1, Alignment):
        seq1 = str(seq1[0])
    if isinstance(seq2, Alignment):
        seq2 = str(seq2[0])
    translated_m = np.chararray(matrix.shape)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            translated_m[i, j] = traceback_decoding[matrix[i, j]]

    return format_dynamic_programming_matrix(seq1, seq2, translated_m,
                                             cell_width)


def format_dynamic_programming_matrix_subset(seq1,seq2,matrix, cell_width=6, num_positions=10):
    """ return first num_positions x num_positions of dynamic programming matrix
    """
    lines = []
    if isinstance(seq1, Alignment):
        seq1 = str(seq1[0])
    if isinstance(seq2, Alignment):
        seq2 = str(seq2[0])
    cell_format = "%" + str(cell_width) + "s"
    line_format = cell_format * (len(seq1[:num_positions]) + 1)
    # print seq1 (start the line with two empty strings)
    lines.append(line_format % tuple([' '] + map(str,list(seq1[:num_positions]))))

    # iterate over the rows and print each (starting with the
    # corresponding base in sequence2)
    for row, base in zip(matrix[:num_positions],' ' + seq2[:num_positions]):
        lines.append(line_format % tuple([base] + map(str,row[:num_positions])))

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
            best_score = max_score_tuple(diag_score, up_score, left_score)
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
            raise ValueError("Invalid value in traceback matrix: %s" % current_value)

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
            best_score = max_score_tuple(diag_score, up_score, left_score,
                                         new_alignment_score)
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
            raise ValueError("Invalid value in traceback matrix: %s" % current_value)

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
                raise ValueError("Invalid value in traceback matrix: %s" % current_value)
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
            best_score = max_score_tuple(diag_score, up_score, left_score,
                                         new_alignment_score)
            current_row.append(best_score[0])
            current_traceback_matrix_row.append(best_score[1])
        # append the current row to the matrix
        sw_matrix.append(current_row)
        traceback_matrix.append(current_traceback_matrix_row)
    return sw_matrix, traceback_matrix

def max_score_tuple(*scores):
    return max(scores, key=lambda x: x[0])

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
       >>> print(a1)
       CGTGCCTAAGGTATGCAAG
       >>> print(a2)
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

def align(sequence1, sequence2, gap_penalty, substitution_matrix, local):
    if local:
        return sw_align(sequence1, sequence2, gap_penalty, substitution_matrix)
    else:
        return nw_align(sequence1, sequence2, gap_penalty, substitution_matrix)

##
# DB searching notebook
##

def local_alignment_search(query, reference_db, aligner=local_pairwise_align_ssw):
    best_score = 0.0
    best_match = None
    best_a1 = None
    best_a2 = None
    for seq_id, seq in reference_db:
        alignment = aligner(query, seq)
        score = alignment.score()
        if score > best_score:
            best_score = score
            best_match = seq_id
            best_a1 = str(alignment[0])
            best_a2 = str(alignment[1])
    return best_a1, best_a2, best_score, best_match

def approximated_local_alignment_search_random(
        query, reference_db, p=0.10, aligner=local_pairwise_align_ssw):
    best_score = 0.0
    best_match = None
    best_a1 = None
    best_a2 = None
    for seq_id, seq in reference_db:
        if random() < p:
            alignment = aligner(query, seq)
            score = alignment.score()
            if score > best_score:
                best_score = score
                best_match = seq_id
                best_a1 = str(alignment[0])
                best_a2 = str(alignment[1])
    return best_a1, best_a2, best_score, best_match

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq)

def approximated_local_alignment_search_gc(
        query, reference_db, reference_db_gc_contents, p=0.05,
        aligner=local_pairwise_align_ssw):
    query_gc = gc_content(query)
    best_score = 0.0
    best_match = None
    best_a1 = None
    best_a2 = None
    for seq_id, seq in reference_db:
        ref_gc = reference_db_gc_contents[seq_id]
        if ref_gc - p < query_gc < ref_gc + p:
            alignment = aligner(query, seq)
            score = alignment.score()
            if score > best_score:
                best_score = score
                best_match = seq_id
                best_a1 = str(alignment[0])
                best_a2 = str(alignment[1])
    return best_a1, best_a2, best_score, best_match

def generate_random_score_distribution(query_sequence,
                                       subject_sequence,
                                       n=99,
                                       aligner=local_pairwise_align_ssw):
    result = []
    random_sequence = list(query_sequence)
    for i in range(n):
        shuffle(random_sequence)
        alignment = aligner(random_sequence,subject_sequence)
        result.append(alignment.score())
    return result

def fraction_better_or_equivalent_alignments(query_sequence,
                                             subject_sequence,
                                             n = 99,
                                             aligner=local_pairwise_align_ssw):
    random_scores = generate_random_score_distribution(query_sequence,
                                                       subject_sequence,
                                                       n,
                                                       aligner=aligner)
    alignment = aligner(query_sequence, subject_sequence)

    count_better = 0
    for e in random_scores:
        if e >= alignment.score():
            count_better += 1

    return count_better / (n + 1)

## MSA notebook

def kmer_distance(sequence1, sequence2, k=3, overlapping=True):
    """Compute the kmer distance between a pair of sequences

    Parameters
    ----------
    sequence1 : BiologicalSequence
    sequence2 : BiologicalSequence
    k : int, optional
        The word length.
    overlapping : bool, optional
        Defines whether the k-words should be overlapping or not
        overlapping.

    Returns
    -------
    float
        Fraction of the set of k-mers from both sequence1 and
        sequence2 that are unique to either sequence1 or
        sequence2.

    Raises
    ------
    ValueError
        If k < 1.

    Notes
    -----
    k-mer counts are not incorporated in this distance metric.

    """
    sequence1_kmers = set(sequence1.k_words(k, overlapping))
    sequence2_kmers = set(sequence2.k_words(k, overlapping))
    all_kmers = sequence1_kmers | sequence2_kmers
    shared_kmers = sequence1_kmers & sequence2_kmers
    number_unique = len(all_kmers) - len(shared_kmers)
    fraction_unique = number_unique / len(all_kmers)
    return fraction_unique

def guide_tree_from_sequences(sequences,
                              distance_fn=kmer_distance,
                              display_tree = False):
    """ Build a UPGMA tree by applying distance_fn to sequences

    Parameters
    ----------
    sequences : skbio.SequenceCollection
      The sequences to be represented in the resulting guide tree.
    sequence_distance_fn : function
      Function that returns and skbio.DistanceMatrix given an
      skbio.SequenceCollection.
    display_tree : bool, optional
      Print the tree before returning.

    Returns
    -------
    skbio.TreeNode

    """
    guide_dm = sequences.distances(distance_fn)
    guide_lm = average(guide_dm.condensed_form())
    guide_tree = to_tree(guide_lm)
    if display_tree:
        guide_d = dendrogram(guide_lm, labels=guide_dm.ids, orientation='right',
               link_color_func=lambda x: 'black')
    return guide_tree

def progressive_msa(sequences, guide_tree, pairwise_aligner):
    """ Perform progressive msa of sequences

    Parameters
    ----------
    sequences : skbio.SequenceCollection
        The sequences to be aligned.
    guide_tree : skbio.TreeNode
        The tree that should be used to guide the alignment process.
    pairwise_aligner : function
        Function that should be used to perform the pairwise alignments,
        for example skbio.alignment.global_pairwise_align_nucleotide. Must
        support skbio.BiologicalSequence objects or skbio.Alignment objects
        as input.

    Returns
    -------
    skbio.Alignment

    """
    c1, c2 = guide_tree.children
    if c1.is_tip():
        c1_aln = sequences[c1.name]
    else:
        c1_aln = progressive_msa(sequences, c1, pairwise_aligner)

    if c2.is_tip():
        c2_aln = sequences[c2.name]
    else:
        c2_aln = progressive_msa(sequences, c2, pairwise_aligner)

    return pairwise_aligner(c1_aln, c2_aln)

def progressive_msa_and_tree(sequences,
                             pairwise_aligner,
                             sequence_distance_fn=kmer_distance,
                             guide_tree=None,
                             display_aln=False,
                             display_tree=False):
    """ Perform progressive msa of sequences and build a UPGMA tree
    Parameters
    ----------
    sequences : skbio.SequenceCollection
        The sequences to be aligned.
    pairwise_aligner : function
        Function that should be used to perform the pairwise alignments,
        for example skbio.Alignment.global_pairwise_align_nucleotide. Must
        support skbio.BiologicalSequence objects or skbio.Alignment objects
        as input.
    sequence_distance_fn : function
        Function that returns and skbio.DistanceMatrix given an
        skbio.SequenceCollection. This will be used to build a guide tree if
        one is not provided.
    guide_tree : skbio.TreeNode, optional
        The tree that should be used to guide the alignment process.
    display_aln : bool, optional
        Print the alignment before returning.
    display_tree : bool, optional
        Print the tree before returning.

    Returns
    -------
    skbio.alignment
    skbio.TreeNode

    """
    if guide_tree is None:
        guide_dm = sequences.distances(sequence_distance_fn)
        guide_lm = average(guide_dm.condensed_form())
        guide_tree = TreeNode.from_linkage_matrix(guide_lm, guide_dm.ids)

    msa = progressive_msa(sequences, guide_tree,
                          pairwise_aligner=pairwise_aligner)
    if display_aln:
        print(msa)

    msa_dm = msa.distances()
    msa_lm = average(msa_dm.condensed_form())
    msa_tree = TreeNode.from_linkage_matrix(msa_lm, msa_dm.ids)
    if display_tree:
        print("\nOutput tree:")
        d = dendrogram(msa_lm, labels=msa_dm.ids, orientation='right',
                   link_color_func=lambda x: 'black', leaf_font_size=24)
    return msa, msa_tree

def iterative_msa_and_tree(sequences,
                           num_iterations,
                           pairwise_aligner,
                           sequence_distance_fn=kmer_distance,
                           display_aln=False,
                           display_tree=False):
    """ Perform progressive msa of sequences and build a UPGMA tree
    Parameters
    ----------
    sequences : skbio.SequenceCollection
       The sequences to be aligned.
    num_iterations : int
       The number of iterations of progressive multiple sequence alignment to
       perform. Must be greater than zero and less than five.
    pairwise_aligner : function
       Function that should be used to perform the pairwise alignments,
       for example skbio.Alignment.global_pairwise_align_nucleotide. Must
       support skbio.BiologicalSequence objects or skbio.Alignment objects
       as input.
    sequence_distance_fn : function
       Function that returns and skbio.DistanceMatrix given an
       skbio.SequenceCollection. This will be used to build a guide tree.
    display_aln : bool, optional
       Print the alignment before returning.
    display_tree : bool, optional
       Print the tree before returning.

    Returns
    -------
    skbio.alignment
    skbio.TreeNode

   """
    if num_iterations > 5:
        raise ValueError("A maximum of five iterations is allowed."
                         "You requested %d." % num_iterations)
    previous_iter_tree = None
    for i in range(num_iterations):
        if i == (num_iterations - 1):
            # only display the last iteration
            display = True
        else:
            display = False
        previous_iter_msa, previous_iter_tree = \
            progressive_msa_and_tree(sequences,
             pairwise_aligner=pairwise_aligner,
             sequence_distance_fn=sequence_distance_fn,
             guide_tree=previous_iter_tree,
             display_aln=display_aln and display,
             display_tree=display_tree and display)

    return previous_iter_msa, previous_iter_tree
