from numpy import full, log,  argmax
from scipy.misc import logsumexp
from functools import reduce

EPSILON = 1E-12
LOG_EPSILON = log(EPSILON)
# amino acids emissions for the viterbi algorithm

AA_DICT = {"A" : 0,"R" : 1,"N" :2 ,"D" :3 ,"C" : 4,"Q" : 5,"E" : 6,"G" : 7,"H" : 8,
           "I" : 9,"L" : 10,"K" : 11,"M" : 12,
           "F" : 13,"P" : 14,"S" : 15,"T" : 16,"W" : 17,"Y" : 18,"V" : 19}

#state for the viterbi algorithm
STATES_DICT = {0 :"O", 1: "H", 2:"H", 3:"H", 4:"H", 5:"E", 6:"E", 7:"E"}


def backward(seq, tau, e):
    """
    The function runs the backward algorithm
    :param seq: the protein sequence to analyze
    :param tau: the transitions matrix
    :param e: the emissions matrix
    :return: the backward table filled with the results
    """
    rows, columns = e.shape[0], len(seq) + 1
    backward_run_mat = full((rows, columns), fill_value=LOG_EPSILON)
    backward_run_mat[0, columns - 1] = 0
    for i in range(columns - 2, -1, -1):
        for j in range(rows):
            temp = backward_run_mat[:, i + 1] + tau[j, :] + e[j, AA_DICT[seq[i]]]
            backward_run_mat[j][i] = logsumexp(temp)
    return backward_run_mat


def forward(seq, tau, e):
    """
    The function runs the forward algorithm
    :param seq: the protein sequence to analyze
    :param tau: the transitions matrix
    :param e: the emissions matrix
    :return: the forward table filled with the results
    """
    rows, columns = e.shape[0], len(seq) + 1
    forward_run_mat = full((rows, columns), fill_value=LOG_EPSILON, dtype=float)
    forward_run_mat[0][0] = 0
    for i in range(1, columns):
        for j in range(rows):
            temp = forward_run_mat[:, i - 1] + tau[:, j] + e[j][
                AA_DICT[seq[i - 1]]]
            forward_run_mat[j][i] = logsumexp(temp)
    return forward_run_mat


def backtrack_sec_track(forward_table, backward_table):
    """
    backtracking the secondary structure from the posterior, calculated using
    the forward and backward algorithm
    :param forward_table: the forward table of the protein
    :param backward_table: the backward table of the protein
    :return: the sequence of the structures
    """
    backtrack = ""
    cols = forward_table.shape[1]
    rows = forward_table.shape[0]
    for j in range(1, cols):  # from 1 - without state 0
        # for each column, choose the highest posterior of that column (iterate through the rows)
        s = argmax([forward_table[i, j] + backward_table[i, j]  for i
                    in range(rows)])
        backtrack += STATES_DICT[s]
    return backtrack


def fb_round(e, tau, seq):
    f = forward(seq, tau, e)
    b = backward(seq, tau, e)
    # print (f)
    # print (b)
    return f, b, backtrack_sec_track(f, b)

# example main
# if __name__ == '__main__':
#     n = full((3,20), 0.05)
#     n[1,4] = 0.09
#     n[1,7] = 0.01
#
#     n[2,9] = 0.075
#     n[2,10] = 0.025
#
#     m = full((3,3), 0.333333)
#     m[0,1] = 0.5
#     m[0,0] = 0.45
#     m[0,2] = 0.05
#     m[1,1] = 0.3333
#     m[1,2] = 0.33333
#     m[1,0] = 0.333333333
#
#     m[2,2] = 0.6666667
#     m[2,1] = 0.33333333
#     m[2,0] = EPSILON
#     print m
#     print n
#     seq = calc_state_seq(log(n), log(m), "AQCHHHVCCCCCCCCCCCCCCCCCCCCVVYWWWWWTTTVAQQQCHH")
#     print (seq)
#
#
