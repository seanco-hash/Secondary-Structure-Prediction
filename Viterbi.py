import numpy as np

EPSILON = 0.0000000000001
LOG_EPSILON = np.log(EPSILON)

# amino acids emissions for the viterbi algorithm
AA_DICT = {"A": 0,"R" : 1,"N" :2 ,"D" :3 ,"C" : 4,"Q" : 5,"E" : 6,"G" : 7,
           "H" : 8,"I" : 9,"L" : 10,"K" : 11,"M" : 12,"F" : 13,"P" : 14,
           "S" : 15,"T" : 16,"W" : 17,"Y" : 18,"V" : 19}


#state for the viterbi algorithm
# STATES_DICT = {"O": 0, "H": 1, "E": 2 }
STATES_DICT = {"O": 0, "H1": 1, "H2": 2, "H3": 3, "H4": 4, "E1": 5, "E2": 6, "E3": 7}


# def viterbi_round( e, tau, seq):
#     """
#     this is a function that runs the viterbi algorithm for a sequence, given numpy arrays
#     :param e:  emission table
#     :param tau:  transition table
#     :param seq: sequence
#     :return: the retraced sequence
#     """
#     rows, cols = e.shape[0], len(seq) + 1  # todo option
#
#     table = np.full((rows, cols), LOG_EPSILON) # initialize viterbi matrix
#     trace = np.zeros((rows, cols), dtype=int) # initialize viterbi trace
#
#     # initialize starting states
#     table[0][0] = np.log(1)
#
#     for i in range(1, cols):
#         for j in range(rows):
#             # calculate the transition and save the trace for backtracking
#             trans_vec = np.add(table[:, i-1], tau[:, j])
#             table[j][i] = np.max(trans_vec) + e[j, AA_DICT[seq[i-1]]]
#             trace[j][i] = np.argmax(trans_vec)
#     score_index = np.argmax(table[:, - 1])
#     retraced = retrace_seq(trace, seq, score_index)
#
#     # #  print the table and trace
#     # for i in range (len(table[1])):
#     #     print table[:,i]
#     # for i in range (len(trace[1])):
#     #     print trace[:,i]
#
#     return retraced

def viterbi_round( e, tau, seq):
    """
    this is a function that runs the viterbi algorithm for a sequence, given numpy arrays
    :param e:  emission table
    :param tau:  transition table
    :param seq: sequence
    :return: the retraced sequence
    """
    rows, cols = e.shape[0], len(seq) # todo option

    table = np.zeros([rows, cols], dtype=float) # initialize viterbi matrix
    trace = np.zeros([rows, cols], dtype=int) # initialize viterbi trace

    for i in range(1, rows):
        table[i, 0] = LOG_EPSILON

    for i in range(1, cols):
        for j in range(0, rows):
            # calculate the transition and save the trace for backtracking
            lastCol = table[:, i-1]
            maxIndex, maxValue = max(enumerate(lastCol + tau[:, j]), key=lambda t:t[1])
            table[j, i] = maxValue + e[j, AA_DICT[seq[i]]]
            trace[j, i] = maxIndex

    i_row, j_col = np.argmax(table[:, - 1]), cols-1
    result = ["O"]
    while j_col > 0:
        if trace[i_row, j_col] == 2:
            result.append('E')
        elif trace[i_row, j_col] == 1:
            result.append('H')
        else:
            result.append('O')
        i_row = trace[i_row, j_col]
        j_col -= 1
    result.reverse()
    result = "".join(result)

    return result


def retrace_seq( trace, seq, score_index ):
    """
    this is a function that retraces a sequence from the trace of the viterbi algorithm
    :param trace: viterbi trace table
    :param seq: string sequence
    :param score_index: index of the maximal score in last column
    :return: string sequence
    """
    #
    classification = []
    prev, counter, last_col = score_index, 0, len(seq)
    rev_dict = reverse_dict(STATES_DICT)
    for i in range(last_col, 0, -1):
        x = str(trace[prev, i])
        classification.append(rev_dict[x])
        prev = trace[prev][i-1]
    return "".join(classification[::-1])


def reverse_dict(dict):
    """
    replaces the keys by the values and vice versa
    :param dict: the priginal dictionary
    :return: the reversed one
    """
    new_dict = {}
    for key, val in dict.items():
        # todo: string number?
        new_dict[str(val)] = key
    return new_dict

# # examplle main
# #
# if __name__ == '__main__':
#     # for testing
#     n = np.full((3,20), 0.05)
#     n[1,4] = 0.09
#     n[1,7] = 0.01
#
#     n[2,9] = 0.075
#     n[2,10] = 0.025
#
#     m = np.full((3,3), 0.333333)
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
#
#
#     print(viterbi_round(np.log(n), np.log(m), "AQCHHHVCCCCCCCCCCCCCCCCCCCCVVYWWWWWTTTVAQQQCHH"))
#
#
