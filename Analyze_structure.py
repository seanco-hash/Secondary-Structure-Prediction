import Viterbi,createEandT, fb
import matplotlib.pyplot as plt
from numpy import logaddexp, log, mean
from math import exp
from operator import add, sub
import sys

SEQ = 0
SOL = 1
STATES_DICT = {"O": 0, "H": 1, "E": 2}


def posterior_table(f, b):
    s1 = list(map(add, f[0], b[0]))
    s2 = list(map(add, f[1], b[1]))
    s3 = list(map(add, f[2], b[2]))
    #print(s1, s2, s3)
    sums = list(map(logaddexp, s1, s2))
    sums = list(map(logaddexp, sums, s3))
    #print(sums)
    s1 = list(map(sub, s1, sums))
    s2 = list(map(sub, s2, sums))
    s3 = list(map(sub, s3, sums))
    #print(s1, s2, s3)
    s1 = list(map(exp, s1))
    s2 = list(map(exp, s2))
    s3 = list(map(exp, s3))
    #print(s1, s2, s3)
    return s1, s2, s3


def plot_structures_probability(s1, s2, s3, size):
    x_axis = list(range(1, size + 1))
    #print( "ken s1", len(s1))
    plt.plot(x_axis, s1, label='Other')
    plt.plot(x_axis, s2, label='Alpha Helix')
    plt.plot(x_axis, s3, label='Beta Sheet')
    plt.legend()
    plt.xlabel('Position in Protein Sequence')
    plt.ylabel('Probability')
    plt.title('Structures probabilities per position in protein')
    plt.show()


def plot_error_analysis(viterbi_errors, posterior_errors,v_wins, fb_wins):
    trials = list(range(1, len(viterbi_errors) + 1))
    #print(trials)
    plt.xlabel('prediction rounds')
    plt.ylabel('error')
    avg1 = '%.3f' % mean(viterbi_errors)
    v_wins = '%.3f' % v_wins
    avg2 = '%.3f' % mean(posterior_errors)
    fb_wins = '%.3f' % fb_wins
    print("v wins", v_wins, "fb wins", fb_wins)
    v = 'viterbi error graph. average error: ' + str(avg1) + ". total wins:", v_wins
    fb = 'f-b error graph. average error: ' + str(avg2)+". total wins:", fb_wins
    plt.plot(trials, viterbi_errors, 'pink', label=v)
    plt.plot(trials, posterior_errors, 'navy', label=fb)
    plt.legend()
    plt.title('Comparison between the viterbi errors and f-b errors')
    plt.show()


def calc_fb_penalty(fb, result):
    penalty = 0
    for i in range(len(result)):
        if result[i] == "H" or result[i] == "E": # H or E in the i'th column
            logaddexp(penalty, fb[STATES_DICT[result[i]]][i])
        else:
            logaddexp(penalty, fb[0][i])  # Other in the i'th column
    return penalty


def calc_error(trial, result):
    length = len(trial)
    errors = 0
    # compare all the positions
    for i in range(length):
        if (trial[i] != result[i]):
            # check if it is O and another ones
            if trial[i] == "O" and result[i] != "H" and result[i] != "E":
                continue
            else:
                errors += 1
    return errors/length


def perform_prediction(data_to_predict, e, tau, test):
    final_predictions, final_errors = [], []
    viterbi_errors, fb_errors, fb_penalties = [], [], []
    viterbi_wins = 0
    fb_wins = 0
    for data in data_to_predict:
        seq = data[SEQ]
        # train: predict the secondary structure using viterbi and posterior
        viterbi_structure = Viterbi.viterbi_round(e, tau, seq)
        f, b, fb_structure = fb.fb_round(e, tau, seq)

        # find error in v and p
        cur_viterbi_error = calc_error(viterbi_structure, data[SOL])
        viterbi_errors.append(cur_viterbi_error)
        cur_fb_error = calc_error(fb_structure, data[SOL])
        fb_errors.append(cur_fb_error)

        # choose the path with the smallest error
        if (cur_viterbi_error < cur_fb_error):
            final_predictions.append((seq, viterbi_structure))
            final_errors.append(cur_viterbi_error)
            viterbi_wins += 1
        else:
            final_predictions.append((seq, fb_structure))
            fb_wins += 1
            final_errors.append(cur_fb_error)
        # calc the posterior
        s1, s2, s3 = posterior_table(f, b)
        # calc the penalties on the posterior (for the graph later that tomi wanted)
        fb_penalties.append(calc_fb_penalty([s1[1:], s2[1:], s3[1:]], data[SOL]))
        if test:
            # plot the probability for each state to emit E, H, ot O
            plot_structures_probability(s1[1:], s2[1:], s3[1:], len(seq))

    # plot viterbi vs. posterior by error percentage
    data_size = len(data_to_predict)
    print(data_size)
    v_total_wins = float(viterbi_wins)/float(data_size)
    fb_total_wins = 1 - v_total_wins
    plot_error_analysis(viterbi_errors, fb_errors, v_total_wins, fb_total_wins)

    # plot f-b penalties
    #todo
    # final predictions is for yulia
    #todo send to yulia


def analyze(data, test):
    test_data, train_data, e, tau = createEandT.getParametrs(data, test)
    # prediction of the train
    perform_prediction(train_data[:10], e, tau, False)
    # prediction of the test
    perform_prediction([test_data[0]], e, tau, True)
    return e,tau

#print(calc_error("HHEEOO", "HHERHH"))
#print(plot_error_analysis([0,0.1, 0.2, 0.4],[0.001, 0.12, 0.21, 0.5],0.15,0.20))
#LOG = log(0.1)
#f = [[LOG, LOG, LOG], [LOG, LOG, LOG], [LOG, LOG, LOG]]
#plot_structures_probability(f[0], f[1], f[2], len(f))

