import sys
import glob
import os
import numpy as np
import fb as fb
import Viterbi as vt

DATA_POSTFIX = "pdbnew.txt"
O_LOC = 0
H_LOC = 1
E_LOC = 2
INDEX_LOC = 3
EPSILON = 0.0000000000001
STATES_NAME = {"O": 0, "H1": 1, "H2": 2, "H3": 3, "H4": 4, "E1": 5, "E2": 6, "E3": 7}
HELP_DICT = {"O": 0, "H": 1, "E": 2}
# key as aa and values the num of the aa in [O, H, E] secondary structure
AA_DICT = {"A": [0, 0, 0, 0], "R": [0, 0, 0, 1], "N": [0, 0, 0, 2],
           "D": [0, 0, 0, 3], "C": [0, 0, 0, 4], "Q": [0, 0, 0, 5],
           "E": [0, 0, 0, 6], "G": [0, 0, 0, 7], "H": [0, 0, 0, 8],
           "I": [0, 0, 0, 9], "L": [0, 0, 0, 10], "K": [0, 0, 0, 11],
           "M": [0, 0, 0, 12], "F": [0, 0, 0, 13], "P": [0, 0, 0, 14],
           "S": [0, 0, 0, 15], "T": [0, 0, 0, 16], "W": [0, 0, 0, 17],
           "Y": [0, 0, 0, 18], "V": [0, 0, 0, 19]}


def getData(dirPath):

    files = glob.glob(os.path.join(dirPath, "*" + DATA_POSTFIX))
    seqList = []

    for file in files:
        # print(file)
        seq = ""
        secStruct = ""
        fileRead = open(file, "r")
        line = fileRead.readline()

        while line:
            seq += line[:1]
            secStruct += line[2:3]
            line = fileRead.readline()

        seqList.append((seq, secStruct))

    return seqList


def getParssedData(dirPath):
    """
    Get file with data and parse it into 2 sequences: The original sequence
     and the secondary structure when each ss different from H or E labeled as
    O(=other)
    :param dirPath: A path for the data directory
    :return:
    """

    files = glob.glob(os.path.join(dirPath, "*" + DATA_POSTFIX))
    seqList = []

    for file in files:
        # print(file)
        seq = ""
        secStruct = ""
        fileRead = open(file, "r")
        line = fileRead.readline()
        print (line)
        while line:
            aa = line[:1]
            ss = line[2:3]

            if aa not in AA_DICT.keys():
                line = fileRead.readline()
                continue
            if ss != "H" and ss != "E":
                ss = "O"
            seq += aa
            secStruct += ss
            line = fileRead.readline()

        seqList.append((seq, secStruct))
        fileRead.close()

    return seqList



def createEmission(seqList):
    # number of H, E, D(=everything else) appearance
    numOfH, numOfE, numOfO  = 0, 0, 0

    # update how many times each aa appear in each secondary structure
    #  and how many times each ss appears
    for seq in seqList:
        for i in range(len(seq[0])):
            if seq[0][i] not in AA_DICT.keys():
                continue
            if seq[1][i] == 'H':
                AA_DICT[seq[0][i]][H_LOC] += 1
                numOfH += 1
            elif seq[1][i] == 'E':
                AA_DICT[seq[0][i]][E_LOC] += 1
                numOfE += 1
            else:
                AA_DICT[seq[0][i]][O_LOC] += 1
                numOfO += 1
    e = np.full([8, 20], fill_value= EPSILON)
    for key in AA_DICT.keys():
        index = AA_DICT[key][INDEX_LOC]

        e[STATES_NAME["O"]][index] = float(AA_DICT[key][O_LOC]) / float(numOfO)

        x = float(AA_DICT[key][H_LOC]) / float(numOfH)
        e[STATES_NAME["H1"]][index] = x
        e[STATES_NAME["H2"]][index] = x
        e[STATES_NAME["H3"]][index] = x
        e[STATES_NAME["H4"]][index] =x
        y = float(AA_DICT[key][E_LOC]) / float(numOfE)
        e[STATES_NAME["E1"]][index] = y
        e[STATES_NAME["E2"]][index] = y
        e[STATES_NAME["E3"]][index] = y




    seqLen = numOfO + numOfH + numOfE
    print(e)
    return np.log(e), seqLen
#
# def createEmission(seqList):
#     # number of H, E, D(=everything else) appearance
#     numOfH, numOfE, numOfO  = 0, 0, 0
#
#     # update how many times each aa appear in each secondary structure
#     #  and how many times each ss appears
#     for seq in seqList:
#         for i in range(len(seq[0])):
#             if seq[0][i] not in AA_DICT.keys():
#                 continue
#             if seq[1][i] == 'H':
#                 AA_DICT[seq[0][i]][H_LOC] += 1
#                 numOfH += 1
#             elif seq[1][i] == 'E':
#                 AA_DICT[seq[0][i]][E_LOC] += 1
#                 numOfE += 1
#             else:
#                 AA_DICT[seq[0][i]][O_LOC] += 1
#                 numOfO += 1
#     e = np.full([3, 20], fill_value= EPSILON)
#     for key in AA_DICT.keys():
#         index = AA_DICT[key][INDEX_LOC]
#         e[O_LOC][index] = float(AA_DICT[key][O_LOC]) / float(numOfO)
#         e[H_LOC][index] = float(AA_DICT[key][H_LOC]) / float(numOfH)
#         e[E_LOC][index] = float(AA_DICT[key][E_LOC]) / float(numOfE)
#
#     seqLen = numOfO + numOfH + numOfE
#     print(e)
#     return np.log(e), seqLen

#
# def createTransition(seqList):
#
#     statesDict = {"O": [0, 0, 0], "H": [0, 0, 0], "E": [0, 0, 0]}
#
#     # update the transition from each state to each state when the list is
#     #  in the order [O, H, E]
#     countTrans(seqList, statesDict)
#
#     t = np.full([3, 3], fill_value= EPSILON)
#
#     for state in statesDict.keys():
#         totalSum = float(sum(statesDict[state]))
#         print(state, totalSum)
#         t[HELP_DICT[state], O_LOC] = float(statesDict[state][O_LOC]) / totalSum
#         t[HELP_DICT[state], H_LOC] = float(statesDict[state][H_LOC]) / totalSum
#         t[HELP_DICT[state], E_LOC] = float(statesDict[state][E_LOC]) / totalSum
#     print(t)
#     return np.log(t)



def createTransition(seqList, pAlpha, pBeta):

    statesDict = {"O": [0, 0, 0], "H": [0, 0, 0], "E": [0, 0, 0]}

    # update the transition from each state to each state when the list is
    #  in the order [O, H, E]
    countTrans(seqList, statesDict)

    t = np.full([8, 8], fill_value= EPSILON)

    # sets O transitions
    totalSum = float(sum(statesDict["O"]))
    t[0, STATES_NAME["O"]] = float(statesDict["O"][O_LOC]) / totalSum
    t[0, STATES_NAME["H1"]] = float(statesDict["O"][H_LOC]) / totalSum
    t[0, STATES_NAME["E1"]] = float(statesDict["O"][E_LOC]) / totalSum

    # sets Alpha transitions:
    t[STATES_NAME["H1"],STATES_NAME["H2"]] =  1
    t[STATES_NAME["H2"],STATES_NAME["H3"]] = 1 - pAlpha[0]
    t[STATES_NAME["H2"], STATES_NAME["H4"]] =  pAlpha[0]
    t[STATES_NAME["H3"], STATES_NAME["H4"]] = 1
    p2 = pAlpha[1]
    #relative transition back to O and E
    totalSum = float(sum(statesDict["H"])) - statesDict["H"][H_LOC]
    oPart = (statesDict["H"][O_LOC] /  totalSum ) * (1-p2)
    ePart =  (statesDict["H"][E_LOC] /  totalSum ) * (1-p2)
    t[STATES_NAME["H4"], STATES_NAME["H4"]] = p2
    t[STATES_NAME["H4"], STATES_NAME["O"]] = oPart
    t[STATES_NAME["H4"], STATES_NAME["E1"]] = ePart

    # set Beta transitions
    t[STATES_NAME["E1"], STATES_NAME["E2"]] = 1-pBeta[0]
    t[STATES_NAME["E1"], STATES_NAME["E3"]] = pBeta[0]
    t[STATES_NAME["E2"], STATES_NAME["E3"]] = 1

    # relative transition back to O and E
    p2 = pBeta[1]
    totalSum = float(sum(statesDict["E"])) - statesDict["E"][E_LOC]
    oPart = (statesDict["E"][O_LOC] / totalSum) * (1 - p2)
    hPart = (statesDict["O"][H_LOC] / totalSum) * (1 - p2)
    t[STATES_NAME["E3"], STATES_NAME["E3"]] = p2
    t[STATES_NAME["E3"], STATES_NAME["O"]] = oPart
    t[STATES_NAME["E3"], STATES_NAME["H1"]] = hPart

    print(t)
    return np.log(t)



def buildHistogrem(seqList):

    HLength = {0: 0, 1: 0, 2: 0, 4: 0, 5: 0, 6: 0}
    ELength = {0: 0, 1: 0, 2: 0, 4: 0, 5: 0, 6: 0}

    currLength = {"H": 0, "E": 0}

    for seq in seqList:
        for i in range(len(seq[0])):
            curr = seq[1][i]
            if curr not in currLength.keys():
                if currLength["E"] != 0:
                    updateEstructure(ELength, currLength)
                elif currLength["H"] != 0:
                    updateHstructure(HLength, currLength)
                continue
            currLength[curr] += 1

            if curr == "H":
                if currLength["E"] != 0:
                    updateEstructure(ELength, currLength)
            elif curr == "E":
                if currLength["H"] != 0:
                    updateHstructure(HLength, currLength)

                    # print(("Hlength" , HLength))
                    # print(("Elength" , ELength))


def updateHstructure(HLength, currLength):
    hCurrLen = currLength["H"]
    if hCurrLen not in HLength.keys():
        HLength[hCurrLen] = 0
    HLength[currLength["H"]] += 1
    currLength["H"] = 0


def updateEstructure(ELength, currLength):
    eCurrLen = currLength["E"]
    if eCurrLen not in ELength.keys():
        ELength[eCurrLen] = 0
    ELength[currLength["E"]] += 1
    currLength["E"] = 0


def countTrans(seqList, statesDict):
    for seq in seqList:
        curr = seq[1][0]
        if curr != "H" and curr != "E":
            curr = "O"
        for i in range(len(seq[1]) - 1):
            next = seq[1][i + 1]
            if next != "H" and next != "E":
                next = "O"
            statesDict[curr][HELP_DICT[next]] += 1
            curr = next


def getParametrs(dirPath, test):
    seqList = getParssedData(dirPath)

    e, segLen = createEmission(seqList)
    pAlpha = [0.86, 0.89]
    pBeta = [0.89, 0.78]
    tau = createTransition(seqList, pAlpha, pBeta) # todo change to real vals

    buildHistogrem(seqList)
    testseqList = getParssedData(test)
    print testseqList[0][0]
    # newTau = tau
    # newTau[1,1] -= 0.1
    # newTau[1,0] += 0.1

    #
    for i in range(100):
        cur_seq = testseqList[i][0]

        seq_viterbi = vt.viterbi_round(e, tau, cur_seq)
        seq_fb = fb.fb_round(e, tau, cur_seq)

        print("viterbi\t" + seq_viterbi)
        print("fd  \t" + seq_fb[2])
        print("real\t" + testseqList[i][1])
     #   print("seq\t\t" + testseqList[i][0])

    return testseqList, seqList, e, tau





    # for i in range(91, 300):
    #     print(i)
    #     cur_seq = seqList[i][0]
    #     seq_viterbi = vt.viterbi_round(e, tau, cur_seq)
    #     seq_fb = fb.fb_round(e, tau, cur_seq)
    #
    #     print("viterbi \n", seq_viterbi)
    #     print("forward backward\n", seq_fb)
    #     print("real\n", seqList[i][1])


    # print(segLen)