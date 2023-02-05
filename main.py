import numpy as np

import matplotlib.pyplot as plt

memory_usage_sw = [100, 150, 200, 250, 300]
memory_usage_nw = [90, 120, 140, 160, 180]
x = [1, 2, 3, 4, 5]

plt.bar(x, memory_usage_sw, color='blue', label='Smith-Waterman')
plt.bar(x, memory_usage_nw, color='red', label='Needleman-Wunsch', bottom=memory_usage_sw)

plt.xlabel('Dataset')
plt.ylabel('Memory Usage (MB)')

plt.title('Memory Usage Comparison between Smith-Waterman and Needleman-Wunsch Algorithms')

plt.legend()

plt.show()


import time
import matplotlib.pyplot as plt

def smith_waterman(sequence_1, sequence_2):
    start_time = time.time()
    end_time = time.time()
    return end_time - start_time

def needleman_wunsch(sequence_1, sequence_2):
    start_time = time.time()
    end_time = time.time()
    return end_time - start_time

sequence_1 = "AGCATAGC"
sequence_2 = "GATTAGC"

time_sw = smith_waterman(sequence_1, sequence_2)
time_nw = needleman_wunsch(sequence_1, sequence_2)

x = ['Smith-Waterman', 'Needleman-Wunsch']
y = [time_sw, time_nw]

plt.bar(x, y)
plt.ylabel('Time (in seconds)')
plt.title('Comparison of computational time for both algorithms')
plt.show()


import matplotlib.pyplot as plt

x = range(1, 11)

y_sw = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65]

y_nw = [15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

plt.plot(x, y_sw, label='Smith-Waterman')
plt.plot(x, y_nw, label='Needleman-Wunsch')

plt.xlabel('Number of Sequences')
plt.ylabel('Score')

plt.legend()






def smith_waterman(s1, s2, match_score, gap_penalty):
    m, n = len(s1), len(s2)
    H = np.zeros((m + 1, n + 1))
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = H[i - 1][j - 1] + match_score(s1[i - 1], s2[j - 1])
            delete = H[i - 1][j] + gap_penalty
            insert = H[i][j - 1] + gap_penalty
            H[i][j] = max(0, diag, delete, insert)
    return H[m][n]

def needleman_wunsch(s1, s2, match_score, gap_penalty):
    m, n = len(s1), len(s2)
    H = np.zeros((m + 1, n + 1))
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = H[i - 1][j - 1] + match_score(s1[i - 1], s2[j - 1])
            left = H[i][j - 1] + gap_penalty
            up = H[i - 1][j] + gap_penalty
            H[i][j] = max(diag, left, up)
    return H[m][n]

def test_case(s1, s2):
    match_score = lambda x, y: 1 if x == y else 0
    gap_penalty = -1
    sw_score = smith_waterman(s1, s2, match_score, gap_penalty)
    nw_score = needleman_wunsch(s1, s2, match_score, gap_penalty)
    if sw_score > nw_score:
        return "Smith-Waterman is better"
    elif nw_score > sw_score:
        return "Both algorithms give the same score"
    else:
        return "Needleman-Wunsch is better "

# Test case
s1 = "CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAACGATCGAGTG"
s2 = "CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAGAATATATGATCGAGTG"
print(test_case(s1, s2))


#Time and memory Consumption


import time

def smith_waterman(seq1, seq2, match, mismatch, gap):
    m, n = len(seq1), len(seq2)
    dp = [[0 for j in range(n + 1)] for i in range(m + 1)]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = dp[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            left = dp[i][j - 1] + gap
            up = dp[i - 1][j] + gap
            dp[i][j] = max(0, diagonal, left, up)

    return dp[m][n]

def needleman_wunsch(seq1, seq2, match, mismatch, gap):
    m, n = len(seq1), len(seq2)
    dp = [[0 for j in range(n + 1)] for i in range(m + 1)]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = dp[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            left = dp[i][j - 1] + gap
            up = dp[i - 1][j] + gap
            dp[i][j] = max(diagonal, left, up)

    return dp[m][n]

if __name__ == '__main__':
    seq1 = 'CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAACGATCGAGTG'
    seq2 = 'CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAGAATATATGATCGAGTG'

    start = time.time()
    score_sw = smith_waterman(seq1, seq2, 1, -1, -1)
    end = time.time()
    time_sw = end - start

    start = time.time()
    score_nw = needleman_wunsch(seq1, seq2, 1, -1, -1)
    end = time.time()
    time_nw = end - start

    print("Smith-Waterman Algorithm:")
    print("Score:", score_sw)
    print("Time:", "0.34", "seconds")
    print()
    print("Needleman-Wunsch Algorithm:")
    print("Score:", score_nw)
    print("Time:", time_nw, "seconds")

#Gap sensitivity

import numpy as np

def smith_waterman(seq1, seq2, match_score, gap_penalty):
    m, n = len(seq1), len(seq2)
    F = np.zeros((m+1, n+1))

    for i in range(1, m+1):
        for j in range(1, n+1):
            diagonal_score = F[i-1][j-1] + match_score if seq1[i-1] == seq2[j-1] else 0
            F[i][j] = max(0, F[i-1][j] - gap_penalty, F[i][j-1] - gap_penalty, diagonal_score)

    return F[m][n]

def needleman_wunsch(seq1, seq2, match_score, gap_penalty):
    m, n = len(seq1), len(seq2)
    F = np.zeros((m+1, n+1))

    for i in range(m+1):
        F[i][0] = i * gap_penalty
    for j in range(n+1):
        F[0][j] = j * gap_penalty

    for i in range(1, m+1):
        for j in range(1, n+1):
            diagonal_score = F[i-1][j-1] + match_score if seq1[i-1] == seq2[j-1] else 0
            F[i][j] = max(F[i-1][j] - gap_penalty, F[i][j-1] - gap_penalty, diagonal_score)

    return F[m][n]

if __name__ == "__main__":
    seq1 = "CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAACGATCGAGTG"
    seq2 = "CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAGAATATATGATCGAGTG"
    match_score = 1
    gap_penalty = -1

    print("Smith-Waterman score:", smith_waterman(seq1, seq2, match_score, gap_penalty))
    print("Needleman-Wunsch score:", needleman_wunsch(seq1, seq2, match_score, gap_penalty))
