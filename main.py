import numpy as np

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
        return "Needleman-Wunsch is better"
    else:
        return "Both algorithms give the same score"

# Test case
s1 = "AGC"
s2 = "AAC"
print(test_case(s1, s2))
