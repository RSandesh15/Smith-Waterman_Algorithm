def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1):
    m, n = len(seq1), len(seq2)
    score = [[0 for j in range(n + 1)] for i in range(m + 1)]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            left = score[i][j - 1] + gap
            up = score[i - 1][j] + gap
            score[i][j] = max(0, diagonal, up, left)

    alignment = []
    i, j = m, n
    while i > 0 and j > 0:
        if score[i][j] == score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
            alignment.append((seq1[i - 1], seq2[j - 1]))
            i -= 1
            j -= 1
        elif score[i][j] == score[i][j - 1] + gap:
            alignment.append(('-', seq2[j - 1]))
            j -= 1
        elif score[i][j] == score[i - 1][j] + gap:
            alignment.append((seq1[i - 1], '-'))
            i -= 1
    alignment.reverse()

    print(score[m][n], alignment)


smith_waterman('ATCG', 'TGCA')

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_score=-1):
    m, n = len(seq1), len(seq2)
    F = [[0 for j in range(n + 1)] for i in range(m + 1)]
    for i in range(m + 1):
        F[i][0] = gap_score * i
    for j in range(n + 1):
        F[0][j] = gap_score * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = F[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = F[i - 1][j] + gap_score
            insert = F[i][j - 1] + gap_score
            F[i][j] = max(match, delete, insert)
    align1, align2 = '', ''
    i, j = m, n
    while i > 0 and j > 0:
        score = F[i][j]
        diag_score = F[i - 1][j - 1]
        up_score = F[i][j - 1]
        left_score = F[i - 1][j]
        if score == diag_score + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif score == up_score + gap_score:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1
        elif score == left_score + gap_score:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
    while i > 0:
        align1 += seq1[i - 1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j - 1]
        j -= 1
    align1 = align1[::-1]
    align2 = align2[::-1]
    print(align1, align2)

needleman_wunsch('ATGCTA', 'GCTACG')
