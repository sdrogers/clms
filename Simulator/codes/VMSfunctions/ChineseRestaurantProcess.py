from random import random
import numpy as np
import math


def discrete_draw(p):
    probs = [float(z) / sum(p) for z in p]
    rv = np.random.multinomial(1, probs)
    return int(np.where(np.random.multinomial(1, probs) == 1)[0])


def Restricted_Crp(alpha, previous_ms2, len_current_ms2):
    n = len(previous_ms2)
    if previous_ms2 == []:
        return 0
    if alpha != math.inf:
        counts = []
        for i in range(0, max(previous_ms2) + 1):
            counts.append(0)
            for j in range(0, len(previous_ms2)):
                if i == previous_ms2[j]:
                    counts[i] += 1
        assign_probs = [None] * (len(counts) + 1)
        index_to_zero = previous_ms2[-(len_current_ms2):]
        for i in range(len(counts)):
            if i in index_to_zero:
                assign_probs[i] = 0
            else:
                assign_probs[i] = counts[i] / (n - 1 + alpha)
        assign_probs[-1] = alpha / (n - 1 + alpha)
        return discrete_draw(assign_probs)
    else:
        return max(previous_ms2) + 1
