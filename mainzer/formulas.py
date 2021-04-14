import collections
import re


def formula2counter(formula):
    parsed = re.findall('([A-Z][a-z]*)(-?[0-9]*)', formula)
    cnt = collections.Counter()
    for e, n in parsed:
        n = int(n) if n else 1
        cnt[e] += n
    return cnt

def counter2formula(cnt):
    return "".join(f"{el}{c}" for el,c in cnt.items())
