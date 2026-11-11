# scratch: gene regions = exons ∪ 100-nt flanks per transcript, merged per gene/strand
import sys, collections
tx = collections.defaultdict(list); info = {}
for l in open(sys.argv[1]):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('exon', 'pseudogenic_exon'): continue
    for p in f[8].split(';'):
        if p.startswith('Parent='):
            for t in p[7:].split(','):
                tx[t].append((int(f[3]) - 1, int(f[4]))); info[t] = (f[0], f[6], t.split('.')[0])
iv = collections.defaultdict(list)
for t, ex in tx.items():
    ex.sort(); c, s, g = info[t]
    iv[(c, s, g)] += ex + [(max(ex[0][0] - 100, 0), ex[0][0]), (ex[-1][1], ex[-1][1] + 100)]
for (c, s, g), xs in sorted(iv.items()):
    xs.sort(); m = []
    for a, b in xs:
        if m and a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
        else: m.append([a, b])
    for a, b in m: print(f"{c}\t{a}\t{b}\t{g}\t0\t{s}")
