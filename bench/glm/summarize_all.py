import re, statistics as st, math, sys
from collections import defaultdict

FN = sys.argv[1] if len(sys.argv) > 1 else 'results_all_074516.tsv'
rows = []; rnd = None
for line in open(FN):
    line = line.rstrip('\n')
    m = re.match(r'### round=(\d+) build=(\S+) N=(\d+)', line)
    if m: rnd = (int(m.group(1)), m.group(2), int(m.group(3))); continue
    if 'ratio=' not in line or not rnd or line.startswith('#'): continue
    p = line.split('\t'); f = {}
    for tok in p:
        if '=' in tok: k, v = tok.split('=', 1); f[k] = v
    rows.append(dict(rnd=rnd[0], build=rnd[1], N=rnd[2], grp=p[1].strip(), op=p[2].strip(),
                     ratio=float(f['ratio'])))

BLD = ['all_O2', 'all_O3', 'all_LTO']
TAG = {'all_O2': 'gcc12 -O2', 'all_O3': 'gcc12 -O3 -march=native', 'all_LTO': 'gcc12 -O3 -march=native -flto'}
SIZES = [8192, 262144, 2097152]

aa = defaultdict(list)
for r in rows:
    if r['grp'] == 'AA': aa[(r['build'], r['N'])].append(r['ratio'])
band = {k: (min(v), max(v)) for k, v in aa.items()}

agg = defaultdict(list)
order = []
for r in rows:
    if r['grp'] != 'AA':
        agg[(r['build'], r['N'], r['grp'], r['op'])].append(r['ratio'])
        if (r['grp'], r['op']) not in order: order.append((r['grp'], r['op']))

def cell(b, N, grp, op):
    v = agg.get((b, N, grp, op))
    if not v: return None, None
    med = st.median(v); lo, hi = band.get((b, N), (0.97, 1.03))
    tag = 'FV' if med > hi else ('GL' if med < lo else '=')
    return med, tag

print('### A/A 噪声带（同一实现互比，中位数层面的本机噪声）')
for b in BLD:
    print(f'- {TAG[b]}: ' + ' | '.join(f'N={N}: {band[(b,N)][0]:.2f}–{band[(b,N)][1]:.2f}' for N in SIZES))

for grp, label in [('P', 'P 组：标量 API 逐元素'), ('B', 'B 组：batch / SoA / stream / fast'), ('S', 'S 组：小 n 批量调用开销')]:
    print(f'\n### {label}')
    print('| 运算 | ' + ' | '.join(f'{TAG[b]} N={N}' for b in ['all_O2','all_O3','all_LTO'] for N in SIZES) + ' |')
    print('|---|' + '---|' * 9)
    for g, op in order:
        if g != grp: continue
        cells = []
        for b in BLD:
            for N in SIZES:
                med, tag = cell(b, N, grp, op)
                cells.append('-' if med is None else (f'**{med:.2f}**' if tag == 'FV' else (f'_{med:.2f}_' if tag == 'GL' else f'{med:.2f}')))
        print(f'| {op} | ' + ' | '.join(cells) + ' |')
    print('（**粗体**=flashvec 快，_斜体_=GLM 快，普通=噪声内打平；比值=t_GLM/t_fv）')

# 统计
print('\n### 胜负统计（按口径）')
for grp in ['P', 'B', 'S']:
    for b in BLD:
        fv = gl = tie = tot = 0
        for g, op in order:
            if g != grp: continue
            for N in SIZES:
                med, tag = cell(b, N, grp, op)
                if med is None: continue
                tot += 1
                if tag == 'FV': fv += 1
                elif tag == 'GL': gl += 1
                else: tie += 1
        print(f'- {grp} / {TAG[b]}: flashvec 胜 {fv} · 打平 {tie} · GLM 胜 {gl} （共 {tot} 格）')

# 每组极值与几何均值
print('\n### 各组几何均值（全部格子，剔除打平后按胜出方分档）')
for grp in ['P', 'B']:
    for b in ['all_O2', 'all_O3', 'all_LTO']:
        vals = [m for (g, op) in order if g == grp for N in SIZES for m in [cell(b, N, grp, op)[0]] if m]
        if not vals: continue
        gm = math.exp(sum(math.log(max(v, 1e-6)) for v in vals) / len(vals))
        print(f'- {grp} / {TAG[b]}: 几何均值 {gm:.2f}x · 最好 {max(vals):.2f}x · 最差 {min(vals):.2f}x')
