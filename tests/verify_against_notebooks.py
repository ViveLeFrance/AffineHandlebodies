"""Diff every example against the values stored in the original notebooks.
This is the migration check: it only runs while the .ipynb files are present.
Once they are gone, tests/snapshot.py is the regression test."""
import json, re, sys, time
import sympy as sp
import handlebodies as hb
from examples import ALL

def sage_expr(s):
    """Parse a Sage-printed polynomial into sympy for comparison."""
    return sp.expand(sp.sympify(s.replace('^', '**'), locals={'I': sp.I}))

def stored(nb):
    txt = ''.join(''.join(''.join(o.get('text','')) for o in c.get('outputs',[]))
                  for c in json.load(open(nb))['cells'] if c['cell_type']=='code')
    g = {}
    m = re.search(r'model fibre over .*? is \{(.*?) == 0\}', txt, re.S)
    g['fibre'] = m.group(1) if m else None
    m = re.search(r'reference fibre of rho_0 over .*? is\s*\{(.*?) == 0\}', txt, re.S)
    g['rho_fibre'] = m.group(1) if m else None
    n = re.findall(r'Critical values: (\d+)\s*\nCritical points: (\d+)', txt)
    if len(n) >= 2:
        g['pi_values'], g['pi_points'] = int(n[0][0]), int(n[0][1])
        g['rho_values'], g['rho_points'] = int(n[1][0]), int(n[1][1])
    m = re.search(r'hence has (\d+) punctures', txt)
    g['punctures'] = int(m.group(1)) if m else None
    g['ranks'] = [int(r) for r in re.findall(r'rank of M at .*? is (\d+)', txt)]
    g['smooth'] = 'The variety is smooth.' in txt
    return g

rows, failures = [], []
for name in ALL:
    t0 = time.time()
    s = hb.load(name)
    got = s.summary()
    t = time.time() - t0
    want = stored(f'{name}.ipynb')
    checks = []
    for k in ['smooth','pi_values','pi_points','rho_values','rho_points','punctures','ranks']:
        if want.get(k) is not None:
            checks.append((k, got[k] == want[k], got[k], want[k]))
    for k in ['fibre','rho_fibre']:
        if want.get(k):
            ok = sp.simplify(sage_expr(want[k]) - sp.sympify(got[k])) == 0
            checks.append((k, ok, got[k], want[k]))
    bad = [c for c in checks if not c[1]]
    rows.append((name, len(checks), len(bad), t))
    for c in bad:
        failures.append((name, c))
    print(f"{name:6} {'ok' if not bad else 'FAIL':>5}  {len(checks)-len(bad)}/{len(checks)} checks  {t:5.1f}s")

print()
if failures:
    print("FAILURES:")
    for n, (k, _, got, want) in failures:
        print(f"  {n}.{k}\n     got  {got}\n     want {want}")
else:
    print("all examples reproduce their stored notebook outputs")
print(f"\ntotal {sum(r[3] for r in rows):.1f}s for {len(rows)} examples")
