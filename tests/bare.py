"""Check that the package still runs on sympy alone.

Singular is optional and detected at runtime, so the examples compute either
way -- and that is what makes the browser build possible, since neither
Singular nor any other binary can be there.  A heavy backend added later keeps
that property only if it is imported where it is used rather than at the top
of a module, so both halves of it are checked here:

    python -m tests.bare

the import stays light, and the answers do not change without Singular.
"""

import json
import os
import pathlib
import subprocess
import sys

from . import snapshot

ROOT = pathlib.Path(__file__).resolve().parent.parent

#: importing these pulls in a compiler, a runtime or a solver, and would make
#: the package unusable in the browser.  Add to it as backends are added.
HEAVY = ('juliacall', 'juliapkg', 'julia', 'homotopycontinuation', 'phcpy',
         'sage', 'torch')


def import_is_light():
    """What the package drags in, measured in a process of its own."""
    code = ('import json, sys\n'
            'import handlebodies\n'
            f'print(json.dumps(sorted(set(sys.modules) & set({HEAVY!r}))))\n')
    out = subprocess.run([sys.executable, '-c', code], capture_output=True,
                         text=True, cwd=ROOT)
    if out.returncode:
        return [f'importing handlebodies failed: {out.stderr.strip()[-400:]}']
    found = json.loads(out.stdout)
    return [f'{name} is imported at module scope' for name in found]


def answers_unchanged():
    """The recorded results, recomputed with the Singular path shut off."""
    os.environ['HANDLEBODIES_NO_SINGULAR'] = '1'
    from handlebodies import algebra
    if algebra.use_singular():                        # the switch is the point
        return ['HANDLEBODIES_NO_SINGULAR did not take effect']
    return snapshot.compare(snapshot.compute(),
                            json.loads(snapshot.BASELINE.read_text()))


def main():
    problems = import_is_light()
    print(f'  import  {"ok" if not problems else "FAIL"}')
    heavy = list(problems)
    problems += answers_unchanged()
    print(f'  sympy   {"ok" if len(problems) == len(heavy) else "FAIL"}')
    print()
    for line in problems:
        print('FAIL', line)
    print(f'{len(problems)} problems' if problems
          else 'the package is light and sympy alone gives the same answers')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
