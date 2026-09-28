"""Carry regenerated Java upstream-test expectations into the C++ copies.

Behaves like a regeneration of the Java test: afterwards every value and tolerance is
the Java's, and hand edits are re-applied by hand.

Each assertion is matched by its label. Before replacing a value, the C++ value is
checked against the Java one as it was *before* the regeneration, so a label whose
meaning drifted is reported and nothing is written. A tolerance that differs from the
old Java is a hand widening; it is reset to the Java's and listed, and its
"Manually widened" comment is removed. The Java's own comments that precede an
assertion -- its manually widened tolerances -- are then carried across.

usage: sync_upstream_tests.py --old-rev <rev> [--beam43 <path>] [--tests A,B] [--dry]

  --old-rev  Beam43 revision holding the Java expectations as the C++ currently
             matches them, i.e. the commit before the regeneration (e.g. abc1234^).
             Required: it is what makes the check meaningful.
  --beam43   defaults to ../../Beam43, relative to this repository
  --tests    comma separated test class names; by default every
             tests/upstream/*UpstreamTest.cpp in this repository
  --dry      report, write nothing

// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
"""
import re
import subprocess
import sys
from pathlib import Path

CPP = Path(__file__).resolve().parent.parent
JAVA_DIR = "rayoptics/src/test/java/org/redukti/rayoptics/upstream/"

args = sys.argv[1:]


def option(flag, default=None):
    return args[args.index(flag) + 1] if flag in args else default


DRY = "--dry" in args
OLD_REV = option("--old-rev")
if not OLD_REV:
    sys.exit(next(l for l in __doc__.splitlines() if l.startswith("usage:")))
BEAM43 = Path(option("--beam43", CPP.parent / "Beam43"))
if not (BEAM43 / ".git").exists():
    sys.exit(f"Beam43 checkout not found: {BEAM43}\nPass --beam43 <path>.")

if option("--tests"):
    NAMES = option("--tests").split(",")
else:
    NAMES = sorted(p.stem for p in (CPP / "tests" / "upstream").glob("*UpstreamTest.cpp"))
if not NAMES:
    sys.exit("no upstream tests found")

JAVA_RE = re.compile(r'assertClose\("([^"]+)",\s*([^,]+?)d,\s*[^;]*?,\s*([0-9.eE+-]+)\);')
CPP_RE = re.compile(r'assertClose\("([^"]+)",(\s*)([^,]+?),(\s*[^;]*?),(\s*)([0-9.eE+-]+)\);')


def git_show(rev, path):
    return subprocess.run(["git", "-C", str(BEAM43), "show", f"{rev}:{path}"],
                          capture_output=True, text=True, check=True).stdout


def java_map(text, with_comments=False):
    out = {}
    lines = text.splitlines()
    for i, line in enumerate(lines):
        m = JAVA_RE.search(line)
        if not m:
            continue
        label, val, tol = m.group(1), m.group(2).strip(), m.group(3)
        assert label not in out, f"duplicate label {label}"
        comments = []
        if with_comments:
            j = i - 1
            while j >= 0 and lines[j].strip().startswith("//"):
                comments.insert(0, lines[j].strip())
                j -= 1
        out[label] = (val, tol, comments)
    return out


def same_double(a, b):
    return float(a) == float(b) and str(a).startswith("-") == str(b).startswith("-")


WIDENED = "Manually widened"


def strip_widened_comments(text):
    """Drops each comment run, starting at a "Manually widened" line, that sits
    directly above an assertClose."""
    lines = text.split("\n")
    out, removed, i = [], 0, 0
    while i < len(lines):
        if lines[i].lstrip().startswith("//") and WIDENED in lines[i]:
            j = i
            while j < len(lines) and lines[j].lstrip().startswith("//"):
                j += 1
            if j < len(lines) and "assertClose(" in lines[j]:
                removed += 1
                i = j
                continue
        out.append(lines[i])
        i += 1
    return "\n".join(out), removed


total_changed = 0
for name in NAMES:
    jpath = JAVA_DIR + name + ".java"
    old = java_map(git_show(OLD_REV, jpath))
    new = java_map(git_show("HEAD", jpath), with_comments=True)
    assert old.keys() == new.keys(), f"{name}: label sets differ between Java revisions"

    cpath = CPP / "tests" / "upstream" / f"{name}.cpp"
    with open(cpath, encoding="utf-8", newline="") as f:
        src = f.read()

    seen = set()
    problems = []
    resets = []
    changed = 0

    def repl(m):
        global changed
        label, ws1, val, expr, ws2, tol = m.groups()
        if label not in old:
            problems.append(f"{label}: not in Java")
            return m.group(0)
        seen.add(label)
        oval, otol, _ = old[label]
        if not same_double(val.strip(), oval):
            problems.append(f"{label}: C++ value {val} != old Java {oval}")
            return m.group(0)
        nval, ntol, _ = new[label]
        if float(tol) != float(otol):
            resets.append(f"{label}: {tol} -> {ntol}")
        if nval != val.strip() or ntol != tol:
            changed += 1
        return f'assertClose("{label}",{ws1}{nval},{expr},{ws2}{ntol});'

    widened_before = src.count(WIDENED)
    src, _ = strip_widened_comments(src)
    out = CPP_RE.sub(repl, src)
    missing = set(old) - seen
    if missing:
        problems.append(f"{len(missing)} Java labels absent from C++, e.g. {sorted(missing)[:3]}")
    if problems:
        print(f"{name}: {len(problems)} PROBLEMS")
        for p in problems[:20]:
            print("   ", p)
        sys.exit(1)

    # Carry across the explanatory comments that precede an assertion in the Java.
    for label, (_, _, comments) in new.items():
        if not comments:
            continue
        anchor = f'assertClose("{label}",'
        assert out.count(anchor) == 1, label
        pos = out.index(anchor)
        line_start = out.rindex("\n", 0, pos) + 1
        indent = out[line_start:pos]
        # The Java's Documentation/ paths are Beam42's; this repository has no such folder.
        block = "".join(f"{indent}{c.replace('See Documentation/', 'See Beam42 Documentation/')}\n"
                        for c in comments)
        if block not in out[:line_start][-len(block) - 1:]:
            out = out[:line_start] + block + out[line_start:]

    # Net of the Java's own comments, which were stripped and then carried back.
    widened_removed = widened_before - out.count(WIDENED)
    print(f"{name}: {len(seen)} assertions matched, {changed} updated")
    if resets or widened_removed:
        print(f"    {len(resets)} hand-widened tolerance(s) reset, "
              f"{widened_removed} 'Manually widened' comment(s) removed; re-widen by hand if still needed:")
        for r in resets:
            print("     ", r)
    total_changed += changed
    if not DRY:
        with open(cpath, "w", encoding="utf-8", newline="") as f:
            f.write(out)

print(f"total updated: {total_changed}{' (dry run)' if DRY else ''}")
