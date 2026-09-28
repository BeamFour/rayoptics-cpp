"""Compare this repository's licence headers with the Java ones they port.

Read-only. Prints, per C++ file, the Java classes it says it ports, their header
categories, and the C++ file's own category, flagging disagreements. AGENTS.md in
Beam43 requires the two sets of attributions to stay in step; this is the check.

usage: licence_audit.py [--all] [--java <Beam43 java source root>]

  --all   list every file, not only the disagreements
  --java  defaults to ../../Beam43/rayoptics/src/main/java, relative to this repo

// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
"""
import os
import re
import sys
from pathlib import Path

CPP_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_JAVA_ROOT = CPP_ROOT.parent / "Beam43" / "rayoptics" / "src" / "main" / "java"

args = sys.argv[1:]
show_all = "--all" in args
java_root = Path(args[args.index("--java") + 1]) if "--java" in args else DEFAULT_JAVA_ROOT
if not java_root.is_dir():
    sys.exit(f"Java source root not found: {java_root}\nPass --java <path>.")


def java_category(text):
    head = text[:1500]
    cats = []
    if "Portions Copyright 2017-2025 Michael J. Hayford" in head:
        cats.append("hayford-portions")
    elif "Copyright 2017-2025 Michael J. Hayford" in head:
        cats.append("hayford")
    if "ported from Goptical" in head:
        cats.append("goptical")
    if "This code is part of Beam42 project" in head:
        cats.append("beam42")
    if "Ported from Minpack" in head:
        cats.append("minpack")
    before_pkg = head.split("package ", 1)[0]
    if not cats:
        cats.append("none" if not before_pkg.strip() else "other:" + before_pkg.strip().splitlines()[0][:60])
    return "+".join(cats)


def cpp_category(text):
    lines = []
    for l in text.splitlines():
        if l.startswith("//") or not l.strip():
            lines.append(l)
            if len(lines) > 12:
                break
        else:
            break
    head = "\n".join(lines)
    cats = []
    if "Portions Copyright 2017-2025 Michael J. Hayford" in head:
        cats.append("hayford-portions")
    elif "Copyright 2017-2025 Michael J. Hayford" in head:
        cats.append("hayford")
    if "ported from Goptical" in head:
        cats.append("goptical")
    if "Portions derived from Goptical" in head:
        cats.append("goptical-portions")
    if "This code is part of Beam42 project" in head:
        cats.append("beam42")
    if "Minpack" in head or "minpack" in head:
        cats.append("minpack")
    return "+".join(cats) if cats else "none"


# Java classes: simple name -> [(fully qualified name, category)]; top-level types only.
java = {}
dupes = set()
for dp, _, fns in os.walk(java_root):
    for fn in fns:
        if not fn.endswith(".java"):
            continue
        p = Path(dp) / fn
        t = p.read_text(encoding="utf-8", errors="replace")
        simple = fn[:-5]
        fq = str(p.relative_to(java_root))[:-5].replace(os.sep, ".")
        if simple in java:
            dupes.add(simple)
        java.setdefault(simple, []).append((fq, java_category(t)))

PORT_RE = re.compile(r"(?:C\+\+ )?[Pp]ort(?:ed)? of (.+)")

cpp_files = []
for sub in ("include", "src"):
    for dp, _, fns in os.walk(CPP_ROOT / sub):
        if "third_party" in dp:
            continue
        for fn in fns:
            if fn.endswith((".h", ".cpp")):
                cpp_files.append(Path(dp) / fn)

rows = []
for p in sorted(cpp_files):
    t = p.read_text(encoding="utf-8", errors="replace")
    rel = str(p.relative_to(CPP_ROOT)).replace(os.sep, "/")
    header_lines = []
    for l in t.splitlines()[:15]:
        if l.startswith("//"):
            header_lines.append(l)
        elif l.strip():
            break
    classes = []
    for l in header_lines:
        m = PORT_RE.search(l)
        if m:
            for tok in re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", m.group(1)):
                simple = tok.split(".")[-1]
                if simple in java and simple not in classes:
                    classes.append(simple)
    base = p.stem
    if not classes and base in java:
        classes.append(base)
    jcats = [(c, fq, cat) for c in classes for fq, cat in java[c]]
    rows.append((rel, cpp_category(t), jcats, header_lines[:1]))

flagged = 0
for rel, ccat, jcats, first in rows:
    want = sorted({cat for _, _, cat in jcats})
    mismatch = not jcats or set(want) != {ccat}
    if mismatch:
        flagged += 1
    if show_all or mismatch:
        js = ", ".join(f"{c}[{cat}]" for c, fq, cat in jcats) or "-- no Java class found"
        print(f"{rel}\n    cpp={ccat}   java: {js}\n    first: {first[0] if first else ''}")

print(f"\n{len(rows)} C++ files, {flagged} to look at.")
print("A file with no Java class of the same name is usually C++-only, and fine.")
if dupes:
    print("Duplicate Java simple names:", sorted(dupes))
