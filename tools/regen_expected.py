"""Replace one EXPECTED_* array in a tests/*Expected.h with lines from a Java dump.

usage: regen_expected.py <header.h> <ARRAY_NAME> <dump.txt> [--dry]

The dump is what the matching Java program in Beam43 writes -- see tools/README.md.
Lines are split as the C++ tests split them (trailing empty lines dropped).
Only backslash and double quote are escaped; tabs stay raw, as in the originals.
Reports how many lines changed and shows the first few differences.

// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
"""
import difflib
import sys

if len(sys.argv) < 4:
    sys.exit(__doc__.strip().splitlines()[2])

hdr_path, name, dump_path = sys.argv[1:4]
dry = "--dry" in sys.argv

with open(dump_path, encoding="utf-8", newline="") as f:
    text = f.read()
# Java writes the platform line separator; the C++ writes \n.
text = text.replace("\r\n", "\n")
new_lines = text.split("\n")
while new_lines and new_lines[-1] == "":
    new_lines.pop()
assert not any("\r" in l for l in new_lines)

with open(hdr_path, encoding="utf-8", newline="") as f:
    src = f.read()
nl = "\r\n" if "\r\n" in src else "\n"
opener = f"inline const char *const {name}[] = {{{nl}"
assert src.count(opener) == 1, f"{name} not found once"
start = src.index(opener) + len(opener)
end = src.index(f"}};{nl}", start)
body = src[start:end]


def unescape(lit):
    out, i = [], 0
    while i < len(lit):
        c = lit[i]
        if c == "\\":
            nxt = lit[i + 1]
            assert nxt in '\\"', f"unexpected escape \\{nxt}"
            out.append(nxt)
            i += 2
        else:
            out.append(c)
            i += 1
    return "".join(out)


old_lines = []
for row in body.split(nl):
    if not row:
        continue
    assert row.startswith('    "') and row.endswith('",'), row[:80]
    old_lines.append(unescape(row[5:-2]))


def escape(s):
    return s.replace("\\", "\\\\").replace('"', '\\"')


changed = sum(1 for a, b in zip(old_lines, new_lines) if a != b) + abs(len(old_lines) - len(new_lines))
print(f"{name}: old={len(old_lines)} new={len(new_lines)} lines differing={changed}")
shown = 0
for d in difflib.unified_diff(old_lines, new_lines, lineterm="", n=0):
    if d.startswith(("---", "+++")):
        continue
    print("   ", d[:160])
    shown += 1
    if shown >= 12:
        print("    ...")
        break

if not dry:
    new_body = "".join(f'    "{escape(l)}",{nl}' for l in new_lines)
    with open(hdr_path, "w", encoding="utf-8", newline="") as f:
        f.write(src[:start] + new_body + src[end:])
    print("written" if changed else "written (no change)")
