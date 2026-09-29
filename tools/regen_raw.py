"""Replace one raw-string EXPECTED_* block in a tests/*Expected.h from a Java dump.

usage: regen_raw.py <header.h> <NAME> <dump.txt> [--dry]

The companion of regen_expected.py, for the blocks written as

    inline const char *NAME = R"JV(...)JV";

rather than as an array of lines. The dump is what the matching Java program in
Beam43 writes -- see tools/README.md. Reports how many lines changed and shows
the first few differences.

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
    new_text = f.read().replace("\r\n", "\n")

with open(hdr_path, encoding="utf-8", newline="") as f:
    src = f.read()
nl = "\r\n" if "\r\n" in src else "\n"

opener = f'inline const char *{name} = R"JV('
if opener not in src:
    opener = f'inline const char *const {name} = R"JV('
assert src.count(opener) == 1, f"{name} not found once"
start = src.index(opener) + len(opener)
end = src.index(')JV";', start)
old_text = src[start:end]

# The header stores the block with the file's own line endings; compare on \n.
old_lines = old_text.replace("\r\n", "\n").split("\n")
new_lines = new_text.split("\n")
changed = sum(1 for a, b in zip(old_lines, new_lines) if a != b) + abs(
    len(old_lines) - len(new_lines))
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
    body = new_text.replace("\n", nl) if nl != "\n" else new_text
    with open(hdr_path, "w", encoding="utf-8", newline="") as f:
        f.write(src[:start] + body + src[end:])
    print("written" if changed else "written (no change)")
