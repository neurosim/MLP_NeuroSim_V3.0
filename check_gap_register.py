#!/usr/bin/env python3
"""Keep the gap register's status tags honest.

The 24 August R&D audit is recorded in two artifacts that drifted apart:

  rd-audit-2026-08-24/GAP_REGISTER.md   one entry per finding, with a status tag
  RD_AUDIT_AND_ROADMAP_2026-08-24.md    the implementation log, table-per-session

The register is the document a reader opens first, and its tags went stale the
moment work started landing: findings closed in code and recorded in the log kept
their [UNVERIFIED] tag for days. That is not cosmetic. The roadmap says so in its
own words, at "I reported the CVC items wrongly for several turns": a whole
reporting pass was spent re-deriving closures that were already recorded, because
a status field in a register is not the same artifact as an implementation log.

This gate makes the two agree by construction. For every finding it collects
independent closure evidence:

  log    a closure-table row in the implementation log whose verdict cell says
         Fixed / Partly fixed / Documented / Refuted / Not a defect
  code   the finding id appearing in the source tree, which is the convention
         every fix in this campaign follows (a `CHO-16:` comment at the fix site
         and the id in the regression test's name)

and fails when the register's tag contradicts that evidence in either direction:
a finding tagged open with evidence of closure, or tagged closed with none.

Run standalone, or as one gate of tools/check_ntn_standards.py.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REGISTER = ROOT / "rd-audit-2026-08-24" / "GAP_REGISTER.md"
LOG = ROOT / "RD_AUDIT_AND_ROADMAP_2026-08-24.md"

# The finding-id prefixes the audit actually used. Kept explicit rather than a
# generic [A-Z]+-\d+ so that a stray "TS 38-331" or "P.676-13" in prose cannot
# be mistaken for a finding.
PREFIXES = (
    "NT", "OBS", "V2X", "SAGIN", "CHO", "CON", "SIONNA", "THZ", "BOTH",
    "AI", "TWIN", "RRC", "SLICE", "FAPI", "TEST", "WF", "ORAN", "CVC",
)
ID_RE = re.compile(r"\b(?:%s)-\d+\b" % "|".join(PREFIXES))

# "#### CHO-16 [MINOR] [CONFIRMED] [CLOSED 2026-08-27] CondEventT1 is ..."
#
# Three bracketed fields, and the third is the one this gate exists for.
# Severity is [CRITICAL|MAJOR|MINOR]. The second field is the VERIFICATION axis
# from the register's own status key -- CONFIRMED means an adversarial
# refutation pass failed to knock the finding down, REFUTED means it succeeded,
# UNVERIFIED means that pass never ran. None of those three say anything about
# whether the defect was subsequently REPAIRED, and conflating the two axes is
# precisely how the register went stale: a finding could be CONFIRMED (real) and
# fixed a week later, and the tag would not move. The third field is the
# remediation axis and is required.
ENTRY_RE = re.compile(
    r"^####\s+((?:%s)-\d+)\s+\[([A-Z]+)\]\s+\[([^\]]+)\](?:\s+\[([^\]]+)\])?\s+(.*)$"
    % "|".join(PREFIXES)
)

# Remediation tags that mean the finding is still outstanding.
OPEN_TAGS = {"OPEN", "DEFERRED", "WONTFIX"}

CLOSED_VERDICT_RE = re.compile(
    r"\b(fixed|partly fixed|documented|refuted|not a defect|no separate work|superseded)\b",
    re.IGNORECASE,
)

# Where a fix marker is allowed to live. Deliberately excludes the two audit
# documents themselves, so the register cannot certify itself.
CODE_ROOTS = ("contrib", "tools", "utils", "examples", "scratch", "distribution",
              ".github")
CODE_FILES = ("README.md", "SCOPE_AND_LIMITATIONS.md")
# Some findings (the OBS parser/replay ones) are defects in the control-center
# backend, which lives beside the ns-3 tree rather than inside it. Scanned as an
# evidence root so a fix landed there is not read as a missing fix here.
EXTRA_ROOTS = (ROOT.parent / "ntn-control-center" / "backend",
               ROOT.parent / "ntn-control-center" / "frontend" / "src")


def parse_register(text: str):
    """Yield (line_no, id, severity, tag, title) for every finding entry."""
    out = []
    for n, line in enumerate(text.splitlines(), 1):
        m = ENTRY_RE.match(line)
        if m:
            remediation = (m.group(4) or "").strip()
            out.append((n, m.group(1), m.group(2), m.group(3).strip(),
                        remediation, m.group(5)))
    return out


def log_closures(text: str):
    """Ids the implementation log records a closure verdict for.

    Only markdown table rows are read, and only their first two cells: the id
    cell and the verdict cell. Prose that merely names a finding is not a
    closure -- an id can appear in a sentence explaining why something was NOT
    done, and counting that would defeat the point of the gate.
    """
    closed = {}
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith("|"):
            continue
        cells = [c.strip() for c in s.strip("|").split("|")]
        if len(cells) < 2:
            continue
        ids = ID_RE.findall(cells[0])
        if not ids:
            continue
        if CLOSED_VERDICT_RE.search(cells[1]):
            verdict = re.sub(r"\*+", "", cells[1]).strip()
            for fid in ids:
                closed.setdefault(fid, verdict)
    return closed


def code_markers(ids):
    """Ids that appear somewhere in the source tree (a fix or test marker)."""
    paths = [str(ROOT / p) for p in CODE_ROOTS if (ROOT / p).exists()]
    paths += [str(ROOT / f) for f in CODE_FILES if (ROOT / f).exists()]
    paths += [str(p) for p in EXTRA_ROOTS if p.exists()]
    if not paths:
        return {}
    pattern = "|".join(re.escape(i) for i in sorted(ids))
    try:
        res = subprocess.run(
            ["grep", "-rIoh", "--exclude-dir=build", "--exclude-dir=.git",
             "--exclude-dir=node_modules", "--exclude-dir=__pycache__",
             "-E", pattern, *paths],
            capture_output=True, text=True, timeout=600,
        )
    except (OSError, subprocess.SubprocessError) as exc:  # pragma: no cover
        print(f"[gap-register] WARN  could not scan the tree for fix markers: {exc}")
        return {}
    found = {}
    for tok in res.stdout.split():
        if tok in ids:
            found[tok] = found.get(tok, 0) + 1
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true",
                    help="print only the verdict line")
    args = ap.parse_args()

    if not REGISTER.exists() or not LOG.exists():
        print("[gap-register] SKIP  audit artifacts not present in this tree")
        return 0

    entries = parse_register(REGISTER.read_text(encoding="utf-8"))
    if not entries:
        print("[gap-register] FAIL  no finding entries parsed from the register")
        return 1

    closures = log_closures(LOG.read_text(encoding="utf-8"))
    ids = {e[1] for e in entries}
    markers = code_markers(ids)

    stale, unsupported, untagged = [], [], []
    for line_no, fid, sev, verif, remed, title in entries:
        evidence = []
        if fid in closures:
            evidence.append("log")
        if fid in markers:
            evidence.append("code")

        if not remed:
            untagged.append((line_no, fid, sev, verif))
            continue

        head = remed.split()[0].upper()
        is_open = head in OPEN_TAGS

        if is_open and evidence:
            stale.append((line_no, fid, sev, remed, "+".join(evidence)))
        elif not is_open and not evidence:
            unsupported.append((line_no, fid, sev, remed))

    if not args.quiet:
        print(f"[gap-register] {len(entries)} findings, "
              f"{len(closures)} with a log closure row, "
              f"{len(markers)} with a fix marker in the tree")

    ok = True
    if untagged:
        ok = False
        print(f"[gap-register] FAIL  {len(untagged)} finding(s) carry no "
              f"remediation tag. Severity and verification do not say whether "
              f"the defect was repaired; add a third field, e.g. "
              f"[CLOSED <date>] or [OPEN].")
        for line_no, fid, sev, verif in untagged:
            print(f"    GAP_REGISTER.md:{line_no}  {fid} [{sev}] [{verif}]")
    if stale:
        ok = False
        print(f"[gap-register] FAIL  {len(stale)} finding(s) tagged OPEN in the "
              f"register that the implementation log or a fix marker in the tree "
              f"records as closed. Retag them; a reader trusts the register first.")
        for line_no, fid, sev, tag, ev in stale:
            print(f"    GAP_REGISTER.md:{line_no}  {fid} [{sev}] [{tag}]  "
                  f"evidence={ev}")
    if unsupported:
        ok = False
        print(f"[gap-register] FAIL  {len(unsupported)} finding(s) tagged closed "
              f"with no closure row in the log and no fix marker in the tree. "
              f"Either the evidence is missing or the tag is wrong.")
        for line_no, fid, sev, tag in unsupported:
            print(f"    GAP_REGISTER.md:{line_no}  {fid} [{sev}] [{tag}]")

    if ok:
        print("[gap-register] PASS  every register status tag agrees with the "
              "implementation log and the fix markers in the tree")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
