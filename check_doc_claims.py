#!/usr/bin/env python3
"""CVC-13: fail the build when user-facing docs claim capabilities the code denies.

This gap was not one bad sentence. The claim "Real ASN.1-compiled E2AP / KPM /
RC over SCTP via FlexRIC" appeared in the docs site, the rendered site HTML, the
Docker Hub description and the main README, while three places IN THE CODE said
the opposite in as many words:

  * contrib/oran-ntn/asn1/asn1-per-codec.h    "NOT bit-conformant Aligned-PER
                                               ... a test-only / in-sim codec"
  * contrib/oran-ntn/flexric-bridge/e2-message-types.h
                                              "toolkit-internal ... don't map
                                               onto an official E2AP codepoint
                                               registry"
  * oran-ntn-repo/flexric-bridge/README.md    "no asn1c/cffi codec, no SCTP
                                               association, and no FlexRIC
                                               interop ship"

The README even contradicted itself within one page: a capability table row
claiming a live FlexRIC wire, and a paragraph fifty lines later calling that
wire a roadmap item.

A prose claim has no test, so it drifts back. This is that test. It is
deliberately narrow: it looks for specific overclaiming phrases rather than
trying to understand the docs, so a false positive is a phrase someone should
reword anyway.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

# (regex, why it is wrong, what to say instead)
BANNED = [
    (r"[Rr]eal ASN\.1[- ]compiled",
     "the shipped codec is explicitly NOT bit-conformant Aligned-PER "
     "(asn1-per-codec.h says so in an honesty note)",
     "call it a PER-style in-simulator codec"),
    (r"live FlexRIC SCTP/E2AP",
     "no SCTP association and no FlexRIC interop ship in the tree",
     "say the live wire is scaffolded, not demonstrated"),
    (r"E2AP[- ]over[- ]SCTP\b(?![^.]*\b(roadmap|scaffold|planned|not (yet )?built|no SCTP)\b)",
     "an unqualified SCTP wire claim",
     "qualify it as roadmapped or scaffolded in the same sentence"),
    # AI-09: the "multi-agent" trainers are single-agent PPO/SAC over N
    # independent HandoverEnv/BeamMgmtEnv copies. There is no centralized
    # critic, no joint observation and no inter-agent coupling, so the N agents
    # cannot see or affect each other. The module docstring conceded this
    # ("rather than re-implement MAPPO from scratch"); the docs did not.

]

# Files a reader takes capability claims from. Source comments and the audit
# register are exempt: the register QUOTES the bad phrases as evidence, and
# banning a phrase from the document that records it would be self-defeating.
DOC_GLOBS = [
    "docs/index.md",
    "site/index.html",
    "ns-3-dev/README.md",
    "ns-3-dev/distribution/docker/HUBDESCRIPTION.md",
    "ns-3-dev/distribution/docker/README.md",
]


# CVC-15: numbers in the README that MUST match committed data, and where the
# truth lives. The originals drifted apart: the README claimed 85,074 actions
# from 13 xApps while xapp_metrics.csv summed to 71,967 from 5; it put the A3
# ping-pong rate at 57% when mc_table.csv gives a3=50.23 and location=57.07,
# i.e. the two were transposed; and it described the 10-seed campaign as
# "Walker-Star, 1200 km, 53 deg" against the manuscript's 780 km / 86.4 deg,
# which is also self-contradictory since a Walker-Star shell is near-polar.
#
# A number with a source of truth can be checked. These are.
NUMERIC_CLAIMS = [
    # (file, must-contain, why)
    ("ns-3-dev/README.md", "71,967",
     "xapp_metrics.csv sums to 71,967 routed actions across 5 active xApps"),
    ("ns-3-dev/README.md", "50.2",
     "mc_table.csv gives a3 ping-pong 50.23%; 57.07% is the LOCATION baseline"),
    ("ns-3-dev/README.md", "780 km",
     "the 10-seed campaign shell, matching tab:kpi in the manuscript"),
]
STALE_CLAIMS = [
    ("ns-3-dev/README.md", "85,074",
     "superseded action count; the committed CSV sums to 71,967"),
    ("ns-3-dev/README.md", "12/12 standards gates",
     "the standards checker runs more gates than that; count them rather than "
     "quoting a stale total"),
]


# AI-09: the "multi-agent" trainers are single-agent PPO/SAC over N independent
# HandoverEnv / BeamMgmtEnv copies. No centralized critic, no joint observation,
# no inter-agent coupling: the N agents cannot see or affect each other. The
# module docstring conceded it ("rather than re-implement MAPPO from scratch");
# the docs named MAPPO/MASAC without the caveat.
#
# A forward-only lookahead cannot express this, because the natural correction
# puts the caveat BEFORE the acronym ("not MAPPO/MASAC"). So the rule is
# line-level: wherever MAPPO or MASAC appears, a caveat word must appear on the
# same line.
CAVEATED_TERMS = [
    ("MAPPO", ("not ", "single-agent", "independent")),
    ("MASAC", ("not ", "single-agent", "independent")),
]


_GROUP_SEP = re.compile(r"(?<=\d)[\u202f\u00a0\u2009 ](?=\d\d\d\b)")


def _normalize_digit_groups(text: str) -> str:
    """Treat 71,967 / 71 967 / 71\u202f967 as the same number.

    A thousands separator is typography, not data. Matching on the literal
    bytes made the gate fail a README that stated the right value with a comma,
    which trains a reader to edit the gate rather than the claim.
    """
    return _GROUP_SEP.sub(",", text)


def check_caveated_terms() -> list:
    out = []
    for rel in DOC_GLOBS:
        p = ROOT / rel
        if not p.exists():
            continue
        for lineno, line in enumerate(
            p.read_text(encoding="utf-8", errors="replace").splitlines(), 1
        ):
            for term, caveats in CAVEATED_TERMS:
                if term in line and not any(c in line for c in caveats):
                    out.append(
                        f"{rel}:{lineno}: '{term}' named without the caveat that the trainer "
                        f"is single-agent PPO/SAC over N independent environment copies\n"
                        f"    line: {line.strip()[:150]}\n"
                        f"    fix : name it alongside the limitation, on the same line")
    return out


def check_numeric_claims() -> list:
    out = []
    for rel, needle, why in NUMERIC_CLAIMS:
        p = ROOT / rel
        if not p.exists():
            continue
        if needle not in _normalize_digit_groups(
                p.read_text(encoding="utf-8", errors="replace")):
            out.append(f"{rel}: expected to state '{needle}' ({why}) and does not")
    for rel, needle, why in STALE_CLAIMS:
        p = ROOT / rel
        if not p.exists():
            continue
        text = _normalize_digit_groups(p.read_text(encoding="utf-8", errors="replace"))
        for lineno, line in enumerate(text.splitlines(), 1):
            # A line that quotes the stale value while explaining that it IS
            # stale is fine; that is how the correction is documented.
            if needle in line and "superseded" not in line and "contradicted" not in line:
                out.append(f"{rel}:{lineno}: stale claim '{needle}' ({why})")
    return out


def main() -> int:
    failures = []
    checked = 0
    for rel in DOC_GLOBS:
        p = ROOT / rel
        if not p.exists():
            continue
        checked += 1
        text = _normalize_digit_groups(p.read_text(encoding="utf-8", errors="replace"))
        for lineno, line in enumerate(text.splitlines(), 1):
            for pattern, why, instead in BANNED:
                if re.search(pattern, line):
                    failures.append(
                        f"{rel}:{lineno}: {why}\n"
                        f"    line: {line.strip()[:160]}\n"
                        f"    fix : {instead}")

    if checked == 0:
        print("[doc-claims] FAIL  no documents were checked; the path list is stale")
        return 1

    failures.extend(check_numeric_claims())
    failures.extend(check_caveated_terms())

    if failures:
        print(f"[doc-claims] FAIL  {len(failures)} claim defect(s) across {checked} documents\n")
        for f in failures:
            print("  " + f.replace("\n", "\n  "))
        return 1

    print(f"[doc-claims] PASS  {checked} user-facing documents: no capability claim the "
          f"code contradicts, and the checked numbers match the committed data")
    return 0


if __name__ == "__main__":
    sys.exit(main())
